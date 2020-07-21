#include <pcl/common/transforms.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <boost/optional.hpp>
#include <algorithm>
#include<eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <iostream>
#include <pcl/point_types.h>
#include <pcl/filters/passthrough.h>
#include <pcl/point_cloud.h>
#include <sophus/sim3.hpp>
#include <sophus/se3.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> JacobianMatrix;

class NullStream
{
    public:
    NullStream() { }
    template<typename T> NullStream& operator<<(T const&) { return *this; }
};
#define LOG(type) NullStream() 

struct IcpParameters
{
    //! Maximum number of allowed iterations
    unsigned int max_iter;
    //! Stopping condition
    /*! ICP stops when the error variation between two iteration is under
        min_variation.
        TODO: Add better convergence criteria */
    double min_variation;
    //! Maximum search distance for correspondances
    /*! Do not look further than this for the kdtree search */
    double max_correspondance_distance;

    //! Initial guess for the registration
    Eigen::Matrix4d initial_guess;

    //! Use MEstimators?
    bool mestimator;

    IcpParameters() : max_iter(10), min_variation(10e-5),
        max_correspondance_distance(std::numeric_limits<double>::max()), mestimator(false) {
        initial_guess = Eigen::Matrix<double, 4, 4>::Identity();
    }

};

// std::ostream &operator<<(std::ostream &s, const IcpParameters &p) {
//     s << "MEstimator: " << std::boolalpha << p.mestimator
//         << "\nMax iterations: " << p.max_iter
//         << "\nMin variation: " << p.min_variation
//         << "\nInitial guess (twist):\n" << p.initial_guess;
//     return s;
//     }


struct IcpResults
{
    //! History of previous registration errors
    /*!
        - First value is the initial error before ICP,
        - Last value is the final error after ICP. */
    std::vector<double> registrationError;

    //! Transformation (SE3) of the final registration transformation
    Eigen::Matrix<double, 4, 4> transformation;
    Eigen::Matrix<double, 4, 4> relativeTransformation;

    // Scale for Sim3 icp
    double scale;

    // True if ICP has converged
    bool has_converged;

    IcpResults() : transformation(Eigen::Matrix<double, 4, 4>::Identity()),
        relativeTransformation(Eigen::Matrix<double, 4, 4>::Identity()),
        scale(1.),
        has_converged(false) {
    }

    boost::optional<double> getLastErrorVariation() const {
        if (registrationError.size() >= 2) {
        return registrationError[registrationError.size() - 1] - registrationError[registrationError.size() - 2];
        } else {
        return boost::none;
        }
    }

    boost::optional<double> getLastError() const {
        if (registrationError.size() > 0) {
        return registrationError[registrationError.size() - 1];
        } else {
        return boost::none;
        }
    }

    void clear() {
        registrationError.clear();
        transformation = Eigen::Matrix<double, 4, 4>::Identity();
    }
};

// std::ostream &operator<<(std::ostream &s, const IcpResults &r) {
//     if (!r.registrationError.empty()) {
//         s << "Initial error: " << r.registrationError[0]
//         << "\nFinal error: " << r.registrationError[r.registrationError.size() - 1]
//         << "\nFinal transformation: \n"
//         << r.transformation
//         << "\nRelative transformation: \n"
//         << r.relativeTransformation
//         << "\nScale factor: " << r.scale
//         << "\nError history: ";
//         for (int i = 0; i < r.registrationError.size(); ++i) {
//         s << r.registrationError[i]  << ", ";
//         }
//     } else {
//         s << "Icp: No Results!";
//     }
//     return s;
//     }


class FixTranslationConstraint
{
protected:
typedef boost::array<bool, 3> FixedAxes;
FixedAxes fixedAxes_;

public:
    FixTranslationConstraint() {
    setFixedAxes(false, false, false);
    }

    FixTranslationConstraint(bool x, bool y, bool z)
    {
    setFixedAxes(x, y, z);
    }

    void setFixedAxes(bool x, bool y, bool z) {
    fixedAxes_[0] = x;
    fixedAxes_[1] = y;
    fixedAxes_[2] = z;
    }

    int numFixedAxes() const{
    return std::count(std::begin(fixedAxes_), std::end(fixedAxes_), true);
    }

    FixedAxes getFixedAxes() const {
    return fixedAxes_;
    }
};

class Constraints
{
    typedef typename Eigen::Matrix<double, 6, 1> Twist;

protected:
    FixTranslationConstraint translationConstraint_;


public:
    Constraints () {
    }

    void setTranslationConstraint(const FixTranslationConstraint &translationConstraint) {
    LOG(INFO) << translationConstraint.getFixedAxes()[0] << ", " << translationConstraint.getFixedAxes()[1] << ", " <<
                translationConstraint.getFixedAxes()[2];
    translationConstraint_ = translationConstraint;
    LOG(INFO) << translationConstraint_.getFixedAxes()[0] << ", " << translationConstraint_.getFixedAxes()[1] << ", " <<
                translationConstraint_.getFixedAxes()[2];
    }

    FixTranslationConstraint getTranslationConstraint() const {
    return translationConstraint_;
    }

    bool hasConstraints() const {
    return translationConstraint_.numFixedAxes() != 0;
    }

    void processJacobian(const JacobianMatrix &J, JacobianMatrix &Jconstrained) {
        Jconstrained = J;
    }


    Twist getTwist(const Eigen::Matrix<double, Eigen::Dynamic, 1> &twist) {
    return twist;
    }
};


class ErrorPointToPoint 
{
public:
    typedef typename pcl::PointCloud<pcl::PointXYZ> Pr;
    typedef typename pcl::PointCloud<pcl::PointXYZ> Pc;
    typedef typename Pc::Ptr PcPtr;
    typedef typename Pr::Ptr PrPtr;
    
    void computeWeights();
    Eigen::Matrix<double, 4, 4> update();

    double getErrorNorm() const {
        return errorVector_.norm();
        }

    //! Compute the error
    /*! \f[ e = P^* - P \f]
    *
    *  Stack the error in vectors of form
    *
    * \f[ eg = [ex_0; ey_0; ez_0; ex_1; ey_1; ez_1; ...; ex_n; ey_n; ez_n]; \f]
    */
    void computeError();

    //! Jacobian of \f$ e(x) \f$, eg \f[ J = \frac{de}{dx} \f]
    /*!
        For a 3D point of coordinates \f$ (X, Y, Z) \f$, the jacobian is
        \f[ \left( \begin{array}{cccccc}
        1  & 0  & 0  &  0  &  Z  & -Y \\
        0  & 1  & 0  & -Z  &  0  &  X \\
        0  & 0  & 1  &  Y  & -X  &  0 \\
        \end{array} \right)
        \f]

        Note:

        We update the pose on the left hand side :
        \f[ \widehat{T} \leftarrow e^x*\widehat{T} \f]
        This means that the pose jacobian is computed  at \f$ x=\widehat{x} \f$,
        Eg;
        \f[ \frac{\partial (e^x*\widehat{T}*P)}{\partial P} =
        \frac{\partial e^x*Pe}{\partial P} = [eye(3) ; skew(Pe)];
        \f]

        If the update was computed on the right hand side :
        \f[ \widehat{T} \leftarrow \widehat{T}*e^x \f]
        The pose jacobian has to be estimated at \f$ x = 0 \f$
        eg \f[ \frac{\partial (\widehat{T}*e^x*P)}{\partial x} = \widehat{T}*[eye(3) skew(P)] \f]
        */
    void computeJacobian();

    void setInputCurrent(const PcPtr &in){
        current_ = in;

        // Resize the data structures
        errorVector_.resize(3 * current_->size(), Eigen::NoChange);
        weightsVector_ = VectorX::Ones(current_->size() * 3);
        J_.setZero(3 * current_->size(), 6);
        }
    void setInputReference(const PcPtr &in){
        reference_ = in;
        }

    JacobianMatrix getJacobian() const {
    return J_;
    }
    VectorX getErrorVector() const {
    return errorVector_;
    }

void setConstraints(const boost::shared_ptr<Constraints> constraints) {
    constraints_ = constraints;
    FixTranslationConstraint translationConstraint = constraints_->getTranslationConstraint();
    LOG(INFO) << translationConstraint.getFixedAxes()[0] << ", " << translationConstraint.getFixedAxes()[1] << ", " << translationConstraint.getFixedAxes()[2];
}

protected:
    PcPtr current_;
    PcPtr reference_;

    //! Vector containing the error for each point
    VectorX errorVector_;
    VectorX weightsVector_;

    //! Corresponding Jacobian
    JacobianMatrix J_;

    //! Constraints
    boost::shared_ptr<Constraints> constraints_ = boost::make_shared <Constraints> ();
};



class IcpPointToPoint
{
public:
    typedef typename pcl::PointCloud<pcl::PointXYZ> Pr;
    typedef typename pcl::PointCloud<pcl::PointXYZ> Pc;
    typedef typename Pc::Ptr PcPtr;
    typedef typename Pr::Ptr PrPtr;

protected:
    // Reference (model) point cloud. This is the cloud that we want to register
    PcPtr P_current_;
    // kd-tree of the model point cloud
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree_;
    // Reference cloud, upon which others will be registered
    PrPtr P_ref_;
    PrPtr P_ref_init_inv;

    // Instance of an error kernel used to compute the error vector, Jacobian...
    ErrorPointToPoint err_;

    // Parameters of the algorithm (rate of convergence, stopping condition...)
    IcpParameters param_;

    // Results of the ICP
    IcpResults r_;

    unsigned int iter_;
    Eigen::Matrix<double, 4, 4> T_;

protected:
    void initialize(const PcPtr &model, const PrPtr &data,
                    const IcpParameters &param);


/**
 * @brief Finds the nearest neighbors between the current cloud (src) and the kdtree
 * (buit from the reference cloud)
 *
 * @param src
 *  The current cloud
 * @param max_correspondance_distance
 *  Max distance in which closest point has to be looked for (in meters)
 * @param indices_src
 * @param indices_target
 * @param distances
 */
    void findNearestNeighbors(const pcl::PointCloud<pcl::PointXYZ>::Ptr &src,
                            const double max_correspondance_distance,
                            std::vector<int> &indices_src,
                            std::vector<int> &indices_target,
                            std::vector<double> &distances);

    void convergenceFailed() {
    r_.has_converged = false;
    r_.transformation = Eigen::Matrix<double, 4, 4>::Identity();
    r_.relativeTransformation = Eigen::Matrix<double, 4, 4>::Identity();
    }

public:
    IcpPointToPoint() : P_current_(new Pc()), P_ref_(new Pr()), P_ref_init_inv(new Pr()), T_(Eigen::Matrix<double, 4, 4>::Identity()) {
    }

    /**
     * \brief Runs the ICP algorithm with given parameters.
     *
     * Runs the ICP according to the templated \c ErrorPointToPoint function,
     * and optimisation parameters \c IcpParameters_
     *
     * \retval void You can get a structure containing the results of the ICP (error, registered point cloud...)
     * by using \c getResults()
    **/
    void run();

    /**
     * @brief Run the next iteration of the ICP optimization
     */
    bool step();

    /**
     * @brief Sets the parameters for the optimisation.
     *
     * All parameters are defined within the \c IcpParameters_ structure.
     *
     * @param param
     *  Parameters to the minimisation
     */
    void setParameters(const IcpParameters &param) {
    param_ = param;
    //T_ = param_.initial_guess;
    // XXX Suboptimal, to fix problem with initial transformation
    if (P_ref_->size() != 0) {
        Eigen::Matrix<double, 4, 4> init_T_inv = param_.initial_guess.inverse();
        pcl::transformPointCloud(*P_ref_, *P_ref_init_inv, init_T_inv);
        kdtree_.setInputCloud(P_ref_init_inv);
        }
    }

    IcpParameters getParameters() const {
    return param_;
    }
    /** \brief Provide a pointer to the input target (e.g., the point cloud that we want to align).
    * \param[in] cloud the input point cloud target
    */
    void setInputCurrent(const PcPtr &in) {
    if (in->size() == 0) {
        LOG(WARNING) << "You are using an empty source cloud!";
    }
    P_current_ = in;
    }
    /**
     * @brief Provide a pointer to the input source (e.g., the target pointcloud
     * that we want to align to)
     *
     * @param[in] cloud	the reference point cloud source
     */
    void setInputReference(const PrPtr &in) {
    if (in->size() == 0) {
        LOG(WARNING) << "You are using an empty reference cloud!";
    }
    if (in->size() != 0) {
        P_ref_ = in;
        Eigen::Matrix<double, 4, 4> init_T_inv = param_.initial_guess.inverse();
        pcl::transformPointCloud(*in, *P_ref_init_inv, init_T_inv);
        kdtree_.setInputCloud(P_ref_init_inv);
    }
    }

    void setError(ErrorPointToPoint err) {
    err_ = err;
    }
    /**
     * @brief Gets the result of the ICP.
     *
     * @return
     * Results of the ICP (call \c run() to run the ICP and generate results)
     */
    IcpResults getResults() const {
    return r_;
    }
};

