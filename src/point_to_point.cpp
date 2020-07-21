#include "point_to_point.hpp"

using namespace Eigen;

inline void rodrigues_so3_exp(const Eigen::Matrix<double, 3, 1> &w, const double A, const double B,
                       Eigen::Matrix<double, 3, 3> &R) {
  {
    const double wx2 = (double)w(0) * w(0);
    const double wy2 = (double)w(1) * w(1);
    const double wz2 = (double)w(2) * w(2);

    R(0, 0) = 1.0 - B * (wy2 + wz2);
    R(1, 1) = 1.0 - B * (wx2 + wz2);
    R(2, 2) = 1.0 - B * (wx2 + wy2);
  }
  {
    const double a = A * w(2);
    const double b = B * (w(0) * w(1));
    R(0, 1) = b - a;
    R(1, 0) = b + a;
  }
  {
    const double a = A * w(1);
    const double b = B * (w(0) * w(2));
    R(0, 2) = b + a;
    R(2, 0) = b - a;
  }
  {
    const double a = A * w(0);
    const double b = B * (w(1) * w(2));
    R(1, 2) = b - a;
    R(2, 1) = b + a;
  }
}


inline Eigen::Matrix<double, 3, 3> skew(const Eigen::Matrix<double, 3, 1> x) {
  Eigen::Matrix<double, 3, 3> A;

  A(0, 0) = 0.0; A(0, 1) = -x(2); A(0, 2) = x(1);
  A(1, 0) = x(2); A(1, 1) = 0.0; A(1, 2) = -x(0);
  A(2, 0) = -x(1); A(2, 1) = x(0); A(2, 2) = 0.0;

  return A;
}

inline Eigen::Matrix<double, 3, 3>
expSO3(const Eigen::Matrix<double, 3, 1> w) {
  using std::sqrt;
  using std::sin;
  using std::cos;
  static const double one_6th = 1.0 / 6.0;
  static const double one_20th = 1.0 / 20.0;

  Eigen::Matrix<double, 3, 3> R;

  const double theta_sq = w.dot(w);
  const double theta = sqrt(theta_sq);
  double A, B;
  //Use a Taylor series expansion near zero. This is required for
  //accuracy, since sin t / t and (1-cos t)/t^2 are both 0/0.
  if (theta_sq < 1e-8) {
    A = 1.0 - one_6th * theta_sq;
    B = 0.5;
  } else {
    if (theta_sq < 1e-6) {
      B = 0.5 - 0.25 * one_6th * theta_sq;
      A = 1.0 - theta_sq * one_6th * (1.0 - one_20th * theta_sq);
    } else {
      const double inv_theta = 1.0 / theta;
      A = sin(theta) * inv_theta;
      B = (1 - cos(theta)) * (inv_theta * inv_theta);
    }
  }
  rodrigues_so3_exp(w, A, B, R);
  return R;
}

inline Eigen::Matrix<double, 4, 4> expSE3(const Eigen::Matrix<double, 6, 1> x) {
  Eigen::Matrix<double, 4, 4> P = Eigen::Matrix<double, 4, 4>::Identity();

  double one_6th = 1.0 / 6.0;
  double one_20th = 1.0 / 20.0;

  Eigen::Matrix<double, 3, 1> w = x.block(3, 0, 3, 1);
  Eigen::Matrix<double, 3, 1> t = x.block(0, 0, 3, 1);

  double theta_sq = w.dot(w);
  double theta = sqrt(theta_sq);
  Eigen::Matrix<double, 3, 1> cross = skew(w) * t;

  Eigen::Matrix<double, 3, 3> R_out;
  Eigen::Matrix<double, 3, 1> t_out;

  double A, B;

  if (theta_sq < 1e-8)
  {
    A = 1.0 - one_6th * theta_sq;
    B = 0.5;
    t_out = t + 0.5 * cross;
  }
  else
  {
    double C;
    if (theta_sq < 1e-6)
    {
      C = one_6th * (1.0 - one_20th * theta_sq);
      A = 1.0 - theta_sq * C;
      B = 0.5 - 0.25 * one_6th * theta_sq;
    }
    else
    {
      double inv_theta = 1.0 / theta;
      A = sin(theta) * inv_theta;
      B = (1 - cos(theta)) * (inv_theta * inv_theta);
      C = (1 - A) * (inv_theta * inv_theta);
    }
    t_out = t +  B * cross + C * skew(w) * cross;
  }
  R_out = expSO3(w);

  P.block(0, 0, 3, 3) = R_out;
  P.block(0, 3, 3, 1) = t_out;

  return P;
}




pcl::PointXYZ substract(const pcl::PointXYZ &p1, const pcl::PointXYZ &p2) {
  pcl::PointXYZ result;
  result.x = p2.x - p1.x;
  result.y = p2.y - p1.y;
  result.z = p2.z - p1.z;
  return result;
}


pcl::PointCloud<pcl::PointXYZ>::Ptr substractPointcloud(
  const pcl::PointCloud<pcl::PointXYZ>::Ptr pc1,
  const pcl::PointCloud<pcl::PointXYZ>::Ptr pc2) {

  if (pc1->size() != pc2->size()) throw
    std::runtime_error("pcltools::substract - Error the point clouds must have the same size!");
  if (pc1->size() == 0 || pc2->size() == 0) throw
    std::runtime_error("pcltools::substract - Error the point clouds must not be empty!");

  typename pcl::PointCloud<pcl::PointXYZ>::Ptr result(new pcl::PointCloud<pcl::PointXYZ>());
  result->reserve(pc1->size());
  for (unsigned int i = 0; i < pc1->size(); ++i) {
    const pcl::PointXYZ &p1 = (*pc1)[i];
    const pcl::PointXYZ &p2 = (*pc2)[i];
    const pcl::PointXYZ &r = substract(p1, p2);
    result->push_back(r);
  }
  return result;
}

void subPointCloud(const  pcl::PointCloud<pcl::PointXYZ>::Ptr &src,
                   const std::vector<int> &indices,
                    pcl::PointCloud<pcl::PointXYZ>::Ptr &dst) {
  dst->clear();
  dst->reserve(indices.size());
  for (unsigned int i = 0; i < indices.size(); i++) {
    dst->push_back((*src)[indices[i]]);
  }
}


void IcpPointToPoint::findNearestNeighbors(
  const pcl::PointCloud<pcl::PointXYZ>::Ptr &src,
  double max_correspondance_distance,
  std::vector<int> &indices_ref,
  std::vector<int> &indices_current,
  std::vector<double> &distances) {
  // We're only interrested in the nearest point
  const int K = 1;
  indices_ref.clear();
  indices_current.clear();
  indices_ref.reserve(src->size());
  indices_current.reserve(src->size());
  distances.clear();
  distances.reserve(src->size());
  std::vector<int> pointIdxNKNSearch(K);
  std::vector<float> pointNKNSquaredDistance(K);

  pcl::PointXYZ pt;
  for (unsigned int i = 0; i < src->size(); i++) {
    // Copy only coordinates from the point (for genericity)
    pt.x = src->at(i).x;
    pt.y = src->at(i).y;
    pt.z = src->at(i).z;

    // Look for the nearest neighbor
    if ( kdtree_.nearestKSearch(pt, K, pointIdxNKNSearch,
                                pointNKNSquaredDistance) > 0 ) {
      double distance = pointNKNSquaredDistance[0];
      if (distance <= max_correspondance_distance * max_correspondance_distance) {
        indices_ref.push_back(i);
        indices_current.push_back(pointIdxNKNSearch[0]);
        distances.push_back(distance);
      } else {
        //LOG(INFO) << "Ignoring, distance too big << " << distance;
      }
    } else {
      LOG(WARNING) << "Could not find a nearest neighbor for point " << i;
    }
  }
}


void IcpPointToPoint::run() {
  // Cleanup
  r_.clear();
  iter_ = 0;
  boost::optional<double> error_variation;

  // Stopping condition. ICP will stop when one of two things
  // happens
  // - The error variation drops below a small threshold min_variation
  // - The number of iteration reaches the maximum max_iter allowed
  bool converged = false;
  while ((converged = step()) && (!error_variation || (error_variation && *error_variation < 0 &&
                                         -*error_variation > param_.min_variation &&
                                         iter_ <= param_.max_iter)))  {
    error_variation = r_.getLastErrorVariation();

    if (error_variation) {
      LOG(INFO) << "Iteration " << iter_ << "/" << param_.max_iter <<
                std::setprecision(8) << ", E=" << r_.getLastError() <<
                ", error_variation=" << *error_variation;
    } else {
      LOG(INFO) << "Iteration " << iter_ << "/" << param_.max_iter <<
                std::setprecision(8) << ", E=" << r_.getLastError() <<
                ", error_variation=none";
    }
  }
  r_.has_converged = converged && (iter_ <= param_.max_iter);
}

bool IcpPointToPoint::step() {
  /**
   * Notations:
   * - P_ref_: reference point cloud \f[ P^* \f]
   * - P_current_: \f[ P \f], current point cloud (CAO model, cloud extracted from one sensor
   * view....)
   * - xk: pose twist to be optimized \f[ \xi \f]
   * - T_(xk): pose in SE3
   * - hat_T: previous pose
   **/

  ++iter_;
  if (P_current_->size() == 0) {
    convergenceFailed();
    return false;
  }

  std::vector<int> indices_ref;
  std::vector<int> indices_current;
  std::vector<double> distances;
  PcPtr P_current_transformed(new Pc());
  PcPtr P_current_phi(new Pc());
  PrPtr P_ref_phi(new Pr());

  pcl::transformPointCloud(*P_current_, *P_current_transformed, T_);
  // XXX only convert if needed!
  pcl::PointCloud<pcl::PointXYZ>::Ptr P_current_transformed_xyz(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::copyPointCloud(*P_current_transformed, *P_current_transformed_xyz);

  if (P_current_transformed_xyz->size() == 0) {
    LOG(ERROR) << "Error: ICP can't run on empty pointclouds!";
    convergenceFailed();
    return false;
  }

  try {
    findNearestNeighbors(P_current_transformed_xyz, param_.max_correspondance_distance,
                         indices_ref, indices_current, distances);
  } catch (...) {
    LOG(WARNING) << "Could not find the nearest neighbors in the KD-Tree, impossible to run ICP without them!";
    return false;
  }

  if (indices_ref.size() == 0) {
    LOG(ERROR) << "Error: No nearest neightbors found";
    convergenceFailed();
    return false;
  }


  // Generate new current point cloud with only the matches in it
  // XXX: Speed improvement possible by using the indices directly instead of
  // generating a new pointcloud. Maybe PCL has stuff to do it.
  subPointCloud(P_current_transformed, indices_ref, P_current_phi);
  subPointCloud(P_ref_init_inv, indices_current, P_ref_phi);

  // Update the reference point cloud to use the previously estimated one
  err_.setInputReference(P_ref_phi);
  err_.setInputCurrent(P_current_phi);
  // Computes the Jacobian
  err_.computeJacobian();

  // Computes the error
  err_.computeError();

  // Initialize mestimator weights from point cloud
  if (param_.mestimator) {
    err_.computeWeights();
  }

  // Transforms the reference point cloud according to new twist
  // Computes the Gauss-Newton update-step
  T_ = err_.update() * T_;

  double E = err_.getErrorNorm();
  r_.registrationError.push_back(E);
  r_.transformation = param_.initial_guess * T_ ;
  r_.relativeTransformation = T_;
  try {
    r_.scale = Sophus::Sim3d(T_).scale();
  } catch (...) {
    LOG(WARNING) << "Invalid icp scale factor, setting to 1!";
    r_.scale = 1;
  }
  if (std::isinf(E)) {
    LOG(WARNING) << "Error is infinite!";
  }
  return true;
}

void ErrorPointToPoint::computeJacobian() {
  const unsigned int n = reference_->size();
  JacobianMatrix J;
  J.setZero(3 * n, 6);
  pcl::PointXYZ p;
  for (unsigned int i = 0; i < n; ++i)
  {
    const pcl::PointXYZ &p_t =  (*reference_)[i];
    p.x = p_t.x;
    p.y = p_t.y;
    p.z = p_t.z;
    J.row(i * 3)     <<  -1,     0,    0,    0,   -p.z,   p.y;
    J.row(i * 3 + 1) <<  0,    -1,    0,  p.z,      0,  -p.x;
    J.row(i * 3 + 2) <<  0,     0,   -1, -p.y,    p.x,     0;
  }
  constraints_->processJacobian(J, J_);
}

void ErrorPointToPoint::computeError() {
  // XXX: Does not make use of eigen's map, possible optimization for floats

  PcPtr pc_e = substractPointcloud(current_, reference_);
  //Eigen::MatrixXf matrixMap = current_->getMatrixXfMap(3, 4, 0) - reference_->getMatrixXfMap(3, 4, 0);

  pcl::PointXYZ p;
  for (unsigned int i = 0; i < pc_e->size(); ++i)
  {
    const pcl::PointXYZ &p_t = (*pc_e)[i];
    p.x = p_t.x;
    p.y = p_t.y;
    p.z = p_t.z;
    errorVector_[i * 3] = p.x;
    errorVector_[i * 3 + 1] =  p.y;
    errorVector_[i * 3 + 2] =  p.z;
  }
  if (!errorVector_.allFinite()) {
    LOG(WARNING) << "Error Vector has NaN values\n!" << errorVector_;
  }
}

Eigen::Matrix<double, 4, 4> ErrorPointToPoint::update() {
  auto Jt = J_.transpose();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> W = weightsVector_.asDiagonal().inverse();
  Eigen::Matrix<double, 6, 1> x = -constraints_->getTwist((Jt*W*J_).ldlt().solve(Jt * W * errorVector_));
  // return update step transformation matrix
  return  expSE3(x);
}

double median(const Eigen::Matrix<double, Eigen::Dynamic, 1> &M) {
  // Work on a copy
  Eigen::Matrix<double, Eigen::Dynamic, 1> copy = M;
  const int len = copy.derived().size();
  if (len % 2 == 0) {
    // Even number of elements,
    // the median is the average of the two central values
    // Sort half the elements
    std::nth_element( copy.data(),
                      copy.data() + len / 2 - 1,
                      copy.data() + len);
    const double n1 = copy(len / 2 - 1);
    std::nth_element( copy.data(),
                      copy.data() + len / 2,
                      copy.data() + len);
    const double n2 = copy(len / 2);
    return (n1 + n2) / double(2);
  } else {
    std::nth_element( copy.data(),
                      copy.data() + len / 2,
                      copy.data() + len);
    // midpoint is the median
    return copy(len / 2 );
  }
}

double median_absolute_deviation(const VectorX& v)
{
  double median_ = median(v);
  // Median centered residual error
  VectorX r = v;
  r = (r-median_*VectorX::Ones(v.size())).cwiseAbs();
//   VectorX r = (v - median_).cwiseAbs();

  // median absolute deviation deviation (MAD)
  double mad = median(r);
  return mad;
}

double hubert_weight(const double z, const double c = 1.345)
{
  double abs_z = std::abs(z);
  if(abs_z < c)
  {
    return 1;
  }
  else
  {
    return c/abs_z;
  }
}

void hubert_weight(const VectorX& r, VectorX& result, double scale, double c = 1.345)
{
  assert(r.size() == result.size());
  for (int i = 0; i < r.size(); ++i) {
    result[i] = hubert_weight(r[i] / scale, c);
  }
}

void ErrorPointToPoint::computeWeights()
{
  double mad = median_absolute_deviation(errorVector_);
  double scale = mad / 0.6745;
  hubert_weight(errorVector_, weightsVector_, scale);
}