#ifndef ICP_DATA_TYPE_HPP
#define ICP_DATA_TYPE_HPP
#include <eigen3/Eigen/Dense>
#include <vector>
using namespace Eigen;
using namespace std;
class PCD_CLOUD{
public:
    vector<Vector3d> pointsXYZ_;
    vector<Vector3d> normals_; 
    bool hasNormal_ = false;
};
#endif