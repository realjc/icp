#ifndef ICP_ICP_P2P_HPP
#define ICP_ICP_P2P_HPP
#include<eigen3/Eigen/Geometry>
#include <vector>
#include "load_data.hpp"
using namespace Eigen;
using namespace std;

void pointToPoint(Matrix3Xd& src, Matrix3Xd& dst, Isometry3d& trans);

#endif