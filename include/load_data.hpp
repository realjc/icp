#ifndef ICP_LOAD_DATA_HPP
#define ICP_LOAD_DATA_HPP
#include <string>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <numeric>
#include <random>
#include <sophus/so3.hpp>
//#include <glog/logging.h>
#include "data_type.hpp"

using namespace std;
using namespace Eigen;

static std::mt19937 generator;

void loadData(std::string path, PCD_CLOUD& cloud);
Matrix3Xd vec2mat(vector<Vector3d>& vec);
Isometry3d  addNoise(Isometry3d& pose, double mean, double std);

#endif