#include "load_data.hpp"

using namespace std;
using namespace Eigen;

void loadData(std::string path, PCD_CLOUD& cloud){
    std::ifstream file(path.c_str(), std::ifstream::in);
    if(file.fail()==true){
       // LOG(ERROR) << "点云读取失败 ..."<<endl;
    }
    int i = 0;
    while (file)
    {
        Vector3d pt, norm;
        file >> pt.x()>>pt.y()>>pt.z()>>norm.x()>>norm.y()>>norm.z();
        cloud.pointsXYZ_.push_back(pt);
        cloud.normals_.push_back(norm);
    }
    return;
}

Matrix3Xd vec2mat(vector<Vector3d>& vec){
    Map<Matrix<double,3,Dynamic,ColMajor> > mat(&vec[0].x(), 3, vec.size());
    return mat;
}

Isometry3d addNoise(Isometry3d& pose, double mean, double std){
    std::normal_distribution<double> nd(mean, std);
    Vector3d noise(nd(generator), nd(generator), nd(generator));
    auto foo = Sophus::SO3d::exp(noise);
    return pose*foo.unit_quaternion();
}
