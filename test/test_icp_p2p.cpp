#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include<iostream>
#include "load_data.hpp"
#include "data_type.hpp"
#include "icp_p2p.hpp"
using namespace Eigen;

int main()
{
    // 在输入点云上添加
    Quaterniond q = Quaterniond::Identity();
    q = q*AngleAxisd(0.1,Vector3d::UnitX());
    q = q*AngleAxisd(0.1,Vector3d::UnitY());
    q = q*AngleAxisd(0.1,Vector3d::UnitZ());
    Vector4d rot(q.x(),q.y(),q.z(),q.w());
    Vector3d trans(0.1,0.1,0.1);
    Isometry3d pose = Translation3d(trans)*Quaterniond(rot);

    PCD_CLOUD pcd;
    std::string path = "/home/bluerov/Downloads/icp/data/cloudXYZ_0.xyz";
    loadData(path, pcd);
    Matrix3Xd src = vec2mat(pcd.pointsXYZ_);
    Matrix3Xd dst = pose*src;
    //在旋转上增加噪声
    Isometry3d init_pose = addNoise(pose,0.1,0.1);

    std::cout<<"pose: "<<endl;
    std::cout<< pose.matrix()<<endl;
    std::cout<<"init_pose: "<<endl;
    std::cout<< init_pose.matrix()<<endl;
    for(int i = 0;i<3;i++){
        pointToPoint(src, dst, init_pose);
        std::cout<<"第 "<<i<<" 次结果"<<endl;
        std::cout<<init_pose.matrix()<<endl;
    }

}