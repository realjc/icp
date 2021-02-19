#include "icp_p2p.hpp"

void pointToPoint(Matrix3Xd& src, Matrix3Xd& dst, Isometry3d& trans){
    int N = src.cols(),M = dst.cols();
    assert(M==N);
    Vector3d ps_mean = src.rowwise().mean();
    Vector3d qs_mean = dst.rowwise().mean();
    Matrix3Xd ps_centered = src.colwise() - ps_mean;
    Matrix3Xd qs_centered = dst.colwise() - qs_mean;
    Matrix3d k = qs_centered*ps_centered.transpose();
    JacobiSVD<Matrix3d> svd(k,ComputeFullU|ComputeFullV);
    Matrix3d r = svd.matrixU()*svd.matrixV().transpose();
    if(r.determinant()<0) r.col(2)*=-1;
    trans.linear() = r;
    trans.translation() = qs_mean-r*ps_mean;
    return;
}