#include <pcl/point_cloud.h>
#include <glog/logging.h>


template<typename PointT>
bool operator==(const PointT &p1, const PointT &p2) {
  const float thr = 0.00001;
  return (fabs(p1.x - p2.x) < thr) &&
         (fabs(p1.y - p2.y) < thr) &&
         (fabs(p1.z - p2.z) < thr);
}


namespace pcltools
{

template<typename PointT>
PointT substract(const PointT &p1, const PointT &p2) {
  PointT result;
  result.x = p2.x - p1.x;
  result.y = p2.y - p1.y;
  result.z = p2.z - p1.z;
  return result;
}

template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr substractPointcloud(const typename
    pcl::PointCloud<PointT>::Ptr pc1,
    const typename pcl::PointCloud<PointT>::Ptr pc2) {

  if (pc1->size() != pc2->size()) throw new
    std::runtime_error("pcltools::substract - Error the point clouds must have the same size!");

  typename pcl::PointCloud<PointT>::Ptr result(new pcl::PointCloud<PointT>());
  result->reserve(pc1->size());
  for (int i = 0; i < pc1->size(); ++i) {
    const PointT &p1 = (*pc1)[i];
    const PointT &p2 = (*pc2)[i];
    const PointT &r = substract(p1, p2);
    result->push_back(r);
  }
  return result;
}

} /* pcltools */
