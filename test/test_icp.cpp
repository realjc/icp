#include <pcl/common/transforms.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/point_types.h>
#include "point_to_point.hpp"


typedef pcl::PointCloud<pcl::PointXYZ> PointCloudXYZ;
typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudXYZRGB;

void pp_callback (const pcl::visualization::PointPickingEvent &event);

int main()
{
    std::cout << "ICP Example" << std::endl;
    /*
    Loads test point cloud
    */
    std::cout << "Loading Model pointcloud" << std::endl;
    PointCloudXYZ::Ptr modelCloud(new PointCloudXYZ());
    PointCloudXYZ::Ptr dataCloud(new PointCloudXYZ());
    std::string model = "../models/valve.pcd";
    if (pcl::io::loadPCDFile<pcl::PointXYZ> (model.c_str(), *modelCloud) == -1) {
    LOG(FATAL) << "Could't read file " << model;
    return (-1);
    }
    std::cout << "Model Point cloud has " << modelCloud->points.size()
        << " points" << std::endl;
    std::string data = "../models/valve_simulation_clean.pcd";
    if (pcl::io::loadPCDFile<pcl::PointXYZ> (data.c_str(), *dataCloud) == -1) {
    LOG(FATAL) << "Could't read file " << model;
    return (-1);
    }
    std::cout << "Model Point cloud has " << dataCloud->points.size()
        << " points" << std::endl;
    Eigen::Matrix4d transformation;
    PointCloudXYZ::Ptr resultCloud(new PointCloudXYZ());

      /*
     Define parameters for the ICP
     */
    IcpParameters icp_param;
    icp_param.mestimator = true;
    icp_param.max_iter = 20;
    icp_param.min_variation = 10e-5;
    icp_param.initial_guess = Eigen::Matrix4d::Identity();
    // Far
    //icp_param.initial_guess(0, 3) = 1.6;
    //icp_param.initial_guess(1, 3) = 0.6;
    //icp_param.initial_guess(2, 3) = 1.6;
    // Less far
    icp_param.initial_guess(0, 3) = 2;
    icp_param.initial_guess(1, 3) = 0.7;
    icp_param.initial_guess(2, 3) = 1;
    icp_param.initial_guess.block<3,3>(0,0) =
        Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitX()).matrix()
        * Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitY()).matrix()
        * Eigen::AngleAxisd(M_PI+0.3, Eigen::Vector3d::UnitZ()).matrix();
    // Almost registered
    //icp_param.initial_guess(0, 3) = 2.176;
    //icp_param.initial_guess(1, 3) = 0.868;
    //icp_param.initial_guess(2, 3) = 1;
    //

    // big valve close
    icp_param.initial_guess(0, 3) = 2.1;
    icp_param.initial_guess(1, 3) = 0.;
    icp_param.initial_guess(2, 3) = 1;
    icp_param.initial_guess.block<3,3>(0,0) =
        1.5 * Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitX()).matrix()
        * Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitY()).matrix()
        * Eigen::AngleAxisd(M_PI+0.3, Eigen::Vector3d::UnitZ()).matrix();
    // std::cout << "ICP Parameters:\n" << icp_param << std::endl;

    std::cout << "Computing PointToPoint ICP" << std::endl;
    /**
     * Point to point
     **/
    IcpPointToPoint icp_algorithm;
    icp_algorithm.setParameters(icp_param);
    icp_algorithm.setInputCurrent(modelCloud);
    icp_algorithm.setInputReference(dataCloud);
    icp_algorithm.run();

    IcpResults icp_results = icp_algorithm.getResults();
    std::cout << "ICP Results:\n" << icp_results.transformation << std::endl;
    pcl::transformPointCloud(*modelCloud, *resultCloud, icp_results.transformation);


    /**
     * Visualize
     **/
    std::cout << "Point cloud colors:\n" <<
                "\twhite  = original\n"
                "\tred    = transformed\n"
                "\tgreen  = registered";
    pcl::visualization::PCLVisualizer viewer("Matrix transformation example");
    viewer.registerPointPickingCallback (pp_callback);

    PointCloudXYZ::Ptr initialRegistrationCloud(new PointCloudXYZ());
    pcl::transformPointCloud(*modelCloud, * initialRegistrationCloud, icp_param.initial_guess);

    // Define R,G,B colors for the point cloud
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>
    source_cloud_color_handler(modelCloud, 255, 0, 0);
    // We add the point cloud to the viewer and pass the color handler
    viewer.addPointCloud(initialRegistrationCloud,
                        source_cloud_color_handler, "original_cloud");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>
    transformed_cloud_color_handler(dataCloud, 100, 100, 100);  // Red
    viewer.addPointCloud(dataCloud, transformed_cloud_color_handler,
                        "transformed_cloud");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>
    registered_cloud_color_handler(dataCloud, 20, 230, 20);  // Green
    viewer.addPointCloud(resultCloud,
                        registered_cloud_color_handler,
                        "registered cloud");

    viewer.addCoordinateSystem(1.0, "cloud", 0);
    // Setting background to a dark grey
    //viewer.setBackgroundColor(0.05, 0.05, 0.05, 0);
    viewer.setBackgroundColor(1., 1., 1., 0);
    viewer.setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "original_cloud");
    viewer.setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "registered cloud");
    viewer.setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "transformed_cloud");

    // std::stringstream r;
    // r << "White: Origial, Red: Transformed, Green: Registered\n";
    // r << icp_results;
    // viewer.addText(r.str(), 0, 0);
    viewer.setShowFPS(false);

    // Display the visualiser until 'q' key is pressed
    while (!viewer.wasStopped()) {
        viewer.spinOnce();
    }
    //}

    return 0;
}


void pp_callback (const pcl::visualization::PointPickingEvent &event)
{
  if (event.getPointIndex () == -1)
    return;
  float x, y, z;
  event.getPoint(x, y, z);
  std::cout << "Point Selected: \n"
            << "\tIndex: " << event.getPointIndex ()
            << "\tCoord: (" <<  x << ", " << y << ", " << z << ")";
}
