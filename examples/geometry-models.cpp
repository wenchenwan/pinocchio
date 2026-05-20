// ============================================================
// 几何模型的加载与位姿更新（C++ 版）
// 演示：
//   1. urdf::buildGeom → 分别加载碰撞/视觉几何模型
//   2. GeometryData → 几何数据结构（存储位姿计算结果）
//   3. forwardKinematics → 正向运动学
//   4. updateGeometryPlacements → 将关节位姿传播到几何体
//
// 与 geometry-models.py 的等价 C++ 实现
// ============================================================

#include "pinocchio/multibody/fcl.hpp"         // GeometryModel, GeometryData
#include "pinocchio/parsers/urdf.hpp"           // urdf::buildModel, urdf::buildGeom

#include "pinocchio/algorithm/joint-configuration.hpp"  // randomConfiguration
#include "pinocchio/algorithm/kinematics.hpp"   // forwardKinematics
#include "pinocchio/algorithm/geometry.hpp"     // updateGeometryPlacements

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  // 支持命令行参数指定模型路径
  const std::string model_path =
    (argc <= 1) ? PINOCCHIO_MODEL_DIR + std::string("/example-robot-data/robots") : argv[1];
  const std::string mesh_dir      = (argc <= 1) ? PINOCCHIO_MODEL_DIR : argv[1];
  const std::string urdf_filename = model_path + "/ur_description/urdf/ur5_robot.urdf";

  // ---- 加载运动学模型 ----
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);

  // ---- 分别加载碰撞/视觉几何模型 ----
  // urdf::buildGeom(model, urdf, geom_type, geom_model, mesh_dir)
  // COLLISION：解析 URDF 中 <collision> 标签（通常是简化形状）
  // VISUAL：解析 URDF 中 <visual> 标签（通常是精细 STL/DAE 网格）
  GeometryModel collision_model;
  pinocchio::urdf::buildGeom(model, urdf_filename, COLLISION, collision_model, mesh_dir);
  GeometryModel visual_model;
  pinocchio::urdf::buildGeom(model, urdf_filename, VISUAL, visual_model, mesh_dir);
  std::cout << "model name: " << model.name << std::endl;

  // ---- 创建数据结构 ----
  Data data(model);
  // GeometryData：存储每个几何体在世界坐标系中的 SE(3) 位姿（oMg[]）
  GeometryData collision_data(collision_model);
  GeometryData visual_data(visual_model);

  // Sample a random configuration
  Eigen::VectorXd q = randomConfiguration(model);
  std::cout << "q: " << q.transpose() << std::endl;

  // ---- 正向运动学 ----
  // 更新 data.oMi[]：每个关节在世界坐标系中的 SE(3) 位姿
  forwardKinematics(model, data, q);

  // ---- 传播到几何体 ----
  // updateGeometryPlacements：
  //   对每个几何体 k，所属关节为 j：
  //   geom_data.oMg[k] = data.oMi[j] * geom_model.geometryObjects[k].placement
  // 必须在 forwardKinematics 之后调用
  updateGeometryPlacements(model, data, collision_model, collision_data);
  updateGeometryPlacements(model, data, visual_model,    visual_data);

  // Print out the placement of each joint of the kinematic tree
  std::cout << "\nJoint placements:" << std::endl;
  for (JointIndex joint_id = 0; joint_id < (JointIndex)model.njoints; ++joint_id)
    std::cout << std::setw(24) << std::left << model.names[joint_id] << ": " << std::fixed
              << std::setprecision(2) << data.oMi[joint_id].translation().transpose() << std::endl;

  // Print out the placement of each collision geometry object
  // collision_data.oMg[k]：第 k 个碰撞几何体在世界坐标系中的位置
  std::cout << "\nCollision object placements:" << std::endl;
  for (GeomIndex geom_id = 0; geom_id < (GeomIndex)collision_model.ngeoms; ++geom_id)
    std::cout << geom_id << ": " << std::fixed << std::setprecision(2)
              << collision_data.oMg[geom_id].translation().transpose() << std::endl;

  // Print out the placement of each visual geometry object
  // visual_data.oMg[k]：第 k 个视觉几何体在世界坐标系中的位置
  std::cout << "\nVisual object placements:" << std::endl;
  for (GeomIndex geom_id = 0; geom_id < (GeomIndex)visual_model.ngeoms; ++geom_id)
    std::cout << geom_id << ": " << std::fixed << std::setprecision(2)
              << visual_data.oMg[geom_id].translation().transpose() << std::endl;
}
