// ============================================================
// 从 URDF 加载模型，计算正向运动学并打印关节位姿
// 演示：
//   1. urdf::buildModel → 解析 URDF 文件，构建运动学树
//   2. randomConfiguration → 在关节限位内随机采样配置
//   3. forwardKinematics → 正向运动学（更新 data.oMi）
//   4. data.oMi[joint_id] → 关节在世界坐标系中的 SE(3) 位姿
//
// PINOCCHIO_MODEL_DIR 由 CMake 定义（指向 pinocchio/models 目录）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"               // urdf::buildModel

#include "pinocchio/algorithm/joint-configuration.hpp"  // randomConfiguration
#include "pinocchio/algorithm/kinematics.hpp"        // forwardKinematics

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  // 支持命令行参数指定 URDF 路径，否则使用 CMake 内置路径
  const std::string urdf_filename =
    (argc <= 1) ? PINOCCHIO_MODEL_DIR
                    + std::string("/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf")
                : argv[1];

  // ---- urdf::buildModel：解析 URDF，构建运动学模型 ----
  // 固定基版本（不添加浮动基关节）
  // 若需要浮动基：urdf::buildModel(urdf_filename, JointModelFreeFlyer(), model)
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);
  std::cout << "model name: " << model.name << std::endl;

  // Create data required by the algorithms
  Data data(model);

  // ---- randomConfiguration：在关节限位内均匀随机采样 ----
  // 对于旋转关节：在 [lowerPositionLimit, upperPositionLimit] 内均匀采样角度
  // 对于球形关节/浮动基（四元数）：在 SO(3)/SE(3) 上均匀采样（保证归一化）
  Eigen::VectorXd q = randomConfiguration(model);
  std::cout << "q: " << q.transpose() << std::endl;

  // ---- forwardKinematics：正向运动学 ----
  // 从根关节开始，递归计算每个关节在世界坐标系中的 SE(3) 变换
  // 结果存储在 data.oMi[joint_id]（world_T_joint）
  // 时间复杂度：O(n)，n = 关节数
  forwardKinematics(model, data, q);

  // ---- 打印每个关节的世界坐标系位置 ----
  for (JointIndex joint_id = 0; joint_id < (JointIndex)model.njoints; ++joint_id)
    std::cout << std::setw(24) << std::left << model.names[joint_id] << ": " << std::fixed
              << std::setprecision(2)
              // oMi[joint_id].translation()：关节原点在世界坐标系中的位置（3D 向量）
              // oMi 命名：o = origin（世界系）, M = SE(3) 变换, i = joint i
              << data.oMi[joint_id].translation().transpose() << std::endl;
}
