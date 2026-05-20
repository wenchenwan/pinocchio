// ============================================================
// 运动学解析导数（Kinematics Derivatives）
// 计算：∂v_joint/∂q, ∂a_joint/∂q, ∂a_joint/∂v, ∂a_joint/∂a ∈ R^{6×nv}
//
// 物理含义：
//   v_partial_dq = ∂v_joint/∂q（速度对配置的偏导 = 几何 Jacobian）
//   a_partial_dq = ∂a_joint/∂q（加速度对配置的偏导，包含离心/Coriolis 项）
//   a_partial_dv = ∂a_joint/∂v（加速度对速度的偏导 = 速度相关的 Jacobian 导数）
//   a_partial_da = ∂a_joint/∂a（加速度对加速度的偏导 = 几何 Jacobian，与 v_partial_dv 相等）
//
// 注意：v_partial_dv = a_partial_da（两者相等，均等于几何 Jacobian）
//   因此 getJointAccelerationDerivatives 不单独返回 v_partial_dv
//
// 参考系选项：
//   LOCAL：在关节局部坐标系中表达（原点=关节原点，轴=关节轴）
//   WORLD：在世界坐标系对齐的坐标系中表达（原点=关节原点，轴=世界轴）
//   LOCAL_WORLD_ALIGNED：与 WORLD 相同（等价）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"  // computeForwardKinematicsDerivatives, getJointAccelerationDerivatives

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  const std::string urdf_filename =
    (argc <= 1) ? PINOCCHIO_MODEL_DIR
                    + std::string("/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf")
                : argv[1];

  // Load the URDF model
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);

  // Build a data related to model
  Data data(model);

  // Sample a random joint configuration as well as random joint velocity and acceleration
  Eigen::VectorXd q = randomConfiguration(model);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

  // ---- computeForwardKinematicsDerivatives：预计算所有关节的运动学导数 ----
  // 等价于先调用 forwardKinematics，再存储用于后续 get 调用的中间量
  // 必须在 getJointVelocityDerivatives 或 getJointAccelerationDerivatives 之前调用
  computeForwardKinematicsDerivatives(model, data, q, v, a);

  // ---- 获取末端关节（最后一个关节）的运动学导数 ----
  // 从索引 1 开始（0 = 宇宙关节），最后一个为末端执行器
  JointIndex joint_id = (JointIndex)(model.njoints - 1);

  // 分配结果矩阵（6×nv，6D 空间速度/加速度导数）
  Data::Matrix6x v_partial_dq(6, model.nv), a_partial_dq(6, model.nv),
    a_partial_dv(6, model.nv), a_partial_da(6, model.nv);
  v_partial_dq.setZero();
  a_partial_dq.setZero();
  a_partial_dv.setZero();
  a_partial_da.setZero();

  // ---- LOCAL 参考系：导数在关节局部坐标系中表达 ----
  // LOCAL：原点=关节原点，轴=关节轴（随关节旋转而变化）
  // 返回：
  //   v_partial_dq = ∂v_joint/∂q（= 局部系下的几何 Jacobian J_loc）
  //   a_partial_dq = ∂a_joint/∂q
  //   a_partial_dv = ∂a_joint/∂v
  //   a_partial_da = ∂a_joint/∂a = J_loc（与 v_partial_dq 相等）
  getJointAccelerationDerivatives(
    model, data, joint_id, LOCAL, v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da);

  // Remark: we are not directly computing the quantity v_partial_dv as it is also equal to
  // a_partial_da.

  // ---- WORLD 参考系：导数在世界坐标系对齐的框架中表达 ----
  // WORLD：原点=关节原点，轴与世界坐标系对齐（不随关节旋转）
  // 适合：与世界坐标系中的任务（如末端执行器位置控制）直接对接
  // 用法相同，只改变 reference_frame 参数
  getJointAccelerationDerivatives(
    model, data, joint_id, WORLD, v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da);

  // 可选：打印 Jacobian（= a_partial_da 在 LOCAL 参考系）
  // std::cout << "Geometric Jacobian (LOCAL):\n" << a_partial_da << std::endl;
}
