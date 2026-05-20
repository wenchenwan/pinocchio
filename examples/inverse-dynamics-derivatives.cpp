// ============================================================
// RNEA 解析导数（Inverse Dynamics Derivatives）
// 计算：∂τ/∂q, ∂τ/∂v, ∂τ/∂a ∈ R^{nv×nv}
//
// 应用场景（DDP/iLQR 轨迹优化）：
//   - ∂τ/∂q：重力梯度 + Coriolis 对位置的偏导（影响 DDP 的 A 矩阵）
//   - ∂τ/∂v：Coriolis 矩阵的偏导（影响 DDP 的 A 矩阵速度块）
//   - ∂τ/∂a = M(q)：质量矩阵（由 RNEA 线性关系决定）
//
// 优势：相比有限差分，解析导数精确且快速（同等精度下快 2~5×）
// 实现：内部使用两次完整的 RNEA 传播（正向+逆向各一次）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"  // computeRNEADerivatives

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

  // ---- 分配导数矩阵（nv×nv）----
  // djoint_torque_dq：∂τ/∂q（重力梯度项 + Coriolis 对 q 的偏导）
  // djoint_torque_dv：∂τ/∂v（Coriolis 矩阵）
  // djoint_torque_da：∂τ/∂a = M(q)（质量矩阵，与 crba 结果一致）
  Eigen::MatrixXd djoint_torque_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd djoint_torque_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd djoint_torque_da = Eigen::MatrixXd::Zero(model.nv, model.nv);

  // ---- computeRNEADerivatives：一次调用同时计算三个导数矩阵 ----
  // 内部顺序：
  //   1. forwardKinematics + RNEA（计算 τ，存储中间量）
  //   2. 反向传播变分（计算 ∂τ/∂(q,v,a)）
  // 注意：调用后 data.tau 同时被更新（包含 τ 本身）
  computeRNEADerivatives(
    model, data, q, v, a, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);

  // data.tau 在 computeRNEADerivatives 内部被更新（等同于调用 rnea）
  std::cout << "Joint torque: " << data.tau.transpose() << std::endl;

  // 可选：访问导数矩阵
  // std::cout << "dtau/dq:\n" << djoint_torque_dq << std::endl;
  // std::cout << "dtau/dv:\n" << djoint_torque_dv << std::endl;
  // std::cout << "dtau/da (= M(q)):\n" << djoint_torque_da << std::endl;
}
