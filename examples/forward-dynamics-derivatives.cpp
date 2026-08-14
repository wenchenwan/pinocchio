// ============================================================
// ABA 解析导数（Forward Dynamics Derivatives）
// 计算：∂q̈/∂q, ∂q̈/∂v, ∂q̈/∂τ ∈ R^{nv×nv}
//
// ABA（Articulated Body Algorithm）正向动力学：
//   q̈ = ABA(q, v, τ) = M(q)⁻¹ · [τ - C(q,v)·v - g(q)]
//
// 导数的物理含义（DDP/iLQR 轨迹优化中的 A, B 矩阵）：
//   ∂q̈/∂q  → 系统矩阵 A 的下半块（位置对加速度的影响，包含重力梯度）
//   ∂q̈/∂v  → 系统矩阵 A 的速度块（Coriolis 项对加速度的影响）
//   ∂q̈/∂τ  = M(q)⁻¹（控制矩阵 B，力矩到加速度的映射）
//
// 优势：相比有限差分，解析导数无截断误差，且计算量只比单次 ABA 略高
// ============================================================

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"  // computeABADerivatives

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

  // Sample a random joint configuration as well as random joint velocity and torque
  Eigen::VectorXd q   = randomConfiguration(model);
  Eigen::VectorXd v   = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);

  // ---- 分配导数矩阵（nv×nv）----
  // djoint_acc_dq：∂q̈/∂q（系统矩阵 A 下半块 q 部分）
  // djoint_acc_dv：∂q̈/∂v（系统矩阵 A 下半块 v 部分）
  // djoint_acc_dtau：∂q̈/∂τ = M(q)⁻¹（控制矩阵 B）
  Eigen::MatrixXd djoint_acc_dq   = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd djoint_acc_dv   = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd djoint_acc_dtau = Eigen::MatrixXd::Zero(model.nv, model.nv);

  // ---- computeABADerivatives：一次调用计算 q̈ + 三个导数矩阵 ----
  // 内部：
  //   1. 调用 ABA 计算 data.ddq（正向动力学结果）
  //   2. 通过 ABA 的变分传播计算导数（无需有限差分）
  // 调用后：
  //   data.ddq = q̈（正向动力学结果）
  //   djoint_acc_dq, djoint_acc_dv, djoint_acc_dtau 被填充
  computeABADerivatives(model, data, q, v, tau, djoint_acc_dq, djoint_acc_dv, djoint_acc_dtau);

  // data.ddq：ABA 计算得到的关节加速度（nv 维向量）
  std::cout << "Joint acceleration: " << data.ddq.transpose() << std::endl;

  // 可选：访问导数矩阵用于 DDP/iLQR
  // djoint_acc_dtau = M(q)⁻¹
  //
  // 与 RNEA 导数不同，本例这个【显式输出参数】版本的 djoint_acc_dtau
  // 是【完整对称】的，可直接使用（实测严格下三角非零且数值正确）。
  //
  // ⚠️ 但同一次调用中的 data.Minv 只填了上三角（它被当作中间缓冲）。
  //   两个"看似等价"的东西三角性不同：
  //     computeABADerivatives(m,d,q,v,tau)            → data.Minv 完整对称
  //     computeABADerivatives(m,d,q,v,tau,dq,dv,dtau) → dtau 完整、data.Minv 仅上三角
  //   所以要用 M⁻¹ 时，请用【你自己传进去的 dtau】，而不是顺手去读 data.Minv。
  //
  //   对比：RNEA 导数的 ∂τ/∂a 在两种形式下【都】只有上三角，
  //   必须 .selfadjointView<Eigen::Upper>() 才能当 M(q) 用
  //   （见 examples/inverse-dynamics-derivatives.cpp 的说明）。
  //
  // std::cout << "M_inv (=ddq/dtau):\n" << djoint_acc_dtau << std::endl;
}
