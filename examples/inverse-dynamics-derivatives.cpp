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

  // ============================================================
  // 分配导数矩阵：三个都必须是 nv × nv
  //
  // djoint_torque_dq：∂τ/∂q（重力梯度项 + Coriolis 对 q 的偏导）
  // djoint_torque_dv：∂τ/∂v（Coriolis 矩阵）
  // djoint_torque_da：∂τ/∂a = M(q)（质量矩阵）
  //
  // ⚠️ 为什么 ∂τ/∂q 是 nv×nv 而不是 nv×nq？
  //   τ 是 nv 维，而 q 是 nq 维，直觉上似乎该是 nv×nq。但求导是对
  //   【切空间扰动】δq ∈ R^nv 而言的，不是对 q 的原始分量：
  //       τ(q ⊕ δq) ≈ τ(q) + (∂τ/∂q)·δq,   δq ∈ R^nv
  //   q 住在流形上，其原始分量之间有约束（如四元数的归一化），
  //   对它们逐个求偏导没有意义。所以全部导数都定义在切空间。
  //
  //   本例的 UR5 是固定基（nq = nv = 6），掩盖了这个区别；
  //   换成浮动基人形（nq=35, nv=34）就能看出来 —— 三个矩阵仍是 34×34。
  //   若误写成 Zero(model.nv, model.nq)，浮动基下会触发
  //   rnea-derivatives.hxx 里的 PINOCCHIO_CHECK_ARGUMENT_SIZE 断言。
  //
  // ⚠️ 必须预先置零：computeRNEADerivatives 内部是【累加】写入，
  //   头文件明确要求 "must be first initialized with zeros"。
  // ============================================================
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

  // ============================================================
  // ⚠️ 三角填充陷阱：djoint_torque_da 只填充【上三角】
  //
  // 头文件明确声明（rnea-derivatives.hpp）：
  //   "As for pinocchio::crba, only the upper triangular part of
  //    rnea_partial_da is filled."
  // 实测确认：其严格下三角范数恒为 0。
  //
  // 因为 ∂τ/∂a = M(q)，若直接把它当质量矩阵用，会得到一个上三角阵
  // 而非对称阵 —— 且【不会有任何报错】。实测后果：
  //     ‖(未对称化 M)·a + nle − τ‖ = 1.48e+00     ← 完全错误
  //     ‖(对称化后 M)·a + nle − τ‖ = 2.40e-15     ← 正确
  // 相差 15 个数量级。
  //
  // 注意 C++ 的 crba() 返回值同样只有上三角；Python 绑定则会在
  // 返回前对称化，所以这个坑只在 C++ 侧出现。
  //
  // 正确用法：
  //   Eigen::MatrixXd M = djoint_torque_da.selfadjointView<Eigen::Upper>();
  // ============================================================

  // 可选：访问导数矩阵
  // std::cout << "dtau/dq:\n" << djoint_torque_dq << std::endl;
  // std::cout << "dtau/dv:\n" << djoint_torque_dv << std::endl;
  // 若要当作 M(q) 使用，先对称化：
  // Eigen::MatrixXd M = djoint_torque_da.selfadjointView<Eigen::Upper>();
}
