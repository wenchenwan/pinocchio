// ============================================================
// 3D 逆运动学（位置 IK）：只控制末端位置，不控制姿态（C++ 版）
// 与 inverse-kinematics.cpp 的区别：
//   - 误差空间：R³（3D 位置差）而非 R⁶（6D 位姿差）
//   - Jacobian：只取 6×nv 矩阵的前 3 行（平移 Jacobian）
//   - 无需 Jlog6 修正（位置误差是线性的，不需要 SE(3) 链式法则）
//   - 阻尼系数更小（1e-12 vs 1e-6），因为问题维度更低，数值更稳定
//
// 数学：
//   err = iMd.translation() ∈ R³（当前→目标的位置差，在关节局部系中）
//   J = -joint_jacobian.topRows<3>() ∈ R^{3×nv}（平移 Jacobian）
//   v = -Jᵀ · (JJᵀ + λ²I)⁻¹ · err（阻尼最小二乘）
//   q ← integrate(model, q, v·DT)
// ============================================================

#include <iostream>

#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

int main(int /* argc */, char ** /* argv */)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);
  pinocchio::Data data(model);

  const int JOINT_ID = 6;
  // 目标位姿（旋转部分在 3D IK 中被忽略，只用位置）
  const pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), Eigen::Vector3d(1., 0., 1.));

  Eigen::VectorXd q = pinocchio::neutral(model);
  const double eps    = 1e-4;
  const int    IT_MAX = 1000;
  const double DT     = 1e-1;
  // 3D IK 阻尼可以更小（问题本身更欠定，数值更稳定）
  const double damp   = 1e-12;

  // 完整 6×nv Jacobian（计算后只取前 3 行）
  pinocchio::Data::Matrix6x joint_jacobian(6, model.nv);
  joint_jacobian.setZero();

  bool success = false;
  Eigen::Vector3d err;   // 3D 位置误差（而非 6D）
  Eigen::VectorXd v(model.nv);

  for (int i = 0;; i++)
  {
    pinocchio::forwardKinematics(model, data, q);

    // iMd = 当前关节位姿的逆 × 目标位姿（变换到关节局部系）
    const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

    // 3D 位置误差：iMd.translation()（在关节局部系中的位置差）
    // 无需 log6，因为位置差是欧氏空间的线性量
    err = iMd.translation();  // in joint frame

    if (err.norm() < eps)
    {
      success = true;
      break;
    }
    if (i >= IT_MAX)
    {
      success = false;
      break;
    }

    // 计算完整的 6×nv Jacobian
    pinocchio::computeJointJacobian(
      model, data, q, JOINT_ID, joint_jacobian);  // joint_jacobian expressed in the joint frame

    // 只取前 3 行（平移部分的 Jacobian J_lin ∈ R^{3×nv}）
    // 负号：与 iMd = oMi⁻¹·oMdes 的推导一致（方向对齐）
    const auto J = -joint_jacobian.topRows<3>();  // Jacobian associated with the error

    // JJᵀ ∈ R^{3×3}（比 6D 版本小，计算更快）
    const Eigen::Matrix3d JJt = J * J.transpose() + damp * Eigen::Matrix3d::Identity();
    v.noalias() = -J.transpose() * JJt.ldlt().solve(err);
    q = pinocchio::integrate(model, q, v * DT);

    if (!(i % 10))
      std::cout << i << ": error = " << err.transpose() << std::endl;
  }

  if (success)
  {
    std::cout << "Convergence achieved!" << std::endl;
  }
  else
  {
    std::cout
      << "\nWarning: the iterative algorithm has not reached convergence to the desired precision"
      << std::endl;
  }

  std::cout << "\nresult: " << q.transpose() << std::endl;
  std::cout << "\nfinal error: " << err.transpose() << std::endl;
}
