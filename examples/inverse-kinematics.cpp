// ============================================================
// 逆运动学（IK）：阻尼最小二乘迭代法（C++ 版）
// 目标：找关节角 q，使第 6 个关节的 SE(3) 位姿 = oMdes
//
// 数学：
//   iMd = oMi⁻¹ · oMdes（当前→目标的相对变换，在关节局部系中）
//   err = log6(iMd) ∈ R⁶（SE(3) 上的对数映射，得到 6D 误差）
//   J_corrected = -Jlog6(iMd⁻¹) · J（链式法则修正，适配 log6 的 Jacobian）
//   v = -Jᵀ · (J·Jᵀ + λ²I)⁻¹ · err（阻尼最小二乘）
//   q ← integrate(model, q, v·DT)（流形积分，维持四元数归一化）
//
// Jlog6 的必要性：
//   err = log6(iMd) 是 iMd 的非线性函数，直接用 J 会有偏差
//   Jlog6 = d(log6)/d(iMd) 是 SE(3) 对数映射的 Jacobian
//   正确的线性化：J_corrected = Jlog6 · d(iMd)/dq = Jlog6 · J（带符号）
// ============================================================

#include <iostream>

#include "pinocchio/multibody/sample-models.hpp"    // buildModels::manipulator
#include "pinocchio/spatial/explog.hpp"             // log6(), Jlog6()
#include "pinocchio/algorithm/kinematics.hpp"       // forwardKinematics
#include "pinocchio/algorithm/jacobian.hpp"         // computeJointJacobian
#include "pinocchio/algorithm/joint-configuration.hpp"  // neutral, integrate

int main(int /* argc */, char ** /* argv */)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);
  pinocchio::Data data(model);

  // 目标关节 ID（第 6 个关节 = 末端执行器）
  const int JOINT_ID = 6;
  // 目标位姿：旋转 = 单位矩阵（无旋转），位置 = (1, 0, 1)
  const pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), Eigen::Vector3d(1., 0., 1.));

  Eigen::VectorXd q = pinocchio::neutral(model);  // 初始配置：零位

  const double eps    = 1e-4;   // 收敛阈值：‖err‖ < eps 时停止
  const int    IT_MAX = 1000;   // 最大迭代次数
  const double DT     = 1e-1;   // 积分步长（相当于速度的时间增量）
  const double damp   = 1e-6;   // 阻尼系数 λ（Levenberg-Marquardt 正则化项）

  // J ∈ R^{6×nv}：关节空间 Jacobian（在关节局部坐标系中）
  pinocchio::Data::Matrix6x J(6, model.nv);
  J.setZero();

  bool success = false;
  typedef Eigen::Matrix<double, 6, 1> Vector6d;
  Vector6d err;
  Eigen::VectorXd v(model.nv);

  for (int i = 0;; i++)
  {
    // 步骤 1：正向运动学，更新所有关节位姿（data.oMi）
    pinocchio::forwardKinematics(model, data, q);

    // 步骤 2：计算当前→目标的相对变换（在当前关节局部系中）
    // actInv：iMd = oMi⁻¹ · oMdes（将 oMdes 变换到关节 i 的局部坐标系）
    const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

    // 步骤 3：计算 6D 误差（SE(3) 对数映射）
    // log6(iMd) ∈ se(3)：线速度误差（3D）+ 角速度误差（3D）
    // 当 iMd = Identity 时，err = 0（收敛）
    err = pinocchio::log6(iMd).toVector();  // in joint frame

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

    // 步骤 4：计算关节 Jacobian（在关节局部系 LOCAL 中）
    pinocchio::computeJointJacobian(model, data, q, JOINT_ID, J);  // J in joint frame

    // 步骤 5：Jlog6 修正——SE(3) 对数映射的 Jacobian
    // Jlog6(M) ∈ R^{6×6}：d(log6(M))/dM 在 M 处的导数
    // 注意：传入 iMd.inverse() 而非 iMd（与 Python 版的 -Jlog 对应）
    pinocchio::Data::Matrix6 Jlog;

    // iMd描述的是局部关节坐标系的位姿误差，对其取逆进行对数变换就会得到左雅可比（左扰动雅可比）
    // Jlog6实际计算的右扰动雅可比（右扰动雅可比），所以需要对其取逆来得到左扰动雅可比
    pinocchio::Jlog6(iMd.inverse(), Jlog);
    // 修正后的 Jacobian（负号来自 iMd 对 q 的偏导方向）
    J = -Jlog * J;

    // 步骤 6：阻尼最小二乘求解速度
    // JJᵀ ∈ R^{6×6}（操作空间惯量逆矩阵的正规方程）
    pinocchio::Data::Matrix6 JJt;
    JJt.noalias() = J * J.transpose();
    // 对角加阻尼：(JJᵀ + λ²I)，防止 Jacobian 奇异时求逆爆炸
    JJt.diagonal().array() += damp;
    // LDLT 分解求解：v = -Jᵀ · (JJᵀ + λ²I)⁻¹ · err（阻尼最小二乘）
    v.noalias() = -J.transpose() * JJt.ldlt().solve(err);

    // 步骤 7：在关节配置流形上积分
    // integrate：q ← q ⊕ (v·DT)，自动处理旋转关节的指数映射
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
