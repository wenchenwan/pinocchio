#include <iostream>
#include "pinocchio/multibody/liegroup/liegroup.hpp"

using namespace pinocchio;

int main()
{
  typedef double Scalar;
  enum { Options = 0 };

  // SE(2)：2D 特殊欧氏群，描述平面刚体的位姿（平面移动机器人、足底地面投影等）
  // 配置向量格式（4 维）：[x, y, cos(θ), sin(θ)]
  //   前 2 维：平面位置；后 2 维：旋转角的余弦/正弦表示（避免角度的不连续性）
  // 切空间格式（3 维）：[Δx, Δy, Δθ]
  typedef SpecialEuclideanOperationTpl<2, Scalar, Options> SE2Operation;
  SE2Operation aSE2;
  SE2Operation::ConfigVector_t pose_s, pose_g;  // 4 维
  SE2Operation::TangentVector_t delta_u;         // 3 维
  delta_u.setZero();

  // 起始配置：位置 [1, 1]，朝向角 45°（π/4）
  pose_s(0) = 1.0;
  pose_s(1) = 1.0;
  pose_s(2) = cos(M_PI / 4.0);  // cos(45°) ≈ 0.707
  pose_s(3) = sin(M_PI / 4.0);  // sin(45°) ≈ 0.707

  // 目标配置：位置 [3, -1]，朝向角 -90°（-π/2）
  pose_g(0) = 3.0;
  pose_g(1) = -1.0;
  pose_g(2) = cos(-M_PI / 2.0);  // cos(-90°) = 0
  pose_g(3) = sin(-M_PI / 2.0);  // sin(-90°) = -1

  // difference(a, b, u)：在切空间中计算从 a 到 b 的差分
  // u = log(a^{-1} * b)，维度为 3（Δx, Δy, Δθ），是流形上最短路径方向
  // 人形应用：步态规划中计算当前步态相位和目标步态之间的差异
  aSE2.difference(pose_s, pose_g, delta_u);
  std::cout << "difference: " << delta_u.transpose() << std::endl;

  // integrate(a, u, b)：从 a 出发沿切向量 u 积分，得到 b = a * exp(u)
  // 验证：以 pose_s 为起点，沿 delta_u 走一步，应恢复 pose_g
  SE2Operation::ConfigVector_t pose_check;
  aSE2.integrate(pose_s, delta_u, pose_check);
  std::cout << "goal configuration (from composition): " << pose_check.transpose() << std::endl;
  std::cout << "goal configuration: " << pose_g.transpose() << std::endl;
}
