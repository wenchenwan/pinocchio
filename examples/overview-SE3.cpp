#include <iostream>
// SE(3) 及其他李群的统一接口头文件
// 提供 difference / integrate / normalize / interpolate 等流形操作
#include "pinocchio/multibody/liegroup/liegroup.hpp"

using namespace pinocchio;

int main()
{
  typedef double Scalar;

  // SE(3)：3D 特殊欧氏群，描述刚体在三维空间中的完整位姿（位置 + 姿态）
  // 配置向量格式（7 维）：[x, y, z, qx, qy, qz, qw]
  //   前 3 维：平移；后 4 维：单位四元数表示旋转
  // 切空间（速度/误差）格式（6 维）：[vx, vy, vz, wx, wy, wz]
  typedef SpecialEuclideanOperationTpl<3, Scalar> SE3Operation;

  SE3Operation aSE3;
  SE3Operation::ConfigVector_t pose_s, pose_g;  // 7 维配置向量
  SE3Operation::TangentVector_t delta_u;         // 6 维切向量（李代数元素）

  // 起始位姿：位置 [1, 1, 1]，四元数 [qx, qy, qz, qw]
  pose_s(0) = 1.0;   // x
  pose_s(1) = 1.0;   // y
  pose_s(2) = 1;     // z
  pose_s(3) = -0.13795;  // qx
  pose_s(4) = 0.13795;   // qy
  pose_s(5) = 0.69352;   // qz
  pose_s(6) = 0.69352;   // qw

  // 目标位姿：位置 [4, 3, 3]
  pose_g(0) = 4;
  pose_g(1) = 3;
  pose_g(2) = 3;
  pose_g(3) = -0.46194;
  pose_g(4) = 0.331414;
  pose_g(5) = 0.800103;
  pose_g(6) = 0.191342;

  // 归一化四元数，确保满足 ||q||=1 的约束（输入可能存在数值误差）
  aSE3.normalize(pose_s);
  std::cout << "pose_s: " << pose_s.transpose() << std::endl;
  aSE3.normalize(pose_g);
  std::cout << "pose_g: " << pose_g.transpose() << std::endl;

  // difference(a, b, u)：计算从 a 到 b 的"最短路径"切向量
  // 数学含义：u = log(a^{-1} * b)，结果 u 在 a 的切空间（李代数）中
  // 物理含义：以 a 为起点，沿 u 走"一步"可到达 b
  // 人形机器人应用：计算当前姿态与目标姿态之间的 6D 误差
  aSE3.difference(pose_s, pose_g, delta_u);
  std::cout << "delta_u: " << delta_u.transpose() << std::endl;

  // integrate(a, u, b)：从 a 出发沿切向量 u 做李群积分，结果写入 b
  // 数学含义：b = a * exp(u)
  // 验证：integrate(pose_s, delta_u) 应该恢复 pose_g
  SE3Operation::ConfigVector_t pose_check;
  aSE3.integrate(pose_s, delta_u, pose_check);
  std::cout << "pose_check: " << pose_check.transpose() << std::endl;

  return 0;
}
