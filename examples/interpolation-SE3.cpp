// ============================================================
// SE(3) 插值（Geodesic Interpolation on Lie Groups）
// 演示：在两个 SE(3) 位姿之间进行测地线插值
//
// 数学原理：
//   SE(3) 是黎曼流形（Lie 群），直接线性插值会破坏矩阵约束
//   正确方法：测地线插值（Geodesic / SLERP 扩展到 SE(3)）
//     interpolate(a, b, t) = a ⊕ (t · log(a⁻¹·b))
//     = a · exp(t · log(a⁻¹·b))
//   当 t=0 时返回 a，t=1 时返回 b，t=0.5 时返回中点
//
// 应用场景：
//   - 轨迹规划：在关键帧之间插值末端执行器位姿
//   - 动画：关节姿态的平滑过渡（SLERP 的 SE(3) 推广）
//   - 位姿估计：融合多个传感器的位姿观测
//
// 配置向量格式（SE(3) 的 7D 参数化）：
//   [x, y, z, qx, qy, qz, qw]（位置 3D + 四元数 4D）
//   注意：Pinocchio 四元数顺序为 (x, y, z, w)（与 Eigen 不同，Eigen 为 (w, x, y, z)）
// ============================================================

#include <iostream>
#include "pinocchio/multibody/liegroup/liegroup.hpp"  // SpecialEuclideanOperationTpl

using namespace pinocchio;

int main()
{
  typedef double Scalar;

  // SpecialEuclideanOperationTpl<3, Scalar>：SE(3) Lie 群操作类模板
  // 提供：integrate(), difference(), interpolate(), normalize() 等流形操作
  typedef SpecialEuclideanOperationTpl<3, Scalar> SE3Operation;
  SE3Operation aSE3;

  // ConfigVector_t：SE(3) 的配置向量类型（7D：位置 3 + 四元数 4）
  // TangentVector_t：SE(3) 的切向量类型（6D：速度 3 + 角速度 3）
  SE3Operation::ConfigVector_t pose_s, pose_g;
  SE3Operation::TangentVector_t delta_u;

  // ---- 起始位姿（Starting configuration）----
  // 位置：(1, 1, 1)
  // 四元数：(-0.13795, 0.13795, 0.69352, 0.69352)（格式：qx, qy, qz, qw）
  pose_s(0) = 1.0;    // x
  pose_s(1) = 1.0;    // y
  pose_s(2) = 1;      // z
  pose_s(3) = -0.13795;  // qx
  pose_s(4) = 0.13795;   // qy
  pose_s(5) = 0.69352;   // qz
  pose_s(6) = 0.69352;   // qw
  // normalize：确保四元数单位化（|q| = 1）
  // 手工输入的四元数可能因数值误差不满足归一化约束
  aSE3.normalize(pose_s);

  // ---- 目标位姿（Goal configuration）----
  // 位置：(4, 3, 3)
  // 四元数：(-0.46194, 0.331414, 0.800103, 0.191342)
  pose_g(0) = 4;
  pose_g(1) = 3;
  pose_g(2) = 3;
  pose_g(3) = -0.46194;
  pose_g(4) = 0.331414;
  pose_g(5) = 0.800103;
  pose_g(6) = 0.191342;
  aSE3.normalize(pose_g);

  // ---- SE(3) 测地线插值 ----
  // interpolate(a, b, t, result)：
  //   t=0.5 → 起始和目标之间的中点位姿
  // 内部：result = a ⊕ (t · (b ⊖ a)) = a · exp(t · log(a⁻¹·b))
  // 保证：result 始终是合法的 SE(3) 元素（旋转矩阵正交，四元数归一化）
  SE3Operation::ConfigVector_t pole_u;
  aSE3.interpolate(pose_s, pose_g, 0.5, pole_u);
  std::cout << "Interpolated configuration: " << pole_u.transpose() << std::endl;

  return 0;
}
