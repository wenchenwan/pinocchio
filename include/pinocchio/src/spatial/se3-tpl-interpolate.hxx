//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  // This function is defined here because it need se3-tpl and explog
  // headers
  // ============================================================
  // SE(3) 测地线插值（SLERP 在刚体位姿上的推广）
  //
  //   Interpolate(A, B, α) = A · exp₆( α · log₆(A⁻¹B) )
  //
  // 三步几何解释：
  //   ① A⁻¹B          把"从 A 到 B 的相对运动"提取出来（在 A 的局部系中）
  //   ② log₆(·)       映射到李代数 se(3)，得到一个 6D 旋量 dv（螺旋运动）
  //   ③ A·exp₆(α·dv)  沿该螺旋走 α 比例的路程，再变换回世界系
  //
  // 性质：
  //   α=0 → A，α=1 → B，中间沿【测地线】（流形上的最短路径）平滑过渡；
  //   结果始终是合法的 SE(3) 元素（R 保持正交），这是逐元素线性插值
  //   做不到的 —— 线性插值会离开 SO(3)，得到非旋转矩阵。
  //
  // 之所以单独放在这个文件：Interpolate 同时依赖 se3-tpl（SE3Tpl 的定义）
  // 和 explog（log6/exp6），必须等两者都完整定义后才能实现，
  // 故在 se3-tpl.hxx 中只作声明。
  // ============================================================
  template<typename Scalar, int Options>
  template<typename OtherScalar>
  SE3Tpl<Scalar, Options> SE3Tpl<Scalar, Options>::Interpolate(
    const SE3Tpl & A, const SE3Tpl & B, const OtherScalar & alpha)
  {
    typedef SE3Tpl<Scalar, Options> ReturnType;
    typedef MotionTpl<Scalar, Options> Motion;

    Motion dv = log6(A.actInv(B));      // ①②：相对变换 → 李代数旋量
    ReturnType res = A * exp6(alpha * dv); // ③：沿旋量走 α 比例后回到群
    return res;
  }
} // namespace pinocchio
