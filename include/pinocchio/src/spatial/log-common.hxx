//
// Copyright (c) 2015-2021 CNRS INRIA

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  template<typename Scalar>
  // ============================================================
  // 前向声明：log 系列的标量分派器
  //
  // 为什么要拆成 xxx_impl 结构体而非直接写函数：
  //   log3/log6 在 θ→0、θ→π 处有奇异，处理方式依赖【标量类型】——
  //   普通 double 可以用 if 分支，而自动微分标量必须用可微的选择算子。
  //   做成模板结构体后即可按 Scalar 偏特化，见 log.hxx 及
  //   autodiff/ 目录下针对 CppAD/CasADi 的特化。
  // ============================================================
  struct log3_impl;
  template<typename Scalar>
  struct Jlog3_impl;

  template<typename Scalar>
  struct log6_impl;
  template<typename Scalar>
  struct Jlog6_impl;

  template<typename Matrix3>
  inline typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3)
    renormalize_rotation_matrix(const Eigen::MatrixBase<Matrix3> & R);

} // namespace pinocchio
