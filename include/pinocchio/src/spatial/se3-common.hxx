//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

// ============================================================
// SE3 家族的公共类型定义宏
//
// 作用：任何 SE3 类（SE3Tpl、以及未来的其他派生类）都需要从
//   traits<Derived> 中导入同一批类型别名。手写十几行 typedef 既冗长又
//   容易漏，故用宏统一注入 —— 这是 Eigen/Pinocchio 处理 CRTP 类型系统
//   的标准做法（Motion/Force/Inertia 各有对应的同类宏）。
//
// 两个入口的区别（TYPENAME 参数）：
//   PINOCCHIO_SE3_TYPEDEF_TPL(D) → 传入 typename，用于【模板类】内部
//       （此时 traits<D>::Scalar 是待决类型，语法上必须写 typename）
//   PINOCCHIO_SE3_TYPEDEF(D)     → 传入空，用于【非模板/已特化类】内部
//       （类型已确定，加 typename 反而是语法错误）
// ============================================================
#define PINOCCHIO_SE3_TYPEDEF_GENERIC(Derived, TYPENAME)                                           \
  typedef TYPENAME traits<Derived>::Scalar Scalar;                                                 \
  typedef TYPENAME traits<Derived>::AngularType AngularType;                                       \
  typedef TYPENAME traits<Derived>::LinearType LinearType;                                         \
  typedef TYPENAME traits<Derived>::AngularRef AngularRef;                                         \
  typedef TYPENAME traits<Derived>::LinearRef LinearRef;                                           \
  typedef TYPENAME traits<Derived>::ConstAngularRef ConstAngularRef;                               \
  typedef TYPENAME traits<Derived>::ConstLinearRef ConstLinearRef;                                 \
  typedef TYPENAME traits<Derived>::ActionMatrixType ActionMatrixType;                             \
  typedef TYPENAME traits<Derived>::HomogeneousMatrixType HomogeneousMatrixType;                   \
  typedef TYPENAME traits<Derived>::PlainType PlainType;                                           \
  static constexpr int Options = traits<Derived>::Options;                                         \
  static constexpr int LINEAR = traits<Derived>::LINEAR;                                           \
  static constexpr int ANGULAR = traits<Derived>::ANGULAR

#define PINOCCHIO_SE3_TYPEDEF_TPL(Derived) PINOCCHIO_SE3_TYPEDEF_GENERIC(Derived, typename)

#define PINOCCHIO_SE3_TYPEDEF(Derived)                                                             \
  PINOCCHIO_SE3_TYPEDEF_GENERIC(Derived, PINOCCHIO_MACRO_EMPTY_ARG)

namespace pinocchio
{
  /* Type returned by the "se3Action" and "se3ActionInverse" functions. */
  // ---- 群作用的返回类型萃取 ----
  // 回答："SE3 作用在类型 D 上，结果是什么类型？"
  //
  // 默认（本主模板）：返回 D 自身 —— 变换一个 Motion 仍得 Motion，
  //   变换一个 Force 仍得 Force。
  // 但有些类型会【特化】此模板来改变返回类型，例如：
  //   MotionZero 被变换后仍是零运动 → 返回 MotionZero（编译期就知道是 0，
  //   后续运算可被完全优化掉）；
  //   各种 Ref/表达式类型被变换后需返回其对应的 Plain 类型。
  //
  // 用法见 se3-tpl.hxx 的 act_impl：
  //   typename SE3GroupAction<D>::ReturnType act_impl(const D & d) const
  template<typename D>
  struct SE3GroupAction
  {
    typedef D ReturnType;
  };
} // namespace pinocchio
