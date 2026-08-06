//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

// ============================================================
// Motion（空间速度 / Twist）家族的公共类型定义宏
//
// 与 SE3 家族的 PINOCCHIO_SE3_TYPEDEF_* 同理：把 traits<Derived> 中的
// 一批类型别名统一注入到各派生类中，避免每个类手写十几行 typedef。
//
// TYPENAME 参数的两个入口：
//   MOTION_TYPEDEF_TPL(D) → 传 typename，用于模板类内部（类型待决）
//   MOTION_TYPEDEF(D)     → 传空，用于非模板/已特化类内部
// ============================================================
#define MOTION_TYPEDEF_GENERIC(Derived, TYPENAME)                                                  \
  typedef TYPENAME traits<Derived>::Scalar Scalar;                                                 \
  typedef TYPENAME traits<Derived>::Vector3 Vector3;                                               \
  typedef TYPENAME traits<Derived>::Vector6 Vector6;                                               \
  typedef TYPENAME traits<Derived>::Matrix4 Matrix4;                                               \
  typedef TYPENAME traits<Derived>::Matrix6 Matrix6;                                               \
  typedef TYPENAME traits<Derived>::ToVectorReturnType ToVectorReturnType;                         \
  typedef TYPENAME traits<Derived>::ToVectorConstReturnType ToVectorConstReturnType;               \
  typedef TYPENAME traits<Derived>::AngularType AngularType;                                       \
  typedef TYPENAME traits<Derived>::LinearType LinearType;                                         \
  typedef TYPENAME traits<Derived>::ConstAngularType ConstAngularType;                             \
  typedef TYPENAME traits<Derived>::ConstLinearType ConstLinearType;                               \
  typedef TYPENAME traits<Derived>::ActionMatrixType ActionMatrixType;                             \
  typedef TYPENAME traits<Derived>::HomogeneousMatrixType HomogeneousMatrixType;                   \
  typedef TYPENAME traits<Derived>::MotionPlain MotionPlain;                                       \
  typedef TYPENAME traits<Derived>::PlainReturnType PlainReturnType;                               \
  static constexpr int LINEAR = traits<Derived>::LINEAR;                                           \
  static constexpr int ANGULAR = traits<Derived>::ANGULAR

#define MOTION_TYPEDEF_TPL(Derived) MOTION_TYPEDEF_GENERIC(Derived, typename)

#define MOTION_TYPEDEF(Derived) MOTION_TYPEDEF_GENERIC(Derived, PINOCCHIO_MACRO_EMPTY_ARG)

namespace pinocchio
{
  ///
  /// \brief Return type of the ation of a Motion onto an object of type D
  ///
  // ---- 李代数作用（叉乘）的返回类型萃取 ----
  // 回答："空间速度 ν 叉乘作用在类型 D 上，结果是什么类型？"
  //
  // 空间代数中有两种叉乘（见 PINOCCHIO_GUIDE.md §2.2.6）：
  //   ν × μ   （Motion × Motion → Motion）  运动叉乘
  //   ν ×* φ  （Motion × Force  → Force ）  力叉乘（对偶版本）
  // 二者返回类型不同，故需要这个萃取在编译期选出正确类型。
  //
  // 默认返回 D 自身；特殊类型（如 MotionZero、各种 Ref/表达式）会特化它。
  // 用于 MotionBase::cross()：
  //   typename MotionAlgebraAction<OtherSpatialType, Derived>::ReturnType
  template<typename D, typename MotionDerived>
  struct MotionAlgebraAction
  {
    typedef D ReturnType;
  };

} // namespace pinocchio
