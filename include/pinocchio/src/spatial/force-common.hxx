//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

// ============================================================
// Force（空间力 / Wrench）家族的公共类型定义宏
//
// 与 MOTION_TYPEDEF_* / PINOCCHIO_SE3_TYPEDEF_* 同一套路：
// 把 traits<Derived> 中的类型别名统一注入各派生类。
//   FORCE_TYPEDEF_TPL(D) → 模板类内部用（传 typename）
//   FORCE_TYPEDEF(D)     → 非模板类内部用（传空）
// ============================================================
#define FORCE_TYPEDEF_GENERIC(Derived, TYPENAME)                                                   \
  typedef TYPENAME traits<Derived>::Scalar Scalar;                                                 \
  typedef TYPENAME traits<Derived>::Vector3 Vector3;                                               \
  typedef TYPENAME traits<Derived>::Vector6 Vector6;                                               \
  typedef TYPENAME traits<Derived>::Matrix6 Matrix6;                                               \
  typedef TYPENAME traits<Derived>::ToVectorReturnType ToVectorReturnType;                         \
  typedef TYPENAME traits<Derived>::ToVectorConstReturnType ToVectorConstReturnType;               \
  typedef TYPENAME traits<Derived>::AngularType AngularType;                                       \
  typedef TYPENAME traits<Derived>::LinearType LinearType;                                         \
  typedef TYPENAME traits<Derived>::ConstAngularType ConstAngularType;                             \
  typedef TYPENAME traits<Derived>::ConstLinearType ConstLinearType;                               \
  typedef TYPENAME traits<Derived>::ForcePlain ForcePlain;                                         \
  static constexpr int LINEAR = traits<Derived>::LINEAR;                                           \
  static constexpr int ANGULAR = traits<Derived>::ANGULAR

#define FORCE_TYPEDEF_TPL(Derived) FORCE_TYPEDEF_GENERIC(Derived, typename)

#define FORCE_TYPEDEF(Derived) FORCE_TYPEDEF_GENERIC(Derived, PINOCCHIO_MACRO_EMPTY_ARG)
