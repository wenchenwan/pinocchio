//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: private, include "pinocchio/spatial/fwd.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial/fwd.hpp"
#endif // PINOCCHIO_LSP

// ============================================================
// spatial 模块的前向声明汇总
//
// 作用：先声明 SE3Tpl / MotionTpl / ForceTpl / InertiaTpl 等模板的
// 存在与默认模板参数，让各头文件之间的【循环依赖】得以解开
//   （例如 SE3 的 act() 要提到 Motion，而 Motion 的变换又要提到 SE3）。
//
// 这里同时给出面向用户的别名 SE3 / Motion / Force / Inertia，
// 它们都是对应 XxxTpl<context::Scalar, context::Options> 的 typedef
// （生成机制见 PINOCCHIO_GUIDE.md §11.2）。
// ============================================================
namespace pinocchio
{
  /// \internal
  namespace internal
  {
    ///  \brief Default return type for the operation: Type*Scalar
    template<typename Type, typename Scalar>
    struct RHSScalarMultiplication
    {
      typedef Type ReturnType;
    };

    ///  \brief Default return type for the operation: Scalar*Type
    template<typename Type, typename Scalar>
    struct LHSScalarMultiplication
    {
      typedef Type ReturnType;
    };

    // for certain Scalar type, it might be needed to proceed to call some normalization procedure
    // in when performing a cast. This struct is an helper to support such modality.
    template<typename Class, typename NewScalar, typename Scalar>
    struct cast_call_normalize_method;
  } // namespace internal
  /// \endinternal

  template<typename Derived>
  class MotionBase;
  template<typename Derived>
  class MotionDense;
  template<typename Vector6ArgType>
  class MotionRef;
  template<typename Scalar, int Options = context::Options>
  class MotionTpl;
  template<typename Scalar, int Options = context::Options>
  struct MotionZeroTpl;

  template<typename Derived>
  class ForceBase;
  template<typename Derived>
  class ForceDense;
  template<typename Vector6ArgType>
  class ForceRef;
  template<typename Scalar, int Options = context::Options>
  class ForceTpl;

  template<class Derived>
  struct SE3Base;
  template<typename _Scalar, int _Options = context::Options>
  struct SE3Tpl;

  template<typename Scalar, int Options = context::Options>
  class Symmetric3Tpl;

  template<class Derived>
  struct InertiaBase;
  template<typename _Scalar, int _Options = context::Options>
  struct InertiaTpl;
  template<typename Scalar, int Options = context::Options>
  struct PseudoInertiaTpl;
  template<typename Scalar, int Options = context::Options>
  struct LogCholeskyParametersTpl;

  using Force = ForceTpl<context::Scalar, context::Options>;
  using Motion = MotionTpl<context::Scalar, context::Options>;
  using MotionZero = MotionZeroTpl<context::Scalar, context::Options>;
  using SE3 = SE3Tpl<context::Scalar, context::Options>;
  using Symmetric3 = Symmetric3Tpl<context::Scalar, context::Options>;
  using Inertia = InertiaTpl<context::Scalar, context::Options>;
  using PseudoInertia = PseudoInertiaTpl<context::Scalar, context::Options>;
  using LogCholeskyParameters = LogCholeskyParametersTpl<context::Scalar, context::Options>;

} // namespace pinocchio
