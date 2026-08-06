//
// Copyright (c) 2015-2019 CNRS INRIA
// Copyright (c) 2015-2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  template<typename Scalar, int Options>
  struct SE3GroupAction<MotionZeroTpl<Scalar, Options>>
  {
    typedef MotionZeroTpl<Scalar, Options> ReturnType;
  };

  template<typename Scalar, int Options, typename MotionDerived>
  struct MotionAlgebraAction<MotionZeroTpl<Scalar, Options>, MotionDerived>
  {
    typedef MotionZeroTpl<Scalar, Options> ReturnType;
  };

  template<typename _Scalar, int _Options>
  struct traits<MotionZeroTpl<_Scalar, _Options>>
  {
    static constexpr int Options = _Options;
    static constexpr int LINEAR = 0;
    static constexpr int ANGULAR = 3;
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
    typedef Eigen::Matrix<Scalar, 4, 4, Options> Matrix4;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Vector6) ToVectorConstReturnType;
    typedef typename PINOCCHIO_EIGEN_REF_TYPE(Vector6) ToVectorReturnType;
    typedef Matrix6 ActionMatrixType;
    typedef Matrix4 HomogeneousMatrixType;
    typedef Vector3 AngularType;
    typedef const Vector3 ConstAngularType;
    typedef Vector3 LinearType;
    typedef const Vector3 ConstLinearType;
    typedef MotionTpl<Scalar, _Options> MotionPlain;
    typedef MotionPlain PlainReturnType;

  }; // traits MotionZeroTpl

  // ============================================================
  // MotionZeroTpl：编译期已知为【零】的空间速度
  //
  // 核心思想：把"这个速度恒为 0"这一事实编码进【类型】里，
  //   使编译器能把相关运算彻底消除，而不是运行时算一堆 0。
  //
  // 典型用途：
  //   · 固定关节 / 世界根节点的速度恒为 0
  //   · RNEA 中根节点的初始速度 ν₀ = 0
  //   · 静力学计算（v = a = 0）时把整条速度传播链优化掉
  //
  // 注意本类【不存储任何数据】（sizeof 约为 1 字节），所有成员函数
  // 要么是 static，要么直接返回零对象 —— 零内存、零运算。
  // ============================================================
  template<typename Scalar, int Options>
  struct MotionZeroTpl : public MotionBase<MotionZeroTpl<Scalar, Options>>
  {
    typedef typename traits<MotionZeroTpl>::MotionPlain MotionPlain;
    typedef typename traits<MotionZeroTpl>::PlainReturnType PlainReturnType;

    // 需要真实数据时才"物化"成一个全零的 MotionTpl
    static PlainReturnType plain()
    {
      return MotionPlain::Zero();
    }

    // ---- 与稠密 Motion 比较：退化为判断对方是否全零 ----
    template<typename D2>
    static bool isEqual_impl(const MotionDense<D2> & other)
    {
      return other.linear().isZero(0) && other.angular().isZero(0);
    }

    // 两个零运动必然相等 —— 编译期即可确定，不做任何运行时比较
    static bool isEqual_impl(const MotionZeroTpl &)
    {
      return true;
    }

    // ---- 加到别的对象上：加 0 等于什么都不做，整个调用会被优化掉 ----
    template<typename D2>
    static void addTo(const MotionBase<D2> &)
    {
    }

    // ---- 赋值给别的对象：直接置零，无需逐分量拷贝 ----
    template<typename D2>
    static void setTo(MotionBase<D2> & other)
    {
      other.setZero();
    }

    // ---- 叉乘：0 × ν = 0，直接返回零对象，不做任何计算 ----
    template<typename M1>
    MotionZeroTpl motionAction(const MotionBase<M1> &) const
    {
      return MotionZeroTpl();
    }

    // ---- 坐标变换：Ad_M · 0 = 0，与 M 无关，故参数名都省略了 ----
    // 这正是本类的价值所在：一次 6×6 变换（36 次乘加）被完全消除
    template<typename S2, int O2, typename D2>
    void se3Action_impl(const SE3Tpl<S2, O2> &, MotionDense<D2> & v) const
    {
      v.setZero();
    }

    template<typename S2, int O2>
    MotionZeroTpl se3Action_impl(const SE3Tpl<S2, O2> &) const
    {
      return MotionZeroTpl(); // 返回类型仍是 MotionZero，"零"的信息不丢失
    }

    template<typename S2, int O2, typename D2>
    void se3ActionInverse_impl(const SE3Tpl<S2, O2> &, MotionDense<D2> & v) const
    {
      v.setZero();
    }

    template<typename S2, int O2>
    MotionZeroTpl se3ActionInverse_impl(const SE3Tpl<S2, O2> &) const
    {
      return MotionZeroTpl();
    }

  }; // struct MotionZeroTpl

  // ---- ν + 0 = ν：直接返回原对象的【引用】----
  // 既不做加法，也不产生任何临时对象；两个重载覆盖左右两侧
  template<typename M1, typename Scalar, int Options>
  inline const M1 & operator+(const MotionBase<M1> & v, const MotionZeroTpl<Scalar, Options> &)
  {
    return v.derived();
  }

  template<typename Scalar, int Options, typename M1>
  inline const M1 & operator+(const MotionZeroTpl<Scalar, Options> &, const MotionBase<M1> & v)
  {
    return v.derived();
  }

} // namespace pinocchio
