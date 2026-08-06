//
// Copyright (c) 2017-2019 CNRS INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<int axis>
  struct SpatialAxis;

  template<int axis, typename MotionDerived>
  struct MotionAlgebraAction<SpatialAxis<axis>, MotionDerived>
  {
    typedef typename MotionDerived::MotionPlain ReturnType;
  };

  // ============================================================
  // SpatialAxis<axis>：编译期已知的【6D 单位旋量轴】
  //
  //   axis 0,1,2 → 平移轴 (e₀,0)、(e₁,0)、(e₂,0)
  //   axis 3,4,5 → 旋转轴 (0,e₀)、(0,e₁)、(0,e₂)
  // 内部用 CartesianAxis<axis % 3> 复用 3D 轴的零成本运算。
  //
  // 用途：单自由度关节的【运动子空间 S】。
  //   例如 JointModelRZ 的 S = SpatialAxis<5>，
  //   于是 ν = S·q̇ 退化成"把 q̇ 写进第 5 个分量"，
  //   τ = Sᵀf 退化成"取出 f 的第 5 个分量" —— 全部无乘法。
  // 这正是 Pinocchio 中 RNEA/ABA 递推能极快的微观原因之一。
  // ============================================================
  template<int _axis>
  struct SpatialAxis //: MotionBase< SpatialAxis<_axis> >
  {
    static constexpr int axis = _axis;
    static constexpr int dim = 6;
    typedef CartesianAxis<_axis % 3> CartesianAxis3;
    static constexpr int LINEAR = 0;
    static constexpr int ANGULAR = 3;

    template<typename Derived1, typename Derived2>
    inline static void cross(const MotionDense<Derived1> & min, const MotionDense<Derived2> & mout);

    template<typename Derived>
    static typename traits<Derived>::MotionPlain cross(const MotionDense<Derived> & min)
    {
      typename MotionDense<Derived>::MotionPlain res;
      cross(min, res);
      return res;
    }

    template<typename Derived1, typename Derived2>
    inline static void cross(const ForceDense<Derived1> & fin, const ForceDense<Derived2> & fout);

    template<typename Derived>
    static typename traits<Derived>::ForcePlain cross(const ForceDense<Derived> & fin)
    {
      typename ForceDense<Derived>::ForcePlain fout;
      cross(fin, fout);
      return fout;
    }

    template<typename Scalar>
    MotionTpl<Scalar> operator*(const Scalar & s) const
    {
      typedef MotionTpl<Scalar> ReturnType;
      ReturnType res;
      for (Eigen::Index i = 0; i < dim; ++i)
        res.toVector()[i] = i == axis ? s : Scalar(0);

      return res;
    }

    template<typename Scalar>
    friend inline MotionTpl<Scalar> operator*(const Scalar & s, const SpatialAxis &)
    {
      return SpatialAxis() * s;
    }

    template<typename Derived>
    friend Derived & operator<<(MotionDense<Derived> & min, const SpatialAxis &)
    {
      typedef typename traits<Derived>::Scalar Scalar;
      min.setZero();
      min.toVector()[axis] = Scalar(1);
      return min.derived();
    }

    template<typename MotionDerived>
    typename MotionDerived::MotionPlain motionAction(const MotionDense<MotionDerived> & m) const
    {
      typename MotionDerived::MotionPlain res;
      if ((LINEAR == 0 && axis < 3) || (LINEAR == 3 && axis >= 3))
      {
        res.angular().setZero();
        CartesianAxis3::cross(-m.angular(), res.linear());
      }
      else
      {
        CartesianAxis3::cross(-m.linear(), res.linear());
        CartesianAxis3::cross(-m.angular(), res.angular());
      }

      return res;
    }

  }; // struct SpatialAxis

  template<int axis>
  template<typename Derived1, typename Derived2>
  inline void
  SpatialAxis<axis>::cross(const MotionDense<Derived1> & min, const MotionDense<Derived2> & mout)
  {
    Derived2 & mout_ = PINOCCHIO_EIGEN_CONST_CAST(Derived2, mout);

    if ((LINEAR == 0 && axis < 3) || (LINEAR == 3 && axis >= 3))
    {
      mout_.angular().setZero();
      CartesianAxis3::cross(min.angular(), mout_.linear());
    }
    else
    {
      CartesianAxis3::cross(min.linear(), mout_.linear());
      CartesianAxis3::cross(min.angular(), mout_.angular());
    }
  }

  template<int axis>
  template<typename Derived1, typename Derived2>
  inline void
  SpatialAxis<axis>::cross(const ForceDense<Derived1> & fin, const ForceDense<Derived2> & fout)
  {
    Derived2 & fout_ = PINOCCHIO_EIGEN_CONST_CAST(Derived2, fout);

    if ((LINEAR == 0 && axis < 3) || (LINEAR == 3 && axis >= 3))
    {
      fout_.linear().setZero();
      CartesianAxis3::cross(fin.linear(), fout_.angular());
    }
    else
    {
      CartesianAxis3::cross(fin.linear(), fout_.linear());
      CartesianAxis3::cross(fin.angular(), fout_.angular());
    }
  }

  typedef SpatialAxis<0> AxisVX;
  typedef SpatialAxis<1> AxisVY;
  typedef SpatialAxis<2> AxisVZ;

  typedef SpatialAxis<3> AxisWX;
  typedef SpatialAxis<4> AxisWY;
  typedef SpatialAxis<5> AxisWZ;
} // namespace pinocchio
