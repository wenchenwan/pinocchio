//
// Copyright (c) 2017-2020 CNRS
// Copyright (c) 2018-2024 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  // ---- 返回类型萃取：MotionDense 的作用结果 = 其具体派生类的结果 ----
  // 这两个特化把"包装层"透明化：对 MotionDense<D> 做变换/叉乘，
  // 返回类型直接沿用 D 的结果类型，不会退化成 MotionDense
  template<typename Derived>
  struct SE3GroupAction<MotionDense<Derived>>
  {
    typedef typename SE3GroupAction<Derived>::ReturnType ReturnType;
  };

  template<typename Derived, typename MotionDerived>
  struct MotionAlgebraAction<MotionDense<Derived>, MotionDerived>
  {
    typedef typename MotionAlgebraAction<Derived, MotionDerived>::ReturnType ReturnType;
  };

  // ============================================================
  // MotionDense：稠密存储 Motion 的【通用实现层】
  //
  // 继承链中的位置：
  //   MotionBase（纯接口）→ MotionDense（本类，实现所有 _impl）
  //                       → MotionTpl（自有内存）/ MotionRef（引用外部内存）
  //
  // 本类只依赖派生类提供 linear()/angular() 两个访问器，
  // 其余全部算法（叉乘、变换、加减）都在这里统一实现 ——
  // 因此 MotionTpl 和 MotionRef 无需各自重复这些代码。
  // ============================================================
  template<typename Derived>
  class MotionDense : public MotionBase<Derived>
  {
  public:
    typedef MotionBase<Derived> Base;
    MOTION_TYPEDEF_TPL(Derived);
    typedef typename traits<Derived>::MotionRefType MotionRefType;

    // 把基类的重载集引入本作用域，避免被本类同名函数"名字隐藏"
    using Base::angular;
    using Base::derived;
    using Base::isApprox;
    using Base::isZero;
    using Base::linear;

    Derived & setZero()
    {
      linear().setZero();
      angular().setZero();
      return derived();
    }
    Derived & setRandom()
    {
      linear().setRandom();
      angular().setRandom();
      return derived();
    }

    // ---- 6×6 伴随作用矩阵 ad_ν：把"ν × ·"写成矩阵乘法 ----
    //
    //   ad_ν = [ ω̂   v̂ ]      按 [v; ω] 顺序
    //          [ 0    ω̂ ]
    //
    // 满足 ad_ν · μ = ν × μ（实测残差 0）。
    // 展开即：(ν×μ).v = v×μ_ω + ω×μ_v，(ν×μ).ω = ω×μ_ω
    // 这里和论文定义的负号相反，原因是论文中速度的定义是【w ; v】,这里是【v ; w】。
    ActionMatrixType toActionMatrix_impl() const
    {
      ActionMatrixType X;
      X.template block<3, 3>(ANGULAR, ANGULAR) = X.template block<3, 3>(LINEAR, LINEAR) =
        skew(angular());                                          // 对角两块都是 ω̂
      X.template block<3, 3>(LINEAR, ANGULAR) = skew(linear());   // 右上 = v̂
      X.template block<3, 3>(ANGULAR, LINEAR).setZero();          // 左下 = 0

      return X;
    }

    // ---- 对偶伴随 ad_ν*：用于力的叉乘 ν ×* φ ----
    // 与上面的区别：v̂ 从【右上】移到【左下】（转置的结果，体现力-速度对偶）
    // 满足 ad*_ν · φ = ν ×* φ（实测残差 5.6e-17）。
    // RNEA 中的陀螺/科氏项 ν ×* (Iν) 正是用它（见 §4.3.3）
    ActionMatrixType toDualActionMatrix_impl() const
    {
      ActionMatrixType X;
      X.template block<3, 3>(ANGULAR, ANGULAR) = X.template block<3, 3>(LINEAR, LINEAR) =
        skew(angular());
      X.template block<3, 3>(ANGULAR, LINEAR) = skew(linear());   // 左下 = v̂
      X.template block<3, 3>(LINEAR, ANGULAR).setZero();          // 右上 = 0

      return X;
    }

    // ---- 4×4 齐次表示 ξ̂ = [ω̂ v; 0 0]，满足 Ṁ = ξ̂·M ----
    HomogeneousMatrixType toHomogeneousMatrix_impl() const
    {
      HomogeneousMatrixType M;
      M.template block<3, 3>(0, 0) = skew(angular());
      M.template block<3, 1>(0, 3) = linear();
      M.template block<1, 4>(3, 0).setZero();   // 最后一行恒为 0（不是 [0 0 0 1]）
      return M;
    }

    // ---- 精确相等：稠密-稠密之间逐分量比较 ----
    template<typename D2>
    bool isEqual_impl(const MotionDense<D2> & other) const
    {
      return linear() == other.linear() && angular() == other.angular();
    }

    // 与非稠密类型（如 MotionZero）比较时反转调用方向，
    // 让对方用自己的特化实现来比（零运动可在编译期判定，无需展开成 6 个数）
    template<typename D2>
    bool isEqual_impl(const MotionBase<D2> & other) const
    {
      return other.derived() == derived();
    }

    // Arithmetic operators
    // ---- 赋值：稠密 ← 稠密 ----
    template<typename D2>
    Derived & operator=(const MotionDense<D2> & other)
    {
      return derived().set(other.derived());
    }

    Derived & operator=(const MotionDense & other)
    {
      return derived().set(other.derived());
    }

    template<typename D2>
    Derived & set(const MotionDense<D2> & other)
    {
      linear() = other.linear();
      angular() = other.angular();
      return derived();
    }

    // ---- 赋值：稠密 ← 任意 Motion（含 MotionZero 等特殊类型）----
    // 用 setTo 反转分派：由来源类型决定如何高效写入（零运动直接 setZero）
    template<typename D2>
    Derived & operator=(const MotionBase<D2> & other)
    {
      other.derived().setTo(derived());
      return derived();
    }

    // ---- 赋值：从 6D 原始向量 [v; ω] 构造 ----
    // 注意分段顺序遵循 LINEAR=0 / ANGULAR=3，即线性在前
    template<typename V6>
    Derived & operator=(const Eigen::MatrixBase<V6> & v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(V6);
      assert(v.size() == 6);
      linear() = v.template segment<3>(LINEAR);
      angular() = v.template segment<3>(ANGULAR);
      return derived();
    }

    MotionPlain operator-() const
    {
      return derived().__opposite__();
    }
    template<typename M1>
    MotionPlain operator+(const MotionDense<M1> & v) const
    {
      return derived().__plus__(v.derived());
    }
    template<typename M1>
    MotionPlain operator-(const MotionDense<M1> & v) const
    {
      return derived().__minus__(v.derived());
    }

    template<typename M1>
    Derived & operator+=(const MotionDense<M1> & v)
    {
      return derived().__pequ__(v.derived());
    }
    template<typename M1>
    Derived & operator+=(const MotionBase<M1> & v)
    {
      v.derived().addTo(derived());
      return derived();
    }

    template<typename M1>
    Derived & operator-=(const MotionDense<M1> & v)
    {
      return derived().__mequ__(v.derived());
    }

    // ============================================================
    // 向量空间运算的实现体（se(3) 是线性空间，可逐分量加减数乘）
    // 双下划线命名是 Pinocchio 内部约定，由 MotionBase 的运算符转发过来
    // ============================================================
    MotionPlain __opposite__() const
    {
      return MotionPlain(-linear(), -angular());
    }

    template<typename M1>
    MotionPlain __plus__(const MotionDense<M1> & v) const
    {
      return MotionPlain(linear() + v.linear(), angular() + v.angular());
    }

    template<typename M1>
    MotionPlain __minus__(const MotionDense<M1> & v) const
    {
      return MotionPlain(linear() - v.linear(), angular() - v.angular());
    }

    // 就地版本：直接改写自身，不产生临时对象（热循环中优先用 += / -=）
    template<typename M1>
    Derived & __pequ__(const MotionDense<M1> & v)
    {
      linear() += v.linear();
      angular() += v.angular();
      return derived();
    }

    template<typename M1>
    Derived & __mequ__(const MotionDense<M1> & v)
    {
      linear() -= v.linear();
      angular() -= v.angular();
      return derived();
    }

    template<typename OtherScalar>
    MotionPlain __mult__(const OtherScalar & alpha) const
    {
      return MotionPlain(alpha * linear(), alpha * angular());
    }

    // 除法转乘倒数：一次除法 + 六次乘法，快于六次除法
    template<typename OtherScalar>
    MotionPlain __div__(const OtherScalar & alpha) const
    {
      return derived().__mult__((OtherScalar)(1) / alpha);
    }

    // ---- 与力的对偶配对：功率 P = φᵀν = f·v + τ·ω ----
    // 标量且与坐标系无关（功率不变性，见 §2.2.3）
    template<typename F1>
    Scalar dot(const ForceBase<F1> & phi) const
    {
      return phi.linear().dot(linear()) + phi.angular().dot(angular());
    }

    // ---- 叉乘入口：再次使用双分派 ----
    // 由对方（Motion 或 Force 或 Inertia）实现 motionAction，
    // 从而自动选中"运动叉乘"还是"力叉乘"，无需在此处分类讨论
    template<typename D>
    typename MotionAlgebraAction<D, Derived>::ReturnType cross_impl(const D & d) const
    {
      return d.motionAction(derived());
    }

    // ---- 运动叉乘 ν₁ × ν₂（李括号）----
    // 注意调用方向：cross_impl 中是 d.motionAction(*this)，
    // 故此处 this = ν₂（被作用者），参数 v = ν₁（施加者），结果 mout = ν₁ × ν₂：
    //
    //   (ν₁ × ν₂).v = v₁ × ω₂ + ω₁ × v₂
    //   (ν₁ × ν₂).ω = ω₁ × ω₂
    //
    // 实测与该公式吻合（残差 0），且满足反对称性 ν₁×ν₂ = −ν₂×ν₁。
    // 这是 RNEA 第一趟中速度积项 ν × Sq̇ 的来源（见 §4.3.2）
    template<typename M1, typename M2>
    void motionAction(const MotionDense<M1> & v, MotionDense<M2> & mout) const
    {
      mout.linear() = v.linear().cross(angular()) + v.angular().cross(linear());
      mout.angular() = v.angular().cross(angular());
    }

    template<typename M1>
    MotionPlain motionAction(const MotionDense<M1> & v) const
    {
      MotionPlain res;
      motionAction(v, res);
      return res;
    }

    template<typename M2>
    bool isApprox(
      const MotionDense<M2> & m2,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isApprox_impl(m2, prec);
    }

    template<typename D2>
    bool isApprox_impl(
      const MotionDense<D2> & m2,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return isApproxOrZero(linear(), m2.linear(), prec)
             && isApproxOrZero(angular(), m2.angular(), prec);
    }

    bool isZero_impl(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return linear().isZero(prec) && angular().isZero(prec);
    }

    // ---- 坐标变换（伴随 Ad_M）：ᵃν = ᵃX_b · ᵇν ----
    //   ᵃω = R·ᵇω                ← 角速度只旋转
    //   ᵃv = R·ᵇv + t × (R·ᵇω)   ← 线速度加牵连项
    //
    // 实现要点：先算 v.angular()，再在算 linear() 时【复用】它，
    //   避免第二次做 R·ω 乘法 —— 这就是先写 angular 行的原因。
    // noalias() 告诉 Eigen 左右无内存重叠，省去中间临时对象。
    template<typename S2, int O2, typename D2>
    void se3Action_impl(const SE3Tpl<S2, O2> & m, MotionDense<D2> & v) const
    {
      v.angular().noalias() = m.rotation() * angular();
      v.linear().noalias() = m.rotation() * linear() + m.translation().cross(v.angular());
    }

    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType se3Action_impl(const SE3Tpl<S2, O2> & m) const
    {
      typename SE3GroupAction<Derived>::ReturnType res;
      se3Action_impl(m, res);
      return res;
    }

    // ---- 逆变换 Ad_{M⁻¹}：ᵇν = (ᵃX_b)⁻¹ · ᵃν ----
    //   ᵇv = Rᵀ(ᵃv − t × ᵃω)
    //   ᵇω = Rᵀ·ᵃω
    // 注意必须先用【原始】的 angular() 算 linear()（此时 v 可能与 *this 是
    // 同一对象），故这里的书写顺序与正向变换相反
    template<typename S2, int O2, typename D2>
    void se3ActionInverse_impl(const SE3Tpl<S2, O2> & m, MotionDense<D2> & v) const
    {
      v.linear().noalias() =
        m.rotation().transpose() * (linear() - m.translation().cross(angular()));
      v.angular().noalias() = m.rotation().transpose() * angular();
    }

    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType
    se3ActionInverse_impl(const SE3Tpl<S2, O2> & m) const
    {
      typename SE3GroupAction<Derived>::ReturnType res;
      se3ActionInverse_impl(m, res);
      return res;
    }

    void disp_impl(std::ostream & os) const
    {
      os << "  v = " << linear().transpose() << std::endl
         << "  w = " << angular().transpose() << std::endl;
    }

    /// \returns a MotionRef on this.
    // 返回指向自身内存的轻量引用视图（不拷贝数据）
    MotionRefType ref()
    {
      return derived().ref();
    }

  protected:
    // 构造函数设为 protected：本类是抽象实现层，只能被派生类继承使用，
    // 不允许直接实例化；拷贝构造被删除以防止切片（slicing）
    MotionDense() {};

    MotionDense(const MotionDense &) = delete;

  }; // class MotionDense

  /// Basic operations specialization
  // ---- 运算符 ^ ：叉乘的简写形式（Featherstone 书中的记法）----
  // v1 ^ v2 = v1 × v2 （Motion × Motion → Motion）

  // 这里使用motionPlain作为返回类型，而不是用M1，避免M1是ref类型时返回ref类型，导致返回值是局部变量的引用，出现悬空引用
  template<typename M1, typename M2>
  typename traits<M1>::MotionPlain operator^(const MotionDense<M1> & v1, const MotionDense<M2> & v2)
  {
    return v1.derived().cross(v2.derived());
  }

  // v ^ f = v ×* f （Motion × Force → Force，对偶叉乘）
  template<typename M1, typename F1>
  typename traits<F1>::ForcePlain operator^(const MotionDense<M1> & v, const ForceBase<F1> & f)
  {
    return v.derived().cross(f.derived());
  }

  // 左乘标量 α·ν
  template<typename M1>
  typename traits<M1>::MotionPlain
  operator*(const typename traits<M1>::Scalar alpha, const MotionDense<M1> & v)
  {
    return v * alpha;
  }

} // namespace pinocchio
