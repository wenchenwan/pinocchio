//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2018-2025 INRIA
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

  // ============================================================
  // MotionBase：空间速度（Twist）的【CRTP 接口层】
  //
  //   ν = [v; ω] ∈ R⁶ ≅ se(3)   ★ Pinocchio 采用【线性在前】
  //     v = linear()  ：坐标系【原点处】的线速度（不是质心速度！）
  //     ω = angular() ：刚体角速度（与参考点无关）
  //
  // 刚体上任意点 P 的速度由速度场给出：v_P = v_O + ω × OP
  // 因此 linear() 是"刚体延拓后经过原点的虚拟质点"的速度。
  // 详见 PINOCCHIO_GUIDE.md §2.2.1。
  //
  // 同样是纯转发层：实现在 MotionDense / MotionTpl / MotionRef / MotionZero。
  // Motion 的继承链比 SE3 深一层：
  //   MotionBase → MotionDense（稠密存储的通用实现）→ MotionTpl / MotionRef
  //              → MotionZero （编译期已知为零的特化，运算可被完全优化掉）
  // ============================================================
  template<class Derived>
  class MotionBase : NumericalBase<Derived>
  {
  public:
    MOTION_TYPEDEF_TPL(Derived);

    // ---- CRTP 基础设施 ----
    Derived & derived()
    {
      return *static_cast<Derived *>(this);
    }
    const Derived & derived() const
    {
      return *static_cast<const Derived *>(this);
    }

    Derived & const_cast_derived() const
    {
      return *const_cast<Derived *>(&derived());
    }

    // ============================================================
    // 分量访问
    // ⚠️ linear() 是【原点处】线速度，不是质心速度：
    //      v_质心 = linear() + angular() × c   （c 为质心相对原点的位置）
    // ============================================================
    ConstAngularType angular() const // 只读 ω
    {
      return derived().angular_impl();
    }
    ConstLinearType linear() const // 只读 v
    {
      return derived().linear_impl();
    }
    AngularType angular() // 可写引用 ω
    {
      return derived().angular_impl();
    }
    LinearType linear() // 可写引用 v
    {
      return derived().linear_impl();
    }

    // 赋值版本：模板化以接受任意 3D Eigen 表达式（切片、运算结果等，零拷贝）
    template<typename V3Like>
    void angular(const Eigen::MatrixBase<V3Like> & w)
    {
      derived().angular_impl(w.derived());
    }

    template<typename V3Like>
    void linear(const Eigen::MatrixBase<V3Like> & v)
    {
      derived().linear_impl(v.derived());
    }

    // ---- 求值为"朴素"类型：把表达式/引用类型固化成 MotionTpl ----
    // 用于打断表达式模板链，避免悬空引用
    operator PlainReturnType() const
    {
      return derived().plain();
    }
    PlainReturnType plain() const
    {
      return derived().plain();
    }

    // ---- 转成 6D 列向量 [v; ω]，便于与稠密矩阵运算对接 ----
    ToVectorConstReturnType toVector() const
    {
      return derived().toVector_impl();
    }
    ToVectorReturnType toVector()
    {
      return derived().toVector_impl();
    }
    operator Vector6() const
    {
      return toVector();
    }

    // ---- 6×6 伴随作用矩阵 ad_ν：把叉乘写成矩阵形式 ----
    //   ad_ν · μ = ν × μ   （Motion × Motion）
    //   ad_ν = [ ω̂  v̂ ]
    //          [ 0   ω̂ ]
    ActionMatrixType toActionMatrix() const
    {
      return derived().toActionMatrix_impl();
    }
    // ---- 对偶版 ad_ν* = −ad_νᵀ：用于力的叉乘 ν ×* φ ----
    // RNEA 中的科氏力/陀螺项 ν ×* (Iν) 就用它（见 §4.3.3）
    ActionMatrixType toDualActionMatrix() const
    {
      return derived().toDualActionMatrix_impl();
    }
    operator Matrix6() const
    {
      return toActionMatrix();
    }

    /**
     * @brief The homogeneous representation of the motion vector \f$ \xi \f$.
     *
     * With \f$ \hat{\xi} = \left( \begin{array}{cc} \omega & v \\ 0 & 0 \\ \end{array} \right) \f$,
     * \f[
     * {}^a\dot{M}_b = \hat{\xi} {}^aM_b
     * \f]
     *
     * @note This function is provided for completeness, but it is not the best
     * way to use Motion quantities in terms of sparsity exploitation and
     * general efficiency. For integration, the recommended way is to use
     * Motion vectors along with the \ref integrate function.
     */
    // ---- 4×4 齐次表示 ξ̂ = [ω̂ v; 0 0]，满足 ᵃṀ_b = ξ̂·ᵃM_b ----
    // 即"位姿对时间的导数"。注释已说明：仅为完备性提供，
    // 实际积分应使用 integrate()（在流形上做指数映射，见 §2.1），
    // 而不是用这个矩阵做数值积分 —— 后者会离开 SE(3)。
    HomogeneousMatrixType toHomogeneousMatrix() const
    {
      return derived().toHomogeneousMatrix_impl();
    }

    void setZero()
    {
      derived().setZero();
    }

    // ---- 比较：精确相等（浮点逐元素）----
    template<typename M2>
    bool operator==(const MotionBase<M2> & other) const
    {
      return derived().isEqual_impl(other.derived());
    }

    template<typename M2>
    bool operator!=(const MotionBase<M2> & other) const
    {
      return !(derived() == other.derived());
    }

    // ============================================================
    // 向量空间运算：se(3) 是【线性空间】，可自由加减和数乘
    // （这正是"速度在切空间"的含义 —— 与位姿 SE(3) 只能做群乘法形成对比）
    // 双下划线命名 __plus__ 等是 Pinocchio 内部约定，避免与用户符号冲突
    // ============================================================
    Derived operator-() const // −ν
    {
      return derived().__opposite__();
    }
    Derived operator+(const MotionBase<Derived> & v) const // ν₁ + ν₂
    {
      return derived().__plus__(v.derived());
    }
    Derived operator-(const MotionBase<Derived> & v) const // ν₁ − ν₂
    {
      return derived().__minus__(v.derived());
    }
    Derived & operator+=(const MotionBase<Derived> & v) // 就地累加（省临时对象）
    {
      return derived().__pequ__(v.derived());
    }
    Derived & operator-=(const MotionBase<Derived> & v)
    {
      return derived().__mequ__(v.derived());
    }

    // 数乘 ν·α。返回类型经 RHSScalarMultiplication 萃取：
    // 支持标量类型提升（如 Motion<float> * double），亦兼容自动微分标量
    template<typename OtherScalar>
    typename internal::RHSScalarMultiplication<Derived, OtherScalar>::ReturnType
    operator*(const OtherScalar & alpha) const
    {
      return derived().__mult__(alpha);
    }

    template<typename OtherScalar>
    Derived operator/(const OtherScalar & alpha) const
    {
      return derived().__div__(alpha);
    }

    // ---- 李括号 / 空间叉乘 ----
    //   d 为 Motion → ν × μ ，结果为 Motion（运动叉乘）
    //   d 为 Force  → ν ×* φ，结果为 Force （力叉乘，对偶版本）
    // 返回类型由 MotionAlgebraAction 在编译期选出。
    // 这是 RNEA 中一切科氏力/离心力/陀螺力矩的来源（§4.3.2、§4.3.3）
    template<typename OtherSpatialType>
    typename MotionAlgebraAction<OtherSpatialType, Derived>::ReturnType
    cross(const OtherSpatialType & d) const
    {
      return derived().cross_impl(d);
    }

    // ---- 带容差比较（推荐用于浮点判等）----
    bool isApprox(
      const Derived & other,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isApprox_impl(other, prec);
    }

    bool isZero(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isZero_impl(prec);
    }

    // ---- 坐标变换：速度用【伴随】Ad_M ----
    //   ᵃν = ᵃX_b · ᵇν ⇒ ᵃv = R·ᵇv + t×(R·ᵇω)，ᵃω = R·ᵇω
    // 被 SE3::act(Motion) 经双分派调用（见 se3-tpl.hxx 的 act_impl）
    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType se3Action(const SE3Tpl<S2, O2> & m) const
    {
      return derived().se3Action_impl(m);
    }

    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType se3ActionInverse(const SE3Tpl<S2, O2> & m) const
    {
      return derived().se3ActionInverse_impl(m);
    }

    // ---- 与力的对偶配对：功率 P = φᵀν = f·v + τ·ω ----
    // 标量，且【与坐标系选择无关】—— 这正是 Motion/Force 互为对偶空间的定义
    // （功率不变性的推导见 §2.2.3）
    template<typename ForceDerived>
    Scalar dot(const ForceDense<ForceDerived> & f) const
    {
      return derived().dot(f.derived());
    }

    void disp(std::ostream & os) const
    {
      derived().disp_impl(os);
    }
    friend std::ostream & operator<<(std::ostream & os, const MotionBase<Derived> & v)
    {
      v.disp(os);
      return os;
    }

    ///
    /// \brief Returns the size of the Derived object in bytes.
    ///
    std::size_t sizeInBytes() const
    {
      return derived().sizeInBytes();
    }

  }; // class MotionBase

  // ---- 左乘标量：α·ν（成员 operator* 只能处理 ν·α，故补一个自由函数）----
  template<typename MotionDerived>
  typename internal::RHSScalarMultiplication<
    MotionDerived,
    typename MotionDerived::Scalar>::ReturnType
  operator*(const typename MotionDerived::Scalar & alpha, const MotionBase<MotionDerived> & motion)
  {
    return motion * alpha;
  }

} // namespace pinocchio
