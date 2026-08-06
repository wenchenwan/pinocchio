//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2018-2025 INRIA
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  /**
   * @brief      Base interface for forces representation.
   * @details    The Class implements all
   *
   *  This class hierarchy represents a spatial force, e.g. a spatial impulse or force associated to
   * a body. The spatial force is the mathematical representation of \f$ se^{*}(3) \f$, the dual of
   * \f$ se(3) \f$.
   *
   * @tparam     Derived  { description }
   */
  // ============================================================
  // ForceBase：空间力（Wrench）的【CRTP 接口层】
  //
  //   φ = [f; τ_O] ∈ R⁶ ≅ se*(3)   ★ Pinocchio 采用【线性在前】
  //     f   = linear()  ：合力（与约化点无关）
  //     τ_O = angular() ：对【坐标系原点】的合力矩（随原点变化！）
  //
  // 力系约化：一组力 {f_k}（作用点 P_k）等效到原点 O：
  //     f = Σ f_k ,  τ_O = Σ OP_k × f_k
  //
  // ⚠️ 常见误解：angular() 不是"纯力偶"，它包含
  //    真正的力偶 + 合力对原点的力矩臂效应。换原点时它会变：
  //        τ_C = τ_O − c × f     （c 为新原点相对旧原点的位置）
  //
  // 与 Motion 的对偶关系（本文件存在的根本原因）：
  //   功率 P = φᵀν = f·v + τ·ω 是标量，且与坐标系选择【无关】。
  //   正因功率必须不变，当速度按伴随 Ad 变换时，力必须按
  //   余伴随 Ad⁻ᵀ 变换 —— 这就是"力生活在速度的对偶空间 se*(3)"的含义。
  //   推导见 PINOCCHIO_GUIDE.md §2.2.2 / §2.2.3。
  //
  // 类层次与 Motion 完全平行：
  //   ForceBase → ForceDense → ForceTpl / ForceRef
  // ============================================================
  template<class Derived>
  class ForceBase : NumericalBase<Derived>
  {
  public:
    FORCE_TYPEDEF_TPL(Derived);

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

    /**
     * @brief      Return the angular part of the force vector
     *
     * @return     The 3D vector associated to the angular part of the 6D force vector
     */
    // 力矩 τ_O：对【坐标系原点】而言。⚠️ 换原点时该值会变（见类头说明）
    ConstAngularType angular() const
    {
      return derived().angular_impl();
    }

    /**
     * @brief      Return the linear part of the force vector
     *
     * @return     The 3D vector associated to the linear part of the 6D force vector
     */
    // 合力 f：与约化点无关，是刚体受力的固有属性
    ConstLinearType linear() const
    {
      return derived().linear_impl();
    }

    /// \copydoc ForceBase::angular
    AngularType angular()
    {
      return derived().angular_impl();
    }

    /// \copydoc ForceBase::linear
    LinearType linear()
    {
      return derived().linear_impl();
    }

    /**
     * @brief      Set the angular part of the force vector
     *
     * @tparam V3Like A vector 3 like type.
     *
     * @param[in]  n
     */
    template<typename V3Like>
    void angular(const Eigen::MatrixBase<V3Like> & n)
    {
      derived().angular_impl(n.derived());
    }

    /**
     * @brief      Set the linear part of the force vector
     *
     * @tparam V3Like A vector 3 like type.
     *
     * @param[in]  f
     */
    template<typename V3Like>
    void linear(const Eigen::MatrixBase<V3Like> & f)
    {
      derived().linear_impl(f.derived());
    }

    /**
     * @brief      Return the force as an Eigen vector.
     *
     * @return     The 6D vector \f$ \phi \f$ such that
     * \f{equation*}
     * {}^{A}\phi = \begin{bmatrix} {}^{A}f \\  {}^{A}\tau \end{bmatrix}
     * \f}
     */
    ToVectorConstReturnType toVector() const
    {
      return derived().toVector_impl();
    }

    /// \copydoc ForceBase::toVector
    ToVectorReturnType toVector()
    {
      return derived().toVector_impl();
    }

    /*
     * @brief C-style cast operator
     * \copydoc ForceBase::toVector
     */
    operator Vector6() const
    {
      return toVector();
    }

    /** \returns true if each coefficients of \c *this and \a other are all exactly equal.
     * \warning When using floating point scalar values you probably should rather use a
     *          fuzzy comparison such as isApprox()
     */
    template<typename F2>
    bool operator==(const ForceBase<F2> & other) const
    {
      return derived().isEqual_impl(other.derived());
    }

    /** \returns true if at least one coefficient of \c *this and \a other does not match.
     */
    template<typename F2>
    bool operator!=(const ForceBase<F2> & other) const
    {
      return !(derived() == other.derived());
    }

    /** \returns true if *this is approximately equal to other, within the precision given by prec.
     */
    bool isApprox(
      const Derived & other,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isApprox_impl(other, prec);
    }

    /** \returns true if the component of the linear and angular part of the Spatial Force are
     * approximately equal to zero, within the precision given by prec.
     */
    bool isZero(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isZero_impl(prec);
    }

    /** \brief Copies the Derived Force into *this
     *  \return a reference to *this
     */
    Derived & operator=(const ForceBase<Derived> & other)
    {
      return derived().setFrom(other.derived());
    }

    /**
     * \brief Replaces *this by *this + other.
     * \return a reference to *this
     */
    Derived & operator+=(const ForceBase<Derived> & phi)
    {
      return derived().__pequ__(phi.derived());
    }

    /**
     * \brief Replaces *this by *this - other.
     * \return a reference to *this
     */
    Derived & operator-=(const ForceBase<Derived> & phi)
    {
      return derived().__mequ__(phi.derived());
    }

    /** \return an expression of the sum of *this and other
     */
    Derived operator+(const ForceBase<Derived> & phi) const
    {
      return derived().__plus__(phi.derived());
    }

    /** \return an expression of *this scaled by the factor alpha
     */
    template<typename OtherScalar>
    ForcePlain operator*(const OtherScalar & alpha) const
    {
      return derived().__mult__(alpha);
    }

    /** \return an expression of *this divided by the factor alpha
     */
    template<typename OtherScalar>
    ForcePlain operator/(const OtherScalar & alpha) const
    {
      return derived().__div__(alpha);
    }

    /** \return an expression of the opposite of *this
     */
    Derived operator-() const
    {
      return derived().__opposite__();
    }

    /** \return an expression of the difference of *this and phi
     */
    Derived operator-(const ForceBase<Derived> & phi) const
    {
      return derived().__minus__(phi.derived());
    }

    /** \return the dot product of *this with m     *
     */
    // ---- 与速度的对偶配对：功率 P = φᵀν = f·v + τ·ω ----
    // 标量，且【与坐标系选择无关】（实测变换前后残差 8.3e-17）。
    // 这正是 Force 与 Motion 互为对偶空间的定义性质，
    // 也是虚功原理 τ_joint = Sᵀf（RNEA 投影步骤）的依据
    template<typename MotionDerived>
    Scalar dot(const MotionDense<MotionDerived> & m) const
    {
      return derived().dot(m.derived());
    }

    /**
     * @brief      Transform from A to B coordinates the Force represented by *this such that
     *             \f{equation*}
     *             {}^{B}f  =  {}^{B}X_A^* * {}^{A}f
     *             \f}
     *
     * @param[in]  m     The rigid transformation \f$ {}^{B}m_A \f$ whose coordinates transform for
     * forces is
     *                   {}^{B}X_A^*
     *
     * @return     an expression of the force expressed in the new coordinates
     */
    // ---- 坐标变换：力用【余伴随】Ad_M⁻ᵀ（注意不是伴随 Ad_M！）----
    //   ᴮf = R·ᴬf                ← 合力只旋转
    //   ᴮτ = R·ᴬτ + t × (R·ᴬf)   ← 力矩加力矩臂效应（经典的"力的平移定理"）
    // 与 Motion::se3Action 对比：t× 项作用在【力】上而非角速度上，
    // 这是伴随矩阵转置的直接结果（见 §2.2.4）
    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType se3Action(const SE3Tpl<S2, O2> & m) const
    {
      return derived().se3Action_impl(m);
    }

    /**
     * @brief      Transform from B to A coordinates the Force represented by *this such that
     *             \f{equation*}
     *             {}^{A}f  =  {}^{A}X_B^* * {}^{A}f
     *             \f}
     *
     * @param[in]  m     The rigid transformation \f$ {}^{B}m_A \f$ whose coordinates transform for
     * forces is
     *                   {}^{B}X_A^*
     *
     * @return     an expression of the force expressed in the new coordinates
     */
    template<typename S2, int O2>
    typename SE3GroupAction<Derived>::ReturnType se3ActionInverse(const SE3Tpl<S2, O2> & m) const
    {
      return derived().se3ActionInverse_impl(m);
    }

    // ---- 被速度叉乘作用：ν ×* φ（对偶叉乘）----
    // 由 Motion::cross(Force) 经双分派调用到这里。
    // RNEA 中连杆受力递推的 ν ×* (Iν) 项即由此计算（见 §4.3.3）
    template<typename M1>
    typename MotionAlgebraAction<Derived, M1>::ReturnType
    motionAction(const MotionDense<M1> & v) const
    {
      return derived().motionAction(v.derived());
    }

    void disp(std::ostream & os) const
    {
      derived().disp_impl(os);
    }
    friend std::ostream & operator<<(std::ostream & os, const ForceBase<Derived> & X)
    {
      X.disp(os);
      return os;
    }

    ///
    /// \brief Returns the size of the Derived object in bytes.
    ///
    std::size_t sizeInBytes() const
    {
      return derived().sizeInBytes();
    }

  }; // class ForceBase

} // namespace pinocchio
