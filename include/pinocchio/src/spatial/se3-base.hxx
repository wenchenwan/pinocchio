//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2018-2025 INRIA
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial/se3.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial/se3.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  /** \brief Base class for rigid transformation.
   *
   * The rigid transform aMb can be seen in two ways:
   *
   * - given a point p expressed in frame B by its coordinate vector \f$ ^bp \f$, \f$ ^aM_b \f$
   * computes its coordinates in frame A by \f$ ^ap = {}^aM_b {}^bp \f$.
   * - \f$ ^aM_b \f$ displaces a solid S centered at frame A into the solid centered in
   * B. In particular, the origin of A is displaced at the origin of B:
   * \f$^aM_b {}^aA = {}^aB \f$.

   * The rigid displacement is stored as a rotation matrix and translation vector by:
   * \f$ ^aM_b x = {}^aR_b x + {}^aAB \f$
   * where \f$^aAB\f$ is the vector from origin A to origin B expressed in coordinates A.
   *
   * \cheatsheet \f$ {}^aM_c = {}^aM_b {}^bM_c \f$
   *
   */
  // ============================================================
  // SE3Base：刚体位姿的【CRTP 接口层】
  //
  // 本文件几乎不含计算逻辑 —— 每个函数都是一行转发：
  //     返回类型 foo() const { return derived().foo_impl(); }
  // 真正的实现在派生类 SE3Tpl（见 se3-tpl.hxx）。
  //
  // 为什么用 CRTP 而非虚函数：编译期静态分派，可内联、零运行时开销。
  // 代价是模板报错冗长（详见 PINOCCHIO_GUIDE.md §11.3）。
  //
  // 阅读建议：本文件当作"接口清单"看；想知道某函数怎么算的，
  //           去 se3-tpl.hxx 找同名的 xxx_impl。
  // ============================================================
  template<class Derived>
  struct SE3Base : NumericalBase<Derived>
  {
    PINOCCHIO_SE3_TYPEDEF_TPL(Derived);

    // ---- CRTP 基础设施：向下转型到具体类型，所有转发的基础 ----
    Derived & derived()
    {
      return *static_cast<Derived *>(this);
    }
    const Derived & derived() const
    {
      return *static_cast<const Derived *>(this);
    }

    // 在 const 成员函数里取得非 const 引用。
    // 用于"逻辑上 const 但需写入输出参数"的场景（如 toActionMatrix 的输出版）
    Derived & const_cast_derived() const
    {
      return *const_cast<Derived *>(&derived());
    }

    // ============================================================
    // 数据访问：位姿 M = (R, t)，R ∈ SO(3)，t ∈ R³
    // 三组重载：只读 / 可写引用 / 赋值
    // 返回的是 Ref（引用而非拷贝），故 M.rotation() *= 1.05 可直接改原对象
    // ============================================================
    ConstAngularRef rotation() const // 只读 R
    {
      return derived().rotation_impl();
    }
    ConstLinearRef translation() const // 只读 t
    {
      return derived().translation_impl();
    }
    AngularRef rotation() // 可写引用 R
    {
      return derived().rotation_impl();
    }
    LinearRef translation() // 可写引用 t
    {
      return derived().translation_impl();
    }
    void rotation(const AngularType & R) // 设置 R
    {
      derived().rotation_impl(R);
    }
    void translation(const LinearType & t) // 设置 t
    {
      derived().translation_impl(t);
    }

    // ============================================================
    // 四种矩阵表示：同一个 M，按作用对象不同转成不同矩阵
    //   4×4 齐次   → 作用于【点】
    //   6×6 Ad     → 作用于【空间速度 Motion】
    //   6×6 Ad⁻¹   → 速度的逆变换
    //   6×6 Ad⁻ᵀ   → 作用于【空间力 Force】
    // 数学推导见 PINOCCHIO_GUIDE.md §2.2.4
    // ============================================================

    // ---- 4×4 齐次矩阵 M = [R t; 0 1] ----
    HomogeneousMatrixType toHomogeneousMatrix() const
    {
      return derived().toHomogeneousMatrix_impl();
    }
    // 隐式转换：允许 Eigen::Matrix4d H = M;
    // ⚠️ 会静默产生矩阵拷贝，热循环中注意开销
    operator HomogeneousMatrixType() const
    {
      return toHomogeneousMatrix();
    }

    /**
     * @brief The action matrix \f$ {}^aX_b \f$ of \f$ {}^aM_b \f$.
     *
     * With \f$ {}^aM_b = \left( \begin{array}{cc} R & t \\ 0 & 1 \\ \end{array} \right) \f$,
     * \f[
     * {}^aX_b = \left( \begin{array}{cc} R & \hat{t} R \\ 0 & R \\ \end{array} \right)
     * \f]
     *
     * \cheatsheet \f$ {}^a\nu_c = {}^aX_b {}^b\nu_c \f$
     */
    // ---- 6×6 伴随矩阵 Ad_M：作用于空间速度 ᵃν = ᵃX_b·ᵇν ----
    ActionMatrixType toActionMatrix() const
    {
      return derived().toActionMatrix_impl();
    }
    // 隐式转换为 6×6（同样注意拷贝开销）
    operator ActionMatrixType() const
    {
      return toActionMatrix();
    }

    // 带输出参数的版本：写入调用者提供的缓冲区，避免返回临时 6×6 矩阵。
    // 形参是 const 引用却要写入 —— 靠 const_cast_derived() 实现，
    // 这是 Eigen 表达"输出参数"的惯用手法。
    template<typename Matrix6Like>
    void toActionMatrix(const Eigen::MatrixBase<Matrix6Like> & action_matrix) const
    {
      derived().toActionMatrix_impl(action_matrix);
    }

    /**
     * @brief The action matrix \f$ {}^bX_a \f$ of \f$ {}^aM_b \f$.
     * \sa toActionMatrix()
     */
    // ---- 6×6 逆伴随 Ad_{M⁻¹}：利用 SE(3) 结构解析求得，无需 6×6 求逆 ----
    ActionMatrixType toActionMatrixInverse() const
    {
      return derived().toActionMatrixInverse_impl();
    }

    template<typename Matrix6Like>
    void toActionMatrixInverse(const Eigen::MatrixBase<Matrix6Like> & action_matrix_inverse) const
    {
      derived().toActionMatrixInverse_impl(action_matrix_inverse.const_cast_derived());
    }

    // ---- 6×6 余伴随 Ad_M⁻ᵀ：作用于空间力（与速度互为对偶）----
    ActionMatrixType toDualActionMatrix() const
    {
      return derived().toDualActionMatrix_impl();
    }

    // ---- 打印：转发到派生类的 disp_impl ----
    void disp(std::ostream & os) const
    {
      static_cast<const Derived *>(this)->disp_impl(os);
    }

    template<typename Matrix6Like>
    void toDualActionMatrix(const Eigen::MatrixBase<Matrix6Like> & dual_action_matrix) const
    {
      derived().toDualActionMatrix_impl(dual_action_matrix);
    }

    // ============================================================
    // 群运算
    // ============================================================

    // ---- 变换复合：ᵃM_c = ᵃM_b · ᵇM_c ----
    template<typename OtherDerived>
    typename SE3GroupAction<Derived>::ReturnType operator*(const SE3Base<OtherDerived> & m2) const
    {
      return derived().act(m2.derived());
    }

    /// ay = aXb.act(by)
    // ---- 泛型群作用：d 可以是 SE3 / Motion / Force / Inertia / 3D 点 ----
    // 返回类型由 SE3GroupAction<D>::ReturnType 在编译期决定；
    // 内部经由"双分派"落到 d.se3Action(*this)（见 se3-tpl.hxx）
    template<typename D>
    typename SE3GroupAction<D>::ReturnType act(const D & d) const
    {
      return derived().act_impl(d);
    }

    /// by = aXb.actInv(ay)
    // ---- 逆作用：比 inverse().act() 快，不构造中间的逆变换对象 ----
    // IK 中的 data.oMi[JOINT_ID].actInv(oMdes) 走的就是这里（见 §4.2）
    template<typename D>
    typename SE3GroupAction<D>::ReturnType actInv(const D & d) const
    {
      return derived().actInv_impl(d);
    }

    // ============================================================
    // 比较运算
    // ============================================================

    // ⚠️ 精确相等（浮点逐元素 ==）。经过任何运算后几乎不可能成立，
    //    实际比较请用 isApprox()
    bool operator==(const Derived & other) const
    {
      return derived().isEqual(other);
    }

    bool operator!=(const Derived & other) const
    {
      return !(*this == other);
    }

    // 带容差的近似比较（推荐）。默认容差：double≈1e-12，float≈1e-5
    bool isApprox(
      const Derived & other,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isApprox_impl(other, prec);
    }

    // friend 函数：让 std::cout << M 可用，内部转发到 disp()
    friend std::ostream & operator<<(std::ostream & os, const SE3Base<Derived> & X)
    {
      X.disp(os);
      return os;
    }

    ///
    /// \returns true if *this is approximately equal to the identity placement, within the
    /// precision given by prec.
    ///
    // ---- 是否近似为单位变换（R≈I 且 t≈0）----
    // IK 收敛判据的底层：ᵢM_des → I 等价于 log6(ᵢM_des) → 0
    bool isIdentity(
      const typename traits<Derived>::Scalar & prec =
        Eigen::NumTraits<typename traits<Derived>::Scalar>::dummy_precision()) const
    {
      return derived().isIdentity(prec);
    }

    ///
    /// \returns true if the rotational part of *this is a rotation matrix (normalized columns),
    /// within the precision given by prec.
    ///
    // ---- 旋转部分是否仍在 SO(3) 上 ----
    // 反复矩阵乘法会让 R 因浮点误差缓慢偏离 SO(3)，需定期检查
    bool isNormalized(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return derived().isNormalized(prec);
    }

    ///
    /// \brief Normalize *this in such a way the rotation part of *this lies on SO(3).
    ///
    // ---- 就地把 R 投影回 SO(3)（取 Frobenius 范数下最接近的旋转矩阵）----
    void normalize()
    {
      derived().normalize();
    }

    ///
    /// \returns a Normalized version of *this, in such a way the rotation part of the returned
    /// transformation lies on SO(3).
    ///
    // ⚠️ 已知缺陷：此处【漏写 return】，函数声明返回 PlainType 却未返回任何值
    //    （对比上方 toHomogeneousMatrix() 等正确写法即可看出差异）。
    //    落到这里会得到未定义的返回值 —— 实测返回未归一化的垃圾数据。
    //
    //    平时不易暴露：SE3Tpl 自己声明了 normalized()，会"名字隐藏"掉本函数，
    //    因此 SE3 对象直接调用是正常的；只有通过 SE3Base<SE3>& 引用调用
    //    （例如写泛型模板代码时）才会踩到。
    //    修复方式：在 derived().normalized() 前补上 return。
    PlainType normalized() const
    {
      derived().normalized();
    }

    ///
    /// \brief Returns the size of the Derived object in bytes.
    ///
    // ---- 对象占用字节数（内存分析用）----
    std::size_t sizeInBytes() const
    {
      return derived().sizeInBytes();
    }

  }; // struct SE3Base

} // namespace pinocchio
