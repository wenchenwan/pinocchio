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
  // ============================================================
  // SE3Tpl：刚体位姿 M = (R, t) 的【实现层】
  //
  // 与 se3-base.hxx 的关系：
  //   SE3Base<Derived>  = CRTP 接口层，每个函数都是 return derived().xxx_impl()
  //   SE3Tpl            = 本文件，提供所有 _impl 的真正实现
  //   → 想看"有哪些接口"读 se3-base.hxx；想看"怎么算的"读本文件
  //
  // 数学约定（详见 PINOCCHIO_GUIDE.md §2.2）：
  //   位姿   M = (R, t)，R ∈ SO(3)，t ∈ R³
  //   点变换 ᵃp = R·ᵇp + t
  //   复合   ᵃM_c = ᵃM_b · ᵇM_c
  // ============================================================

  template<typename _Scalar, int _Options>
  struct traits<SE3Tpl<_Scalar, _Options>>
  {
    static constexpr int Options = _Options;
    // 6D 空间向量中线性/角度分量的起始下标：Pinocchio 采用【线性在前】
    //   Motion ν = [v(0:3); ω(3:6)]，Force φ = [f(0:3); τ(3:6)]
    // 注意：这两个常量在本文件里还被复用为"0 和 3 的通用偏移量"
    //   （如 4×4 齐次矩阵中 block(LINEAR, ANGULAR) 其实是 block(0,3) = 平移列）
    static constexpr int LINEAR = 0;
    static constexpr int ANGULAR = 3;
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 4, 1, Options> Vector4;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
    typedef Eigen::Matrix<Scalar, 4, 4, Options> Matrix4;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef Matrix3 AngularType;
    typedef typename PINOCCHIO_EIGEN_REF_TYPE(Matrix3) AngularRef;
    typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Matrix3) ConstAngularRef;
    typedef Vector3 LinearType;
    typedef typename PINOCCHIO_EIGEN_REF_TYPE(Vector3) LinearRef;
    typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Vector3) ConstLinearRef;
    typedef Matrix6 ActionMatrixType;
    typedef Matrix4 HomogeneousMatrixType;
    typedef SE3Tpl<Scalar, Options> PlainType;
  }; // traits SE3Tpl

  template<typename _Scalar, int _Options>
  struct SE3Tpl : public SE3Base<SE3Tpl<_Scalar, _Options>>
  {

    PINOCCHIO_SE3_TYPEDEF_TPL(SE3Tpl);
    typedef SE3Base<SE3Tpl<_Scalar, _Options>> Base;
    typedef Eigen::Quaternion<Scalar, Options> Quaternion;
    typedef typename traits<SE3Tpl>::Vector3 Vector3;
    typedef typename traits<SE3Tpl>::Matrix3 Matrix3;
    typedef typename traits<SE3Tpl>::Matrix4 Matrix4;
    typedef typename traits<SE3Tpl>::Vector4 Vector4;
    typedef typename traits<SE3Tpl>::Matrix6 Matrix6;

    // 把基类的 rotation()/translation() 重载集引入本类作用域
    // （否则本类自己的 rotation_impl 等会因"名字隐藏"遮蔽基类同名函数）
    using Base::rotation;
    using Base::translation;

    // ---- 默认构造：注意 rot/trans 均【未初始化】（Eigen 默认行为，不清零）----
    // 想要单位变换请用 SE3Tpl::Identity()
    SE3Tpl()
    : rot()
    , trans() {};

    // ---- 由四元数 + 平移构造（最常用：URDF/消息里的位姿通常是四元数）----
    // quat.matrix() 内部完成四元数→旋转矩阵的转换
    template<typename QuaternionLike, typename Vector3Like>
    SE3Tpl(
      const Eigen::QuaternionBase<QuaternionLike> & quat,
      const Eigen::MatrixBase<Vector3Like> & trans)
    : rot(quat.matrix())
    , trans(trans)
    {
      // 编译期断言：第二个参数必须是 3 维向量（尺寸错误在编译期就报出来）
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(Vector3Like, 3);
    }

    // ---- 由旋转矩阵 + 平移构造 ----
    // 不校验 R 是否真的属于 SO(3)；若来源可疑，构造后应调用 normalize()
    template<typename Matrix3Like, typename Vector3Like>
    SE3Tpl(const Eigen::MatrixBase<Matrix3Like> & R, const Eigen::MatrixBase<Vector3Like> & trans)
    : rot(R)
    , trans(trans)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(Vector3Like, 3);
      PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(Matrix3Like, 3, 3);
    }

    // ---- 拷贝构造 ----
    SE3Tpl(const SE3Tpl & other)
    {
      *this = other;
    }

    // ---- 跨标量类型构造（如 SE3Tpl<float> → SE3Tpl<double>）----
    // explicit：禁止隐式转换，避免精度悄悄丢失
    // 转换经由 cast<>()，其中可能触发旋转部分的重新归一化（见下方 cast()）
    template<typename S2, int O2>
    explicit SE3Tpl(const SE3Tpl<S2, O2> & other)
    {
      *this = other.template cast<Scalar>();
    }

    // ---- 由 4×4 齐次矩阵构造 ----
    // block(LINEAR,LINEAR)=block(0,0) → 左上 3×3 旋转块
    // block(LINEAR,ANGULAR)=block(0,3) → 右上 3×1 平移列
    // （此处 LINEAR/ANGULAR 被当作 0/3 的通用偏移量复用，与 6D 分量含义无关）
    template<typename Matrix4Like>
    explicit SE3Tpl(const Eigen::MatrixBase<Matrix4Like> & m)
    : rot(m.template block<3, 3>(LINEAR, LINEAR))
    , trans(m.template block<3, 1>(LINEAR, ANGULAR))
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(Matrix4Like, 4, 4);
    }

    // ---- 构造单位变换的内部技巧 ----
    // 参数 int 只是"标签"，不使用其值；用来与默认构造函数区分开
    // 外部请用 SE3Tpl::Identity()（它就是 return SE3Tpl(1)）
    explicit SE3Tpl(int)
    : rot(AngularType::Identity())
    , trans(LinearType::Zero())
    {
    }

    // ---- 同标量、不同 Eigen 内存对齐选项（Options）之间的转换 ----
    template<int O2>
    SE3Tpl(const SE3Tpl<Scalar, O2> & clone)
    : rot(clone.rotation())
    , trans(clone.translation())
    {
    }

    // ---- 赋值：同标量、跨对齐选项 ----
    template<int O2>
    SE3Tpl & operator=(const SE3Tpl<Scalar, O2> & other)
    {
      rot = other.rotation();
      trans = other.translation();
      return *this;
    }

    ///
    /// \brief Copy assignment operator.
    ///
    /// \param[in] other SE3 to copy
    ///
    SE3Tpl & operator=(const SE3Tpl & other)
    {
      rot = other.rotation();
      trans = other.translation();
      return *this;
    }

    // ---- 单位变换 I = (I₃, 0)：群的幺元 ----
    static SE3Tpl Identity()
    {
      return SE3Tpl(1); // 调用上面那个"int 标签"构造函数
    }

    // 就地置为单位变换（复用已有对象，避免重新分配）
    SE3Tpl & setIdentity()
    {
      rot.setIdentity();
      trans.setZero();
      return *this;
    }

    /// aXb = bXa.inverse()
    // ---- 逆变换：M⁻¹ = (Rᵀ, −Rᵀt) ----
    // 关键：SE(3) 的逆有【解析式】，只需一次转置 + 一次矩阵向量乘，
    //       无需真正做 4×4 矩阵求逆（后者既慢又有数值误差）。
    // 推导：M⁻¹ 应满足 M⁻¹(Rp+t) = p ⇒ Rᵀ(Rp+t) − Rᵀt = p ✓
    SE3Tpl inverse() const
    {
      return SE3Tpl(rot.transpose(), -rot.transpose() * trans);
    }

    // ---- 随机位姿（主要用于单元测试/基准测试）----
    static SE3Tpl Random()
    {
      return SE3Tpl().setRandom();
    }

    // 就地随机化。
    // 注意旋转部分用 quaternion::uniformRandom：在 SO(3) 上【均匀】采样，
    // 而不是对三个欧拉角各自均匀取值（后者在球面上分布不均，会偏向两极）。
    SE3Tpl & setRandom()
    {
      Quaternion q;
      quaternion::uniformRandom(q);
      rot = q.matrix();
      trans.setRandom(); // 平移在 [-1,1]³ 内均匀
      return *this;
    }

    // ============================================================
    // 四种矩阵表示（对应 se3-base.hxx 中的四个 to*Matrix 接口）
    // 同一个 M，按作用对象不同转成不同矩阵：
    //   4×4 齐次矩阵   → 作用于【点】
    //   6×6 伴随 Ad    → 作用于【空间速度 Motion】
    //   6×6 Ad⁻¹       → 逆向变换速度
    //   6×6 余伴随 Ad⁻ᵀ → 作用于【空间力 Force】
    // ============================================================

    // ---- 4×4 齐次矩阵：M = [R t; 0 1]，作用于齐次坐标点 ----
    HomogeneousMatrixType toHomogeneousMatrix_impl() const
    {
      HomogeneousMatrixType M;
      M.template block<3, 3>(LINEAR, LINEAR) = rot;     // 左上 3×3 = R
      M.template block<3, 1>(LINEAR, ANGULAR) = trans;  // 右上 3×1 = t
      M.template block<1, 3>(ANGULAR, LINEAR).setZero();// 左下 1×3 = 0
      M(3, 3) = 1;
      return M;
    }

    /// Vb.toVector() = bXa.toMatrix() * Va.toVector()
    // ---- 6×6 伴随矩阵 Ad_M：作用于空间速度（Motion / Twist）----
    //
    //   ᵃX_b = [ R   t̂R ]      按 Pinocchio 的 [v; ω] 顺序
    //          [ 0   R  ]
    //
    // 展开即 §2.2.4 的两行：
    //   ᵃv = R·ᵇv + t × (R·ᵇω)   ← 线速度多了牵连项 t×ω
    //   ᵃω = R·ᵇω                ← 角速度只旋转
    //
    // 实现技巧：不显式构造反对称矩阵 t̂ 再相乘（那要 27 次乘法），
    //   而是逐列做叉乘 t × R.col(i)，等价于 (t̂R).col(i)，更快也更省内存。
    //
    // 形参声明成 const 引用却要写入 → 用 const_cast_derived() 取消 const。
    // 这是 Eigen 表达"输出参数"的惯用手法（见 §11 关于模板推导的说明）。
    template<typename Matrix6Like>
    void toActionMatrix_impl(const Eigen::MatrixBase<Matrix6Like> & action_matrix) const
    {
      typedef Eigen::Block<Matrix6Like, 3, 3> Block3;

      Matrix6Like & M = action_matrix.const_cast_derived();
      // 对角两块都是 R（一次链式赋值同时写两块）
      M.template block<3, 3>(ANGULAR, ANGULAR) = M.template block<3, 3>(LINEAR, LINEAR) = rot;
      M.template block<3, 3>(ANGULAR, LINEAR).setZero();   // 左下块恒为 0
      Block3 B = M.template block<3, 3>(LINEAR, ANGULAR);  // 右上块 = t̂R

      B.col(0) = trans.cross(rot.col(0));
      B.col(1) = trans.cross(rot.col(1));
      B.col(2) = trans.cross(rot.col(2));
    }

    // 返回值版本：内部分配临时对象再调用上面的就地版本
    // 热循环中建议用带输出参数的版本，避免每次构造 6×6 临时矩阵
    ActionMatrixType toActionMatrix_impl() const
    {
      ActionMatrixType res;
      toActionMatrix_impl(res);
      return res;
    }

    // ---- 6×6 逆伴随矩阵 Ad_{M⁻¹} = (Ad_M)⁻¹ ----
    //
    // 由 M⁻¹ = (Rᵀ, −Rᵀt) 代入伴随公式：
    //   Ad_{M⁻¹} = [ Rᵀ   −(R̂ᵀt)Rᵀ ]
    //              [ 0     Rᵀ       ]
    //
    // 为什么单独提供而不去求 Ad_M 的逆：
    //   直接对 6×6 求逆是 O(6³) 且有数值误差；这里利用 SE(3) 结构
    //   只需转置 + 三次叉乘，实测 ‖X·X⁻¹ − I‖ ≈ 5.3e-16。
    //
    // 右上块的推导（用到 R(a×b) = (Ra)×(Rb)，即旋转与叉乘可交换）：
    //   Rᵀ(eᵢ × t) = (Rᵀeᵢ) × (Rᵀt) = −(Rᵀt) × (Rᵀeᵢ)
    //              = 矩阵 −(R̂ᵀt)Rᵀ 的第 i 列  ✓
    // 所以下面每列算的是 Rᵀ·(eᵢ × t)。
    template<typename Matrix6Like>
    void
    toActionMatrixInverse_impl(const Eigen::MatrixBase<Matrix6Like> & action_matrix_inverse) const
    {
      typedef Eigen::Block<Matrix6Like, 3, 3> Block3;

      Matrix6Like & M = action_matrix_inverse.const_cast_derived();
      // 对角两块都是 Rᵀ
      M.template block<3, 3>(ANGULAR, ANGULAR) = M.template block<3, 3>(LINEAR, LINEAR) =
        rot.transpose();
      // 左下块最终应为 0，但这里【先借用作临时缓冲区】存放叉乘中间结果，
      // 算完右上块后再 setZero() 清掉 —— 省下一个 3×3 临时变量
      Block3 C = M.template block<3, 3>(ANGULAR, LINEAR); // used as temporary
      Block3 B = M.template block<3, 3>(LINEAR, ANGULAR);

// 宏展开为两步：
//   ① v3_out = e_{axis_id} × v3_in   （CartesianAxis 特化：与单位轴叉乘只需取负和换位，无乘法）
//   ② res.col(axis_id) = Rᵀ · v3_out （noalias 告诉 Eigen 无别名，省去临时对象）
#define PINOCCHIO_INTERNAL_COMPUTATION(axis_id, v3_in, v3_out, R, res)                             \
  CartesianAxis<axis_id>::cross(v3_in, v3_out);                                                    \
  res.col(axis_id).noalias() = R.transpose() * v3_out;

      // 三次调用都复用 C.col(0) 作临时列 —— 上一次的值已被消费，可安全覆盖
      PINOCCHIO_INTERNAL_COMPUTATION(0, trans, C.col(0), rot, B);
      PINOCCHIO_INTERNAL_COMPUTATION(1, trans, C.col(0), rot, B);
      PINOCCHIO_INTERNAL_COMPUTATION(2, trans, C.col(0), rot, B);

#undef PINOCCHIO_INTERNAL_COMPUTATION

      C.setZero(); // 还原左下块为 0（此前被借用作临时缓冲）
    }

    ActionMatrixType toActionMatrixInverse_impl() const
    {
      ActionMatrixType res;
      toActionMatrixInverse_impl(res);
      return res;
    }

    // ---- 6×6 余伴随矩阵 Ad_M⁻ᵀ：作用于空间力（Force / Wrench）----
    //
    //   ᵃX_b* = [ R    0 ]      按 [f; τ] 顺序
    //           [ t̂R   R ]
    //
    // 展开：
    //   ᵃf = R·ᵇf                ← 合力只旋转
    //   ᵃτ = R·ᵇτ + t × (R·ᵇf)   ← 力矩要加力矩臂效应（经典力学的"力的平移定理"）
    //
    // 与 toActionMatrix 的对称性：t̂R 从【右上】搬到了【左下】
    //   —— 这正是矩阵转置的结果，是"力与速度互为对偶"的直接体现。
    // 实测 ‖Ad* − Ad⁻ᵀ‖ = 0（严格相等），且功率 φᵀν 在变换下不变（残差 8.3e-17）。
    template<typename Matrix6Like>
    void toDualActionMatrix_impl(const Eigen::MatrixBase<Matrix6Like> & dual_action_matrix) const
    {
      typedef Eigen::Block<Matrix6Like, 3, 3> Block3;

      Matrix6Like & M = dual_action_matrix.const_cast_derived();
      M.template block<3, 3>(ANGULAR, ANGULAR) = M.template block<3, 3>(LINEAR, LINEAR) = rot;
      M.template block<3, 3>(LINEAR, ANGULAR).setZero();   // 右上块为 0（与速度版相反）
      Block3 B = M.template block<3, 3>(ANGULAR, LINEAR);  // 左下块 = t̂R

      B.col(0) = trans.cross(rot.col(0));
      B.col(1) = trans.cross(rot.col(1));
      B.col(2) = trans.cross(rot.col(2));
    }

    ActionMatrixType toDualActionMatrix_impl() const
    {
      ActionMatrixType res;
      toDualActionMatrix_impl(res);
      return res;
    }

    // ---- 打印：由基类的 disp() 和 operator<< 转发到这里 ----
    void disp_impl(std::ostream & os) const
    {
      os << "  R =\n" << rot << std::endl << "  p = " << trans.transpose() << std::endl;
    }

    /// --- GROUP ACTIONS ON M6, F6 and I6 ---
    // ============================================================
    // 群作用：把 M 施加到各种几何对象上
    //
    // 关键设计——【双分派】：
    //   SE3 不需要知道如何作用于每一种类型；而是反过来，
    //   每种类型自己实现 se3Action() / se3ActionInverse()。
    //   传 Motion  → 走伴随   Ad
    //   传 Force   → 走余伴随 Ad⁻ᵀ
    //   传 Inertia → 走惯量的合同变换
    //   返回类型由 SE3GroupAction<D>::ReturnType 在编译期决定。
    // 这样新增几何类型时无需改动 SE3Tpl。
    // ============================================================

    /// ay = aXb.act(by)
    // 泛型正向作用：ᵃy = ᵃM_b · ᵇy
    template<typename D>
    typename SE3GroupAction<D>::ReturnType act_impl(const D & d) const
    {
      return d.se3Action(*this); // ← 反向调用：对象自己知道该怎么变换
    }

    /// by = aXb.actInv(ay)
    // 泛型逆向作用：ᵇy = (ᵃM_b)⁻¹ · ᵃy
    // 比先 inverse() 再 act() 快：各类型的 se3ActionInverse 直接用 Rᵀ，
    // 不构造中间的逆变换对象
    template<typename D>
    typename SE3GroupAction<D>::ReturnType actInv_impl(const D & d) const
    {
      return d.se3ActionInverse(*this);
    }

    // ---- 作用于纯 Eigen 对象（3D 点，或 3×N 点集）----
    // 与下面 act_impl(Vector3) 的区别：这里接受任意 Eigen 表达式/矩阵，
    // 可一次变换整列点云；.eval() 强制求值，避免返回悬空的表达式模板
    template<typename EigenDerived>
    typename EigenDerived::PlainObject
    actOnEigenObject(const Eigen::MatrixBase<EigenDerived> & p) const
    {
      return (rotation() * p + translation()).eval(); // p ↦ R·p + t
    }

    // Eigen::Map 版重载：Map 没有 PlainObject，故显式返回 Vector3
    template<typename MapDerived>
    Vector3 actOnEigenObject(const Eigen::MapBase<MapDerived> & p) const
    {
      return Vector3(rotation() * p + translation());
    }

    // 逆作用：p ↦ Rᵀ(p − t)。注意是先减平移再转置旋转，顺序不能反
    template<typename EigenDerived>
    typename EigenDerived::PlainObject
    actInvOnEigenObject(const Eigen::MatrixBase<EigenDerived> & p) const
    {
      return (rotation().transpose() * (p - translation())).eval();
    }

    template<typename MapDerived>
    Vector3 actInvOnEigenObject(const Eigen::MapBase<MapDerived> & p) const
    {
      return Vector3(rotation().transpose() * (p - translation()));
    }

    // ---- 作用于单个 3D 点（Vector3 的非模板重载，优先级高于泛型版本）----
    // 点变换含平移项 t；而向量/速度变换不含 —— 这是点与向量的本质区别
    Vector3 act_impl(const Vector3 & p) const
    {
      return Vector3(rotation() * p + translation());
    }

    Vector3 actInv_impl(const Vector3 & p) const
    {
      return Vector3(rotation().transpose() * (p - translation()));
    }

    // ---- 作用于另一个 SE3：变换复合 ᵃM_c = ᵃM_b · ᵇM_c ----
    // 展开 [R₁ t₁]·[R₂ t₂] = [R₁R₂  R₁t₂+t₁]
    //      [0  1 ] [0  1 ]   [0     1      ]
    // 只做 3×3 乘法 + 一次矩阵向量乘，比 4×4 矩阵相乘省约 40% 运算量
    template<int O2>
    SE3Tpl act_impl(const SE3Tpl<Scalar, O2> & m2) const
    {
      return SE3Tpl(rot * m2.rotation(), translation() + rotation() * m2.translation());
    }

    // 逆复合：M₁⁻¹·M₂ = [R₁ᵀR₂,  R₁ᵀ(t₂−t₁)]
    // IK 中的 data.oMi[JOINT_ID].actInv(oMdes) 走的就是这条路径（见 §4.2）
    template<int O2>
    SE3Tpl actInv_impl(const SE3Tpl<Scalar, O2> & m2) const
    {
      return SE3Tpl(
        rot.transpose() * m2.rotation(), rot.transpose() * (m2.translation() - translation()));
    }

    // ---- 精确相等（逐元素 ==，浮点直接比较）----
    // 由基类 operator== 转发而来。
    // ⚠️ 经过任何浮点运算后几乎不可能精确相等，实际比较请用 isApprox()
    template<int O2>
    bool isEqual(const SE3Tpl<Scalar, O2> & m2) const
    {
      return (rotation() == m2.rotation() && translation() == m2.translation());
    }

    // ---- 带容差的近似相等（推荐使用）----
    // 默认容差 dummy_precision()：double 约 1e-12，float 约 1e-5
    template<int O2>
    bool isApprox_impl(
      const SE3Tpl<Scalar, O2> & m2,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return pinocchio::isApprox(rotation(), m2.rotation(), prec)
             && pinocchio::isApprox(translation(), m2.translation(), prec);
    }

    // ---- 是否近似为单位变换（R≈I 且 t≈0）----
    // IK 收敛判据的底层：ᵢM_des → I 等价于 log6(ᵢM_des) → 0
    bool isIdentity(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return rotation().isIdentity(prec) && translation().isZero(prec);
    }

    // ============================================================
    // 数据访问的 _impl（由基类 rotation()/translation() 转发）
    // 返回 Ref 类型（引用而非拷贝），故 M.rotation() *= 1.05 能改到原对象
    // ============================================================
    ConstAngularRef rotation_impl() const  // 只读取 R
    {
      return rot;
    }
    AngularRef rotation_impl()             // 可写引用 R
    {
      return rot;
    }
    void rotation_impl(const AngularType & R) // 设置 R（不校验是否属于 SO(3)）
    {
      rot = R;
    }
    ConstLinearRef translation_impl() const   // 只读取 t
    {
      return trans;
    }
    LinearRef translation_impl()              // 可写引用 t
    {
      return trans;
    }
    void translation_impl(const LinearType & p) // 设置 t
    {
      trans = p;
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    // ---- 标量类型转换：double ↔ float ↔ long double ↔ 自动微分标量 ----
    // 这是 §11.2 所讲"换标量不改算法"的落地点：整个 Pinocchio 的
    // autodiff / codegen / 多精度能力都建立在各类型的 cast<>() 之上。
    template<typename NewScalar>
    SE3Tpl<NewScalar, Options> cast() const
    {
      typedef SE3Tpl<NewScalar, Options> ReturnType;
      ReturnType res(rot.template cast<NewScalar>(), trans.template cast<NewScalar>());

      // During the cast, it may appear that the matrix is not normalized correctly.
      // Force the normalization of the rotation part of the matrix.
      // 仅当"旧标量的误差在新标量下可见"时才归一化，见文件末尾的特化
      internal::cast_call_normalize_method<SE3Tpl, NewScalar, Scalar>::run(res);
      return res;
    }

    // ---- 旋转部分是否仍在 SO(3) 上（列正交且单位长）----
    // 反复矩阵乘法会让 R 因浮点误差缓慢偏离 SO(3)，需要定期检查/修复
    bool isNormalized(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return isUnitary(rot, prec);
    }

    // ---- 就地把 R 投影回 SO(3) ----
    // orthogonalProjection 取与 R 最接近的旋转矩阵（SVD/极分解意义下：
    //   R = UΣVᵀ ⇒ 投影结果为 UVᵀ），是 Frobenius 范数下的最优正交近似
    void normalize()
    {
      rot = orthogonalProjection(rot);
    }

    // ---- 返回规范化后的副本（不修改自身）----
    // ⚠️ 注意：基类 SE3Base::normalized()（se3-base.hxx:219）转发时
    //    漏写了 return，经由 SE3Base<SE3>& 引用调用会拿到未定义的返回值。
    //    直接用 SE3Tpl 对象调用（走本函数）则正常 —— 因为本函数会
    //    "名字隐藏"掉基类同名版本。
    PlainType normalized() const
    {
      PlainType res(*this);
      res.normalize();
      return res;
    }

    ///
    /// \brief Linear interpolation on the SE3 manifold.
    ///
    /// \param[in] A Initial transformation.
    /// \param[in] B Target transformation.
    /// \param[in] alpha Interpolation factor in [0 ... 1].
    ///
    /// \returns An interpolated transformation between A and B.
    ///
    /// \note This is similar to the SLERP operation which acts initially for rotation but applied
    /// here to rigid transformation.
    ///
    // ---- SE(3) 流形上的测地线插值 ----
    //   Interpolate(A, B, α) = A · exp(α · log(A⁻¹B))
    // α=0 返回 A，α=1 返回 B，中间沿测地线（最短路径）平滑过渡。
    // 不能对 R 做逐元素线性插值 —— 那会离开 SO(3)（结果不再是旋转矩阵）。
    // 这里只有【声明】，定义在 se3-tpl-interpolate.hxx（依赖 explog，需后置）
    template<typename OtherScalar>
    static SE3Tpl Interpolate(const SE3Tpl & A, const SE3Tpl & B, const OtherScalar & alpha);

    ///
    /// \brief Returns the size of the SE3Tpl object in bytes.
    ///
    // constexpr：编译期常量。对 double 为 12 个标量（9+3）× 8 字节 = 96 字节
    static constexpr std::size_t sizeInBytes()
    {
      return sizeof(SE3Tpl);
    }

  protected:
    // ---- 数据成员：只存 R 和 t，不存 4×4 齐次矩阵 ----
    // 存 3×3+3×1（12 个标量）而非 4×4（16 个），省内存也省运算
    // —— 齐次矩阵最后一行恒为 [0 0 0 1]，没有存储价值
    AngularType rot;   // 旋转 R ∈ SO(3)
    LinearType trans;  // 平移 t ∈ R³

  }; // class SE3Tpl

  namespace internal
  {
    // ============================================================
    // cast<>() 后是否需要重新归一化旋转部分？—— 编译期分派
    // ============================================================

    // 特化①：标量类型没变（Scalar → Scalar），什么都不做（零开销）
    template<typename Scalar, int Options>
    struct cast_call_normalize_method<SE3Tpl<Scalar, Options>, Scalar, Scalar>
    {
      template<typename T>
      static void run(T &)
      {
      }
    };

    // 特化②：标量类型改变时，按精度关系决定是否归一化
    //
    // 判据：cast<New>(eps_old) > eps_new ?
    //   含义是"旧类型的固有误差，在新类型的精度下是否仍然可见"。
    //
    //   float → double：eps_old=1.2e-7 > eps_new=2.2e-16 → 【归一化】
    //     float 矩阵积累的 1e-7 级误差，到了 double 下就成了显著的非正交性，
    //     实测：raw float 矩阵在 double 容差下 isUnitary=false，
    //           经 cast<double>() 后 isNormalized=true ✓
    //
    //   double → float：eps_old=2.2e-16 > eps_new=1.2e-7 ? 否 → 【跳过】
    //     旧误差远小于 float 本身的舍入误差，归一化没有意义，白费算力。
    template<typename Scalar, int Options, typename NewScalar>
    struct cast_call_normalize_method<SE3Tpl<Scalar, Options>, NewScalar, Scalar>
    {
      template<typename T>
      static void run(T & self)
      {
        if (
          pinocchio::cast<NewScalar>(Eigen::NumTraits<Scalar>::epsilon())
          > Eigen::NumTraits<NewScalar>::epsilon())
          self.normalize();
      }
    };

  } // namespace internal

} // namespace pinocchio
