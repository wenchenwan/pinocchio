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
  template<typename _Scalar, int _Options>
  struct traits<MotionTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, _Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, _Options> Vector6;
    typedef Eigen::Matrix<Scalar, 4, 4, _Options> Matrix4;
    typedef Eigen::Matrix<Scalar, 6, 6, _Options> Matrix6;
    typedef Matrix6 ActionMatrixType;
    typedef Matrix4 HomogeneousMatrixType;
    typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Vector6) ToVectorConstReturnType;
    typedef typename PINOCCHIO_EIGEN_REF_TYPE(Vector6) ToVectorReturnType;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type LinearType;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type AngularType;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstLinearType;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstAngularType;
    typedef MotionTpl<Scalar, _Options> MotionPlain;
    typedef const MotionPlain & PlainReturnType;
    static constexpr int LINEAR = 0;
    static constexpr int ANGULAR = 3;
    static constexpr int Options = _Options;

    typedef MotionRef<Vector6> MotionRefType;
  }; // traits MotionTpl

  template<typename _Scalar, int _Options>
  class MotionTpl : public MotionDense<MotionTpl<_Scalar, _Options>>
  {
  public:
    typedef MotionDense<MotionTpl> Base;
    MOTION_TYPEDEF_TPL(MotionTpl);
    static constexpr int Options = _Options;

    using Base::operator=;
    using Base::angular;
    using Base::linear;

    using Base::__mequ__;
    using Base::__minus__;
    using Base::__mult__;
    using Base::__opposite__;
    using Base::__pequ__;
    using Base::__plus__;

    // Constructors
    // ---- 默认构造：m_data【未初始化】（Eigen 默认不清零）----
    // 需要零运动请用 MotionTpl::Zero()
    MotionTpl()
    {
    }

    // ---- 由 (线速度 v, 角速度 ω) 两段 3D 向量构造 ----
    // ⚠️ 参数顺序是【线性在前】，与 Featherstone 原书相反，容易写反
    template<typename V1, typename V2>
    MotionTpl(const Eigen::MatrixBase<V1> & v, const Eigen::MatrixBase<V2> & w)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(V1);
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(V2);
      linear() = v;
      angular() = w;
    }

    // ---- 由单个 6D 向量 [v; ω] 构造 ----
    // explicit：避免 Vector6 意外隐式转成 Motion
    template<typename V6>
    explicit MotionTpl(const Eigen::MatrixBase<V6> & v)
    : m_data(v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(V6);
    }

    MotionTpl(const MotionTpl & other)
    {
      *this = other;
    }

    // ---- 跨标量类型构造（如 MotionTpl<float> → MotionTpl<double>）----
    // 与 SE3Tpl 不同：Motion 是线性空间，cast 后无需重新归一化
    template<typename S2, int O2>
    explicit MotionTpl(const MotionTpl<S2, O2> & other)
    {
      *this = other.template cast<Scalar>();
    }

    // ---- 同标量、不同 Eigen 对齐选项之间的转换 ----
    template<int O2>
    explicit MotionTpl(const MotionTpl<Scalar, O2> & clone)
    : m_data(clone.toVector())
    {
    }

    // Same explanation as converting constructor from MotionBase
    template<
      typename M2,
      std::enable_if_t<!std::is_convertible_v<MotionDense<M2>, MotionTpl>, bool> = true>
    explicit MotionTpl(const MotionDense<M2> & clone)
    {
      linear() = clone.linear();
      angular() = clone.angular();
    }

    // MotionBase implement a conversion function to PlainReturnType.
    // Usually, PlainReturnType is defined as MotionTpl.
    // In this case, this converting constructor is redundant and
    // create a warning with -Wconversion
    template<
      typename M2,
      std::enable_if_t<!std::is_convertible_v<MotionBase<M2>, MotionTpl>, bool> = true>
    explicit MotionTpl(const MotionBase<M2> & clone)
    {
      *this = clone;
    }

    ///
    /// \brief Copy assignment operator.
    ///
    /// \param[in] other MotionTpl to copy
    ///
    MotionTpl & operator=(const MotionTpl & clone) // Copy assignment operator
    {
      m_data = clone.toVector();
      return *this;
    }

    // initializers
    // ---- 零运动 ν = 0 ----
    // 注：若在编译期就知道是零，用 MotionZero 更优（运算可被完全优化掉）
    static MotionTpl Zero()
    {
      return MotionTpl(Vector6::Zero());
    }
    static MotionTpl Random()
    {
      return MotionTpl(Vector6::Random());
    }

    // 本类已是"朴素"类型，plain() 直接返回自身
    inline PlainReturnType plain() const
    {
      return *this;
    }

    // ---- 底层 6D 向量：返回引用，零拷贝 ----
    ToVectorConstReturnType toVector_impl() const
    {
      return m_data;
    }
    ToVectorReturnType toVector_impl()
    {
      return m_data;
    }

    // Getters
    // ---- 分量访问：对同一段 m_data 取【视图】(segment)，不复制数据 ----
    // 这就是为什么 v.linear() 可作左值：它是引用而非拷贝
    ConstAngularType angular_impl() const
    {
      return m_data.template segment<3>(ANGULAR); // 下标 3..5
    }
    ConstLinearType linear_impl() const
    {
      return m_data.template segment<3>(LINEAR); // 下标 0..2
    }
    AngularType angular_impl()
    {
      return m_data.template segment<3>(ANGULAR);
    }
    LinearType linear_impl()
    {
      return m_data.template segment<3>(LINEAR);
    }

    template<typename V3>
    void angular_impl(const Eigen::MatrixBase<V3> & w)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(V3, 3);
      angular_impl() = w;
    }
    template<typename V3>
    void linear_impl(const Eigen::MatrixBase<V3> & v)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(V3, 3);
      linear_impl() = v;
    }

    // Specific operators for MotionTpl and MotionRef
    // ============================================================
    // 针对 MotionTpl / MotionRef 的【专用运算重载】
    //
    // 为什么要覆盖 MotionDense 里已有的通用版本：
    //   通用版逐段操作 linear() 和 angular()（两次 3D 运算）；
    //   这里两者内存连续，可直接对整个 6D 向量 m_data 一次性运算 ——
    //   更利于 SIMD 向量化，也少一次函数调用。
    // 每个运算都有 MotionTpl 和 MotionRef 两个重载（后者内存可能不连续）。
    // ============================================================
    template<int O2>
    MotionPlain __plus__(const MotionTpl<Scalar, O2> & v) const
    {
      return MotionPlain(m_data + v.toVector());
    }

    template<typename Vector6ArgType>
    MotionPlain __plus__(const MotionRef<Vector6ArgType> & v) const
    {
      return MotionPlain(m_data + v.toVector());
    }

    template<int O2>
    MotionPlain __minus__(const MotionTpl<Scalar, O2> & v) const
    {
      return MotionPlain(m_data - v.toVector());
    }

    template<typename Vector6ArgType>
    MotionPlain __minus__(const MotionRef<Vector6ArgType> & v) const
    {
      return MotionPlain(m_data - v.toVector());
    }

    template<int O2>
    MotionTpl & __pequ__(const MotionTpl<Scalar, O2> & v)
    {
      m_data += v.toVector();
      return *this;
    }

    template<typename Vector6ArgType>
    MotionTpl & __pequ__(const MotionRef<Vector6ArgType> & v)
    {
      m_data += v.toVector();
      return *this;
    }

    template<int O2>
    MotionTpl & __mequ__(const MotionTpl<Scalar, O2> & v)
    {
      m_data -= v.toVector();
      return *this;
    }

    template<typename Vector6ArgType>
    MotionTpl & __mequ__(const MotionRef<Vector6ArgType> & v)
    {
      m_data -= v.toVector();
      return *this;
    }

    template<typename OtherScalar>
    MotionPlain __mult__(const OtherScalar & alpha) const
    {
      return MotionPlain(alpha * m_data);
    }

    // ---- 返回指向自身数据的轻量引用视图（不拷贝）----
    // 用途：把本对象的内存交给需要 MotionRef 接口的算法就地修改
    MotionRef<Vector6> ref()
    {
      return MotionRef<Vector6>(m_data);
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    // ---- 标量类型转换 ----
    // 与 SE3Tpl::cast 的区别：Motion 属于线性空间 se(3)，没有
    // "必须留在流形上"的约束，因此转换后无需归一化处理
    template<typename NewScalar>
    MotionTpl<NewScalar, Options> cast() const
    {
      typedef MotionTpl<NewScalar, Options> ReturnType;
      ReturnType res(linear().template cast<NewScalar>(), angular().template cast<NewScalar>());
      return res;
    }

    ///
    /// \brief Returns the size of the MotionTpl object in bytes.
    ///
    // double 情形：6 个标量 × 8 字节 = 48 字节
    static constexpr std::size_t sizeInBytes()
    {
      return sizeof(MotionTpl);
    }

  protected:
    // ---- 数据成员：单个连续的 6D 向量 [v(0:3); ω(3:6)] ----
    // 存成一整块（而非两个 Vector3）是为了 toVector() 能零拷贝返回，
    // 且整体运算可被 SIMD 向量化
    Vector6 m_data;

  }; // class MotionTpl

} // namespace pinocchio
