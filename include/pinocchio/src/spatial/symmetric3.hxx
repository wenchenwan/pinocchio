//
// Copyright (c) 2014-2021 CNRS INRIA
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
  struct traits<Symmetric3Tpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options>
  // ============================================================
  // Symmetric3Tpl：3×3 【对称】矩阵的压缩表示
  //
  // 用途：存储转动惯量张量 Ī（物理上必然对称）。
  //
  // 压缩存储：对称阵只有 6 个独立元素，故用一个 Vector6 存
  //     [ a₀  a₁  a₃ ]
  //     [ a₁  a₂  a₄ ]     ← 只存 a₀..a₅ 六个数
  //     [ a₃  a₄  a₅ ]
  // 相比存满 9 个数：省 1/3 内存，且所有运算（加法、合同变换、
  // 与向量相乘）都能利用对称性减少约一半的乘法次数。
  //
  // 这类"利用结构省算力"的做法贯穿整个 spatial 模块 ——
  // SE3 存 (R,t) 而非 4×4、Inertia 存 10 个参数而非 6×6，同一思路。
  // ============================================================
  class Symmetric3Tpl : public NumericalBase<Symmetric3Tpl<_Scalar, _Options>>
  {
  public:
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
    typedef Eigen::Matrix<Scalar, 2, 2, Options> Matrix2;
    typedef Eigen::Matrix<Scalar, 3, 2, Options> Matrix32;

  public:
    Symmetric3Tpl()
    {
    }

    template<typename Sc, int Opt>
    explicit Symmetric3Tpl(const Eigen::Matrix<Sc, 3, 3, Opt> & I)
    {
      assert(check_expression_if_real<Scalar>(pinocchio::isZero((I - I.transpose()))));
      m_data(0) = I(0, 0);
      m_data(1) = I(1, 0);
      m_data(2) = I(1, 1);
      m_data(3) = I(2, 0);
      m_data(4) = I(2, 1);
      m_data(5) = I(2, 2);
    }

    explicit Symmetric3Tpl(const Vector6 & I)
    : m_data(I)
    {
    }

    Symmetric3Tpl(const Symmetric3Tpl & other)
    {
      *this = other;
    }

    template<typename S2, int O2>
    explicit Symmetric3Tpl(const Symmetric3Tpl<S2, O2> & other)
    {
      *this = other.template cast<Scalar>();
    }

    ///
    /// \brief Copy assignment operator.
    ///
    /// \param[in] other Symmetric3 to copy
    ///
    Symmetric3Tpl & operator=(const Symmetric3Tpl & clone) // Copy assignment operator
    {
      m_data = clone.m_data;
      return *this;
    }

    Symmetric3Tpl(
      const Scalar & a0,
      const Scalar & a1,
      const Scalar & a2,
      const Scalar & a3,
      const Scalar & a4,
      const Scalar & a5)
    {
      m_data << a0, a1, a2, a3, a4, a5;
    }

    static Symmetric3Tpl Zero()
    {
      return Symmetric3Tpl(Vector6::Zero());
    }
    void setZero()
    {
      m_data.setZero();
    }

    static Symmetric3Tpl Random()
    {
      return RandomPositive();
    }
    void setRandom()
    {
      Scalar a = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             b = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             c = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             d = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             e = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             f = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0;

      m_data << a, b, c, d, e, f;
    }

    static Symmetric3Tpl Identity()
    {
      return Symmetric3Tpl(Scalar(1), Scalar(0), Scalar(1), Scalar(0), Scalar(0), Scalar(1));
    }
    void setIdentity()
    {
      m_data << Scalar(1), Scalar(0), Scalar(1), Scalar(0), Scalar(0), Scalar(1);
    }

    template<typename Vector3Like>
    void setDiagonal(const Eigen::MatrixBase<Vector3Like> & diag)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(Vector3Like, 3);
      m_data[0] = diag[0];
      m_data[2] = diag[1];
      m_data[5] = diag[2];
    }

    /* Required by Inertia::operator== */
    bool operator==(const Symmetric3Tpl & other) const
    {
      return m_data == other.m_data;
    }

    bool operator!=(const Symmetric3Tpl & other) const
    {
      return !(*this == other);
    }

    bool isApprox(
      const Symmetric3Tpl & other,
      const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return m_data.isApprox(other.m_data, prec);
    }

    bool isZero(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      return m_data.isZero(prec);
    }

    void fill(const Scalar value)
    {
      m_data.fill(value);
    }

    // ---- 3×3 对称矩阵求逆：伴随矩阵(cofactor) / 行列式，闭式解 ----
    // S⁻¹ = adj(S)/det(S)。因对称，只需算 6 个余子式（上三角），
    // 下三角对称填。det 用第一行按余子式展开：det = Σ a_1j·C_1j。
    template<typename Matrix3Like>
    void inverse(const Eigen::MatrixBase<Matrix3Like> & res_) const
    {
      Matrix3Like & res = res_.const_cast_derived();
      const Scalar &a11 = m_data[0], a21 = m_data[1], a22 = m_data[2], a31 = m_data[3],
                   a32 = m_data[4], a33 = m_data[5];

      // 余子式（伴随矩阵元素）：res(i,j) = C_ij，对称故上下三角一次算
      res(0, 0) = a33 * a22 - a32 * a32;
      res(1, 0) = res(0, 1) = -(a33 * a21 - a32 * a31);
      res(2, 0) = res(0, 2) = a32 * a21 - a22 * a31;
      res(1, 1) = a33 * a11 - a31 * a31;
      res(2, 1) = res(1, 2) = -(a32 * a11 - a21 * a31);
      res(2, 2) = a22 * a11 - a21 * a21;

      const Scalar det = a11 * res(0, 0) + a21 * res(0, 1) + a31 * res(0, 2); // 沿第一列展开
      res /= det;                                                            // S⁻¹ = adj/det
    }

    Matrix3 inverse() const
    {
      Matrix3 res;
      inverse(res);
      return res;
    }

    // ---- SkewSquare(v)：表示对称矩阵 [v]×² （叉乘平方，本身对称） ----
    // [v]×² = v vᵀ − ‖v‖²·I，展开各元即下面 6 个数。
    // 用途：平行轴定理项 —— 惯量在原点与质心之间平移时的修正 −m[c]×²。
    // 用轻量代理类（不立即算成 Symmetric3）+ operator- 融合，见下。
    struct SkewSquare
    {
      const Vector3 & v;
      SkewSquare(const Vector3 & v)
      : v(v)
      {
      }
      operator Symmetric3Tpl() const
      {
        const Scalar &x = v[0], &y = v[1], &z = v[2];
        // [v]×² 的 6 个独立元（对角为 −(其余两分量平方和)）
        return Symmetric3Tpl(-y * y - z * z, x * y, -x * x - z * z, x * z, y * z, -x * x - y * y);
      }
    }; // struct SkewSquare

    // S − [v]×²：把"减去叉乘平方"融合成一趟（不先物化 [v]×²），平行轴定理常用
    Symmetric3Tpl operator-(const SkewSquare & v) const
    {
      const Scalar &x = v.v[0], &y = v.v[1], &z = v.v[2];
      return Symmetric3Tpl(
        m_data[0] + y * y + z * z, m_data[1] - x * y, m_data[2] + x * x + z * z, m_data[3] - x * z,
        m_data[4] - y * z, m_data[5] + x * x + y * y);
    }

    Symmetric3Tpl & operator-=(const SkewSquare & v)
    {
      const Scalar &x = v.v[0], &y = v.v[1], &z = v.v[2];
      m_data[0] += y * y + z * z;
      m_data[1] -= x * y;
      m_data[2] += x * x + z * z;
      m_data[3] -= x * z;
      m_data[4] -= y * z;
      m_data[5] += x * x + y * y;
      return *this;
    }

    // ---- AlphaSkewSquare(m,v)：表示 m·[v]×²（带标量系数的叉乘平方）----
    // 平行轴项常带质量系数：−m[c]×²。同样用代理类 + operator- 融合，避免临时对象。
    struct AlphaSkewSquare
    {
      const Scalar & m;
      const Vector3 & v;

      AlphaSkewSquare(const Scalar & m, const SkewSquare & v)
      : m(m)
      , v(v.v)
      {
      }
      AlphaSkewSquare(const Scalar & m, const Vector3 & v)
      : m(m)
      , v(v)
      {
      }

      operator Symmetric3Tpl() const
      {
        const Scalar &x = v[0], &y = v[1], &z = v[2];
        return Symmetric3Tpl(
          -m * (y * y + z * z), m * x * y, -m * (x * x + z * z), m * x * z, m * y * z,
          -m * (x * x + y * y));
      }
    };

    friend AlphaSkewSquare operator*(const Scalar & m, const SkewSquare & sk)
    {
      return AlphaSkewSquare(m, sk);
    }

    Symmetric3Tpl operator-(const AlphaSkewSquare & v) const
    {
      const Scalar &x = v.v[0], &y = v.v[1], &z = v.v[2];
      return Symmetric3Tpl(
        m_data[0] + v.m * (y * y + z * z), m_data[1] - v.m * x * y,
        m_data[2] + v.m * (x * x + z * z), m_data[3] - v.m * x * z, m_data[4] - v.m * y * z,
        m_data[5] + v.m * (x * x + y * y));
    }

    Symmetric3Tpl & operator-=(const AlphaSkewSquare & v)
    {
      const Scalar &x = v.v[0], &y = v.v[1], &z = v.v[2];
      m_data[0] += v.m * (y * y + z * z);
      m_data[1] -= v.m * x * y;
      m_data[2] += v.m * (x * x + z * z);
      m_data[3] -= v.m * x * z;
      m_data[4] -= v.m * y * z;
      m_data[5] += v.m * (x * x + y * y);
      return *this;
    }

    const Vector6 & data() const
    {
      return m_data;
    }
    Vector6 & data()
    {
      return m_data;
    }

    // static Symmetric3Tpl SkewSq( const Vector3 & v )
    // {
    //   const Scalar & x = v[0], & y = v[1], & z = v[2];
    //   return Symmetric3Tpl(-y*y-z*z,
    //          x*y, -x*x-z*z,
    //          x*z, y*z, -x*x-y*y );
    // }

    /* Shoot a positive definite matrix. */
    static Symmetric3Tpl RandomPositive()
    {
      Scalar a = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             b = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             c = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             d = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             e = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0,
             f = Scalar(std::rand()) / RAND_MAX * 2.0 - 1.0;
      return Symmetric3Tpl(
        a * a + b * b + d * d, a * b + b * c + d * e, b * b + c * c + e * e, a * d + b * e + d * f,
        b * d + c * e + e * f, d * d + e * e + f * f);
    }

    // ---- 展开成完整 3×3 对称矩阵（需要显式矩阵时用；operator Matrix3() 隐式转换）----
    Matrix3 matrix() const
    {
      Matrix3 res;
      res(0, 0) = m_data(0);
      res(0, 1) = m_data(1);
      res(0, 2) = m_data(3);
      res(1, 0) = m_data(1);
      res(1, 1) = m_data(2);
      res(1, 2) = m_data(4);
      res(2, 0) = m_data(3);
      res(2, 1) = m_data(4);
      res(2, 2) = m_data(5);
      return res;
    }
    operator Matrix3() const
    {
      return matrix();
    }

    // ---- vtiv(v) = vᵀ S v（二次型，标量）----
    // 对称矩阵二次型 = 对角项 + 2×非对角项：
    //   Sxx x² + Syy y² + Szz z² + 2(Sxy xy + Sxz xz + Syz yz)
    // 直接用 6 个存储元算，不物化 3×3，运算量最小。用途：ωᵀ I_c ω（动能项）。
    Scalar vtiv(const Vector3 & v) const
    {
      const Scalar & x = v[0];
      const Scalar & y = v[1];
      const Scalar & z = v[2];

      const Scalar xx = x * x;
      const Scalar xy = x * y;
      const Scalar xz = x * z;
      const Scalar yy = y * y;
      const Scalar yz = y * z;
      const Scalar zz = z * z;

      return m_data(0) * xx + m_data(2) * yy + m_data(5) * zz            // 对角项
             + 2. * (m_data(1) * xy + m_data(3) * xz + m_data(4) * yz);  // 2×非对角项
    }

    ///
    /// \brief Performs the operation \f$ M = [v]_{\times} S_{3} \f$.
    ///        This operation is equivalent to applying the cross product of v on each column of S.
    ///
    /// \tparam Vector3, Matrix3
    ///
    /// \param[in]  v  a vector of dimension 3.
    /// \param[in]  S3 a symmetric matrix of dimension 3x3.
    /// \param[out] M  an output matrix of dimension 3x3.
    ///
    template<typename Vector3, typename Matrix3>
    static void vxs(
      const Eigen::MatrixBase<Vector3> & v,
      const Symmetric3Tpl & S3,
      const Eigen::MatrixBase<Matrix3> & M)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(Vector3, 3);
      PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(Matrix3, 3, 3);

      const Scalar & a = S3.data()[0];
      const Scalar & b = S3.data()[1];
      const Scalar & c = S3.data()[2];
      const Scalar & d = S3.data()[3];
      const Scalar & e = S3.data()[4];
      const Scalar & f = S3.data()[5];

      const typename Vector3::RealScalar & v0 = v[0];
      const typename Vector3::RealScalar & v1 = v[1];
      const typename Vector3::RealScalar & v2 = v[2];

      // M = [v]× S：对 S 的每一列做 v×（结果一般【不再对称】，故返回满 3×3）
      // 直接用 v 与 S 的 6 个元展开，不构造 [v]× 也不做通用矩阵乘
      Matrix3 & M_ = PINOCCHIO_EIGEN_CONST_CAST(Matrix3, M);
      M_(0, 0) = d * v1 - b * v2;
      M_(1, 0) = a * v2 - d * v0;
      M_(2, 0) = b * v0 - a * v1;

      M_(0, 1) = e * v1 - c * v2;
      M_(1, 1) = b * v2 - e * v0;
      M_(2, 1) = c * v0 - b * v1;

      M_(0, 2) = f * v1 - e * v2;
      M_(1, 2) = d * v2 - f * v0;
      M_(2, 2) = e * v0 - d * v1;
    }

    ///
    /// \brief Performs the operation \f$ [v]_{\times} S \f$.
    ///        This operation is equivalent to applying the cross product of v on each column of S.
    ///
    /// \tparam Vector3
    ///
    /// \param[in]  v  a vector of dimension 3.
    ///
    /// \returns the result \f$ [v]_{\times} S \f$.
    ///
    template<typename Vector3>
    Matrix3 vxs(const Eigen::MatrixBase<Vector3> & v) const
    {
      Matrix3 M;
      vxs(v, *this, M);
      return M;
    }

    ///
    /// \brief Performs the operation \f$ M = S_{3} [v]_{\times} \f$.
    ///
    /// \tparam Vector3, Matrix3
    ///
    /// \param[in]  v  a vector of dimension 3.
    /// \param[in]  S3 a symmetric matrix of dimension 3x3.
    /// \param[out] M  an output matrix of dimension 3x3.
    ///
    template<typename Vector3, typename Matrix3>
    static void svx(
      const Eigen::MatrixBase<Vector3> & v,
      const Symmetric3Tpl & S3,
      const Eigen::MatrixBase<Matrix3> & M)
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(Vector3, 3);
      PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(Matrix3, 3, 3);

      const Scalar & a = S3.data()[0];
      const Scalar & b = S3.data()[1];
      const Scalar & c = S3.data()[2];
      const Scalar & d = S3.data()[3];
      const Scalar & e = S3.data()[4];
      const Scalar & f = S3.data()[5];

      const typename Vector3::RealScalar & v0 = v[0];
      const typename Vector3::RealScalar & v1 = v[1];
      const typename Vector3::RealScalar & v2 = v[2];

      // M = S [v]×：对 S 的每一行做 ×v（= (−[v]× S)ᵀ）。同样返回满 3×3
      Matrix3 & M_ = PINOCCHIO_EIGEN_CONST_CAST(Matrix3, M);
      M_(0, 0) = b * v2 - d * v1;
      M_(1, 0) = c * v2 - e * v1;
      M_(2, 0) = e * v2 - f * v1;

      M_(0, 1) = d * v0 - a * v2;
      M_(1, 1) = e * v0 - b * v2;
      M_(2, 1) = f * v0 - d * v2;

      M_(0, 2) = a * v1 - b * v0;
      M_(1, 2) = b * v1 - c * v0;
      M_(2, 2) = d * v1 - e * v0;
    }

    /// \brief Performs the operation \f$ M = S_{3} [v]_{\times} \f$.
    ///
    /// \tparam Vector3
    ///
    /// \param[in]  v  a vector of dimension 3.
    ///
    /// \returns the result \f$ S [v]_{\times} \f$.
    ///
    template<typename Vector3>
    Matrix3 svx(const Eigen::MatrixBase<Vector3> & v) const
    {
      Matrix3 M;
      svx(v, *this, M);
      return M;
    }

    Symmetric3Tpl operator+(const Symmetric3Tpl & s2) const
    {
      return Symmetric3Tpl(m_data + s2.m_data);
    }

    Symmetric3Tpl operator-(const Symmetric3Tpl & s2) const
    {
      return Symmetric3Tpl(m_data - s2.m_data);
    }

    Symmetric3Tpl & operator+=(const Symmetric3Tpl & s2)
    {
      m_data += s2.m_data;
      return *this;
    }

    Symmetric3Tpl & operator-=(const Symmetric3Tpl & s2)
    {
      m_data -= s2.m_data;
      return *this;
    }

    Symmetric3Tpl & operator*=(const Scalar s)
    {
      m_data *= s;
      return *this;
    }

    // ---- rhsMult: vout = S · vin（对称矩阵 × 向量）----
    // 直接用 6 个存储元展开（3 行各 3 项），不物化 3×3。operator* 转调它。
    template<typename V3in, typename V3out>
    static void rhsMult(
      const Symmetric3Tpl & S3,
      const Eigen::MatrixBase<V3in> & vin,
      const Eigen::MatrixBase<V3out> & vout)
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(V3in, Vector3);
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(V3out, Vector3);

      V3out & vout_ = PINOCCHIO_EIGEN_CONST_CAST(V3out, vout);

      vout_[0] = S3.m_data(0) * vin[0] + S3.m_data(1) * vin[1] + S3.m_data(3) * vin[2]; // 行0
      vout_[1] = S3.m_data(1) * vin[0] + S3.m_data(2) * vin[1] + S3.m_data(4) * vin[2]; // 行1
      vout_[2] = S3.m_data(3) * vin[0] + S3.m_data(4) * vin[1] + S3.m_data(5) * vin[2]; // 行2
    }

    template<typename V3>
    Vector3 operator*(const Eigen::MatrixBase<V3> & v) const
    {
      Vector3 res;
      rhsMult(*this, v, res);
      return res;
    }

    // Matrix3 operator*(const Matrix3 &a) const
    // {
    //   Matrix3 r;
    //   for(unsigned int i=0; i<3; ++i)
    //     {
    //       r(0,i) = m_data(0) * a(0,i) + m_data(1) * a(1,i) + m_data(3) * a(2,i);
    //       r(1,i) = m_data(1) * a(0,i) + m_data(2) * a(1,i) + m_data(4) * a(2,i);
    //       r(2,i) = m_data(3) * a(0,i) + m_data(4) * a(1,i) + m_data(5) * a(2,i);
    //     }
    //   return r;
    // }

    // ---- (i,j) 元素访问：从 6 元紧凑存储反查 ----
    // 列主序上三角打包：索引 0,1,2,3,4,5 ↔ (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)。
    // i,j 都不为 2 时下标 = i+j；含 2 时需 +1 修正（跳过打包间隙）。
    const Scalar & operator()(const int i, const int j) const
    {
      return ((i != 2) && (j != 2)) ? m_data[i + j] : m_data[i + j + 1];
    }

    template<typename Matrix3Like>
    Symmetric3Tpl operator-(const Eigen::MatrixBase<Matrix3Like> & S) const
    {
      assert(check_expression_if_real<Scalar>(pinocchio::isZero(S - S.transpose())));
      return Symmetric3Tpl(
        m_data(0) - S(0, 0), m_data(1) - S(1, 0), m_data(2) - S(1, 1), m_data(3) - S(2, 0),
        m_data(4) - S(2, 1), m_data(5) - S(2, 2));
    }

    template<typename Matrix3Like>
    Symmetric3Tpl operator+(const Eigen::MatrixBase<Matrix3Like> & S) const
    {
      assert(check_expression_if_real<Scalar>(pinocchio::isZero(S - S.transpose())));
      return Symmetric3Tpl(
        m_data(0) + S(0, 0), m_data(1) + S(1, 0), m_data(2) + S(1, 1), m_data(3) + S(2, 0),
        m_data(4) + S(2, 1), m_data(5) + S(2, 2));
    }

    /* --- Symmetric R*S*R' and R'*S*R products --- */
  public: // private:
    /** \brief Computes L for a symmetric matrix A.
     */
    // decomposeltI：rotate 的辅助，把 S 拆出一个 3×2 中间量 L，
    // 使 R S Rᵀ 能用更少乘法算出（见 rotate 里的乘法/加法计数注释）。
    Matrix32 decomposeltI() const
    {
      Matrix32 L;
      L << m_data(0) - m_data(5), m_data(1), m_data(1), m_data(2) - m_data(5), 2 * m_data(3),
        m_data(4) + m_data(4);
      return L;
    }

    // ---- rotate(R): 合同变换 R S Rᵀ（把惯量从一个坐标系旋到另一个）----
    // 结果仍对称，故返回 Symmetric3。用 decomposeltI 分解 + 分块乘，
    // 比朴素 R·S·Rᵀ（两次 3×3 乘）省一半以上乘法（代码里标了各步 m/a 计数）。
    // 用途：Inertia::se3Action 里把 I_c 旋到新系（R I_c Rᵀ）。
    /* R*S*R' */
    template<typename D>
    Symmetric3Tpl rotate(const Eigen::MatrixBase<D> & R) const
    {
      PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(D, 3, 3);
      assert(
        check_expression_if_real<Scalar>(isUnitary(R.transpose() * R))
        && "R is not a Unitary matrix");

      Symmetric3Tpl Sres;

      // 4 a
      const Matrix32 L(decomposeltI());

      // Y = R' L   ===> (12 m + 8 a)
      const Matrix2 Y(R.template block<2, 3>(1, 0) * L);

      // Sres= Y R  ===> (16 m + 8a)
      Sres.m_data(1) = Y(0, 0) * R(0, 0) + Y(0, 1) * R(0, 1);
      Sres.m_data(2) = Y(0, 0) * R(1, 0) + Y(0, 1) * R(1, 1);
      Sres.m_data(3) = Y(1, 0) * R(0, 0) + Y(1, 1) * R(0, 1);
      Sres.m_data(4) = Y(1, 0) * R(1, 0) + Y(1, 1) * R(1, 1);
      Sres.m_data(5) = Y(1, 0) * R(2, 0) + Y(1, 1) * R(2, 1);

      // r=R' v ( 6m + 3a)
      const Vector3 r(
        -R(0, 0) * m_data(4) + R(0, 1) * m_data(3), -R(1, 0) * m_data(4) + R(1, 1) * m_data(3),
        -R(2, 0) * m_data(4) + R(2, 1) * m_data(3));

      // Sres_11 (3a)
      Sres.m_data(0) = L(0, 0) + L(1, 1) - Sres.m_data(2) - Sres.m_data(5);

      // Sres + D + (Ev)x ( 9a)
      Sres.m_data(0) += m_data(5);
      Sres.m_data(1) += r(2);
      Sres.m_data(2) += m_data(5);
      Sres.m_data(3) -= r(1);
      Sres.m_data(4) += r(0);
      Sres.m_data(5) += m_data(5);

      return Sres;
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    Symmetric3Tpl<NewScalar, Options> cast() const
    {
      return Symmetric3Tpl<NewScalar, Options>(m_data.template cast<NewScalar>());
    }

    friend std::ostream & operator<<(std::ostream & os, const Symmetric3Tpl<Scalar, Options> & S3)
    {
      os << "m_data: " << S3.m_data.transpose() << "\n";
      return os;
    }

    // TODO: adjust code
    //    bool isValid() const
    //    {
    //      return
    //         m_data(0) >= Scalar(0)
    //      && m_data(2) >= Scalar(0)
    //      && m_data(5) >= Scalar(0);
    //    }

  protected:
    Vector6 m_data;

  }; // class Symmetric3Tpl

} // namespace pinocchio
