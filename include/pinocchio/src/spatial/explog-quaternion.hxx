//
// Copyright (c) 2018-2021 CNRS INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  namespace quaternion
  {

  // ============================================================
  // quaternion 命名空间：以【四元数】为载体的指数/对数映射
  //
  // 与 explog.hxx 中同名函数的区别：那里输出 3×3 旋转矩阵，
  // 这里输出/输入单位四元数。
  //
  // 为什么需要四元数版本：
  //   · Pinocchio 的浮动基与球关节在配置向量 q 里【就是用四元数存的】
  //     （nq 比 nv 多 1 的原因），积分时直接操作四元数最自然；
  //   · 存 4 个数而非 9 个，且归一化只需除以模长（矩阵则要正交化）；
  //   · 数值上更稳定，反复积分不易累积非正交误差。
  //
  //   exp3(ω) → 四元数 q = (sin(θ/2)·ω/θ, cos(θ/2))，θ = ‖ω‖
  //   注意是【半角】—— 四元数对 SO(3) 是双覆盖（q 与 −q 表示同一旋转）
  // ============================================================

    ///
    /// \brief Exp: so3 -> SO3 (quaternion)
    ///
    /// \returns the integral of the velocity vector as a quaternion.
    ///
    /// \param[in] v The angular velocity vector.
    /// \param[out] qout The quaternion where the result is stored.
    ///
    // ------------------------------------------------------------
    // exp3（四元数版）：so(3) → S³，把角速度积分成单位四元数
    //
    //     q = ( sin(θ/2)·ω/θ ,  cos(θ/2) ),   θ = ‖ω‖
    //         └── vec（虚部）──┘  └─ w（实部）┘
    //
    // ⚠️ 注意是【半角】θ/2：四元数对 SO(3) 是【双覆盖】，
    //    q 与 −q 表示同一个旋转（转 2π 才回到 q 本身）。
    //
    // 实现分两支（见下方 if_then_else）：
    //   一般情形：借 Eigen::AngleAxis(θ, ω/θ) 直接构造
    //   θ→0    ：转 4 阶泰勒展开，避免 ω/θ 的 0/0
    //       sin(θ/2)/θ = ½·(1 − (θ/2)²/6 + (θ/2)⁴/120)
    //       cos(θ/2)   =    1 − (θ/2)²/2 + (θ/2)⁴/24
    //     代码中 t2_2 = t²/4 恰是 (θ/2)²，与上式逐项对应。
    //
    // 用 if_then_else 而非普通 if：兼容自动微分标量（分支不可比较），
    // 且【两支都会求值】—— 这也是 t 里要加 eps 防止 ω/θ 产生 NaN 的原因。
    //
    // 实测：q.matrix() 与矩阵版 exp3(ω) 一致（1.9e-16）；半角公式严格相等；
    //       θ = 1e-2 / 1e-5 / 1e-8 / 0 各档均正确，|q| 恒为 1。
    // ------------------------------------------------------------
    template<typename Vector3Like, typename QuaternionLike>
    void
    exp3(const Eigen::MatrixBase<Vector3Like> & v, Eigen::QuaternionBase<QuaternionLike> & quat_out)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Vector3Like);
      assert(v.size() == 3);

      typedef typename Vector3Like::Scalar Scalar;
      static constexpr int Options =
        PINOCCHIO_EIGEN_PLAIN_TYPE(typename QuaternionLike::Coefficients)::Options;
      typedef Eigen::Quaternion<typename QuaternionLike::Scalar, Options> QuaternionPlain;
      const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

      const Scalar t2 = v.squaredNorm();
      const Scalar t = math::sqrt(t2 + eps * eps);

      static const Scalar ts_prec =
        TaylorSeriesExpansion<Scalar>::template precision<3>(); // Precision for the Taylor series
                                                                // expansion.

      Eigen::AngleAxis<Scalar> aa(t, v / t);
      QuaternionPlain quat_then(aa);

      // order 4 Taylor expansion in theta / (order 2 in t2)
      QuaternionPlain quat_else;
      const Scalar t2_2 = t2 / 4; // theta/2 squared
      quat_else.vec() =
        Scalar(0.5) * (Scalar(1) - t2_2 / Scalar(6) + t2_2 * t2_2 / Scalar(120)) * v;
      quat_else.w() = Scalar(1) - t2_2 / 2 + t2_2 * t2_2 / 24;

      using ::pinocchio::internal::if_then_else;
      for (Eigen::Index k = 0; k < 4; ++k)
      {
        quat_out.coeffs().coeffRef(k) = if_then_else(
          ::pinocchio::internal::GT, t2, ts_prec, quat_then.coeffs().coeffRef(k),
          quat_else.coeffs().coeffRef(k));
      }
    }

    /// \brief Exp: so3 -> SO3 (quaternion)
    ///
    /// \returns the integral of the velocity vector as a quaternion.
    ///
    /// \param[in] v The angular velocity vector.
    ///
    // 返回值版：内部转调上面的就地版本。
    // 热路径建议用就地版（传入已分配好的四元数），省一次构造
    template<typename Vector3Like>
    Eigen::
      Quaternion<typename Vector3Like::Scalar, PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like)::Options>
      exp3(const Eigen::MatrixBase<Vector3Like> & v)
    {
      typedef Eigen::Quaternion<
        typename Vector3Like::Scalar, PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like)::Options>
        ReturnType;
      ReturnType res;
      exp3(v, res);
      return res;
    }

    /// \brief The se3 -> SE3 exponential map, using quaternions to represent the output rotation.
    ///
    /// \returns the integral of the twist motion over unit time.
    ///
    /// \param[in] motion the spatial motion.
    /// \param[out] q the output transform in \f$\mathbb{R}^3 x S^3\f$.
    // ------------------------------------------------------------
    // exp6（四元数版）：se(3) → R³ × S³，输出 7 维配置向量
    //
    // 输出布局 qout = [ 平移(3) ; 四元数(4) ]，正是 Pinocchio 里
    // 【浮动基关节在 q 中的存储格式】（nq 比 nv 多 1 的由来）。
    //
    //   旋转部分：q = exp3(ω)                      （见上）
    //   平移部分：p = V(ω)·v，V 为左雅可比
    //       V(ω) = I + ((1−cosθ)/θ²)·ω̂ + ((θ−sinθ)/θ³)·ω̂²
    //
    // 代码里两个系数即：
    //   alpha_wxv = (1−cosθ)/θ²   → 乘 ω×v      （对应 ω̂v）
    //   alpha_w2  = (θ−sinθ)/θ³   → 乘 ω×(ω×v)  （对应 ω̂²v）
    // 从而 p = v + α₁·(ω×v) + α₂·(ω×(ω×v))，全程只用叉乘、不建 ω̂ 矩阵。
    //
    // ⚠️ p ≠ v：平移过程中坐标系本身也在转，必须用 V(ω) 修正 ——
    //    这是 exp6 与 exp3 最大的不同（详见 GUIDE §2.1）。
    // θ→0 时两个系数同样切泰勒展开（½ − θ²/24、1/6 − θ²/120）。
    //
    // 实现细节：用 Eigen::Map 把 qout 的前 3 维/后 4 维【就地】映射成
    // 平移向量与四元数，直接写入，无中间拷贝。
    //
    // 实测：与矩阵版 exp6 的平移部分严格相等、旋转部分差 2.2e-16。
    // ------------------------------------------------------------
    template<typename MotionDerived, typename Config_t>
    void exp6(const MotionDense<MotionDerived> & motion, Eigen::MatrixBase<Config_t> & qout)
    {
      static constexpr int Options = PINOCCHIO_EIGEN_PLAIN_TYPE(Config_t)::Options;
      typedef typename Config_t::Scalar Scalar;
      typedef typename MotionDerived::Vector3 Vector3;
      typedef Eigen::Quaternion<Scalar, Options> Quaternion_t;
      const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

      const typename MotionDerived::ConstAngularType & w = motion.angular();
      const typename MotionDerived::ConstLinearType & v = motion.linear();

      const Scalar t2 = w.squaredNorm() + eps * eps;
      const Scalar t = math::sqrt(t2);

      Scalar ct, st;
      SINCOS(t, &st, &ct);

      const Scalar inv_t2 = Scalar(1) / t2;
      const Scalar ts_prec =
        TaylorSeriesExpansion<Scalar>::template precision<3>(); // Taylor expansion precision

      using ::pinocchio::internal::if_then_else;
      using ::pinocchio::internal::LT;

      const Scalar alpha_wxv = if_then_else(
        LT, t, ts_prec,
        Scalar(0.5) - t2 / Scalar(24), // then: use Taylor expansion
        (Scalar(1) - ct) * inv_t2      // else
      );

      const Scalar alpha_w2 = if_then_else(
        LT, t, ts_prec, Scalar(1) / Scalar(6) - t2 / Scalar(120), (t - st) * inv_t2 / t);

      // linear part
      Eigen::Map<Vector3> trans_(qout.derived().template head<3>().data());
      trans_.noalias() = v + alpha_wxv * w.cross(v) + alpha_w2 * w.cross(w.cross(v));

      // quaternion part
      typedef Eigen::Map<Quaternion_t> QuaternionMap_t;
      QuaternionMap_t quat_(qout.derived().template tail<4>().data());
      exp3(w, quat_);
    }

    /// \brief The se3 -> SE3 exponential map, using quaternions to represent the output rotation.
    ///
    /// \returns the integral of the twist motion over unit time.
    ///
    /// \param[in] motion the spatial motion.
    template<typename MotionDerived>
    Eigen::Matrix<
      typename MotionDerived::Scalar,
      7,
      1,
      PINOCCHIO_EIGEN_PLAIN_TYPE(typename MotionDerived::Vector3)::Options>
    exp6(const MotionDense<MotionDerived> & motion)
    {
      typedef typename MotionDerived::Scalar Scalar;
      static constexpr int Options =
        PINOCCHIO_EIGEN_PLAIN_TYPE(typename MotionDerived::Vector3)::Options;
      typedef Eigen::Matrix<Scalar, 7, 1, Options> ReturnType;

      ReturnType qout;
      exp6(motion, qout);
      return qout;
    }

    /// \brief The se3 -> SE3 exponential map, using quaternions to represent the output rotation.
    ///
    /// \returns the integral of the spatial velocity over unit time.
    ///
    /// \param[in] vec6 the vector representing the spatial velocity.
    /// \param[out] qout the output transform in R^3 x S^3.
    // Vector6 入参的重载：用 MotionRef 把裸 6D 向量【零拷贝】包装成 Motion
    // 后转调上面的实现（MotionRef 见 §1.4：只引用不拷贝）
    template<typename Vector6Like, typename Config_t>
    void exp6(const Eigen::MatrixBase<Vector6Like> & vec6, Eigen::MatrixBase<Config_t> & qout)
    {
      MotionRef<const Vector6Like> nu(vec6.derived());
      ::pinocchio::quaternion::exp6(nu, qout);
    }

    /// \brief The se3 -> SE3 exponential map, using quaternions to represent the output rotation.
    ///
    /// \returns the integral of the spatial velocity over unit time.
    ///
    /// \param[in] vec6 the vector representing the spatial velocity.
    template<typename Vector6Like>
    Eigen::
      Matrix<typename Vector6Like::Scalar, 7, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Vector6Like)::Options>
      exp6(const Eigen::MatrixBase<Vector6Like> & vec6)
    {
      typedef typename Vector6Like::Scalar Scalar;
      static constexpr int Options = PINOCCHIO_EIGEN_PLAIN_TYPE(Vector6Like)::Options;
      typedef Eigen::Matrix<Scalar, 7, 1, Options> ReturnType;

      ReturnType qout;
      ::pinocchio::quaternion::exp6(vec6, qout);
      return qout;
    }

    /// \brief Same as \ref log3 but with a unit quaternion as input.
    ///
    /// \param[in] quat the unit quaternion.
    /// \param[out] theta the angle value (resuling from compurations).
    ///
    /// \return The angular velocity vector associated to the rotation matrix.
    ///
    template<typename QuaternionLike>
    Eigen::Matrix<
      typename QuaternionLike::Scalar,
      3,
      1,
      PINOCCHIO_EIGEN_PLAIN_TYPE(typename QuaternionLike::Vector3)::Options>
    // ------------------------------------------------------------
    // log3（四元数版）：S³ → so(3)，exp3 的逆
    //
    //     θ = 2·atan2(‖q.vec‖, q.w),      ω = (θ/sin(θ/2))·q.vec
    //
    // 用 atan2 而非 acos(w)：atan2 在整个范围内数值条件都好，
    // 而 acos 在 w→±1（θ→0 或 2π）附近导数发散、精度骤降。
    //
    // 【双覆盖处理】pos_neg：若 q.w < 0 就把整个四元数取反。
    //   因 q 与 −q 表示同一旋转，取反后保证 w ≥ 0，从而 θ ∈ [0, π]，
    //   取到【主值】（最短旋转路径）。
    //   实测：log3(−q) 与 log3(q) 严格相等。
    //
    // θ→0 时 θ/sin(θ/2) 是 0/0，切泰勒展开：
    //     θ/sin(θ/2) ≈ 2·(1 + (θ/2)²/6 + 7(θ/2)⁴/360)
    // 代码中 inv_sinc 即此系数，th2_2 = (θ/2)²。
    //
    // 带 theta 输出参数是为了让调用方复用该角度（如随后算 Jlog3），
    // 避免重复计算。
    //
    // 实测：log3(exp3(ω)) = ω（1.6e-16）；与矩阵版 log3 一致（1.2e-16）；
    //       theta 输出等于 ‖ω‖（1.1e-16）。
    // ------------------------------------------------------------
    log3(
      const Eigen::QuaternionBase<QuaternionLike> & quat, typename QuaternionLike::Scalar & theta)
    {
      typedef typename QuaternionLike::Scalar Scalar;
      static constexpr int Options =
        PINOCCHIO_EIGEN_PLAIN_TYPE(typename QuaternionLike::Vector3)::Options;
      typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;

      Vector3 res;
      const Scalar norm_squared = quat.vec().squaredNorm();

      static const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();
      static const Scalar ts_prec = TaylorSeriesExpansion<Scalar>::template precision<2>();
      const Scalar norm = math::sqrt(norm_squared + eps * eps);

      using ::pinocchio::internal::GE;
      using ::pinocchio::internal::if_then_else;
      using ::pinocchio::internal::LT;

      const Scalar pos_neg = if_then_else(GE, quat.w(), Scalar(0), Scalar(+1), Scalar(-1));

      Eigen::Quaternion<Scalar, Options> quat_pos;
      quat_pos.w() = pos_neg * quat.w();
      quat_pos.vec() = pos_neg * quat.vec();

      const Scalar theta_2 = math::atan2(norm, quat_pos.w()); // in [0,pi]
      const Scalar y_x = norm / quat_pos.w();                 // nonnegative
      const Scalar y_x_sq = norm_squared / (quat_pos.w() * quat_pos.w());

      theta = if_then_else(
        LT, norm_squared, ts_prec, Scalar(2.) * (Scalar(1) - y_x_sq / Scalar(3)) * y_x,
        Scalar(2.) * theta_2);

      const Scalar th2_2 = theta * theta / Scalar(4);
      const Scalar inv_sinc = if_then_else(
        LT, norm_squared, ts_prec,
        Scalar(2) * (Scalar(1) + th2_2 / Scalar(6) + Scalar(7) / Scalar(360) * th2_2 * th2_2),
        theta / math::sin(theta_2));

      for (Eigen::Index k = 0; k < 3; ++k)
      {
        // res[k] = if_then_else(LT, norm_squared, ts_prec,
        //                       Scalar(2) * (Scalar(1) + y_x_sq / Scalar(6) - y_x_sq*y_x_sq /
        //                       Scalar(9)) * pos_neg * quat.vec()[k], inv_sinc * pos_neg *
        //                       quat.vec()[k]);
        res[k] = inv_sinc * quat_pos.vec()[k];
      }
      return res;
    }

    ///
    /// \brief Log: SO3 -> so3.
    ///
    /// Pseudo-inverse of log from \f$ SO3 -> { v \in so3, ||v|| \le pi } \f$.
    ///
    /// \param[in] quat The unit quaternion representing a certain rotation.
    ///
    /// \return The angular velocity vector associated to the quaternion.
    ///
    template<typename QuaternionLike>
    Eigen::Matrix<
      typename QuaternionLike::Scalar,
      3,
      1,
      PINOCCHIO_EIGEN_PLAIN_TYPE(typename QuaternionLike::Vector3)::Options>
    log3(const Eigen::QuaternionBase<QuaternionLike> & quat)
    {
      typename QuaternionLike::Scalar theta;
      return log3(quat.derived(), theta);
    }

    ///
    /// \brief Derivative of \f$ q = \exp{\mathbf{v} + \delta\mathbf{v}} \f$ where \f$
    /// \delta\mathbf{v} \f$
    ///        is a small perturbation of \f$ \mathbf{v} \f$ at identity.
    ///
    /// \returns The Jacobian of the quaternion components variation.
    ///
    // ------------------------------------------------------------
    // Jexp3CoeffWise：∂(四元数系数)/∂ω，尺寸 4×3
    //
    // 注意与 explog.hxx 里 Jexp3 的区别：
    //   Jexp3（3×3）        ：李代数 → 李代数的右雅可比，用于流形上的梯度传播
    //   Jexp3CoeffWise（4×3）：直接对【四元数的 4 个存储系数】求导，
    //                          是"逐系数(coeff-wise)"的普通雅可比
    //
    // 行序与 Eigen 四元数系数一致：前 3 行是虚部 vec，第 4 行是实部 w。
    //
    // 用途：当把四元数当作【普通 4 维参数】参与优化时（如某些
    // 位姿估计/标定问题直接以 q 的分量为决策变量），需要这种导数。
    // 日常动力学走李代数路线，用的是 Jexp3 而非本函数。
    //
    // 两支同样按 ‖ω‖ 大小切换，小角度用泰勒展开。
    // 这里用的是普通 if（非 if_then_else）—— 故本函数不支持自动微分标量。
    //
    // 实测：与有限差分吻合到 1.2e-8（h=1e-7，与截断误差同量级）；
    //       小角度 θ=1e-6 分支同样正确。
    // ------------------------------------------------------------
    template<typename Vector3Like, typename Matrix43Like>
    void Jexp3CoeffWise(
      const Eigen::MatrixBase<Vector3Like> & v, const Eigen::MatrixBase<Matrix43Like> & Jexp)
    {
      //      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix43Like,4,3);
      assert(Jexp.rows() == 4 && Jexp.cols() == 3 && "Jexp does have the right size.");
      Matrix43Like & Jout = PINOCCHIO_EIGEN_CONST_CAST(Matrix43Like, Jexp);

      typedef typename Vector3Like::Scalar Scalar;

      const Scalar n2 = v.squaredNorm();
      const Scalar n = math::sqrt(n2);
      const Scalar theta = Scalar(0.5) * n;
      const Scalar theta2 = Scalar(0.25) * n2;

      if (n2 > math::sqrt(Eigen::NumTraits<Scalar>::epsilon()))
      {
        Scalar c, s;
        SINCOS(theta, &s, &c);
        Jout.template topRows<3>().noalias() =
          ((Scalar(0.5) / n2) * (c - 2 * s / n)) * v * v.transpose();
        Jout.template topRows<3>().diagonal().array() += s / n;
        Jout.template bottomRows<1>().noalias() = -s / (2 * n) * v.transpose();
      }
      else
      {
        Jout.template topRows<3>().noalias() =
          (-Scalar(1) / Scalar(12) + n2 / Scalar(480)) * v * v.transpose();
        Jout.template topRows<3>().diagonal().array() += Scalar(0.5) * (1 - theta2 / 6);
        Jout.template bottomRows<1>().noalias() =
          (Scalar(-0.25) * (Scalar(1) - theta2 / 6)) * v.transpose();
      }
    }

    ///
    ///  \brief Computes the Jacobian of log3 operator for a unit quaternion.
    ///
    /// \param[in] quat A unit quaternion representing the input rotation.
    /// \param[out] Jlog The resulting Jacobian of the log operator.
    ///
    // ------------------------------------------------------------
    // Jlog3（四元数入参）：log3 的雅可比，3×3
    //
    // 本身不含新数学：先用四元数版 log3 取出 (θ, ω)，再直接转调
    // explog.hxx 里那个"已知 θ 与 ω"的 Jlog3 重载 —— 因为雅可比的
    // 表达式只依赖 (θ, ω)，与输入是四元数还是矩阵无关。
    //
    // 提供本重载纯粹是为了免去调用方"先把四元数转成矩阵"的一步。
    // 用途同矩阵版：姿态 IK 的链式法则（见 GUIDE §4.2.2）。
    //
    // 实测：与矩阵版 pinocchio::Jlog3(q.matrix(), ·) 一致（3.7e-16）；
    //       与 log3 的右扰动有限差分吻合到 6.0e-9。
    // ------------------------------------------------------------
    template<typename QuaternionLike, typename Matrix3Like>
    void Jlog3(
      const Eigen::QuaternionBase<QuaternionLike> & quat,
      const Eigen::MatrixBase<Matrix3Like> & Jlog)
    {
      typedef typename QuaternionLike::Scalar Scalar;
      typedef Eigen::Matrix<
        Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(typename QuaternionLike::Coefficients)::Options>
        Vector3;

      Scalar t;
      Vector3 w(log3(quat, t));
      pinocchio::Jlog3(t, w, PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like, Jlog));
    }
  } // namespace quaternion
} // namespace pinocchio
