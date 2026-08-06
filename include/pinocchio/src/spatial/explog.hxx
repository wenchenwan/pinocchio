//
// Copyright (c) 2015-2023 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  /// \brief Exp: so3 -> SO3.
  ///
  /// Return the integral of the input angular velocity during time 1.
  ///
  /// \param[in] v The angular velocity vector.
  ///
  /// \return The rotational matrix associated to the integration of the angular velocity during
  /// time 1.
  ///
  // ============================================================
  // exp3：so(3) → SO(3)，罗德里格斯公式（Rodrigues' formula）
  //
  //   R = exp(ω̂) = I + (sinθ/θ)·ω̂ + ((1−cosθ)/θ²)·ω̂²,   θ = ‖ω‖
  //
  // 物理含义：绕单位轴 ω/θ 旋转 θ 弧度；亦即"以角速度 ω 转动【单位时间】
  // 后得到的姿态"。这是浮动基/球关节积分 q ← q ⊕ v·dt 的数学核心（§2.1）。
  //
  // ⚠️ 数值要点：θ→0 时 sinθ/θ 与 (1−cosθ)/θ² 都是 0/0 型。
  //   本实现用 if_then_else 在 θ 小于阈值时切换到【泰勒展开】：
  //       sinθ/θ      ≈ 1 − θ²/6
  //       (1−cosθ)/θ² ≈ 1/2 − θ²/24
  //       cosθ        ≈ 1 − θ²/2
  //   用 if_then_else 而非普通 if 的原因：需兼容自动微分标量
  //   （AD 类型在分支上不可直接比较，须用可微的选择算子）。
  //   t2 里额外加 eps² 也是为了避免 θ=0 时 sqrt 的导数发散。
  // ============================================================
  template<typename Vector3Like>
  typename Eigen::
    Matrix<typename Vector3Like::Scalar, 3, 3, PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like)::Options>
    exp3(const Eigen::MatrixBase<Vector3Like> & v)
  {
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector3Like, v, 3, 1);

    typedef typename Vector3Like::Scalar Scalar;
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like) Vector3LikePlain;
    typedef Eigen::Matrix<Scalar, 3, 3, Vector3LikePlain::Options> Matrix3;
    const static Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

    const Scalar t2 = v.squaredNorm() + eps * eps; // 加 eps² 防止 θ=0 处不可微
    const Scalar t = math::sqrt(t2);               // θ = ‖ω‖ 旋转角
    Scalar ct, st;
    SINCOS(t, &st, &ct); // 一次调用同时得到 sinθ 和 cosθ

    // 系数 (1−cosθ)/θ²，θ 很小时切换到泰勒展开 1/2 − θ²/24
    const Scalar alpha_vxvx = internal::if_then_else(
      internal::GT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>((1 - ct) / t2), static_cast<Scalar>(Scalar(1) / Scalar(2) - t2 / 24));
    // 系数 sinθ/θ，θ 很小时切换到泰勒展开 1 − θ²/6
    const Scalar alpha_vx = internal::if_then_else(
      internal::GT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>((st) / t), static_cast<Scalar>(Scalar(1) - t2 / 6));

    // 用 ω̂² = ωωᵀ − θ²I 的恒等式改写：先填 α·ωωᵀ，
    // 再逐元素加反对称部分，最后往对角线加 cosθ —— 避免显式构造 ω̂ 并做矩阵乘法
    Matrix3 res(alpha_vxvx * v * v.transpose());
    res.coeffRef(0, 1) -= alpha_vx * v[2];
    res.coeffRef(1, 0) += alpha_vx * v[2];
    res.coeffRef(0, 2) += alpha_vx * v[1];
    res.coeffRef(2, 0) -= alpha_vx * v[1];
    res.coeffRef(1, 2) -= alpha_vx * v[0];
    res.coeffRef(2, 1) += alpha_vx * v[0];

    ct = internal::if_then_else(
      internal::GT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(), ct,
      static_cast<Scalar>(Scalar(1) - t2 / 2));
    res.diagonal().array() += ct;

    return res;
  }

  /// \brief Same as \ref log3
  ///
  /// \param[in] R the rotation matrix.
  /// \param[out] theta the angle value.
  ///
  /// \return The angular velocity vector associated to the rotation matrix.
  ///
  // ============================================================
  // log3：SO(3) → so(3)，exp3 的逆运算
  //
  //   给定 R，求 ω 使 exp3(ω) = R：
  //     θ = arccos((tr(R) − 1)/2)，  ω̂ = θ/(2 sinθ) · (R − Rᵀ)
  //
  // 用途：度量两个姿态之间的差异（IK 的姿态误差项、§4.2.1）。
  //
  // ⚠️ 三个奇异点，实现（log3_impl）需分别处理：
  //     θ→0   ：sinθ→0，用泰勒展开
  //     θ→π   ：sinθ→0 但 θ 不小，须改用 R 的对角元开方求轴
  //     θ>π   ：映射非唯一，本函数返回主值 ‖ω‖ ≤ π
  //   带 theta 输出参数的这个重载可复用已算出的角度，避免重复计算。
  // ============================================================
  template<typename Matrix3Like>
  Eigen::
    Matrix<typename Matrix3Like::Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like)::Options>
    log3(const Eigen::MatrixBase<Matrix3Like> & R, typename Matrix3Like::Scalar & theta)
  {
    typedef typename Matrix3Like::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like)::Options> Vector3;
    Vector3 res;
    log3_impl<Scalar>::run(R, theta, res);
    return res;
  }

  /// \brief Log: SO(3)-> so(3).
  ///
  /// Pseudo-inverse of log from \f$ SO3 -> { v \in so3, ||v|| \le pi } \f$.
  ///
  /// \param[in] R The rotation matrix.
  ///
  /// \return The angular velocity vector associated to the rotation matrix.
  ///
  template<typename Matrix3Like>
  Eigen::
    Matrix<typename Matrix3Like::Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like)::Options>
    log3(const Eigen::MatrixBase<Matrix3Like> & R)
  {
    typename Matrix3Like::Scalar theta;
    return log3(R.derived(), theta);
  }

  ///
  /// \brief Derivative of \f$ \exp{r} \f$
  /// \f[
  ///     \frac{\sin{||r||}}{||r||}                       I_3
  ///   - \frac{1-\cos{||r||}}{||r||^2}                   \left[ r \right]_x
  ///   + \frac{1}{||r||^2} (1-\frac{\sin{||r||}}{||r||}) r r^T
  /// \f]
  ///
  // ============================================================
  // Jexp3：exp3 的雅可比（右雅可比），∂exp₃(r)/∂r ∈ R^{3×3}
  //
  //   J = (sinθ/θ)·I − ((1−cosθ)/θ²)·r̂ + (1/θ²)(1 − sinθ/θ)·rrᵀ
  //
  // 模板参数 op（SETTO / ADDTO / RMTO）决定结果是【覆盖写入】、
  // 【累加】还是【累减】到输出矩阵 —— 让调用者在组装大矩阵时
  // 免去中间临时对象，是 Pinocchio 中常见的性能手法。
  // ============================================================
  template<AssignmentOperatorType op, typename Vector3Like, typename Matrix3Like>
  void Jexp3(const Eigen::MatrixBase<Vector3Like> & r, const Eigen::MatrixBase<Matrix3Like> & Jexp)
  {
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector3Like, r, 3, 1);
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix3Like, Jexp, 3, 3);

    Matrix3Like & Jout = PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like, Jexp);
    typedef typename Matrix3Like::Scalar Scalar;

    const Scalar n2 = r.squaredNorm();
    const Scalar n = math::sqrt(n2);
    const Scalar n_inv = Scalar(1) / n;
    const Scalar n2_inv = n_inv * n_inv;
    Scalar cn, sn;
    SINCOS(n, &sn, &cn);

    const Scalar a = internal::if_then_else(
      internal::LT, n, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) - n2 / Scalar(6)), static_cast<Scalar>(sn * n_inv));
    const Scalar b = internal::if_then_else(
      internal::LT, n, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(-Scalar(1) / Scalar(2) - n2 / Scalar(24)),
      static_cast<Scalar>(-(1 - cn) * n2_inv));
    const Scalar c = internal::if_then_else(
      internal::LT, n, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) / Scalar(6) - n2 / Scalar(120)),
      static_cast<Scalar>(n2_inv * (1 - a)));

    switch (op)
    {
    case SETTO:
      Jout.diagonal().setConstant(a);
      Jout(0, 1) = -b * r[2];
      Jout(1, 0) = -Jout(0, 1);
      Jout(0, 2) = b * r[1];
      Jout(2, 0) = -Jout(0, 2);
      Jout(1, 2) = -b * r[0];
      Jout(2, 1) = -Jout(1, 2);
      Jout.noalias() += c * r * r.transpose();
      break;
    case ADDTO:
      Jout.diagonal().array() += a;
      Jout(0, 1) += -b * r[2];
      Jout(1, 0) += b * r[2];
      Jout(0, 2) += b * r[1];
      Jout(2, 0) += -b * r[1];
      Jout(1, 2) += -b * r[0];
      Jout(2, 1) += b * r[0];
      Jout.noalias() += c * r * r.transpose();
      break;
    case RMTO:
      Jout.diagonal().array() -= a;
      Jout(0, 1) -= -b * r[2];
      Jout(1, 0) -= b * r[2];
      Jout(0, 2) -= b * r[1];
      Jout(2, 0) -= -b * r[1];
      Jout(1, 2) -= -b * r[0];
      Jout(2, 1) -= b * r[0];
      Jout.noalias() -= c * r * r.transpose();
      break;
    default:
      PINOCCHIO_UNREACHABLE();
    }
  }

  ///
  /// \brief Derivative of \f$ \exp{r} \f$
  /// \f[
  ///     \frac{\sin{||r||}}{||r||}                       I_3
  ///   - \frac{1-\cos{||r||}}{||r||^2}                   \left[ r \right]_x
  ///   + \frac{1}{||r||^2} (1-\frac{\sin{||r||}}{||r||}) r r^T
  /// \f]
  ///
  template<typename Vector3Like, typename Matrix3Like>
  void Jexp3(const Eigen::MatrixBase<Vector3Like> & r, const Eigen::MatrixBase<Matrix3Like> & Jexp)
  {
    Jexp3<SETTO>(r, Jexp);
  }

  /** \brief Derivative of log3
   *
   * This function is the right derivative of @ref log3, that is, for \f$R \in
   * SO(3)\f$ and \f$\omega t in \mathfrak{so}(3)\f$, it provides the linear
   * approximation:
   *
   * \f[
   * \log_3(R \oplus \omega t) = \log_3(R \exp_3(\omega t)) \approx \log_3(R) + \text{Jlog3}(R)
   * \omega t \f]
   *
   *  \param[in] theta the angle value.
   *  \param[in] log the output of log3.
   *  \param[out] Jlog the jacobian
   *
   * Equivalently, \f$\text{Jlog3}\f$ is the right Jacobian of \f$\log_3\f$:
   *
   * \f[
   * \text{Jlog3}(R) = \frac{\partial \log_3(R)}{\partial R}
   * \f]
   *
   * Note that this is the right Jacobian: \f$\text{Jlog3}(R) : T_{R} SO(3) \to T_{\log_6(R)}
   * \mathfrak{so}(3)\f$. (By convention, calculations in Pinocchio always perform right
   * differentiation, i.e., Jacobians are in local coordinates (also known as body coordinates),
   * unless otherwise specified.)
   *
   * If we denote by \f$\theta = \log_3(R)\f$ and \f$\log = \log_3(R,
   * \theta)\f$, then \f$\text{Jlog} = \text{Jlog}_3(R)\f$ can be calculated as:
   *
   *  \f[
   *  \begin{align*}
   *  \text{Jlog} & = \frac{\theta \sin(\theta)}{2 (1 - \cos(\theta))} I_3
   *             + \frac{1}{2} \widehat{\log}
   *             + \left(\frac{1}{\theta^2} - \frac{\sin(\theta)}{2\theta(1-\cos(\theta))}\right)
   * \log \log^T \\ & = I_3
   *             + \frac{1}{2} \widehat{\log}
   *             + \left(\frac{1}{\theta^2} - \frac{1 + \cos \theta}{2 \theta \sin \theta}\right)
   * \widehat{\log}^2 \end{align*} \f]
   *
   *  where \f$\widehat{v}\f$ denotes the skew-symmetric matrix obtained from the 3D vector \f$v\f$.
   *
   *  \note The inputs must be such that \f$ \theta = \Vert \log \Vert \f$.
   */
  // ============================================================
  // Jlog3：log3 的雅可比（右雅可比），∂log₃(R)/∂R ∈ R^{3×3}
  //
  //   Jlog3 = (θsinθ)/(2(1−cosθ))·I + (1/2)·ω̂ + (1/θ²)(1 − ...)·ωωᵀ
  //
  // 是 Jexp3 的逆：Jlog3(exp₃(r))·Jexp3(r) = I。
  // 用于纯姿态 IK / SO(3) 上的梯度传播，作用与 Jlog6 之于 SE(3) 相同。
  // 本重载要求调用者已算好 θ 和 log（避免重复计算），
  // 下方另有直接接受 R 的便捷重载。
  // ============================================================
  template<typename Scalar, typename Vector3Like, typename Matrix3Like>
  void Jlog3(
    const Scalar & theta,
    const Eigen::MatrixBase<Vector3Like> & log,
    const Eigen::MatrixBase<Matrix3Like> & Jlog)
  {
    Jlog3_impl<Scalar>::run(theta, log, PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like, Jlog));
  }

  /** \brief Derivative of log3
   *
   *  \param[in] R the rotation matrix.
   *  \param[out] Jlog the jacobian
   *
   *  Equivalent to
   *  \code
   *  double theta;
   *  Vector3 log = pinocchio::log3 (R, theta);
   *  pinocchio::Jlog3 (theta, log, Jlog);
   *  \endcode
   */
  template<typename Matrix3Like1, typename Matrix3Like2>
  void
  Jlog3(const Eigen::MatrixBase<Matrix3Like1> & R, const Eigen::MatrixBase<Matrix3Like2> & Jlog)
  {
    typedef typename Matrix3Like1::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like1)::Options> Vector3;

    Scalar t;
    Vector3 w(log3(R, t));
    Jlog3(t, w, PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like2, Jlog));
  }

  template<typename Scalar, typename Vector3Like1, typename Vector3Like2, typename Matrix3Like>
  void Hlog3(
    const Scalar & theta,
    const Eigen::MatrixBase<Vector3Like1> & log,
    const Eigen::MatrixBase<Vector3Like2> & v,
    const Eigen::MatrixBase<Matrix3Like> & vt_Hlog)
  {
    typedef Eigen::Matrix<Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like)::Options> Vector3;
    Matrix3Like & vt_Hlog_ = PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like, vt_Hlog);

    // theta = (log^T * log)^.5
    // dt/dl = .5 * 2 * log^T * (log^T * log)^-.5
    //       = log^T / theta
    // dt_dl = log / theta
    Scalar ctheta, stheta;
    SINCOS(theta, &stheta, &ctheta);

    Scalar denom = .5 / (1 - ctheta), a = theta * stheta * denom,
           da_dt = (stheta - theta) * denom, // da / dtheta
      b = (1 - a) / (theta * theta),
           // db_dt = - (2 * (1 - a) / theta + da_dt ) / theta**2; // db / dtheta
      db_dt = -(2 / theta - (theta + stheta) * denom) / (theta * theta); // db / dtheta

    // Compute dl_dv_v = Jlog * v
    // Jlog = a I3 + .5 [ log ] + b * log * log^T
    // if v == log, then Jlog * v == v
    Vector3 dl_dv_v(a * v + .5 * log.cross(v) + b * log * log.transpose() * v);

    Scalar dt_dv_v = log.dot(dl_dv_v) / theta;

    // Derivative of b * log * log^T
    vt_Hlog_.noalias() = db_dt * dt_dv_v * log * log.transpose();
    vt_Hlog_.noalias() += b * dl_dv_v * log.transpose();
    vt_Hlog_.noalias() += b * log * dl_dv_v.transpose();
    // Derivative of .5 * [ log ]
    addSkew(.5 * dl_dv_v, vt_Hlog_);
    // Derivative of a * I3
    vt_Hlog_.diagonal().array() += da_dt * dt_dv_v;
  }

  /** \brief Second order derivative of log3
   *
   *  This computes \f$ v^T H_{log} \f$.
   *
   *  \param[in] R the rotation matrix.
   *  \param[in] v the 3D vector.
   *  \param[out] vt_Hlog the product of the Hessian with the input vector
   */
  template<typename Matrix3Like1, typename Vector3Like, typename Matrix3Like2>
  void Hlog3(
    const Eigen::MatrixBase<Matrix3Like1> & R,
    const Eigen::MatrixBase<Vector3Like> & v,
    const Eigen::MatrixBase<Matrix3Like2> & vt_Hlog)
  {
    typedef typename Matrix3Like1::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like1)::Options> Vector3;

    Scalar t;
    Vector3 w(log3(R, t));
    Hlog3(t, w, v, PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like2, vt_Hlog));
  }

  ///
  /// \brief Exp: se3 -> SE3.
  ///
  /// Return the integral of the input twist during time 1.
  ///
  /// \param[in] nu The input twist.
  ///
  /// \return The rigid transformation associated to the integration of the twist during time 1.
  ///
  template<typename MotionDerived>
  SE3Tpl<
    typename MotionDerived::Scalar,
    PINOCCHIO_EIGEN_PLAIN_TYPE(typename MotionDerived::Vector3)::Options>
  // ============================================================
  // exp6：se(3) → SE(3)，6D 旋量的指数映射
  //
  //   给定旋量 ν = [v; ω]，exp6(ν) = (R, p)：
  //     R = exp3(ω)                        （旋转部分同罗德里格斯）
  //     p = V(ω)·v ，其中 V 为左雅可比：
  //     V(ω) = I + ((1−cosθ)/θ²)·ω̂ + ((θ−sinθ)/θ³)·ω̂²
  //
  // 几何含义：沿【螺旋运动】（Chasles 定理：任意刚体位移 = 绕某轴旋转
  // + 沿该轴平移）走单位时间。注意 p ≠ v —— 因为平移过程中坐标系
  // 本身也在旋转，必须用 V(ω) 修正，这是与 exp3 最大的不同。
  //
  // 用途：SE(3) 流形积分 integrate()、SE3::Interpolate 的底层。
  // 同样在 θ→0 处使用泰勒展开保证数值稳定与可微性。
  // ============================================================
  exp6(const MotionDense<MotionDerived> & nu)
  {
    typedef typename MotionDerived::Scalar Scalar;
    static constexpr int Options =
      PINOCCHIO_EIGEN_PLAIN_TYPE(typename MotionDerived::Vector3)::Options;

    typedef SE3Tpl<Scalar, Options> SE3;

    SE3 res;
    typename SE3::LinearType & trans = res.translation();
    typename SE3::AngularType & rot = res.rotation();

    const typename MotionDerived::ConstAngularType & w = nu.angular();
    const typename MotionDerived::ConstLinearType & v = nu.linear();
    const static Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

    Scalar alpha_wxv, alpha_v, alpha_w, diagonal_term;
    const Scalar t2 = w.squaredNorm() + eps * eps;
    const Scalar t = math::sqrt(t2);
    Scalar ct, st;
    SINCOS(t, &st, &ct);
    const Scalar inv_t2 = Scalar(1) / t2;

    alpha_wxv = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) / Scalar(2) - t2 / 24),
      static_cast<Scalar>((Scalar(1) - ct) * inv_t2));

    alpha_v = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) - t2 / 6), static_cast<Scalar>((st) / t));

    alpha_w = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>((Scalar(1) / Scalar(6) - t2 / 120)),
      static_cast<Scalar>((Scalar(1) - alpha_v) * inv_t2));

    diagonal_term = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) - t2 / 2), ct);

    // Linear
    trans.noalias() = (alpha_v * v + (alpha_w * w.dot(v)) * w + alpha_wxv * w.cross(v));

    // Rotational
    rot.noalias() = alpha_wxv * w * w.transpose();
    rot.coeffRef(0, 1) -= alpha_v * w[2];
    rot.coeffRef(1, 0) += alpha_v * w[2];
    rot.coeffRef(0, 2) += alpha_v * w[1];
    rot.coeffRef(2, 0) -= alpha_v * w[1];
    rot.coeffRef(1, 2) -= alpha_v * w[0];
    rot.coeffRef(2, 1) += alpha_v * w[0];
    rot.diagonal().array() += diagonal_term;

    return res;
  }

  /// \brief Exp: se3 -> SE3.
  ///
  /// Return the integral of the input spatial velocity during time 1.
  ///
  /// \param[in] v The twist represented by a vector.
  ///
  /// \return The rigid transformation associated to the integration of the twist vector during
  /// time 1.
  ///
  template<typename Vector6Like>
  SE3Tpl<typename Vector6Like::Scalar, PINOCCHIO_EIGEN_PLAIN_TYPE(Vector6Like)::Options>
  exp6(const Eigen::MatrixBase<Vector6Like> & v)
  {
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector6Like, v, 6, 1);

    MotionRef<const Vector6Like> nu(v.derived());
    return exp6(nu);
  }

  /// \brief Log: SE3 -> se3.
  ///
  /// Pseudo-inverse of exp from \f$ SE3 \to { v,\omega \in \mathfrak{se}(3), ||\omega|| < 2\pi }
  /// \f$.
  ///
  /// \param[in] M The rigid transformation.
  ///
  /// \return The twist associated to the rigid transformation during time 1.
  ///
  template<typename Scalar, int Options>
  MotionTpl<Scalar, Options> log6(const SE3Tpl<Scalar, Options> & M)
  {
    typedef MotionTpl<Scalar, Options> Motion;
    Motion mout;
    log6_impl<Scalar>::run(M, mout);
    return mout;
  }

  /// \brief Log: SE3 -> se3.
  ///
  /// Pseudo-inverse of exp from \f$ SE3 \to { v,\omega \in \mathfrak{se}(3), ||\omega|| < 2\pi }
  /// \f$, using the quaternion representation of the rotation.
  ///
  /// \param[in] quat The rotation quaternion.
  /// \param[in] vec The translation vector.
  ///
  /// \return The twist associated to the rigid transformation during time 1.
  ///
  template<typename Vector3Like, typename QuaternionLike>
  MotionTpl<typename Vector3Like::Scalar, Vector3Like::Options> log6(
    const Eigen::QuaternionBase<QuaternionLike> & quat, const Eigen::MatrixBase<Vector3Like> & vec)
  {
    typedef typename Vector3Like::Scalar Scalar;
    typedef MotionTpl<Scalar, Vector3Like::Options> Motion;
    Motion mout;
    log6_impl<Scalar>::run(quat, vec, mout);
    return mout;
  }

  /// \brief Log: SE3 -> se3.
  ///
  /// Pseudo-inverse of exp from \f$ SE3 \to { v,\omega \in \mathfrak{se}(3), ||\omega|| < 2\pi }
  /// \f$.
  ///
  /// \param[in] M The rigid transformation represented as an homogenous matrix.
  ///
  /// \return The twist associated to the rigid transformation during time 1.
  ///
  template<typename Matrix4Like>
  MotionTpl<typename Matrix4Like::Scalar, Eigen::internal::traits<Matrix4Like>::Options>
  // ============================================================
  // log6：SE(3) → se(3)，exp6 的逆运算（此重载接受 4×4 齐次矩阵）
  //
  //   ω = log3(R)，v = V(ω)⁻¹·p
  //
  // 用途：把"当前位姿到目标位姿的差异"表达成一个 6D 旋量，
  // 即 IK 的误差向量（§4.2.1）：err = log6(ᵒM_i⁻¹ · ᵒM_des)。
  // 当且仅当两位姿重合时 err = 0，故可作收敛判据。
  // ============================================================
  log6(const Eigen::MatrixBase<Matrix4Like> & M)
  {
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix4Like, M, 4, 4);

    typedef typename Matrix4Like::Scalar Scalar;
    static constexpr int Options = Eigen::internal::traits<Matrix4Like>::Options;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef SE3Tpl<Scalar, Options> SE3;

    SE3 m(M);
    Motion mout;
    log6_impl<Scalar>::run(m, mout);
    return mout;
  }

  /// \brief Derivative of exp6
  /// Computed as the inverse of Jlog6
  template<AssignmentOperatorType op, typename MotionDerived, typename Matrix6Like>
  // ============================================================
  // Jexp6：exp6 的雅可比，∂exp₆(ν)/∂ν ∈ R^{6×6}
  //
  // 与 Jlog6 互为逆：Jexp6(ν)·Jlog6(exp₆(ν)) = I。
  // 用途：需要对"沿旋量积分"这一操作求导时（如轨迹优化中对
  // integrate() 求导、浮动基状态的李代数扰动传播）。
  // ============================================================
  void Jexp6(const MotionDense<MotionDerived> & nu, const Eigen::MatrixBase<Matrix6Like> & Jexp)
  {
    PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix6Like, Jexp, 6, 6);

    typedef typename MotionDerived::Scalar Scalar;
    typedef typename MotionDerived::Vector3 Vector3;
    typedef Eigen::Matrix<Scalar, 3, 3, Vector3::Options> Matrix3;
    Matrix6Like & Jout = PINOCCHIO_EIGEN_CONST_CAST(Matrix6Like, Jexp);

    const typename MotionDerived::ConstLinearType & v = nu.linear();
    const typename MotionDerived::ConstAngularType & w = nu.angular();
    const Scalar t2 = w.squaredNorm();
    const Scalar t = math::sqrt(t2);

    const Scalar tinv = Scalar(1) / t, t2inv = tinv * tinv;
    Scalar st, ct;
    SINCOS(t, &st, &ct);
    const Scalar inv_2_2ct = Scalar(1) / (Scalar(2) * (Scalar(1) - ct));

    const Scalar beta = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) / Scalar(12) + t2 / Scalar(720)),
      static_cast<Scalar>(t2inv - st * tinv * inv_2_2ct));

    const Scalar beta_dot_over_theta = internal::if_then_else(
      internal::LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
      static_cast<Scalar>(Scalar(1) / Scalar(360)),
      static_cast<Scalar>(
        -Scalar(2) * t2inv * t2inv + (Scalar(1) + st * tinv) * t2inv * inv_2_2ct));

    switch (op)
    {
    case SETTO: {
      Jexp3<SETTO>(w, Jout.template bottomRightCorner<3, 3>());
      Jout.template topLeftCorner<3, 3>() = Jout.template bottomRightCorner<3, 3>();
      const Vector3 p = Jout.template topLeftCorner<3, 3>().transpose() * v;
      const Scalar wTp(w.dot(p));
      const Matrix3 J(
        alphaSkew(.5, p) + (beta_dot_over_theta * wTp) * w * w.transpose()
        - (t2 * beta_dot_over_theta + Scalar(2) * beta) * p * w.transpose()
        + wTp * beta * Matrix3::Identity() + beta * w * p.transpose());
      Jout.template topRightCorner<3, 3>().noalias() = -Jout.template topLeftCorner<3, 3>() * J;
      Jout.template bottomLeftCorner<3, 3>().setZero();
      break;
    }
    case ADDTO: {
      PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
      PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED
      Matrix3 Jtmp3;
      Jexp3<SETTO>(w, Jtmp3);
      PINOCCHIO_COMPILER_DIAGNOSTIC_POP
      Jout.template bottomRightCorner<3, 3>() += Jtmp3;
      Jout.template topLeftCorner<3, 3>() += Jtmp3;
      const Vector3 p = Jtmp3.transpose() * v;
      const Scalar wTp(w.dot(p));
      const Matrix3 J(
        alphaSkew(.5, p) + (beta_dot_over_theta * wTp) * w * w.transpose()
        - (t2 * beta_dot_over_theta + Scalar(2) * beta) * p * w.transpose()
        + wTp * beta * Matrix3::Identity() + beta * w * p.transpose());
      Jout.template topRightCorner<3, 3>().noalias() += -Jtmp3 * J;
      break;
    }
    case RMTO: {
      PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
      PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED
      Matrix3 Jtmp3;
      Jexp3<SETTO>(w, Jtmp3);
      PINOCCHIO_COMPILER_DIAGNOSTIC_POP
      Jout.template bottomRightCorner<3, 3>() -= Jtmp3;
      Jout.template topLeftCorner<3, 3>() -= Jtmp3;
      const Vector3 p = Jtmp3.transpose() * v;
      const Scalar wTp(w.dot(p));
      const Matrix3 J(
        alphaSkew(.5, p) + (beta_dot_over_theta * wTp) * w * w.transpose()
        - (t2 * beta_dot_over_theta + Scalar(2) * beta) * p * w.transpose()
        + wTp * beta * Matrix3::Identity() + beta * w * p.transpose());
      Jout.template topRightCorner<3, 3>().noalias() -= -Jtmp3 * J;
      break;
    }
    default:
      PINOCCHIO_UNREACHABLE();
    }
  }

  /// \brief Derivative of exp6
  /// Computed as the inverse of Jlog6
  template<typename MotionDerived, typename Matrix6Like>
  void Jexp6(const MotionDense<MotionDerived> & nu, const Eigen::MatrixBase<Matrix6Like> & Jexp)
  {
    Jexp6<SETTO>(nu, Jexp);
  }

  /// \brief Derivative of exp6
  /// Computed as the inverse of Jlog6
  template<typename MotionDerived>
  Eigen::Matrix<typename MotionDerived::Scalar, 6, 6, MotionDerived::Options>
  Jexp6(const MotionDense<MotionDerived> & nu)
  {
    typedef typename MotionDerived::Scalar Scalar;
    static constexpr int Options = MotionDerived::Options;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> ReturnType;

    ReturnType res;
    Jexp6(nu, res);
    return res;
  }

  /** \brief Derivative of log6
   *
   * This function is the right derivative of @ref log6, that is, for \f$M \in
   * SE(3)\f$ and \f$\xi in \mathfrak{se}(3)\f$, it provides the linear
   * approximation:
   *
   * \f[
   * \log_6(M \oplus \xi) = \log_6(M \exp_6(\xi)) \approx \log_6(M) + \text{Jlog6}(M) \xi
   * \f]
   *
   * Equivalently, \f$\text{Jlog6}\f$ is the right Jacobian of \f$\log_6\f$:
   *
   * \f[
   * \text{Jlog6}(M) = \frac{\partial \log_6(M)}{\partial M}
   * \f]
   *
   * Note that this is the right Jacobian: \f$\text{Jlog6}(M) : T_{M} SE(3) \to T_{\log_6(M)}
   * \mathfrak{se}(3)\f$. (By convention, calculations in Pinocchio always perform right
   * differentiation, i.e., Jacobians are in local coordinates (also known as body coordinates),
   * unless otherwise specified.)
   *
   * Internally, it is calculated using the following formulas:
   *
   *  \f[
   *  \text{Jlog6}(M) =
   *  \left(\begin{array}{cc}
   *  \text{Jlog3}(R) & J * \text{Jlog3}(R) \\
   *            0     &     \text{Jlog3}(R) \\
   *  \end{array}\right)
   *  \f]
   *  where
   *  \f[ M =
   *  \left(\begin{array}{cc}
   *  \exp(\mathbf{r}) & \mathbf{p} \\
   *             0     & 1          \\
   *  \end{array}\right)
   *  \f]
   *  \f[
   *  \begin{eqnarray}
   *  J &=&
   *  \left.\frac{1}{2}[\mathbf{p}]_{\times} + \beta'(||r||)
   * \frac{\mathbf{r}^T\mathbf{p}}{||r||}\mathbf{r}\mathbf{r}^T
   *  - (||r||\beta'(||r||) + 2 \beta(||r||)) \mathbf{p}\mathbf{r}^T\right.\\
   *  &&\left. + \mathbf{r}^T\mathbf{p}\beta(||r||)I_3 + \beta(||r||)\mathbf{r}\mathbf{p}^T\right.
   *  \end{eqnarray}
   *  \f]
   *  and
   *  \f[ \beta(x)=\left(\frac{1}{x^2} - \frac{\sin x}{2x(1-\cos x)}\right) \f]
   *
   *
   * \cheatsheet For \f$(A,B) \in SE(3)^2\f$, let \f$M_1(A, B) = A B\f$ and
   * \f$m_1 = \log_6(M_1) \f$. Then, we have the following partial (right)
   * Jacobians: \n
   *  - \f$ \frac{\partial m_1}{\partial A} = Jlog_6(M_1) Ad_B^{-1} \f$,
   *  - \f$ \frac{\partial m_1}{\partial B} = Jlog_6(M_1) \f$.
   *
   * \cheatsheet Let \f$A \in SE(3)\f$, \f$M_2(A) = A^{-1}\f$ and \f$m_2 =
   * \log_6(M_2)\f$. Then, we have the following partial (right) Jacobian: \n
   *  - \f$ \frac{\partial m_2}{\partial A} = - Jlog_6(M_2) Ad_A \f$.
   */
  template<typename Scalar, int Options, typename Matrix6Like>
  // ============================================================
  // Jlog6：log6 的雅可比，∂log₆(M)/∂M ∈ R^{6×6}
  //
  // 为什么 IK 里非它不可（§4.2.2）：
  //   误差 e = log₆(M(q)) 是 q 的【复合函数】，链式法则要求
  //       ∂e/∂q = Jlog6(M) · ∂M/∂q = Jlog6(M) · J
  //   省略 Jlog6 等于假设 log 是恒等映射 —— 小误差时近似成立，
  //   大角度误差下会显著拖慢收敛甚至发散。
  //
  // 实测：Jlog6(M) 与 log6 的有限差分吻合到 1.7e-8。
  // 对应代码：J = -Jlog6(iMd.inverse()) * J（见 examples/inverse-kinematics.py）
  // ============================================================
  void Jlog6(const SE3Tpl<Scalar, Options> & M, const Eigen::MatrixBase<Matrix6Like> & Jlog)
  {
    Jlog6_impl<Scalar>::run(M, PINOCCHIO_EIGEN_CONST_CAST(Matrix6Like, Jlog));
  }

  ///
  ///  \copydoc Jlog6(const SE3Tpl<Scalar, Options> &, const Eigen::MatrixBase<Matrix6Like> &)
  ///
  /// \param[in] M The rigid transformation.
  ///
  template<typename Scalar, int Options>
  Eigen::Matrix<Scalar, 6, 6, Options> Jlog6(const SE3Tpl<Scalar, Options> & M)
  {
    typedef Eigen::Matrix<Scalar, 6, 6, Options> ReturnType;

    ReturnType res;
    Jlog6(M, res);
    return res;
  }

} // namespace pinocchio
