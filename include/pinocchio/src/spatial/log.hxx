//
// Copyright (c) 2015-2021 CNRS INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/spatial.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/spatial.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  namespace internal
  {
    // ============================================================
    // compute_theta_axis<i0>：θ≈π 奇异情形下，从旋转矩阵 R 的【对角】反解转轴
    //
    // 背景：θ≈π 时 sinθ→0，一般公式 ω=θ/(2sinθ)·unSkew(R−Rᵀ) 退化（0/0）。
    //       但此时 R ≈ 2 a·aᵀ − I（a 为单位转轴），于是 R_ii = 2a_i²−1，
    //       可从对角元开方求出 |a_i|。传入的 val_i = 2R_ii − tr + 1 = 4a_i²(θ=π 时)。
    //
    // 本函数以第 i0 个分量为"主轴分量"（其平方最大、开方最稳）来重建整根轴：
    //   · axis[i0] = ±√val/2                （符号由 R(i2,i1) 与 R(i1,i2) 大小定，保证 s≠0）
    //   · axis[i1], axis[i2] 由对称部分 (R_ij+R_ji)/(2s) 得出
    //   · w = (R(i2,i1)−R(i1,i2))/(2s) ~ cos(θ/2) 相关量 → angle = 2·atan2(‖axis‖, w)
    // 调 3 次（i0=0,1,2）取 val 最大者，见 log3_impl。
    // ============================================================
    template<long i0, typename Matrix3, typename Vector3>
    void compute_theta_axis(
      const typename Matrix3::Scalar & val,
      const Eigen::MatrixBase<Matrix3> & R,
      typename Matrix3::Scalar & angle,
      const Eigen::MatrixBase<Vector3> & _axis)
    {
      typedef typename Matrix3::Scalar Scalar;
      static const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

      static const long i1 = (i0 + 1) % 3; // 另两个分量下标（循环）
      static const long i2 = (i0 + 2) % 3;
      Vector3 & axis = _axis.const_cast_derived();

      // s = ±√val ≈ ±2|a_i0|；符号取 sign(R(i2,i1)−R(i1,i2)) 以保证 s≠0 且轴方向一致
      const Scalar s =
        math::sqrt(val + eps + eps * eps)
        * if_then_else(
          GE, R.coeff(i2, i1), R.coeff(i1, i2), Scalar(1.),
          Scalar(-1.)); // Ensure value in sqrt is non negative and that s is non zero
      axis[i0] = s / Scalar(2);                                                     // 主轴分量 ≈ a_i0
      axis[i1] = Scalar(1) / (Scalar(2) * s) * (R.coeff(i1, i0) + R.coeff(i0, i1)); // 由对称部分求
      axis[i2] = Scalar(1) / (Scalar(2) * s) * (R.coeff(i2, i0) + R.coeff(i0, i2));
      const Scalar w = Scalar(1) / (Scalar(2) * s) * (R.coeff(i2, i1) - R.coeff(i1, i2)); // ~cos(θ/2) 相关

      const Scalar axis_norm = axis.norm();
      angle = Scalar(2) * math::atan2(axis_norm, w); // θ = 2·atan2(sin(θ/2), cos(θ/2))
      axis /= axis_norm;                             // 归一化为单位轴
    }
  } // namespace internal

  /// \brief Renormalize a rotation matrix.
  // ============================================================
  // renormalize_rotation_matrix：把可能略微漂移的 R 快速拉回 SO(3)
  //   用 Gram-Schmidt 式正交化（比 §3.5 的 SVD 投影廉价，log3 内部用它兜底）：
  //     col0 ← 归一化                       （定第一根轴）
  //     col1 ← 归一化                       （暂定第二根轴）
  //     col2 ← col0 × col1                  （叉乘得正交第三轴，保右手系）
  //     col0 ← col1 × col2                  （回填 col0，使三轴严格互相正交）
  //   注意：不是最优（Frobenius 最近）投影，但足够 log3 消除小数值误差。
  // ============================================================
  template<typename Matrix3>
  inline typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3)
    renormalize_rotation_matrix(const Eigen::MatrixBase<Matrix3> & R)
  {
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3) Rout;
    Rout.col(0).noalias() = R.col(0) / R.col(0).norm();  // e0 = 归一化 col0
    Rout.col(1).noalias() = R.col(1) / R.col(1).norm();  // e1 = 归一化 col1（暂）
    Rout.col(2).noalias() = Rout.col(0).cross(Rout.col(1)); // e2 = e0 × e1（正交、右手）
    Rout.col(0).noalias() = Rout.col(1).cross(Rout.col(2)); // 回填 e0 = e1 × e2，保证三轴严格正交
    return Rout;
  }

  /// \brief Generic evaluation of log3 function
  template<typename _Scalar>
  // ============================================================
  // log3_impl：log3 的实际计算，按标量类型特化
  //
  // 三种情形分别处理（这是 log3 全部复杂度的来源）：
  //   ① θ 接近 0  ：sinθ→0，改用泰勒展开，避免 0/0
  //   ② θ 接近 π  ：sinθ→0 但 θ 不小，(R−Rᵀ) 退化，
  //                 须改从 R 的【对角元】开方求转轴方向，再定符号
  //   ③ 一般情形  ：ω = θ/(2sinθ)·unSkew(R − Rᵀ)
  // ============================================================
  struct log3_impl
  {
    template<typename Matrix3Like, typename Vector3Out>
    static void run(
      const Eigen::MatrixBase<Matrix3Like> & R,
      typename Matrix3Like::Scalar & theta,
      const Eigen::MatrixBase<Vector3Out> & angle_axis)
    {
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix3Like, R, 3, 3);
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector3Out, angle_axis, 3, 1);
      using namespace internal;

      typedef typename Matrix3Like::Scalar Scalar;
      typedef Eigen::Matrix<Scalar, 3, 1, PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like)::Options> Vector3;
      static const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

      const static Scalar PI_value = PI<Scalar>();
      Vector3Out & angle_axis_ = angle_axis.const_cast_derived();

      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix3Like) Matrix3;
      const Matrix3 Rnormed = renormalize_rotation_matrix(R); // 先消除输入 R 的数值漂移

      // cosθ = (tr R − 1)/2  ——旋转角与迹的标准关系
      const Scalar tr = Rnormed.trace();
      const Scalar cos_value = (tr - Scalar(1)) / Scalar(2);

      const Scalar prec = TaylorSeriesExpansion<Scalar>::template precision<2>();
      // ---- 情形②预算：θ≈π 的奇异分支（无条件先算好，最后用 if_then_else 选） ----
      Vector3 angle_axis_singular;
      Scalar theta_singular;

      {
        // val_i = 2R_ii − tr + 1（θ=π 时 = 4a_i²），三选一取最大分量开方最稳
        Vector3 val_singular;
        val_singular.array() = Scalar(2) * Rnormed.diagonal().array() - tr + Scalar(1);
        Vector3 axis_0, axis_1, axis_2;
        Scalar theta_0, theta_1, theta_2;

        // 以每个对角分量为"主轴"各重建一次（见 compute_theta_axis）
        internal::compute_theta_axis<0>(val_singular[0], Rnormed, theta_0, axis_0);
        internal::compute_theta_axis<1>(val_singular[1], Rnormed, theta_1, axis_1);
        internal::compute_theta_axis<2>(val_singular[2], Rnormed, theta_2, axis_2);

        // 取 val 最大者对应的 (θ, axis)（嵌套 if_then_else = 三选一 argmax）
        theta_singular = if_then_else(
          GE, val_singular[0], val_singular[1],
          if_then_else(GE, val_singular[0], val_singular[2], theta_0, theta_2),
          if_then_else(GE, val_singular[1], val_singular[2], theta_1, theta_2));

        for (int k = 0; k < 3; ++k)
          angle_axis_singular[k] = if_then_else(
            GE, val_singular[0], val_singular[1],
            if_then_else(GE, val_singular[0], val_singular[2], axis_0[k], axis_2[k]),
            if_then_else(GE, val_singular[1], val_singular[2], axis_1[k], axis_2[k]));
      }
      // ---- 一般情形的角度 θ_nominal = acos(cosθ)，两端各切更稳的写法 ----
      // tr→3(θ→0) 或 tr→−1(θ→π) 时 acos 导数发散，改用 √(2(1−cosθ)) 展开
      const Scalar acos_expansion = math::sqrt(Scalar(2) * (Scalar(1) - cos_value) + eps * eps);
      const Scalar theta_nominal = if_then_else(
        LE, tr, static_cast<Scalar>(Scalar(3) - prec),
        if_then_else(
          GE, tr, static_cast<Scalar>(Scalar(-1) + prec),
          math::acos(cos_value),                         // 中间区：直接 acos
          static_cast<Scalar>(PI_value - acos_expansion) // 近 π：π − √(2(1−cos))
          ),
        static_cast<Scalar>(acos_expansion) // 近 0：√(2(1−cos)) ≈ θ
      );
      assert(
        check_expression_if_real<Scalar>(theta_nominal == theta_nominal)
        && "theta contains some NaN"); // theta != NaN

      // unSkew(R) 取 R 的反对称部分对应向量 = sinθ·axis（一般情形的原料）
      Vector3 antisymmetric_R;
      unSkew(Rnormed, antisymmetric_R);
      const Scalar norm_antisymmetric_R_squared = antisymmetric_R.squaredNorm();

      // 系数 t = θ/sinθ；θ→0 时切泰勒 1 + ‖·‖²/6 + 3‖·‖⁴/40（避免 0/0）
      const Scalar t = if_then_else(
        GE, theta_nominal, prec,
        static_cast<Scalar>(theta_nominal / sin(theta_nominal)), // 一般：θ/sinθ
        static_cast<Scalar>(
          Scalar(1.) + norm_antisymmetric_R_squared / Scalar(6)
          + norm_antisymmetric_R_squared * norm_antisymmetric_R_squared * Scalar(3)
              / Scalar(40)) // θ→0 泰勒
      );

      // ---- 最终选择：cosθ 接近 −1(θ≈π) → 用奇异分支；否则用一般公式 ----
      theta = if_then_else(
        GE, cos_value, static_cast<Scalar>(Scalar(-1.) + prec), theta_nominal, theta_singular);

      // 一般情形：ω = (θ/sinθ)·unSkew(R−Rᵀ)/... = t·antisymmetric_R；奇异：θ_π·axis
      for (int k = 0; k < 3; ++k)
        angle_axis_[k] = if_then_else(
          GE, cos_value, static_cast<Scalar>(Scalar(-1.) + prec),
          static_cast<Scalar>(t * antisymmetric_R[k]),               // 一般公式
          static_cast<Scalar>(theta_singular * angle_axis_singular[k])); // θ≈π 奇异
    }
  };

  /// \brief Generic evaluation of Jlog3 function
  template<typename _Scalar>
  // ---- Jlog3_impl：Jlog3 的实际计算 ----
  // 同样需在 θ→0 处切泰勒展开：θsinθ/(2(1−cosθ)) → 1 − θ²/12
  struct Jlog3_impl
  {
    template<typename Scalar, typename Vector3Like, typename Matrix3Like>
    static void run(
      const Scalar & theta,
      const Eigen::MatrixBase<Vector3Like> & log,
      const Eigen::MatrixBase<Matrix3Like> & Jlog)
    {
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector3Like, log, 3, 1);
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix3Like, Jlog, 3, 3);

      // Jlog3(R) = α·ω ωᵀ + diag·I + ½[ω]×  （log 即 ω=log3(R)，θ=‖ω‖）
      using namespace internal;
      Scalar ct, st;
      SINCOS(theta, &st, &ct);
      const Scalar st_1mct = st / (Scalar(1) - ct); // sinθ/(1−cosθ) = cot(θ/2)
      const Scalar prec = TaylorSeriesExpansion<Scalar>::template precision<3>();

      // α = 1/θ² − sinθ/(2θ(1−cosθ))；θ→0 切泰勒 1/12 + θ²/720（外积项系数）
      const Scalar alpha = if_then_else(
        LT, theta, prec,
        static_cast<Scalar>(Scalar(1) / Scalar(12) + theta * theta / Scalar(720)),       // θ→0 泰勒
        static_cast<Scalar>(Scalar(1) / (theta * theta) - st_1mct / (Scalar(2) * theta)) // 一般
      );

      // 对角项 = ½·θ·cot(θ/2)；θ→0 切泰勒 ½(2 − θ²/6)
      const Scalar diag_value = if_then_else(
        LT, theta, prec,
        static_cast<Scalar>(Scalar(0.5) * (Scalar(2) - theta * theta / Scalar(6))), // θ→0 泰勒
        static_cast<Scalar>(Scalar(0.5) * (theta * st_1mct))                        // 一般
      );

      Matrix3Like & Jlog_ = PINOCCHIO_EIGEN_CONST_CAST(Matrix3Like, Jlog);
      Jlog_.noalias() = alpha * log * log.transpose(); // α·ω ωᵀ（外积项）
      Jlog_.diagonal().array() += diag_value;          // + diag·I

      // Jlog += ½[ω]×（反对称项，addSkew 就地累加，见 skew.hxx）
      addSkew(Scalar(0.5) * log, Jlog_);
    }
  };

  /// \brief Generic evaluation of log6 function
  template<typename _Scalar>
  // ---- log6_impl：log6 的实际计算 ----
  // 先由 log3 得到 ω，再求平移部分 v = V(ω)⁻¹·p。
  // V⁻¹ 的系数在 θ→0 处同样奇异，需泰勒展开处理
  struct log6_impl
  {
    template<typename Scalar, int Options, typename MotionDerived>
    static void run(const SE3Tpl<Scalar, Options> & M, MotionDense<MotionDerived> & mout)
    {
      typedef SE3Tpl<Scalar, Options> SE3;
      typedef typename SE3::Vector3 Vector3;

      typename SE3::ConstAngularRef R = M.rotation();
      typename SE3::ConstLinearRef p = M.translation();

      using namespace internal;

      Vector3 antisymmetric_R;
      unSkew(R, antisymmetric_R);
      const Scalar norm_antisymmetric_R_squared = antisymmetric_R.squaredNorm();

      // ω = log3(R)（θ∈[0,π]）；线速度 v = V(ω)⁻¹·p，下面把 V⁻¹ 展开成 α,β 两系数
      Scalar theta;
      const Scalar tr = R.trace();
      const Vector3 w(log3(R, theta)); // t in [0,π]
      const Scalar & t2 = norm_antisymmetric_R_squared;

      Scalar st, ct;
      SINCOS(theta, &st, &ct);
      // α = θsinθ/(2(1−cosθ))（V⁻¹ 的单位阵项系数）；θ→0 切泰勒
      const Scalar alpha = if_then_else(
        GE, tr,
        static_cast<Scalar>(Scalar(3) - TaylorSeriesExpansion<Scalar>::template precision<2>()),
        static_cast<Scalar>(Scalar(1) - t2 / Scalar(12) - t2 * t2 / Scalar(720)), // θ→0 泰勒
        static_cast<Scalar>(theta * st / (Scalar(2) * (Scalar(1) - ct)))          // 一般
      );

      // β = 1/θ² − sinθ/(2θ(1−cosθ))（V⁻¹ 的 ω ωᵀ 项系数）；θ→0 切泰勒
      const Scalar beta = if_then_else(
        GE, tr,
        static_cast<Scalar>(Scalar(3) - TaylorSeriesExpansion<Scalar>::template precision<2>()),
        static_cast<Scalar>(Scalar(1) / Scalar(12) + t2 / Scalar(720)), // θ→0 泰勒
        static_cast<Scalar>(
          Scalar(1) / (theta * theta) - st / (Scalar(2) * theta * (Scalar(1) - ct))) // 一般
      );

      // v = V⁻¹p = α·p − ½ω×p + β(ωᵀp)ω   （V⁻¹ = αI − ½[ω]× + β ω ωᵀ）
      mout.linear().noalias() = alpha * p - Scalar(0.5) * w.cross(p) + (beta * w.dot(p)) * w;
      mout.angular() = w;
    }

    // ---- log6 的四元数重载：从 (四元数 quat, 平移 vec) 直接求 twist ----
    // 优点：避免先把四元数转成旋转矩阵；用半角 ct_2=cos(θ/2)、st_2=sin(θ/2) 更稳。
    // 数学同上：v = V⁻¹·vec，只是这里用四元数版 log3 求 ω、用半角 cot 表达 β。
    template<typename Vector3Like, typename QuaternionLike, typename MotionDerived>
    static void run(
      const Eigen::QuaternionBase<QuaternionLike> & quat,
      const Eigen::MatrixBase<Vector3Like> & vec,
      MotionDense<MotionDerived> & mout)
    {
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Vector3Like, vec, 3, 1);

      typedef typename Vector3Like::Scalar Scalar;
      static constexpr int Options = PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like)::Options;
      typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
      const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();

      using namespace internal;

      // 取 sign(quat.w)：四元数 q 与 −q 表示同一旋转，统一到 w≥0 半球避免歧义
      const Scalar pos_neg = if_then_else(GE, quat.w(), Scalar(0), Scalar(+1), Scalar(-1));

      Scalar theta;
      Vector3 w(quaternion::log3(quat, theta)); // ω，四元数版 log3（构造上无奇异）
      const Scalar t2 = w.squaredNorm();

      // Scalar st,ct; SINCOS(theta,&st,&ct);
      Scalar st_2, ct_2;
      ct_2 = pos_neg * quat.w();
      st_2 = math::sqrt(quat.vec().squaredNorm() + eps * eps);
      const Scalar cot_th_2 = ct_2 / st_2; // cot(θ/2) = cos(θ/2)/sin(θ/2)
      // const Scalar cot_th_2 = ( st / (Scalar(1) - ct) ); // cotan of half angle

      // we use formula (9.26) from
      // https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf for the linear
      // part of the Log map. A Taylor series expansion of cotan can be used up to order 4
      const Scalar th_2_squared = t2 / Scalar(4); // (theta / 2) squared

      //      const Scalar alpha = if_then_else(LE,theta,TaylorSeriesExpansion<Scalar>::template
      //      precision<3>(),
      //                                        static_cast<Scalar>(Scalar(1) - t2/Scalar(12) -
      //                                        t2*t2/Scalar(720)), // then
      //                                        static_cast<Scalar>(theta * cot_th_2 /(Scalar(2)))
      //                                        // else
      //                                        );

      // β = 1/θ² − cot(θ/2)/(2θ)（V⁻¹ 的二次项系数）；θ→0 切泰勒 beta_alt
      const Scalar beta_alt = (Scalar(1) / Scalar(3) - th_2_squared / Scalar(45)) / Scalar(4);
      const Scalar beta = if_then_else(
        LE, theta, TaylorSeriesExpansion<Scalar>::template precision<3>(),
        static_cast<Scalar>(beta_alt),                                       // θ→0 泰勒
        static_cast<Scalar>(Scalar(1) / t2 - cot_th_2 * Scalar(0.5) / theta) // 一般
        // static_cast<Scalar>(Scalar(1) / t2 - st/(Scalar(2)*theta*(Scalar(1)-ct))) // else
      );

      // v = vec − ½ω×vec + β·ω×(ω×vec)  （用二重叉乘等价 V⁻¹·vec，α 项已并入常数 1）
      mout.linear().noalias() = vec - Scalar(0.5) * w.cross(vec) + beta * w.cross(w.cross(vec));
      mout.angular() = w;
    }
  };

  // ============================================================
  // Jlog6_impl：log6 的 6×6 右雅可比 Jlog6(M)。IK 里把 SE(3) 误差反传到关节空间要乘它。
  //
  // 分块结构（旋转/平移各 3 维）：
  //   Jlog6 = [ A  B ；
  //             0  D ]，其中 A = D = Jlog3(R)（旋转块复用），B 是平移对旋转的耦合块，
  //   左下块为 0（角量不受平移影响）。B 由一个临时 C 乘 A 得到（C 借作暂存后清零）。
  // ============================================================
  template<typename _Scalar>
  struct Jlog6_impl
  {
    template<typename Scalar, int Options, typename Matrix6Like>
    static void run(const SE3Tpl<Scalar, Options> & M, const Eigen::MatrixBase<Matrix6Like> & Jlog)
    {
      PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix6Like, Jlog, 6, 6);

      typedef SE3Tpl<Scalar, Options> SE3;
      typedef typename SE3::Vector3 Vector3;
      Matrix6Like & value = PINOCCHIO_EIGEN_CONST_CAST(Matrix6Like, Jlog);

      typename SE3::ConstAngularRef R = M.rotation();
      typename SE3::ConstLinearRef p = M.translation();

      using namespace internal;

      Scalar t;
      Vector3 w(log3(R, t)); // 先取 ω=log3(R)、θ=t

      // value is decomposed as following:
      // value = [ A, B;
      //           C, D ]
      typedef Eigen::Block<Matrix6Like, 3, 3> Block33;
      Block33 A = value.template topLeftCorner<3, 3>();
      Block33 B = value.template topRightCorner<3, 3>();
      Block33 C = value.template bottomLeftCorner<3, 3>();
      Block33 D = value.template bottomRightCorner<3, 3>();

      Jlog3(t, w, A); // A = Jlog3(R)
      D = A;          // D 与 A 相同（旋转块复用）

      const Scalar t2 = t * t;
      const Scalar tinv = Scalar(1) / t, t2inv = tinv * tinv;

      Scalar st, ct;
      SINCOS(t, &st, &ct);
      const Scalar inv_2_2ct = Scalar(1) / (Scalar(2) * (Scalar(1) - ct));

      // β 及其对 θ 的导数相关量（构造耦合块所需的标量系数），θ→0 均切泰勒
      const Scalar beta = if_then_else(
        LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
        static_cast<Scalar>(Scalar(1) / Scalar(12) + t2 / Scalar(720)), // θ→0 泰勒
        static_cast<Scalar>(t2inv - st * tinv * inv_2_2ct)              // 一般
      );

      const Scalar beta_dot_over_theta = if_then_else(
        LT, t, TaylorSeriesExpansion<Scalar>::template precision<3>(),
        static_cast<Scalar>(Scalar(1) / Scalar(360)), // θ→0 泰勒
        static_cast<Scalar>(
          -Scalar(2) * t2inv * t2inv + (Scalar(1) + st * tinv) * t2inv * inv_2_2ct) // 一般
      );

      // ---- 用 C 作临时暂存，拼出耦合项，再 B = C·A ----
      const Scalar wTp = w.dot(p); // ωᵀp
      const Vector3 v3_tmp(
        (beta_dot_over_theta * wTp) * w - (t2 * beta_dot_over_theta + Scalar(2) * beta) * p);
      // C can be treated as a temporary variable
      C.noalias() = v3_tmp * w.transpose();       // 外积项
      C.noalias() += beta * w * p.transpose();    // + β ω pᵀ
      C.diagonal().array() += wTp * beta;         // + (ωᵀp)β·I
      addSkew(Scalar(.5) * p, C);                 // + ½[p]×

      B.noalias() = C * A; // 耦合块 B = C·A（Jlog3）
      C.setZero();         // 左下块清零
    }
  };

} // namespace pinocchio
