//
// Copyright (c) 2015-2024 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: private, include "pinocchio/multibody.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/multibody.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<
    typename NewScalar,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
    // 定义cast之后的返回类型
  struct CastType<NewScalar, ModelTpl<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<NewScalar, Options, JointCollectionTpl> type;
  };

  // traits 特化：把 Model 的关联类型（Scalar / Data / JointCollection）暴露出去。
  // 注意 Data 是从 Model 反查出来的，因此 Model 与 Data 严格配对、不能混用。
  template<typename _Scalar, int _Options, template<typename, int> class JointCollectionTpl>
  struct traits<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    typedef JointCollectionTpl<Scalar, Options> JointCollection;
  };

  // ===========================================================================
  // ModelTpl：运动学树的「静态定义」，即"这台机器人是什么"。
  //   - 只读语义：addJoint / addFrame 建模完成后，算法阶段不再修改它；
  //     所有算法签名都是 (const Model &, Data &)，可被多线程共享。
  //   - 与 DataTpl 的分工：Model 存不随位形变化的量（拓扑、惯量、限位），
  //     Data 存随 q/v 变化的中间结果（oMi、J、M 等）。
  //   - 三个基类：Serializable 提供序列化；NumericalBase 提供 Scalar typedef；
  //     ModelEntity 是 CRTP 实体标记。
  // ===========================================================================
  template<typename _Scalar, int _Options, template<typename, int> class JointCollectionTpl>
  struct ModelTpl
  : serialization::Serializable<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  , NumericalBase<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  , ModelEntity<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  {

    typedef typename traits<ModelTpl>::Scalar Scalar; // 标量类型（double / float / CppAD 等）
    static constexpr int Options = traits<ModelTpl>::Options; // Eigen 对齐/存储选项

    typedef typename traits<ModelTpl>::JointCollection JointCollection; // 关节"菜单"（variant）
    typedef typename traits<ModelTpl>::Data Data;                       // 配对的 Data 类型

    // --- 空间代数基本类型（见 spatial/） ---
    typedef SE3Tpl<Scalar, Options> SE3;         // 刚体位姿（旋转 + 平移）
    typedef MotionTpl<Scalar, Options> Motion;   // 空间速度/加速度（旋量）
    typedef ForceTpl<Scalar, Options> Force;     // 空间力（力旋量：力 + 力矩）
    typedef InertiaTpl<Scalar, Options> Inertia; // 空间惯量（质量 + 质心 + 转动惯量）
    typedef FrameTpl<Scalar, Options> Frame;     // 附着在关节上的操作坐标系

    // --- 索引类型：全部是 std::size_t 的别名，仅用于语义区分 ---
    typedef pinocchio::Index Index;           // 通用下标
    typedef pinocchio::JointIndex JointIndex; // 关节下标（0 = universe）
    typedef pinocchio::GeomIndex GeomIndex;   // 几何体下标（GeometryModel 用）
    typedef pinocchio::FrameIndex FrameIndex; // 坐标系下标
    typedef std::vector<Index> IndexVector;   // 下标列表（父子关系/支撑链/子树等）

    // --- 关节的类型擦除包装：variant，可在同一容器里存放任意关节类型 ---
    typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel; // 关节的静态定义
    typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;   // 关节的计算缓存

    typedef std::vector<JointModel> JointModelVector;
    typedef std::vector<JointData> JointDataVector;

    typedef std::vector<Frame> FrameVector;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;              // 动态列向量
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> MatrixXs; // 动态矩阵
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;                            // 三维向量

    typedef std::vector<Inertia> InertiaVector;
    typedef std::vector<SE3> SE3Vector;

    // --- 稀疏模式相关（记录每个关节影响雅可比的哪些列） ---
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector; // 稠密布尔掩码
    typedef std::vector<Eigen::Index> EigenIndexVector;                    // 非零列的下标列表
    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectorOfEigenIndexVector;

    /// \brief Dense vectorized version of a joint configuration vector.
    typedef VectorXs ConfigVectorType;

    /// \brief Map between a string (key) and a configuration vector
    typedef std::map<std::string, ConfigVectorType> ConfigVectorMap;

    /// \brief Dense vectorized version of a joint tangent vector (e.g. velocity, acceleration,
    /// etc).
    ///        It also handles the notion of co-tangent vector (e.g. torque, etc).
    typedef VectorXs TangentVectorType;

    // =========================================================================
    // 三套维度。因为位形空间是弯曲流形而速度住在切空间，三者一般互不相等：
    //   nq         位形向量 q 的长度      —— 求解器的"状态"
    //   nv         切空间 v 的维数        —— M / tau / qdot 的维数
    //   nvExtended 雅可比列数            —— 展开 mimic 后的维数
    // 无 mimic 关节时 nv == nvExtended。详见下方 idx_* / n*s 索引表的说明。
    // =========================================================================

    /// \brief Dimension of the configuration vector representation.
    int nq; // 位形空间维数：球副占 4（四元数）、浮动基占 7，故可能 > nv

    /// \brief Dimension of the velocity vector space.
    int nv; // 独立自由度个数：广义速度/加速度/力矩的维数

    /// \brief Dimension of the jacobian space.
    int nvExtended; // 雅可比列空间维数：mimic 关节虽不占 nv 但仍占一列，故可能 > nv

    /// \brief Number of joints.
    int njoints; // 关节数，含下标 0 的 universe，故最少为 1

    /// \brief Number of bodies.
    int nbodies; // 刚体数（同样含 universe）

    /// \brief Number of operational frames.
    int nframes; // 操作坐标系数量，等于 frames.size()

    // ---- 以下所有 std::vector 均以 JointIndex 为下标，第 0 项恒为 universe ----

    /// \brief Vector of spatial inertias supported by each joint.
    InertiaVector inertias; // inertias[i]：关节 i 所承载刚体的空间惯量（在关节系下表达）

    /// \brief Vector of joint placements: placement of a joint *i* wrt its parent joint frame.
    SE3Vector jointPlacements; // jointPlacements[i]：关节 i 相对父关节系的固定安装位姿 ^{li}M_i

    /// \brief Vector of joint models.
    JointModelVector joints; // joints[i]：关节 i 的模型（variant，内含 nq/nv/idx_q/idx_v 等）

    // ---- q / v / 雅可比 三套「关节 → 全局向量切片」查找表 ----
    // 用法：q 里属于关节 i 的那一段是 q.segment(idx_qs[i], nqs[i])。
    // 实际代码中不要手写 segment，用 jmodel.jointConfigSelector(q) /
    // jointVelocitySelector(v)，它们内部就是这两组下标。
    // 建表方式见 addJoint：idx_qs[i] = sum_{k<i} nqs[k]（前缀和），
    // 因为 parents[i] < i，所以下标天然按拓扑序递增。

    /// \brief Vector of starting index of the *i*th  joint in the configuration space.
    std::vector<int> idx_qs; // 关节 i 在 q 中的起始下标

    /// \brief Vector of dimension of the  joint configuration subspace.
    std::vector<int> nqs; // 关节 i 在 q 中占的维数（仅自己，不含子树）

    /// \brief Starting index of the *i*th joint in the tangent configuration space.
    std::vector<int> idx_vs; // 关节 i 在 v 中的起始下标（M / tau / J 的行列均按此排布）

    /// \brief Dimension of the *i*th joint tangent subspace.
    std::vector<int> nvs; // 关节 i 在 v 中占的维数

    /// \brief Starting index of the *i*th joint in the jacobian space.
    std::vector<int> idx_vExtendeds; // 关节 i 在 data.J 中的起始列（mimic 关节独占一列）

    /// \brief Dimension of the *i*th joint jacobian subspace.
    std::vector<int> nvExtendeds; // 关节 i 在 data.J 中占的列数

    // ---- 树的拓扑结构 ----

    /// \brief Vector of parent joint indexes. The parent of joint *i*, denoted *li*, corresponds to
    /// li==parents[i].
    std::vector<JointIndex> parents; // parents[i]：关节 i 的父关节。恒有 parents[i] < i，
                                     // 这是反向递推（RNEA/CRBA 的 backward pass）能倒序
                                     // 遍历 i = njoints-1 ... 1 的前提

    /// \brief Vector of children index. Chidren of the *i*th joint, denoted *mu(i)* corresponds to
    /// the set (i==parents[k] for k in mu(i)).
    std::vector<IndexVector> children; // children[i]：关节 i 的全部直接子关节（parents 的逆映射）

    /// \brief Vector of mimicking joints in the tree (with type MimicTpl)
    std::vector<JointIndex> mimicking_joints; // 所有 mimic（从动）关节的 id

    /// \brief Vector of mimicked joints in the tree (can be any joint type)
    /// The i-th element of this vector correspond to the mimicked joint of the i-th mimicking
    /// vector in mimicking_joints
    std::vector<JointIndex> mimicked_joints; // 与上表按下标一一对应的被模仿（主动）关节 id

    /// \brief Name of the joints.
    std::vector<std::string> names; // names[i]：关节名，names[0] == "universe"

    /// \brief Map of reference configurations, indexed by user given names.
    ConfigVectorMap referenceConfigurations; // 具名参考位形（如 URDF/SRDF 里的 "half_sitting"）

    // ---- 电机/传动相关的动力学参数，长度均为 nv，按 idx_v 排布 ----

    /// \brief Vector of armature values expressed at the joint level
    /// This vector may contain the contribution of rotor inertia effects for instance.
    VectorXs armature; // 电枢惯量：直接加到 M 的对角线上（M.diagonal() += armature），
                       // 用于补偿减速器折算后的转子惯量，也能改善 M 的条件数

    /// \brief Vector of rotor inertia parameters
    TangentVectorType rotorInertia; // 转子自身惯量 I_r（未经减速比折算的原始值）

    /// \brief Vector of rotor gear ratio parameters
    TangentVectorType rotorGearRatio; // 减速比 r，与上者合成 armature = r^2 * I_r

    /// \brief Vector of joint friction parameters
    /// Deprecated in favor of lowerDryFrictionLimit and upperDryFrictionLimit
    PINOCCHIO_DEPRECATED TangentVectorType & friction; // 已弃用：是 upperDryFrictionLimit 的引用

    /// \brief Vector of joint friction parameters
    TangentVectorType lowerDryFrictionLimit; // 库仑（干）摩擦力矩下界，通常为负

    /// \brief Vector of joint friction parameters
    TangentVectorType upperDryFrictionLimit; // 库仑（干）摩擦力矩上界

    /// \brief Vector of joint damping parameters
    TangentVectorType damping; // 粘滞阻尼系数 b：阻尼力矩 = -b * v，与速度成正比

    // ---- 限位。effort/velocity/friction 长度为 nv；position 长度为 nq ----

    /// \brief Vector of minimal joint torques
    TangentVectorType lowerEffortLimit; // 关节力矩下界

    /// \brief Vector of maximal joint torques
    TangentVectorType upperEffortLimit; // 关节力矩上界

    /// \brief Vector of maximal joint torques
    /// Deprecated in favor of lowerEffortLimit and upperEffortLimit
    PINOCCHIO_DEPRECATED TangentVectorType & effortLimit; // 已弃用：upperEffortLimit 的引用

    /// \brief Vector of minimal joint velocities
    TangentVectorType lowerVelocityLimit; // 关节速度下界

    /// \brief Vector of maximal joint velocities
    TangentVectorType upperVelocityLimit; // 关节速度上界

    /// \brief Vector of maximal joint velocities
    /// Deprecated in favor of lowerVelocityLimit and upperVelocityLimit
    PINOCCHIO_DEPRECATED TangentVectorType & velocityLimit; // 已弃用：upperVelocityLimit 的引用

    /// \brief Lower joint configuration limit
    ConfigVectorType lowerPositionLimit; // 位形下界，长度 nq，按 idx_q 排布（注意不是 nv！）

    /// \brief Upper joint configuration limit
    ConfigVectorType upperPositionLimit; // 位形上界，长度 nq
                                         // 无界维度（SO(2)/SO(3) 的分量）填 ±inf，
                                         // 可用 hasConfigurationLimit() 判断哪些维度有意义

    /// \brief Joint configuration limit margin
    ConfigVectorType positionLimitMargin; // 位形限位的安全裕度，长度 nq；
                                          // 碰撞/规划时把可行域再向内收缩这么多

    /// \brief Vector of operational frames registered on the model.
    FrameVector frames; // 所有操作坐标系（末端、传感器、固定关节折叠后的残留等）。
                        // 坐标系不是自由度，只是挂在某个关节上的固定偏移：oMf = oMi * placement

    // ---- 预计算的树遍历结果：把 O(n) 的路径搜索变成 O(1) 查表 ----

    /// \brief Vector of joint supports.
    /// supports[j] corresponds to the vector of indices of the joints located on the path between
    /// joint *j*  and "universe".
    /// The first element of supports[j] is "universe", the last one is the index of joint *j*
    /// itself.
    std::vector<IndexVector> supports; // 支撑链：从 universe 到关节 j 的完整路径（含首尾）。
                                       // 关节 j 的雅可比只有这条链上的列非零

    /// \brief Vector of mimic supports joints.
    /// mimic_joint_supports[j] corresponds to the vector of mimic joints indices located on the
    /// path between joint *j*  and "universe". The first element of mimic_joint_supports[j] is
    /// "universe". If *j* is a mimic, the last element is the index of joint *j* itself.
    std::vector<IndexVector> mimic_joint_supports; // supports[j] 中只保留 mimic 关节的子序列，
                                                   // 用于雅可比折叠时的第二趟累加

    /// \brief Vector of joint subtrees.
    /// subtree[j] corresponds to the subtree supported by the joint *j*.
    /// The first element of subtree[j] is the index of the joint *j* itself.
    std::vector<IndexVector> subtrees; // 子树：以关节 j 为根的全部后代（含 j 自己）。
                                       // CRBA 用它确定复合刚体惯量要累加哪些刚体

    /// \brief Sparsity pattern for each joint.
    /// sparsity_pattern_vector[i] is a boolean vector of size nv indicating which columns
    /// of the Jacobian are nonzero for joint i.
    VectorOfBooleanVector sparsity_pattern_vector; // 稠密布尔掩码版：适合做按位与/或

    /// \brief Colwise span indexes for each joints.
    /// span_indexes_vector[i] lists the column indexes of nonzero entries for joint i.
    VectorOfEigenIndexVector span_indexes_vector; // 同一信息的下标列表版：适合直接遍历非零列

    /// \brief Spatial gravity of the model.
    Motion gravity; // 重力的空间加速度表示（linear = g，angular = 0）。
                    // RNEA 通过把基座加速度初始化为 -gravity 来"免费"算出重力项

    /// \brief Default 3D gravity vector (=(0,0,-9.81)).
    static const Vector3 gravity981; // 默认重力向量 (0,0,-9.81)，定义在本文件末尾

    /// \brief Model name.
    std::string name; // 模型名（一般取自 URDF 的 <robot name="...">）

    /// \brief Default constructor. Builds an empty model with no joints.
    // 默认构造：建出一棵「只有 universe 的空树」。
    // 所有以 JointIndex 为下标的容器都被初始化成长度 1、内容为零/单位元的状态，
    // 这就是"下标 0 恒为 universe"这条全局约定的源头；真实关节从 1 开始编号。
    // 三个 PINOCCHIO_DEPRECATED 引用成员（friction/effortLimit/velocityLimit）
    // 必须在初始化列表里绑定，故此处用 DIAGNOSTIC_PUSH/POP 抑制弃用告警。
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    ModelTpl()
    : nq(0)          // 空树没有自由度
    , nv(0)
    , nvExtended(0)
    , njoints(1)     // 但已经有一个"关节"：universe
    , nbodies(1)
    , nframes(0)     // frames 在函数体里通过 addFrame 补上，故这里是 0
    , inertias(1, Inertia::Zero())
    , jointPlacements(1, SE3::Identity())
    , joints(1)
    , idx_qs(1, 0)   // universe 不占 q/v/雅可比的任何一维，六张索引表全填 0
    , nqs(1, 0)
    , idx_vs(1, 0)
    , nvs(1, 0)
    , idx_vExtendeds(1, 0)
    , nvExtendeds(1, 0)
    , parents(1, 0)  // universe 的父亲是它自己，用于终止向上遍历
    , children(1)
    , names(1)
    , friction(upperDryFrictionLimit)   // 三个弃用名只是新名字的别名，不占额外存储
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    , supports(1, IndexVector(1, 0))            // universe 的支撑链就是它自己
    , mimic_joint_supports(1, IndexVector(1, 0))
    , subtrees(1)
    , sparsity_pattern_vector(1)
    , span_indexes_vector(1)
    , gravity(gravity981, Vector3::Zero())      // 线速度部分 = (0,0,-9.81)，角速度部分 = 0
    {
      names[0] = "universe"; // Should be "universe joint (trivial)"
      // FIXME Should the universe joint be a FIXED_JOINT even if it is
      // in the list of joints ? See comment in definition of
      // Model::addJointFrame and Model::addBodyFrame
      addFrame(Frame("universe", 0, 0, SE3::Identity(), FIXED_JOINT));
    }
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    ///
    /// \brief Copy constructor by casting
    ///
    /// \param[in] other model to copy to *this
    ///
    // 换标量类型的拷贝构造：double 模型 -> float / CppAD::AD<double> 等。
    // 自动微分场景常用（把整个模型提升为 AD 类型后再求导）。
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    template<typename S2, int O2>
    explicit ModelTpl(const ModelTpl<S2, O2> & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other.template cast<Scalar>();
    }
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    ///
    /// \brief Copy constructor from another collection
    ///
    /// \param[in] other model to copy to *this
    ///
    // 换关节集合的拷贝构造：标量不变，但目标 JointCollection 可能包含
    // 自定义关节类型。用于扩展了关节"菜单"之后在两种 Model 之间转换。
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    template<template<typename, int> class OtherJointCollectionTpl>
    ModelTpl(const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other;
    }
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    ///
    /// \brief Copy constructor.
    ///
    /// \param[in] other model to copy to *this
    ///
    // 普通拷贝构造。之所以不能用编译器生成的版本，是因为那三个引用成员
    // 必须重新绑定到「本对象」的向量上，否则会悬垂到源对象去。
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    ModelTpl(const ModelTpl & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other;
    }
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    /// \returns A new copy of *this with the Scalar type casted to NewScalar.
    // 标量类型转换：返回一份把所有数值成员都转成 NewScalar 的新模型。
    // 返回类型由文件开头的 CastType 特化给出。典型用途：
    //   Model 转 ModelTpl<CppAD::AD<double>>，再对整条动力学链求解析导数。
    template<typename NewScalar>
    typename CastType<NewScalar, ModelTpl>::type cast() const;

    ///
    /// \brief Equality comparison operator.
    ///
    /// \returns true if *this is equal to other.
    ///
    // 逐字段深比较（维度 / 拓扑 / 惯量 / 限位 / 关节模型全都要一致）。
    // 常用于校验"手工搭的模型"与"URDF 解析出的模型"是否等价。
    bool operator==(const ModelTpl & other) const;

    ///
    /// \brief Assignment operator from another collection.
    ///
    ///
    // 跨 JointCollection 赋值。真正的实现在这一版里，下面的同类型版本转调它。
    template<template<typename, int> class OtherJointCollectionTpl>
    ModelTpl & operator=(const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other);

    ///
    /// \brief Assignment operator.
    ///
    ///
    // 同类型赋值：显式指定模板实参转调上面那个版本。
    // 同样不能用编译器默认版本——引用成员不可重新赋值。
    ModelTpl & operator=(const ModelTpl & other)
    {
      (*this).template operator= <JointCollectionTpl>(other);
      return *this;
    }

    ///
    /// \returns true if *this is NOT equal to other.
    ///
    bool operator!=(const ModelTpl & other) const
    {
      return !(*this == other);
    }

    // =========================================================================
    // addJoint 重载家族（共 6 个）。它是唯一能改变树拓扑的入口，职责是：
    //   1. 分配新的 joint_id = njoints++；
    //   2. 调 jmodel.setIndexes(id, nq, nv, nvExtended) —— 传入的是「累加之前」
    //      的总维数，于是它天然就是本关节的 idx_q/idx_v/idx_vExtended（前缀和）；
    //   3. 累加 nq/nv/nvExtended，并向六张索引表 push_back；
    //   4. resize 所有限位/摩擦/阻尼向量并写入本关节那一段；
    //   5. 维护 parents/children/subtrees/supports/稀疏模式。
    //
    // 【重要约束】必须按深度优先顺序添加，保证 parents[i] < i，
    //            否则反向递推与稀疏 Cholesky 的假设全部失效。
    // 【注意】它不会自动建同名 Frame，需要的话另外调 addJointFrame。
    //
    // 下面各重载只是参数丰俭不同，最终都转调「参数最全」的那一个：
    //   (min/max effort, min/max velocity, min/max config, config_limit_margin,
    //    min/max friction, damping)
    // =========================================================================

    ///
    /// \brief Add a joint to the kinematic tree with infinite bounds.
    ///
    /// \remarks This method does not add a Frame of same name to the vector of frames.
    ///         Use Model::addJointFrame.
    /// \remarks The inertia supported by the joint is set to Zero.
    /// \remark Joints need to be added to the tree in a depth-first order.
    ///
    // 最简版本：不给任何限位，全部取 ±inf；惯量置零（之后用 appendBodyToJoint 补）。
    /// \tparam JointModelDerived The type of the joint model.
    ///
    /// \param[in] parent Index of the parent joint.
    /// \param[in] joint_model The joint model.
    /// \param[in] joint_placement Placement of the joint inside its parent joint.
    /// \param[in] joint_name Name of the joint. If empty, the name is random.
    ///
    /// \return The index of the new joint.
    ///
    /// \sa Model::appendBodyToJoint
    ///
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] max_effort Maximal joint torque.
    /// \param[in] max_velocity Maximal joint velocity.
    /// \param[in] min_config Lower joint configuration.
    /// \param[in] max_config Upper joint configuration.
    ///
    // 对称限位版（已弃用）：只给上界，下界自动取 -max_*。
    // 注意 effort/velocity 长度为 nv，而 min_config/max_config 长度为 nq。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] max_effort Maximal joint torque.
    /// \param[in] max_velocity Maximal joint velocity.
    /// \param[in] min_config Lower joint configuration.
    /// \param[in] max_config Upper joint configuration.
    /// \param[in] config_limit_margin Joint configuration limit margin.
    ///
    // 同上，额外指定位形限位的安全裕度（已弃用）。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & config_limit_margin);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] min_effort Minimal joint torque.
    /// \param[in] min_velocity Minimal joint velocity.
    /// \param[in] min_friction Minimal joint friction parameters.
    /// \param[in] max_friction Maximal joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    // 非对称限位版（已弃用）：上下界独立给出，并带摩擦/阻尼。
    // 实现里把 config_limit_margin 补零后转调最全的那一版。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & min_effort,
      const VectorXs & max_effort,
      const VectorXs & min_velocity,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & min_friction,
      const VectorXs & max_friction,
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] config_limit_margin Joint configuration limit margin.
    /// \param[in] min_effort Minimal joint torque.
    /// \param[in] min_velocity Minimal joint velocity.
    /// \param[in] min_friction Minimal joint friction parameters.
    /// \param[in] max_friction Maximal joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    // ★ 参数最全的版本：其余 5 个重载最终都收敛到这里，真正干活的实现。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & min_effort,
      const VectorXs & max_effort,
      const VectorXs & min_velocity,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & config_limit_margin,
      const VectorXs & min_friction,
      const VectorXs & max_friction,
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    ///
    /// \param[in] config_limit_margin Joint configuration limit margin.
    /// \param[in] friction Joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    // 对称限位 + 对称摩擦 + 裕度：friction 同时作为上下界（取 ±friction）。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & config_limit_margin,
      const VectorXs & friction,
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    ///
    /// \param[in] friction Joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    // 同上但不给裕度（内部补零）。URDF 解析器最常走这一条。
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & friction,
      const VectorXs & damping);

    ///
    /// \brief Add a joint to the frame tree.
    ///
    /// \param[in] jointIndex Index of the joint.
    /// \param[in] frameIndex Index of the parent frame. If negative,
    ///            the parent frame is the frame of the parent joint.
    ///
    /// \return The index of the new frame
    ///
    // 给已存在的关节补一个 JOINT 类型的 Frame，使它能按名字被 getFrameId 找到。
    // addJoint 刻意不做这件事，所以建模流程通常是 addJoint 之后紧跟本函数。
    FrameIndex addJointFrame(const JointIndex joint_index, int previous_frame_index = -1);

    ///
    /// \brief Append a body to a given joint of the kinematic tree.
    ///
    /// \param[in] joint_index Index of the supporting joint.
    /// \param[in] Y Spatial inertia of the body.
    /// \param[in] body_placement The relative placement of the body regarding to the parent joint.
    /// Set default to the Identity placement.
    ///
    /// \sa Model::addJoint
    ///
    // 给关节挂上刚体质量属性：把 Y 经 body_placement 变换到关节系后
    // 【累加】到 inertias[joint_index] 上（是 +=，不是覆盖）。
    // 因此同一关节可以挂多个刚体，等价于把它们复合成一个等效刚体——
    // 这正是 URDF 里"固定关节被折叠"后惯量合并的实现方式。
    void appendBodyToJoint(
      const JointIndex joint_index,
      const Inertia & Y,
      const SE3 & body_placement = SE3::Identity());

    ///
    /// \brief Add a body to the frame tree.
    ///
    /// \param[in] body_name Name of the body.
    /// \param[in] parentJoint Index of the parent joint.
    /// \param[in] body_placement The relative placement of the body regarding to the parent joint.
    /// Set default to the Identity placement. \param[in] parentFrame Index of the parent frame. If
    /// negative,
    ///            the parent frame is the frame of the parent joint.
    ///
    /// \return The index of the new frame
    ///
    // 注册一个 BODY 类型的 Frame（纯几何标记，不携带质量、不产生自由度）。
    // 与 appendBodyToJoint 的区别：那个改惯量，这个只加一个可查询的坐标系。
    FrameIndex addBodyFrame(
      const std::string & body_name,
      const JointIndex & parentJoint,
      const SE3 & body_placement = SE3::Identity(),
      int parentFrame = -1);

    // ---- 名字 → 下标 的查询接口 ----
    // 【共同陷阱】找不到时不抛异常，而是返回"容器当前长度"这个越界哨兵值。
    // 因此结果必须先用对应的 exist* 校验，或至少与 nframes/njoints 比一比，
    // 否则一旦之后又往模型里加了元素，这个哨兵值就会静默变成合法下标。

    ///
    /// \brief Return the index of a body given by its name.
    ///
    /// \warning If no body is found, return the number of elements at time T.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent body)
    ///
    /// \param[in] name Name of the body.
    ///
    /// \return Index of the body.
    ///
    // 注意返回的是 FrameIndex（frames 里的下标），不是 JointIndex。
    FrameIndex getBodyId(const std::string & name) const;

    ///
    /// \brief Check if a body given by its name exists.
    ///
    /// \param[in] name Name of the body.
    ///
    /// \return True if the body exists in the kinematic tree.
    ///
    bool existBodyName(const std::string & name) const;

    ///
    /// \brief Return the index of a joint given by its name.
    ///
    /// \warning If no joint is found, return the number of elements at time T.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent joint)
    /// \param[in] name Name of the joint.
    ///
    /// \return Index of the joint.
    ///
    // 关节名 → JointIndex，可直接用于 data.oMi[id] 等按关节下标的容器。
    JointIndex getJointId(const std::string & name) const;

    ///
    /// \brief Check if a joint given by its name exists.
    ///
    /// \param[in] name Name of the joint.
    ///
    /// \return True if the joint exists in the kinematic tree.
    ///
    bool existJointName(const std::string & name) const;

    ///
    /// \brief Returns the index of a frame given by its name.
    ///        \sa Model::existFrame to check if the frame exists or not.
    ///
    /// \warning If no frame is found, returns the size of the vector of Model::frames.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent frame).
    ///
    /// \param[in] name Name of the frame.
    /// \param[in] type Type of the frame.
    ///
    /// \return Index of the frame.
    ///
    // type 是【位掩码】而非单一枚举值，默认匹配全部五种类型。
    // 同名但不同类型的 Frame 可以共存（例如关节名与连杆名相同时），
    // 此时必须显式传 type 才能选中想要的那一个。
    FrameIndex getFrameId(
      const std::string & name,
      const FrameType & type = (FrameType)(JOINT | FIXED_JOINT | BODY | OP_FRAME | SENSOR)) const;

    ///
    /// \brief Checks if a frame given by its name exists.
    ///
    /// \param[in] name Name of the frame.
    /// \param[in] type Type of the frame.
    ///
    /// \return Returns true if the frame exists.
    ///
    bool existFrame(
      const std::string & name,
      const FrameType & type = (FrameType)(JOINT | FIXED_JOINT | BODY | OP_FRAME | SENSOR)) const;

    ///
    /// \brief Adds a frame to the kinematic tree.
    ///        The inertia stored within the frame will be happened to the inertia supported by the
    ///        joint (frame.parentJoint).
    ///
    /// \param[in] frame The frame to add to the kinematic tree.
    /// \param[in] append_inertia Append the inertia contained in the Frame to the inertia supported
    /// by the joint.
    ///
    /// \return Returns the index of the frame if it has been successfully added or if it already
    /// exists in the kinematic tree.
    ///
    // 添加坐标系，并可选地把 frame 自带的惯量并入其父关节（append_inertia=true）。
    // 幂等：同名同类型已存在时直接返回原下标，不会重复插入。
    // Frame 不增加自由度，只是"挂在某关节上的固定偏移"，因此加多少个都不影响 nq/nv。
    FrameIndex addFrame(const Frame & frame, const bool append_inertia = true);

    ///
    /// \brief Check the validity of the attributes of Model with respect to the specification of
    /// some algorithms.
    ///
    /// The method is a template so that the checkers can be defined in each algorithms.
    /// \param[in] checker a class, typically defined in the algorithm module, that
    /// validates the attributes of model.
    ///
    /// \return true if the Model is valid, false otherwise.
    ///
    // 模板化的自检入口：具体的校验规则由各算法模块提供的 checker 定义
    // （如 ABA 要求父下标递增、CRBA 要求子树连续等），实现了"算法自带前置条件"。
    template<typename D>
    bool check(const AlgorithmCheckerBase<D> & checker) const
    {
      return checker.checkModel(*this);
    }

    ///
    /// \brief Check if joints have configuration limits
    ///
    /// \return Returns list of boolean of size model.nq.
    ///
    // 逐维回答"这一维在位形空间里有没有上下界"，长度 nq，与 lowerPositionLimit 对齐。
    // false 表示该维住在紧致无边流形上（SO(2)/SO(3)），限位这个概念根本没定义，
    // 而不是"暂时没设置"。它是关节【类型】决定的静态属性，与实际填的数值无关，
    // 因此判断可行域时还需再叠加一次 isfinite 检查。
    std::vector<bool> hasConfigurationLimit() const;

    ///
    /// \brief Check if joints have configuration limits
    ///
    /// \return Returns list of boolean of size model.nq.
    ///
    // 同上的切空间版本，长度为【nv】而非 nq（上方 doxygen 注释写成 nq 是笔误）。
    // 差别只出现在 nq != nv 的关节上：FreeFlyer 在 q 里给出 7 项、在 v 里给出 6 项。
    // 做 IK/轨迹优化时决策变量是切空间增量，必须用这个版本才能与下标对齐。
    std::vector<bool> hasConfigurationLimitInTangent() const;

    /// Run check(fusion::list) with DEFAULT_CHECKERS as argument.
    // 用默认检查器集合做一次自检。
    bool check() const;

    ///
    /// \brief Run checkData on data and current model.
    ///
    /// \param[in] data to be checked wrt *this.
    ///
    /// \return true if the data is valid, false otherwise.
    ///
    // 校验一个 Data 是否与本 Model 匹配（各容器尺寸是否符合 nq/nv/njoints）。
    // 算法入口处普遍有 assert(model.check(data))，防止误用别的模型建出的 Data。
    bool check(const Data & data) const;

    ///
    /// \brief Create a Data structure associated with the current model
    ///
    // 工厂方法：按当前模型的维度分配好全部工作区。建模完成后调用一次即可，
    // 之后反复复用同一个 Data 跑算法（实时循环里绝不要每帧新建）。
    // 多线程时应为每个线程各建一份 Data，Model 则可只读共享。
    Data createData() const;

    /// Returns a vector of the children joints of the kinematic tree.
    /// \remark: a child joint is a node without any child joint.
    // 返回所有【叶子】关节（children 为空的节点），即各条运动链的末端。
    std::vector<JointIndex> getChildJoints() const;

  protected:
    ///
    /// \brief Add the joint_id to its parent subtrees.
    ///
    /// \param[in] joint_id The id of the joint to add to the subtrees
    ///
    // addJoint 的内部辅助：沿 parents 一路上溯到 universe，
    // 把新关节的 id 追加进沿途每个祖先的 subtrees 里，从而增量维护子树表。
    void addJointIndexToParentSubtrees(const JointIndex joint_id);
  };

} // namespace pinocchio

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

namespace pinocchio
{
  namespace details
  {
    struct FilterFrame
    {
      const std::string & name;
      const FrameType & typeMask;

      FilterFrame(const std::string & name, const FrameType & typeMask)
      : name(name)
      , typeMask(typeMask)
      {
      }

      template<typename Scalar, int Options>
      bool operator()(const FrameTpl<Scalar, Options> & frame) const
      {
        return (typeMask & frame.type) && (name == frame.name);
      }
    };
  } // namespace details

  // 静态变量类内声明，类外定义。默认重力向量 = (0,0,-9.81)。
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  const typename ModelTpl<Scalar, Options, JointCollectionTpl>::Vector3
    ModelTpl<Scalar, Options, JointCollectionTpl>::gravity981((Scalar)0, (Scalar)0, (Scalar)-9.81);

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline std::ostream &
  operator<<(std::ostream & os, const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
  {
    typedef typename ModelTpl<Scalar, Options, JointCollectionTpl>::Index Index;

    os << "Nb joints = " << model.njoints << " (nq=" << model.nq << ",nv=" << model.nv << ")"
       << std::endl;
    for (Index i = 0; i < (Index)(model.njoints); ++i)
    {
      os << "  Joint " << i << " " << model.names[i] << ": parent=" << model.parents[i]
         << std::endl;
    }

    return os;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & min_effort,
    const VectorXs & max_effort,
    const VectorXs & min_velocity,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & min_joint_friction,
    const VectorXs & max_joint_friction,
    const VectorXs & joint_damping)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));
    return addJoint(
      parent, joint_model, joint_placement, joint_name, min_effort, max_effort, min_velocity,
      max_velocity, min_config, max_config, config_limit_margin, min_joint_friction,
      max_joint_friction, joint_damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & min_effort,
    const VectorXs & max_effort,
    const VectorXs & min_velocity,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin,
    const VectorXs & min_joint_friction,
    const VectorXs & max_joint_friction,
    const VectorXs & joint_damping)
  {
    // ============================================================
    // addJoint 主实现（13 参数版，其余重载最终都汇聚到这里）。
    // 作用：把一个关节追加到运动学树，并同步维护 Model 的全部并行数组与拓扑。
    // 十个阶段：① 校验 → ② 分配 id/拷贝关节/写回索引 → ③ 读回本关节尺寸
    //   → ④ push 树结构字段 → ⑤ 累加全局维度 → ⑥ resize+写入限位/驱动向量
    //   → ⑦ subtrees → ⑧ supports → ⑨ 雅可比稀疏模式 → ⑩ mimic 记账。
    // 前提：关节必须按【深度优先】添加，保证 parents[i] < i（递推算法的不变量）。
    // ============================================================

    // —— 阶段①：一致性断言 + 入参尺寸/范围校验 ——
    //
    // 【为什么全部校验都堆在最前面】
    // 后面阶段②~⑩ 会向十几个并行容器 push_back，中途抛异常会留下一棵半更新的树
    // （njoints 与各容器长度错位），而且没有任何回滚机制。所以这里是"要么全做、
    // 要么不做"的守门人：先把能查的错一次查完，之后的写入过程不再有失败路径。
    //
    // 【两种检查性质完全不同，不要混为一谈】
    //   assert(...)            受 NDEBUG 控制，Release 里整条编译掉（零开销）。
    //                          只断言【库自身的内部不变量】；触发 = Pinocchio 的 bug
    //                          或用户直接篡改了公开字段，不属于正常的输入错误。
    //   PINOCCHIO_CHECK_*(...) 只受 PINOCCHIO_NO_THROW 控制，【Release 下照样生效】，
    //                          失败抛 std::invalid_argument（见 macros.hxx:188-226）。
    //                          校验【用户入参】，是公开 API 契约的一部分。
    //                          addJoint 只在建模期调用，这点开销无所谓。

    // (a) 内部不变量：所有以 JointIndex 为下标的并行数组必须始终等长。
    //     它们在阶段④成组 push_back，长度一旦错位，后续所有按 id 的索引都会读错。
    // 四个并行数组长度必须都等于 njoints（不变量：它们始终同步增长）。
    assert(
      (njoints == (int)joints.size()) && (njoints == (int)inertias.size())
      && (njoints == (int)parents.size()) && (njoints == (int)jointPlacements.size()));
    // (b) 关节自报的三个维度不能为负（防御自定义关节类型把 traits 写错）。
    assert((joint_model.nq() >= 0) && (joint_model.nv() >= 0) && (joint_model.nvExtended() >= 0));
    assert(joint_model.nq() >= joint_model.nv()); // 位形维 ≥ 速度维（弯流形关节 nq>nv）

    // (c) 尺寸校验之一：所有【切空间量】长度必须等于 joint_model.nv()。
    //     力矩、速度、干摩擦都是作用在广义速度上的量，故按 nv 计。
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_effort.size(), joint_model.nv(), "The joint minimal effort vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_joint_friction.size(), joint_model.nv(),
      "The joint minimal dry friction vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_velocity.size(), joint_model.nv(),
      "The joint minimal velocity vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_effort.size(), joint_model.nv(), "The joint maximum effort vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_joint_friction.size(), joint_model.nv(),
      "The joint maximum dry friction vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_velocity.size(), joint_model.nv(),
      "The joint maximum velocity vector is not of right size");

    // (d) 区间合法性：逐分量要求 min <= max，否则可行域为空集，
    //     后续 randomConfiguration / 规划器会静默给出无意义结果而不是报错。
    //     compareAll(..., LE) 是"全部分量都满足"的归约，不是任意一个满足即可。
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      compareAll(min_effort, max_effort, internal::ComparisonOperators::LE),
      "Some components of min_effort are greater than max_effort");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      compareAll(min_joint_friction, max_joint_friction, internal::ComparisonOperators::LE),
      "Some components of min_dry_friction are greater than max_dry_friction");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      compareAll(min_velocity, max_velocity, internal::ComparisonOperators::LE),
      "Some components of min_velocity are greater than max_velocity");

    // (e) 尺寸校验之二：所有【位形量】长度必须等于 joint_model.nq()，而不是 nv！
    //     这是本函数最容易踩的地方——同一个关节，effort 传 nv 维、config 传 nq 维。
    //     以 FreeFlyer 为例：min_effort 是 6 维，min_config 却是 7 维。
    //     两者搞反时不会静默出错，正是靠这几行抛异常挡住。
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_config.size(), joint_model.nq(),
      "The joint lower configuration bound is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_config.size(), joint_model.nq(),
      "The joint upper configuration bound is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      config_limit_margin.size(), joint_model.nq(),
      "The joint config limit margin is not of right size");
    // 注意 damping 又回到 nv（它乘的是速度），别被上面三行带偏。
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      joint_damping.size(), joint_model.nv(), "The joint damping vector is not of right size");

    // (f) 拓扑合法性：父关节必须是【已经存在】的节点。
    //     新关节拿到的 id 就是当前 njoints，所以这条等价于 parent < 新 id，
    //     从而保证 parents[i] < i ——反向递推（RNEA/CRBA 倒序遍历）与稀疏 Cholesky
    //     都依赖这个不变量。这也是"必须按深度优先添加"的强制点。
    //     JointIndex 是无符号类型，故不必再检查 parent >= 0。
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      parent < (JointIndex)njoints, "The index of the parent joint is not valid.");

    // —— 阶段②：分配关节 id、把关节拷贝进 joints、把全局起点写回关节对象 ——
    JointIndex joint_id = (JointIndex)(njoints++); // 新 id = 旧 njoints，然后 njoints 自增

    joints.push_back(JointModel(joint_model.derived())); // 拷贝一份存入树（variant 包装）
    JointModel & jmodel = joints.back();
    // 把"本关节在全局 q/v/vExtended 里的起点"= 当前累计的 nq/nv/nvExtended 写回关节自身，
    // 之后 jmodel.idx_q()/idx_v() 就能从全局向量里定位自己那几维（见 joint-model-base 的 i_q/i_v）。
    jmodel.setIndexes(joint_id, nq, nv, nvExtended);

    // —— 阶段③：从关节对象读回它的尺寸与刚写入的起点 ——
    const int joint_nq = jmodel.nq();
    const int joint_idx_q = jmodel.idx_q();
    const int joint_nv = jmodel.nv();
    const int joint_idx_v = jmodel.idx_v();
    const int joint_nvExtended = jmodel.nvExtended();
    const int joint_idx_vExtended = jmodel.idx_vExtended();

    assert(joint_idx_q >= 0);
    assert(joint_idx_v >= 0);
    assert(joint_idx_vExtended >= 0);

    // —— 阶段④：追加逐关节的静态属性（并行数组同步增长一格） ——
    inertias.push_back(Inertia::Zero());   // 惯量先置零，之后 appendBodyToJoint 再累加刚体
    parents.push_back(parent);             // 记录父关节
    children.push_back(IndexVector());     // 本关节暂无子节点
    children[parent].push_back(joint_id);  // 反向：把自己登记进父的子列表
    jointPlacements.push_back(joint_placement); // 相对父关节系的固定位姿 ^{parent}M_i
    names.push_back(joint_name);

    // —— 阶段⑤：累加全局维度，并 push 本关节的段长/起点索引 ——
    nq += joint_nq;
    nqs.push_back(joint_nq);
    idx_qs.push_back(joint_idx_q);
    nv += joint_nv;
    nvs.push_back(joint_nv);
    idx_vs.push_back(joint_idx_v);
    nvExtended += joint_nvExtended;
    nvExtendeds.push_back(joint_nvExtended);
    idx_vExtendeds.push_back(joint_idx_vExtended);

    // —— 阶段⑥：把各限位/驱动向量扩容到新的 nv/nq，并写入本关节那一段 ——
    // conservativeResize 保留已有元素、只在尾部扩容；再用 jointVelocitySelector（idx_v,nv 段）
    // 或 jointConfigSelector（idx_q,nq 段）精确写入本关节对应的那几维。
    // fixed/静止关节（nq==0）无需写限位，故用 if 跳过。
    if (joint_nq > 0 && joint_nv > 0)
    {
      upperEffortLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperEffortLimit) = max_effort;
      lowerEffortLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(lowerEffortLimit) = min_effort;
      upperVelocityLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperVelocityLimit) = max_velocity;
      lowerVelocityLimit.conservativeResize(nv);
      // ⚠️ 疑似上游 BUG：此处赋 max_velocity，按同类项模式（lowerEffort=min_effort、
      //    lowerPosition=min_config、lowerDryFriction=min_joint_friction）应为 min_velocity。
      jmodel.jointVelocitySelector(lowerVelocityLimit) = max_velocity;
      lowerPositionLimit.conservativeResize(nq);
      jmodel.jointConfigSelector(lowerPositionLimit) = min_config;
      upperPositionLimit.conservativeResize(nq);
      jmodel.jointConfigSelector(upperPositionLimit) = max_config;
      positionLimitMargin.conservativeResize(nq);
      jmodel.jointConfigSelector(positionLimitMargin) = config_limit_margin;

      armature.conservativeResize(nv);
      jmodel.jointVelocitySelector(armature).setZero();
      rotorInertia.conservativeResize(nv);
      jmodel.jointVelocitySelector(rotorInertia).setZero();
      rotorGearRatio.conservativeResize(nv);
      jmodel.jointVelocitySelector(rotorGearRatio).setOnes();
      upperDryFrictionLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperDryFrictionLimit) = max_joint_friction;
      lowerDryFrictionLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(lowerDryFrictionLimit) = min_joint_friction;
      damping.conservativeResize(nv);
      jmodel.jointVelocitySelector(damping) = joint_damping;
    }

    // —— 阶段⑦：subtrees（子树）——本关节子树先只含自己，再把自己加进所有祖先的子树。
    // Init and add joint index to its parent subtrees.
    subtrees.push_back(IndexVector(1));
    subtrees[joint_id][0] = joint_id;
    addJointIndexToParentSubtrees(joint_id); // 沿 parents 上溯，把 joint_id 追加进每个祖先 subtree

    // —— 阶段⑧：supports（支撑路径）——继承父的"根→父"路径，再把自己接在末尾 = 根→本关节。
    // Init and add joint index to the supports
    supports.push_back(supports[parent]);
    supports[joint_id].push_back(joint_id);

    // —— 阶段⑨：雅可比稀疏模式 ——
    // 先把已存在的布尔向量都扩容到新 nv（conservativeResize 不初始化新元素，故手动清零尾部）。
    if (joint_nq > 0 && joint_nv > 0)
    {
      for (auto & sparsity : sparsity_pattern_vector)
      {
        const Eigen::Index old_size = sparsity.size();
        sparsity.conservativeResize(nv);
        sparsity.tail(nv - old_size).setZero();
      }
    }

    // 构建本关节的非零列集合 extended_support：= 支撑路径上所有祖先关节各自的 v 段 + 本关节 v 段。
    // 含义：末端关节 j 的空间速度由"根到 j 路径上所有关节的速度"决定，故这些列在雅可比里非零。
    // Build sparsity pattern and span indexes of the new joint.
    EigenIndexVector extended_support;
    extended_support.reserve(size_t(nv));
    const auto & jsupport = supports[joint_id];
    for (size_t j = 1; j < jsupport.size() - 1; ++j) // 遍历祖先（跳过首元 universe 与末元自己）
    {
      const JointIndex jsupport_id = jsupport[j];
      const int jsupport_nv = nvs[jsupport_id];
      const int jsupport_idx_v = idx_vs[jsupport_id];
      for (int k = 0; k < jsupport_nv; ++k) // 展开该祖先占据的 [idx_v, idx_v+nv) 列
        extended_support.push_back(jsupport_idx_v + k);
    }
    for (int k = 0; k < joint_nv; ++k) // 末尾补上本关节自己的 v 列
    {
      extended_support.push_back(joint_idx_v + k);
    }

    // span_indexes_vector 存"非零列索引列表"，sparsity_pattern_vector 存等价的布尔掩码。
    BooleanVector sparsity_pattern = BooleanVector::Zero(nv);
    for (const auto col_id : extended_support)
      sparsity_pattern[col_id] = true;

    sparsity_pattern_vector.push_back(std::move(sparsity_pattern));
    span_indexes_vector.push_back(std::move(extended_support));

    // —— 阶段⑩：mimic（耦合关节）记账 ——
    // mimic_joint_supports 继承父的 mimic 路径；仅当本关节是 mimic 类型时才追加自己，
    // 并登记 跟随者→被跟随者 的对应（见 §mimic 与 nvExtended）。
    // Update mimicking.
    mimic_joint_supports.push_back(mimic_joint_supports[parent]);
    if (
      const auto & jmodel_ =
        boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(&jmodel))
    {
      mimicking_joints.push_back(jmodel.id());           // 本(跟随)关节 id
      mimicked_joints.push_back(jmodel_->jmodel().id()); // 它所镜像的(被跟随)关节 id
      mimic_joint_supports[joint_id].push_back(joint_id);
    }
    return joint_id; // 返回新关节的全局 id，供 appendBodyToJoint / addJointFrame 引用
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin,
    const VectorXs & friction,
    const VectorXs & damping)
  {

    return addJoint(
      parent, joint_model, joint_placement, joint_name, -max_effort, max_effort, -max_velocity,
      max_velocity, min_config, max_config, config_limit_margin, -friction, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & friction,
    const VectorXs & damping)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, -max_effort, max_effort, -max_velocity,
      max_velocity, min_config, max_config, config_limit_margin, -friction, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));
    const VectorXs friction = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));
    const VectorXs damping = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin)
  {
    const VectorXs friction = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));
    const VectorXs damping = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name)
  {
    const VectorXs max_effort =
      VectorXs::Constant(joint_model.nv(), std::numeric_limits<Scalar>::max());
    const VectorXs max_velocity =
      VectorXs::Constant(joint_model.nv(), std::numeric_limits<Scalar>::max());
    const VectorXs min_config =
      VectorXs::Constant(joint_model.nq(), -std::numeric_limits<Scalar>::max());
    const VectorXs max_config =
      VectorXs::Constant(joint_model.nq(), std::numeric_limits<Scalar>::max());
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  FrameIndex ModelTpl<Scalar, Options, JointCollectionTpl>::addJointFrame(
    const JointIndex joint_index, int previous_frame_index)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      joint_index < joints.size(),
      "The joint index is larger than the number of joints in the model.");
    if (previous_frame_index < 0)
    {
      // FIXED_JOINT is required because the parent can be the universe and its
      // type is FIXED_JOINT
      previous_frame_index =
        (int)getFrameId(names[parents[joint_index]], (FrameType)(JOINT | FIXED_JOINT));
    }
    assert((size_t)previous_frame_index < frames.size() && "Frame index out of bound");

    // Add a the joint frame attached to itself to the frame vector - redundant information but
    // useful.
    return addFrame(Frame(
      names[joint_index], joint_index, (FrameIndex)previous_frame_index, SE3::Identity(), JOINT));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  template<typename NewScalar>
  typename CastType<NewScalar, ModelTpl<Scalar, Options, JointCollectionTpl>>::type
  ModelTpl<Scalar, Options, JointCollectionTpl>::cast() const
  {
    typedef ModelTpl<NewScalar, Options, JointCollectionTpl> ReturnType;

    ReturnType res;

    res.nq = nq;
    res.nv = nv;
    res.nvExtended = nvExtended;
    res.njoints = njoints;
    res.nbodies = nbodies;
    res.nframes = nframes;
    res.parents = parents;
    res.children = children;
    res.names = names;
    res.subtrees = subtrees;
    res.supports = supports;
    res.mimic_joint_supports = mimic_joint_supports;
    res.mimicking_joints = mimicking_joints;
    res.mimicked_joints = mimicked_joints;
    res.gravity = gravity.template cast<NewScalar>();
    res.name = name;
    res.sparsity_pattern_vector = sparsity_pattern_vector;
    res.span_indexes_vector = span_indexes_vector;
    res.idx_qs = idx_qs;
    res.nqs = nqs;
    res.idx_vs = idx_vs;
    res.nvs = nvs;
    res.idx_vExtendeds = idx_vExtendeds;
    res.nvExtendeds = nvExtendeds;
    // Eigen Vectors
    res.armature = armature.template cast<NewScalar>();
    res.damping = damping.template cast<NewScalar>();
    res.rotorInertia = rotorInertia.template cast<NewScalar>();
    res.rotorGearRatio = rotorGearRatio.template cast<NewScalar>();
    res.upperEffortLimit = upperEffortLimit.template cast<NewScalar>();
    res.lowerEffortLimit = lowerEffortLimit.template cast<NewScalar>();
    res.upperDryFrictionLimit = upperDryFrictionLimit.template cast<NewScalar>();
    res.lowerDryFrictionLimit = lowerDryFrictionLimit.template cast<NewScalar>();
    res.lowerVelocityLimit = lowerVelocityLimit.template cast<NewScalar>();
    res.upperVelocityLimit = upperVelocityLimit.template cast<NewScalar>();
    res.lowerPositionLimit = lowerPositionLimit.template cast<NewScalar>();
    res.upperPositionLimit = upperPositionLimit.template cast<NewScalar>();
    res.positionLimitMargin = positionLimitMargin.template cast<NewScalar>();

    typename ConfigVectorMap::const_iterator it;
    for (it = referenceConfigurations.begin(); it != referenceConfigurations.end(); it++)
    {
      res.referenceConfigurations.insert(
        std::make_pair(it->first, it->second.template cast<NewScalar>()));
    }

    // reserve vectors
    res.inertias.resize(inertias.size());
    res.jointPlacements.resize(jointPlacements.size());
    res.joints.resize(joints.size());

    // copy into vectors
    for (size_t k = 0; k < joints.size(); ++k)
    {
      res.inertias[k] = inertias[k].template cast<NewScalar>();
      res.jointPlacements[k] = jointPlacements[k].template cast<NewScalar>();
      res.joints[k] = joints[k].template cast<NewScalar>();
    }

    res.frames.resize(frames.size());
    for (size_t k = 0; k < frames.size(); ++k)
    {
      res.frames[k] = frames[k].template cast<NewScalar>();
    }

    return res;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  template<template<typename, int> class OtherJointCollectionTpl>
  ModelTpl<Scalar, Options, JointCollectionTpl> &
  ModelTpl<Scalar, Options, JointCollectionTpl>::operator=(
    const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other)
  {
    this->nq = other.nq;
    this->nv = other.nv;
    this->nvExtended = other.nvExtended;
    this->njoints = other.njoints;
    this->nbodies = other.nbodies;
    this->nframes = other.nframes;
    this->inertias = other.inertias;
    this->jointPlacements = other.jointPlacements;
    this->joints.clear();
    this->joints.reserve(other.joints.size());
    for (const auto & other_joint : other.joints)
    {
      this->joints.push_back(other_joint);
    }
    this->idx_qs = other.idx_qs;
    this->nqs = other.nqs;
    this->idx_vs = other.idx_vs;
    this->nvs = other.nvs;
    this->idx_vExtendeds = other.idx_vExtendeds;
    this->nvExtendeds = other.nvExtendeds;
    this->parents = other.parents;
    this->children = other.children;
    this->names = other.names;
    this->referenceConfigurations = other.referenceConfigurations;
    this->armature = other.armature;
    this->rotorInertia = other.rotorInertia;
    this->rotorGearRatio = other.rotorGearRatio;
    this->lowerDryFrictionLimit = other.lowerDryFrictionLimit;
    this->upperDryFrictionLimit = other.upperDryFrictionLimit;
    this->damping = other.damping;
    this->lowerEffortLimit = other.lowerEffortLimit;
    this->upperEffortLimit = other.upperEffortLimit;
    this->lowerVelocityLimit = other.lowerVelocityLimit;
    this->upperVelocityLimit = other.upperVelocityLimit;
    this->lowerPositionLimit = other.lowerPositionLimit;
    this->upperPositionLimit = other.upperPositionLimit;
    this->positionLimitMargin = other.positionLimitMargin;
    this->frames = other.frames;
    this->supports = other.supports;
    this->subtrees = other.subtrees;
    this->mimic_joint_supports = other.mimic_joint_supports;
    this->mimicking_joints = other.mimicking_joints;
    this->mimicked_joints = other.mimicked_joints;
    this->gravity = other.gravity;
    this->name = other.name;
    this->sparsity_pattern_vector = other.sparsity_pattern_vector;
    this->span_indexes_vector = other.span_indexes_vector;
    return *this;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool ModelTpl<Scalar, Options, JointCollectionTpl>::operator==(const ModelTpl & other) const
  {
    bool res = other.nq == nq && other.nv == nv && other.nvExtended == nvExtended
               && other.njoints == njoints && other.nbodies == nbodies && other.nframes == nframes
               && other.parents == parents && other.children == children && other.names == names
               && other.subtrees == subtrees && other.mimic_joint_supports == mimic_joint_supports
               && other.sparsity_pattern_vector == sparsity_pattern_vector
               && other.span_indexes_vector == span_indexes_vector
               && other.mimicking_joints == mimicking_joints
               && other.mimicked_joints == mimicked_joints && other.gravity == gravity
               && other.name == name;

    res &= other.idx_qs == idx_qs && other.nqs == nqs && other.idx_vs == idx_vs && other.nvs == nvs
           && other.idx_vExtendeds == idx_vExtendeds && other.nvExtendeds == nvExtendeds;

    if (other.referenceConfigurations.size() != referenceConfigurations.size())
      return false;

    typename ConfigVectorMap::const_iterator it = referenceConfigurations.begin();
    typename ConfigVectorMap::const_iterator it_other = other.referenceConfigurations.begin();
    for (long k = 0; k < (long)referenceConfigurations.size(); ++k)
    {
      if (it->second.size() != it_other->second.size())
        return false;
      if (it->second != it_other->second)
        return false;
      std::advance(it, 1);
      std::advance(it_other, 1);
    }
    if (other.armature.size() != armature.size())
      return false;
    res &= other.armature == armature;
    if (!res)
      return res;

    if (other.damping.size() != damping.size())
      return false;
    res &= other.damping == damping;
    if (!res)
      return res;

    if (other.rotorInertia.size() != rotorInertia.size())
      return false;
    res &= other.rotorInertia == rotorInertia;
    if (!res)
      return res;

    if (other.rotorGearRatio.size() != rotorGearRatio.size())
      return false;
    res &= other.rotorGearRatio == rotorGearRatio;
    if (!res)
      return res;

    if (other.lowerEffortLimit.size() != lowerEffortLimit.size())
      return false;
    res &= other.lowerEffortLimit == lowerEffortLimit;
    if (!res)
      return res;

    if (other.upperEffortLimit.size() != upperEffortLimit.size())
      return false;
    res &= other.upperEffortLimit == upperEffortLimit;
    if (!res)
      return res;

    if (other.lowerDryFrictionLimit.size() != lowerDryFrictionLimit.size())
      return false;
    res &= other.lowerDryFrictionLimit == lowerDryFrictionLimit;
    if (!res)
      return res;

    if (other.upperDryFrictionLimit.size() != upperDryFrictionLimit.size())
      return false;
    res &= other.upperDryFrictionLimit == upperDryFrictionLimit;
    if (!res)
      return res;

    if (other.lowerVelocityLimit.size() != lowerVelocityLimit.size())
      return false;
    res &= other.lowerVelocityLimit == lowerVelocityLimit;
    if (!res)
      return res;

    if (other.upperVelocityLimit.size() != upperVelocityLimit.size())
      return false;
    res &= other.upperVelocityLimit == upperVelocityLimit;
    if (!res)
      return res;

    if (other.lowerPositionLimit.size() != lowerPositionLimit.size())
      return false;
    res &= other.lowerPositionLimit == lowerPositionLimit;
    if (!res)
      return res;

    if (other.upperPositionLimit.size() != upperPositionLimit.size())
      return false;
    res &= other.upperPositionLimit == upperPositionLimit;

    if (other.positionLimitMargin.size() != positionLimitMargin.size())
      return false;
    res &= other.positionLimitMargin == positionLimitMargin;

    if (!res)
      return res;

    for (size_t k = 1; k < inertias.size(); ++k)
    {
      res &= other.inertias[k] == inertias[k];
      if (!res)
        return res;
    }

    for (size_t k = 1; k < other.jointPlacements.size(); ++k)
    {
      res &= other.jointPlacements[k] == jointPlacements[k];
      if (!res)
        return res;
    }

    res &= other.joints == joints && other.frames == frames;

    return res;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void ModelTpl<Scalar, Options, JointCollectionTpl>::appendBodyToJoint(
    const typename ModelTpl::JointIndex joint_index, const Inertia & Y, const SE3 & body_placement)
  {
    const Inertia & iYf = Y.se3Action(body_placement);
    inertias[joint_index] += iYf;
    nbodies++;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addBodyFrame(
    const std::string & body_name,
    const JointIndex & parentJoint,
    const SE3 & body_placement,
    int parentFrame)
  {
    if (parentFrame < 0)
    {
      // FIXED_JOINT is required because the parent can be the universe and its
      // type is FIXED_JOINT
      parentFrame = (int)getFrameId(names[parentJoint], (FrameType)(JOINT | FIXED_JOINT));
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT((size_t)parentFrame < frames.size(), "Frame index out of bound");
    return addFrame(Frame(body_name, parentJoint, (FrameIndex)parentFrame, body_placement, BODY));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getBodyId(const std::string & name) const
  {
    return getFrameId(name, BODY);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool
  ModelTpl<Scalar, Options, JointCollectionTpl>::existBodyName(const std::string & name) const
  {
    return existFrame(name, BODY);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getJointId(const std::string & name) const
  {
    typedef std::vector<std::string>::iterator::difference_type it_diff_t;
    it_diff_t res = std::find(names.begin(), names.end(), name) - names.begin();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (res < INT_MAX), "Id superior to int range. Should never happen.");
    return ModelTpl::JointIndex(res);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool
  ModelTpl<Scalar, Options, JointCollectionTpl>::existJointName(const std::string & name) const
  {
    return (names.end() != std::find(names.begin(), names.end(), name));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getFrameId(
    const std::string & name, const FrameType & type) const
  {
    typename std::vector<Frame>::const_iterator it =
      std::find_if(frames.begin(), frames.end(), details::FilterFrame(name, type));
    std::ostringstream os;
    os << "Several frames match the filter - please specify the FrameType (name=\"" << name
       << "\", type=\"" << type << "\")";
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      ((it == frames.end()
        || (std::find_if(boost::next(it), frames.end(), details::FilterFrame(name, type)) == frames.end()))),
      os.str().c_str());
    return FrameIndex(it - frames.begin());
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool ModelTpl<Scalar, Options, JointCollectionTpl>::existFrame(
    const std::string & name, const FrameType & type) const
  {
    return std::find_if(frames.begin(), frames.end(), details::FilterFrame(name, type))
           != frames.end();
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addFrame(
    const Frame & frame, const bool append_inertia)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      frame.parentJoint < (JointIndex)njoints, "The index of the parent joint is not valid.");

    //    TODO: fix it
    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(frame.inertia.isValid(),
    //                                   "The input inertia is not valid.")

    // Check if the frame.name exists with the same type
    if (existFrame(frame.name, frame.type))
    {
      return getFrameId(frame.name, frame.type);
    }
    // else: we must add a new frames to the current stack
    frames.push_back(frame);
    if (append_inertia)
      inertias[frame.parentJoint] += frame.placement.act(frame.inertia);
    nframes++;
    return FrameIndex(nframes - 1);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void ModelTpl<Scalar, Options, JointCollectionTpl>::addJointIndexToParentSubtrees(
    const JointIndex joint_id)
  {
    for (JointIndex parent = parents[joint_id]; parent > 0; parent = parents[parent])
      subtrees[parent].push_back(joint_id);

    // Also add joint_id to the universe
    subtrees[0].push_back(joint_id);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<bool> ModelTpl<Scalar, Options, JointCollectionTpl>::hasConfigurationLimit() const
  {
    std::vector<bool> vec;
    for (Index i = 1; i < (Index)(njoints); ++i)
    {
      const std::vector<bool> & cf_limits = joints[i].hasConfigurationLimit();
      vec.insert(vec.end(), cf_limits.begin(), cf_limits.end());
    }
    return vec;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<bool>
  ModelTpl<Scalar, Options, JointCollectionTpl>::hasConfigurationLimitInTangent() const
  {
    std::vector<bool> vec;
    for (Index i = 1; i < (Index)(njoints); ++i)
    {
      const std::vector<bool> & cf_limits = joints[i].hasConfigurationLimitInTangent();
      vec.insert(vec.end(), cf_limits.begin(), cf_limits.end());
    }
    return vec;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<JointIndex> ModelTpl<Scalar, Options, JointCollectionTpl>::getChildJoints() const
  {
    std::vector<JointIndex> res;
    for (JointIndex joint_id = 1; joint_id < JointIndex(njoints); ++joint_id)
    {
      if (this->children[joint_id].empty())
        res.push_back(joint_id);
    }
    return res;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool ModelTpl<Scalar, Options, JointCollectionTpl>::check() const
  {
    return this->check(makeDefaultCheckerList);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool ModelTpl<Scalar, Options, JointCollectionTpl>::check(const Data & data) const
  {
    return checkData(*this, data);
  }

} // namespace pinocchio

#ifdef PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

namespace pinocchio
{
  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::ModelTpl();

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI JointIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addJoint(
    const JointIndex, const JointModel &, const SE3 &, const std::string &);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI JointIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addJoint(
    const JointIndex,
    const JointModel &,
    const SE3 &,
    const std::string &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI JointIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addJoint(
    const JointIndex,
    const JointModel &,
    const SE3 &,
    const std::string &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &,
    const context::VectorXs &);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI FrameIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addJointFrame(
    const JointIndex, int);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI void
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::appendBodyToJoint(
    const JointIndex, const Inertia &, const SE3 &);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI FrameIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addBodyFrame(
    const std::string &, const JointIndex &, const SE3 &, int);

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI FrameIndex
  ModelTpl<context::Scalar, context::Options, JointCollectionDefaultTpl>::addFrame(
    const Frame &, const bool);

} // namespace pinocchio

#endif // ifdef PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
