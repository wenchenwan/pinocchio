// ============================================================
// 实验 3：CRTP + Boost.Fusion 访问者（本章最核心的机制）
//
// 三个递进的例子：
//   ① 只读访问者：遍历打印每个关节的类型/维度
//   ② 手写 forwardKinematics：复刻 Pinocchio 内部算法骨架
//   ③ 与官方 forwardKinematics 逐关节对比，验证结果完全一致
//
// 读完你会明白：所有 Pinocchio 算法都是"访问者 + 递推"这一个模式。
//
// 运行：make run-03
// ============================================================
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iostream>
#include <iomanip>

namespace bf = boost::fusion;
using namespace pinocchio;

// ============================================================
// ① 只读访问者：打印关节信息
//
// 三要素（所有 Pinocchio 算法都长这样）：
//   (a) 继承 JointUnaryVisitorBase<自己>          ← CRTP，把"自己"告诉基类
//   (b) typedef ... ArgsType                       ← 声明要透传的参数包
//   (c) static void algo(jmodel, ...)              ← 对每种关节类型各实例化一份
// ============================================================
struct PrintJointVisitor
: public fusion::JointUnaryVisitorBase<PrintJointVisitor>   // (a) CRTP
{
  typedef bf::vector<const Model &> ArgsType;               // (b) 参数包

  template<typename JointModel>                             // (c) 编译期特化
  static void algo(const JointModelBase<JointModel> & jmodel, const Model & model)
  {
    // 进到这里时，jmodel 已经是【具体类型】（如 JointModelRX），
    // 不是基类指针 —— 所以没有虚函数开销，且可以内联。
    std::cout << "   id=" << std::setw(2) << jmodel.id()
              << "  " << std::setw(14) << std::left << model.names[jmodel.id()]
              << "  " << std::setw(20) << jmodel.shortname()
              << "  nq=" << jmodel.nq()
              << "  nv=" << jmodel.nv()
              << "  idx_q=" << std::setw(2) << jmodel.idx_q()
              << "  idx_v=" << jmodel.idx_v()
              << std::endl;
  }
};

// ============================================================
// ② 手写 FK 访问者：复刻 kinematics.hxx 里的 ForwardKinematicZeroStep
//
// 对照第 4.1 节的公式：
//   ᵒM_i = ᵒM_{λ(i)} · (jointPlacement_i · M_J(q_i))
//                       └──── liMi ────┘
// ============================================================
struct MyForwardKinematicsStep
: public fusion::JointUnaryVisitorBase<MyForwardKinematicsStep>
{
  // 这趟需要 model、data、q 三个参数
  typedef bf::vector<const Model &, Data &, const Eigen::VectorXd &> ArgsType;

  template<typename JointModel>
  static void algo(const JointModelBase<JointModel> & jmodel,
                   JointDataBase<typename JointModel::JointDataDerived> & jdata,
                   const Model & model, Data & data, const Eigen::VectorXd & q)
  {
    const JointIndex i      = jmodel.id();
    const JointIndex parent = model.parents[i];

    // (1) 关节自身的运动学：把 q 里属于本关节的分量算成 SE3 变换
    //     每种关节的差异全部封装在这里（旋转关节填 cos/sin，浮动基读 7 维位姿…）
    jmodel.calc(jdata.derived(), q);

    // (2) 关节相对父连杆的完整变换 = 固定安装位姿 × 关节自由度产生的变换
    data.liMi[i] = model.jointPlacements[i] * jdata.M();

    // (3) 递推到世界系（父已算好，因为我们按 id 从小到大遍历）
    if (parent > 0)
      data.oMi[i] = data.oMi[parent] * data.liMi[i];
    else
      data.oMi[i] = data.liMi[i];
  }
};

int main()
{
  std::cout << "=== 实验 3：CRTP + 访问者 ===\n\n";

  Model model;
  buildModels::humanoidRandom(model);          // 带浮动基的人形
  Data data_mine(model), data_ref(model);

  // ---------------------------------------------------------
  // ① 遍历打印关节信息
  // ---------------------------------------------------------
  std::cout << "① 用只读访问者遍历模型（前 8 个关节）：\n";
  for (JointIndex i = 1; i < (JointIndex)std::min<int>(9, model.njoints); ++i)
    PrintJointVisitor::run(model.joints[i], PrintJointVisitor::ArgsType(model));

  std::cout << "\n   注意第 1 个是 JointModelFreeFlyer：nq=7 而 nv=6\n"
            << "   —— 这正是浮动基\"四元数 4 维、角速度 3 维\"造成的 nq≠nv\n"
            << "   模型总计：nq=" << model.nq << ", nv=" << model.nv
            << ", njoints=" << model.njoints << "\n\n";

  // ---------------------------------------------------------
  // ② 用【自己写的】访问者跑一遍 FK
  // ---------------------------------------------------------
  Eigen::VectorXd q = neutral(model);
  for (int i = 7; i < q.size(); ++i)           // 跳过浮动基的 7 维，只动关节角
    q[i] = 0.2 + 0.05 * i;

  std::cout << "② 用手写访问者执行 forwardKinematics ...\n";
  for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    MyForwardKinematicsStep::run(model.joints[i], data_mine.joints[i],
                                 MyForwardKinematicsStep::ArgsType(model, data_mine, q));

  // ---------------------------------------------------------
  // ③ 与官方实现对比
  // ---------------------------------------------------------
  forwardKinematics(model, data_ref, q);       // Pinocchio 官方版本

  double max_err = 0.0;
  for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
  {
    // 用 SE3 的对数映射度量两个位姿的差异（见第 2.1 节）
    const double e = log6(data_mine.oMi[i].actInv(data_ref.oMi[i])).toVector().norm();
    max_err = std::max(max_err, e);
  }

  std::cout << "③ 与官方 forwardKinematics 对比：\n"
            << "   遍历全部 " << model.njoints - 1 << " 个关节，"
            << "max‖log6(ᵒM_mine⁻¹ · ᵒM_ref)‖ = "
            << std::scientific << std::setprecision(3) << max_err << "\n"
            << "   → " << (max_err < 1e-12 ? "完全一致 ✓" : "存在差异 ✗") << "\n\n";

  std::cout << "【结论】\n"
            << "  你刚刚复刻了 Pinocchio 内部算法的完整骨架。所有内置算法都是这个模式：\n"
            << "    · 把 algo() 里换成\"算受力并反向传递\" → 就是 RNEA 的一个 Step\n"
            << "    · RNEA 有 2 个 Step（对应 4.3 节两趟递推）\n"
            << "    · ABA  有 3 个 Step（对应 4.4 节三趟递推）\n"
            << "  读任何 .hxx 的套路：数 struct ...Step 的个数 = 算法有几趟。\n";

  return 0;
}
