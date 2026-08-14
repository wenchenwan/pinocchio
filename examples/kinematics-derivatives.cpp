// ============================================================
// 运动学解析导数（Kinematics Derivatives）
// 计算：∂v_joint/∂q, ∂a_joint/∂q, ∂a_joint/∂v, ∂a_joint/∂a ∈ R^{6×nv}
//
// 物理含义：
//   v_partial_dq = ∂v_joint/∂q（速度对配置的偏导 = 几何 Jacobian）
//   a_partial_dq = ∂a_joint/∂q（加速度对配置的偏导，包含离心/Coriolis 项）
//   a_partial_dv = ∂a_joint/∂v（加速度对速度的偏导 = 速度相关的 Jacobian 导数）
//   a_partial_da = ∂a_joint/∂a（加速度对加速度的偏导 = 几何 Jacobian，与 v_partial_dv 相等）
//
// 注意：v_partial_dv = a_partial_da（两者相等，均等于几何 Jacobian）
//   因此 getJointAccelerationDerivatives 不单独返回 v_partial_dv
//
// 参考系选项：
//   LOCAL：在关节局部坐标系中表达（原点=关节原点，轴=关节轴）
//   WORLD：在世界坐标系对齐的坐标系中表达（原点=关节原点，轴=世界轴）
//   LOCAL_WORLD_ALIGNED：与 WORLD 相同（等价）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"  // computeForwardKinematicsDerivatives, getJointAccelerationDerivatives

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  const std::string urdf_filename =
    (argc <= 1) ? PINOCCHIO_MODEL_DIR
                    + std::string("/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf")
                : argv[1];

  // Load the URDF model
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);

  // Build a data related to model
  Data data(model);

  // Sample a random joint configuration as well as random joint velocity and acceleration
  Eigen::VectorXd q = randomConfiguration(model);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

  // ---- computeForwardKinematicsDerivatives：预计算所有关节的运动学导数 ----
  // 等价于先调用 forwardKinematics，再存储用于后续 get 调用的中间量
  // 必须在 getJointVelocityDerivatives 或 getJointAccelerationDerivatives 之前调用
  computeForwardKinematicsDerivatives(model, data, q, v, a);

  // ---- 获取末端关节（最后一个关节）的运动学导数 ----
  // 从索引 1 开始（0 = 宇宙关节），最后一个为末端执行器
  JointIndex joint_id = (JointIndex)(model.njoints - 1);

  // ============================================================
  // 四个结果矩阵的含义（均为 6×nv）
  //
  // 出发点是这两条运动学关系：
  //     v = J(q)·q̇                        （末端 6D 空间速度）
  //     a = J(q)·q̈ + J̇(q,q̇)·q̇            （末端 6D 空间加速度）
  // 四个矩阵就是 v、a 对状态量 (q, q̇, q̈) 的偏导：
  //
  //   ┌────────────────┬──────────────┬────────────────────────────────┐
  //   │ 变量           │ 数学         │ 直观理解                        │
  //   ├────────────────┼──────────────┼────────────────────────────────┤
  //   │ v_partial_dq   │ ∂v/∂q        │ 姿态变化 → 末端速度变化         │
  //   │ a_partial_dq   │ ∂a/∂q        │ 姿态变化 → 末端加速度变化       │
  //   │ a_partial_dv   │ ∂a/∂q̇        │ 关节速度变化 → 末端加速度变化   │
  //   │ a_partial_da   │ ∂a/∂q̈  = J   │ 关节加速度变化 → 末端加速度变化 │
  //   └────────────────┴──────────────┴────────────────────────────────┘
  //
  // ★ 最重要的一条：a_partial_da 恒等于【几何 Jacobian】
  //   因为 a = J·q̈ + J̇·q̇ 中第二项与 q̈ 无关，对 q̈ 求导只剩 J。
  //   同理 ∂v/∂q̇ = J（因 v = J·q̇）—— 所以这两个量【是同一个矩阵】，
  //   这正是本函数不单独返回 v_partial_dv 的原因（见下方原注释）。
  //   UR5 有 6 个自由度，故 J ∈ R^{6×6}，其第 i 列表示
  //   "只让第 i 个关节加速时，末端 6D 空间加速度的变化"。
  //
  // ⚠️ v_partial_dq【不是】Jacobian —— 这是最容易搞混的一点。
  //   两连杆平面臂实测（q=(0.3,0.5)，LOCAL 系）：
  //     v_partial_dq 第1列 = (0, 0, 0, 0, 0, 0)
  //                  第2列 = (1.755, -0.959, 0, 0, 0, 0)
  //     J            第1列 = (0.479, 0.878, 0, 0, 0, 1)
  //                  第2列 = (0, 0, 0, 0, 0, 1)
  //   完全不同。第 1 列全零的物理含义：在 LOCAL 系下，转动第一个关节
  //   相当于整条臂绕根部旋转，"末端看自己"的速度不变。
  //
  // a_partial_dq / a_partial_dv 刻画的是 J̇·q̇ 这一项的敏感度，
  // 即科氏力、向心加速度等速度相关效应的来源。
  //
  // ⚠️ 必须先 setZero()：这些算法沿运动树【累加】写入各关节的贡献，
  //   不是覆盖写入。不清零 = 往垃圾值上累加，且不会报错。
  // ============================================================
  Data::Matrix6x v_partial_dq(6, model.nv), a_partial_dq(6, model.nv),
    a_partial_dv(6, model.nv), a_partial_da(6, model.nv);
  v_partial_dq.setZero();
  a_partial_dq.setZero();
  a_partial_dv.setZero();
  a_partial_da.setZero();

  // ---- LOCAL 参考系：导数在关节局部坐标系中表达 ----
  // LOCAL：原点=关节原点，轴=关节轴（随关节旋转而变化）
  getJointAccelerationDerivatives(
    model, data, joint_id, LOCAL, v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da);

  // Remark: we are not directly computing the quantity v_partial_dv as it is also equal to
  // a_partial_da.
  // 上面这句原注释即前述要点：∂v/∂q̇ 与 ∂a/∂q̈ 都等于 J，返回两遍是浪费。

  // ---- WORLD 参考系：换个参考系再算一遍 ----
  // 用法相同，只改 reference_frame 参数。但要注意三个参考系的区别：
  //   LOCAL              原点=关节原点，轴=关节轴（随关节转）
  //   WORLD              原点=【世界原点】，轴=世界轴
  //   LOCAL_WORLD_ALIGNED 原点=关节原点，轴=世界轴
  //
  // ⚠️ 末端执行器任务通常要用 LOCAL_WORLD_ALIGNED 而非 WORLD：
  //   WORLD 的线性分量是"世界原点处虚拟点"的量，含 p×ω 牵连项；
  //   LOCAL_WORLD_ALIGNED 的线性分量才是【末端那一点】的真实量。
  //   （详见 PINOCCHIO_GUIDE.md §2.2.7）
  //
  // 换参考系后四个矩阵的数值全变，但上面 a_partial_da == J 的恒等式
  // 在每个参考系内部依然成立（J 也取同一参考系）。
  getJointAccelerationDerivatives(
    model, data, joint_id, WORLD, v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da);

  // 可选：打印几何 Jacobian —— 注意是 a_partial_da 而非 v_partial_dq
  // std::cout << "Geometric Jacobian:\n" << a_partial_da << std::endl;
}
