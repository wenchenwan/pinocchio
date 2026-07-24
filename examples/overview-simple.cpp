// ============================================================
// 最简示例：手工构建机械臂模型，计算重力补偿力矩
// 演示：
//   1. buildModels::manipulator → 内置 6-DOF 串联机械臂
//   2. Data(model) → 模型/数据分离（model 只读，data 可写）
//   3. neutral(model) → 零位配置（关节角全零）
//   4. rnea(model, data, q, v=0, a=0) → 纯重力补偿力矩
//
// RNEA（递归 Newton-Euler 算法）：
//   输入：q（配置）, v（速度）, a（加速度）
//   输出：τ = M(q)·a + C(q,v)·v + g(q)
//   当 v=0, a=0 时：τ = g(q)（重力项，即维持静止所需力矩）
// ============================================================

#include <iostream>

#include "pinocchio/multibody/sample-models.hpp"    // 内置示例模型
#include "pinocchio/algorithm/joint-configuration.hpp"  // neutral()
#include "pinocchio/algorithm/rnea.hpp"             // rnea()

int main()
{
  // 创建空模型，通过 manipulator() 填充：6 个 RZ 关节串联，类似 UR5
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  // Data：算法的可变工作空间（存储中间计算结果，如 oMi、f、v_J 等）
  // model 是只读的，可被多个 Data 共享（多线程安全）
  pinocchio::Data data(model);

  // neutral(model)：返回"零位"配置向量（nq 维）
  // 对于旋转关节（JointModelRZ）：零位 = 0 rad
  // 对于球形关节/浮动基（四元数表示）：零位 = 单位四元数 [0,0,0,1]
  Eigen::VectorXd q = pinocchio::neutral(model);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);  // 速度 = 0
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);  // 加速度 = 0

  // rnea(model, data, q, v, a) → 递归 Newton-Euler 逆动力学
  // 返回引用 data.tau（内部已计算并存储）
  // 在 v=0, a=0 下：tau = g(q)（重力补偿力矩向量）
  const Eigen::VectorXd & tau = pinocchio::rnea(model, data, q, v, a);
  std::cout << "tau = " << tau.transpose() << std::endl;
}
