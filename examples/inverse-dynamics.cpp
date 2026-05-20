// Copyright 2023 Inria
// SPDX-License-Identifier: BSD-2-Clause

// ============================================================
// 逆动力学（Inverse Dynamics）：RNEA 计算关节力矩
//
// 物理含义：给定运动学轨迹 (q, v, a)，计算驱动该运动所需的关节力矩 τ
// 应用场景：
//   - 前馈力矩计算（轨迹跟踪控制器的前馈项）
//   - 系统辨识（用已知运动反求惯量参数）
//   - 重力补偿（v=0, a=0 时 τ = g(q)）
//
// RNEA（Recursive Newton-Euler Algorithm）：
//   正向传播：从根到末端计算各连杆的速度和加速度
//   反向传播：从末端到根计算各连杆的关节力矩
//   时间复杂度：O(n)
//
// 完整动力学方程：τ = M(q)·a + C(q,v)·v + g(q)
//   RNEA 一次调用直接给出总力矩，不分别计算 M, C, g
// ============================================================

#include <iostream>

#include "pinocchio/algorithm/joint-configuration.hpp"  // randomConfiguration
#include "pinocchio/algorithm/rnea.hpp"                 // rnea
#include "pinocchio/parsers/urdf.hpp"                   // urdf::buildModel

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own
// directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  // Change to your own URDF file here, or give a path as command-line argument
  const std::string urdf_filename = (argc <= 1)
                                      ? PINOCCHIO_MODEL_DIR
                                          + std::string("/example-robot-data/robots/"
                                                        "ur_description/urdf/ur5_robot.urdf")
                                      : argv[1];

  // Load the URDF model
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);

  // Build a data frame associated with the model
  Data data(model);

  // Sample a random joint configuration, joint velocities and accelerations
  Eigen::VectorXd q = randomConfiguration(model);      // in rad for the UR5
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv); // in rad/s for the UR5（v=0：静止）
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv); // in rad/s² for the UR5（a=0：重力补偿）

  // ---- RNEA（逆动力学）----
  // 当 v=0, a=0 时：rnea 返回重力补偿力矩 g(q)
  // 结果同时存储在 data.tau（返回值是对 data.tau 的 const 引用）
  Eigen::VectorXd tau = pinocchio::rnea(model, data, q, v, a);

  // Print out to the vector of joint torques (in N.m)
  // data.tau 和 tau 指向同一内存（tau 是 data.tau 的引用拷贝）
  std::cout << "Joint torques: " << data.tau.transpose() << std::endl;
  return 0;
}
