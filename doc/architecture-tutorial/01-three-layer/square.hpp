// ============================================================
// 【声明层】对应 Pinocchio 的 include/pinocchio/algorithm/xxx.hpp
//
// 职责：只写函数声明 + 文档注释。使用者 include 的就是这个文件。
// 关键：末尾必须把 .hxx（实现层）拉进来，否则非默认标量会链接失败。
// ============================================================
#pragma once

/// 计算 x 的平方。
/// 模板化的意义：同一份代码可用于 double / float / 自动微分标量。
template<typename Scalar>
Scalar square(Scalar x);

// ★ 关键的一行：把模板实现拉进来
//   删掉它，main.cpp 里的 square(3.0f) 就会链接报错（见 README 实验一）
//   对照真实文件：kinematics.hpp 末尾的
//   #include "pinocchio/src/algorithm/kinematics.hxx"
#include "square.hxx"
