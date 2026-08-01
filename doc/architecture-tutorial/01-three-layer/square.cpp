// ============================================================
// 【实例化层】对应 Pinocchio 的 src/algorithm/xxx.cpp
//
// 职责：对"默认标量"做一次显式实例化，编进 .so 库。
// 目的：让下游项目直接链接现成的 double 版本，不必每次重新实例化模板
//       —— 这正是 libpinocchio.so 存在的意义（编译加速）。
// 结论：读算法实现永远去 .hxx，不要来这里，这里没有逻辑。
// ============================================================
// 直接 include 实现层：编译这个"库"时必须能看到模板定义体，
// 否则无法生成机器码。真实仓库中，src/algorithm/kinematics.cpp
// 通过 include kinematics.hpp（其末尾又 include 了 .hxx）达到同样效果。
#include "square.hxx"

// 声明也要可见（真实项目由 .hpp 提供；这里单独写出来是为了让
// 实验一"注释掉 .hpp 里的 include"时，本文件依然能独立编译成功，
// 从而演示出"double 可用、float 失败"的对比）
template<typename Scalar>
Scalar square(Scalar x);

// 显式实例化：强制编译器在此生成 double 版本的机器码，编进 .o / .so
// 对照真实文件：src/algorithm/kinematics.cpp 里的
//   template ... void forwardKinematics<context::Scalar, ...>(...);
template double square<double>(double);
