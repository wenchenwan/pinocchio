// ============================================================
// 【实现层】对应 Pinocchio 的 include/pinocchio/src/algorithm/xxx.hxx
//
// 职责：模板函数的真正实现体（算法的循环/递推写在这里）。
// 注意路径：真实仓库里这一层在 include/pinocchio/【src】/algorithm/ 下，
//          很多人只在 algorithm/ 里翻找不到实现，就是漏了中间的 src/。
// ============================================================
#pragma once

// IWYU pragma: private, include "square.hpp"
// ↑ 这行注释是给工具看的：标注"我是私有实现，请 include square.hpp 而不是我"
//   Pinocchio 每个 .hxx 顶部都有同样的标注。

template<typename Scalar>
Scalar square(Scalar x)
{
  return x * x;
}
