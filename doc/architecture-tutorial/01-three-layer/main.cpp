// ============================================================
// 【使用层】对应你自己的项目
//
// 实验目的：亲眼看到"默认标量走预编译、其他标量走编译期实例化"的差异。
// 运行：make run-01
// ============================================================
#include "square.hpp"
#include <iostream>
#include <iomanip>

int main()
{
  std::cout << "=== 实验 1：三层文件布局 ===\n\n";

  // ① double：square.cpp 里已经显式实例化过
  //    即使删掉 square.hpp 末尾的 #include "square.hxx"，这行依然能链接成功
  //    —— 因为 square.o 里已有现成的 double 版本机器码
  std::cout << "square(3.0)   = " << square(3.0) << "   <- double：用 square.cpp 预实例化的版本\n";

  // ② float：square.cpp 没有实例化过它
  //    只能靠 square.hpp 末尾 include 进来的 .hxx 在【编译期】临时生成
  //    删掉那行 include，这里就会报 undefined reference
  std::cout << "square(3.0f)  = " << square(3.0f) << "   <- float：编译期临时实例化（依赖 .hxx 可见）\n";

  // ③ long double：同理，也是编译期临时生成
  std::cout << std::setprecision(20);
  std::cout << "square(3.0L)  = " << square(3.0L) << " <- long double：同上\n";

  std::cout << "\n【结论】\n"
            << "  .hpp 末尾 include .hxx  → 让任意标量都能在编译期实例化\n"
            << "  .cpp 显式实例化 double  → 让最常用的版本预编译进库，省编译时间\n"
            << "\n试试 README 里的\"实验一\"：注释掉 square.hpp 末尾那行 include，\n"
            << "观察 float/long double 如何链接失败，而 double 依然正常。\n";
  return 0;
}
