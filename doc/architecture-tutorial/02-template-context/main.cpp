// ============================================================
// 实验 2：Tpl 模板 + context 默认标量 + cast<>()
//
// 验证三件事：
//   ① pinocchio::Model 只是 ModelTpl<double,0> 的别名（编译期断言）
//   ② 同一份算法代码换个标量就能跑出不同精度
//   ③ cast<> 是"逐字段换标量重造模型"，不是简单的类型转换
//
// 运行：make run-02
// ============================================================
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <type_traits>

int main()
{
  using namespace pinocchio;

  std::cout << "=== 实验 2：模板 + context 默认标量 ===\n\n";

  // ---------------------------------------------------------
  // ① 编译期证明：Model 就是 ModelTpl<context::Scalar, context::Options>
  //    static_assert 在【编译期】检查，能编过就说明别名等价成立
  // ---------------------------------------------------------
  static_assert(std::is_same<Model, ModelTpl<context::Scalar, context::Options>>::value,
                "Model 必须等于 ModelTpl<context::Scalar, context::Options>");
  static_assert(std::is_same<context::Scalar, double>::value,
                "默认 context 标量应为 double");

  std::cout << "① 别名等价（编译期 static_assert 已通过）：\n"
            << "   pinocchio::Model    == ModelTpl<context::Scalar, context::Options>\n"
            << "   context::Scalar     == double\n"
            << "   context::Options    == " << context::Options << "\n\n";

  // ---------------------------------------------------------
  // ② 构建默认（double）模型，跑一次 RNEA
  // ---------------------------------------------------------
  Model model;
  buildModels::manipulator(model);
  Data data(model);

  // 用一个非零配置，否则零位下重力力矩恰好为 0，看不出差异
  Eigen::VectorXd q = neutral(model);
  for (int i = 0; i < q.size(); ++i)
    q[i] = 0.3 + 0.1 * i;                       // 一个确定性的"歪着"的姿态
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

  rnea(model, data, q, v, a);
  std::cout << "② double 版 RNEA（v=a=0，即重力补偿力矩 g(q)）：\n"
            << "   tau = " << std::fixed << std::setprecision(12)
            << data.tau.transpose() << "\n\n";

  // ---------------------------------------------------------
  // ③ cast<long double>()：换标量重造整个模型
  //    注意 ModelTpl 的 cast 会逐字段把 double 转成新标量
  //    —— 惯量、关节位姿、限位……全部重建一遍
  // ---------------------------------------------------------
  typedef ModelTpl<long double> ModelLD;
  typedef DataTpl<long double>  DataLD;

  ModelLD model_ld = model.cast<long double>();
  DataLD  data_ld(model_ld);

  ModelLD::ConfigVectorType  q_ld = q.cast<long double>();  // 配置也要换标量
  ModelLD::TangentVectorType v_ld = ModelLD::TangentVectorType::Zero(model_ld.nv);
  ModelLD::TangentVectorType a_ld = ModelLD::TangentVectorType::Zero(model_ld.nv);

  // ★ 完全同一个函数模板 rnea()，只是标量类型不同
  rnea(model_ld, data_ld, q_ld, v_ld, a_ld);

  std::cout << "③ long double 版 RNEA（同一份 rnea() 源码，换了标量）：\n"
            << "   tau = " << std::setprecision(12)
            << data_ld.tau.transpose() << "\n\n";

  // ---------------------------------------------------------
  // ④ 精度对比：打印两者差异（体现"换标量"的实际意义）
  // ---------------------------------------------------------
  Eigen::VectorXd tau_ld_as_double = data_ld.tau.template cast<double>();
  double max_diff = (data.tau - tau_ld_as_double).cwiseAbs().maxCoeff();

  std::cout << "④ 两种精度的结果差异：\n"
            << "   max|tau_double - tau_longdouble| = "
            << std::scientific << std::setprecision(3) << max_diff << "\n"
            << "   （量级约为 double 的机器精度 " << std::numeric_limits<double>::epsilon()
            << "，说明 double 结果本身是可靠的）\n\n";

  std::cout << "【结论】\n"
            << "  Model/Data/SE3 等都是 XxxTpl<Scalar,...> 的 typedef 别名；\n"
            << "  算法写成模板 → 换标量即可获得：高精度验证、自动微分（CppAD/CasADi）、\n"
            << "  区间分析等能力，而【算法源码一行都不用改】。\n"
            << "  这就是 Pinocchio 能提供解析梯度的底层前提。\n";

  return 0;
}
