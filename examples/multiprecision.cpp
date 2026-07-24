// ============================================================
// 多精度算术（Multiprecision Arithmetic）
// 演示：用 Boost.Multiprecision 实现 100 位十进制精度的 RNEA
//
// 应用场景：
//   - 数值验证：检验双精度算法的精度误差上界
//   - 解析导数的参考值：高精度有限差分（步长可以很小而不损失精度）
//   - 病态问题（ill-conditioned）：当机器人处于奇异配置附近，双精度可能精度不足
//   - 代码验证：对比 double 版本与高精度版本的差异
//
// Pinocchio 的模板设计：
//   ModelTpl<Scalar>：以 Scalar 为数值类型的模板化模型
//   DataTpl<Scalar>：以 Scalar 为数值类型的模板化数据
//   model.cast<float_100>()：将 double 模型转换为 100 位精度模型
//   所有算法（rnea, aba, crba 等）都支持任意精度数值类型
// ============================================================

#include "pinocchio/math/multiprecision.hpp"  // 支持 boost::multiprecision 的 Pinocchio 适配层

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"

// Boost.Multiprecision：任意精度浮点数支持
// cpp_dec_float<100>：100 位十进制精度（≈ 332 位二进制精度，约 100 位有效数字）
// et_off：禁用表达式模板（提高编译速度，简化调试，不影响精度）
#include <boost/multiprecision/cpp_dec_float.hpp>

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

  // ---- 加载标准双精度模型 ----
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);
  Data data(model);

  // ---- 定义 100 位精度数值类型 ----
  typedef ::boost::multiprecision::number<
    ::boost::multiprecision::cpp_dec_float<100>,  // 100 位十进制精度
    ::boost::multiprecision::et_off>              // 禁用表达式模板
    float_100;

  // 模板化的 Model 和 Data（以 float_100 为数值类型）
  typedef ModelTpl<float_100> ModelMulti;
  typedef DataTpl<float_100>  DataMulti;

  // ---- 将 double 模型转换为 100 位精度模型 ----
  // model.cast<float_100>()：逐元素将 double 转为 float_100（无精度损失，只是类型升级）
  ModelMulti model_multi = model.cast<float_100>();
  DataMulti  data_multi(model_multi);

  // ---- 生成测试数据（100 位精度）----
  // randomConfiguration 的多精度版本：在关节限位内均匀采样（高精度随机数）
  ModelMulti::ConfigVectorType  q_multi = randomConfiguration(model_multi);
  // Random(n)：生成 n 维高精度随机向量
  ModelMulti::TangentVectorType v_multi = ModelMulti::TangentVectorType::Random(model.nv);
  ModelMulti::TangentVectorType a_multi = ModelMulti::TangentVectorType::Random(model.nv);

  // ---- 转换为 double（用于对比）----
  // cast<double>()：将 float_100 向量截断回双精度
  Model::ConfigVectorType  q = q_multi.cast<double>();
  Model::TangentVectorType v = v_multi.cast<double>();
  Model::TangentVectorType a = a_multi.cast<double>();

  // ---- 双精度 RNEA ----
  rnea(model, data, q, v, a);

  // ---- 100 位精度 RNEA ----
  // 完全相同的算法，但以 float_100 为数值类型运行
  // 结果存储在 data_multi.tau（100 位精度的关节力矩向量）
  rnea(model_multi, data_multi, q_multi, v_multi, a_multi);

  // ---- 对比输出 ----
  // max_digits10：保证输出精度足以区分 float_100 的每个有效数字
  std::cout << "Joint torque standard arithmetic:\n"
            << std::setprecision(std::numeric_limits<float_100>::max_digits10) << data.tau
            << std::endl;
  std::cout << "Joint torque multiprecision arithmetic:\n"
            << std::setprecision(std::numeric_limits<float_100>::max_digits10) << data_multi.tau
            << std::endl;
  // 对比两个输出的差异，即可评估 double 精度下的舍入误差量级
}
