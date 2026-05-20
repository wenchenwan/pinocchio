// ============================================================
// 碰撞检测（C++ 版）：使用 URDF + SRDF 的完整碰撞检测流程
// 演示：
//   1. urdf::buildGeom → 加载碰撞几何模型
//   2. addAllCollisionPairs + srdf::removeCollisionPairs → 碰撞对管理
//   3. srdf::loadReferenceConfigurations → 加载 half_sitting 等命名配置
//   4. computeCollisions → 批量碰撞检测
//   5. computeCollision → 单对碰撞检测
//   6. updateGeometryPlacements → 更新几何体位姿（正向运动学后）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"   // urdf::buildModel, urdf::buildGeom
#include "pinocchio/parsers/srdf.hpp"   // srdf::removeCollisionPairs, srdf::loadReferenceConfigurations

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/geometry.hpp"     // updateGeometryPlacements
#include "pinocchio/collision/collision.hpp"    // computeCollisions, computeCollision

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own modeldirectory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int /*argc*/, char ** /*argv*/)
{
  using namespace pinocchio;
  const std::string robots_model_path = PINOCCHIO_MODEL_DIR;

  // Talos 降阶版：PAL Robotics 人形机器人（去掉手部精细关节后的版本）
  const std::string urdf_filename =
    robots_model_path
    + std::string("/example-robot-data/robots/talos_data/robots/talos_reduced.urdf");
  // SRDF：语义机器人描述格式，包含自碰撞排除对和命名配置
  const std::string srdf_filename =
    robots_model_path + std::string("/example-robot-data/robots/talos_data/srdf/talos.srdf");

  // Load the URDF model contained in urdf_filename
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);  // 固定基（Talos URDF 默认无浮动基）

  // Build the data associated to the model
  Data data(model);

  // ---- 加载碰撞几何模型 ----
  // urdf::buildGeom：解析 URDF 中 <collision> 标签下的几何体
  // COLLISION：只加载碰撞几何（区别于 VISUAL）
  // robots_model_path：网格文件（STL/DAE）的搜索根目录
  GeometryModel geom_model;
  pinocchio::urdf::buildGeom(
    model, urdf_filename, pinocchio::COLLISION, geom_model, robots_model_path);

  // ---- 碰撞对管理 ----
  // 先添加所有可能的碰撞对（N*(N-1)/2 个）
  geom_model.addAllCollisionPairs();
  // 然后从 SRDF 移除"永不碰撞"的对（相邻连杆、对称关节等）
  // 通常能去掉 70~90% 的碰撞对，大幅减少检测计算量
  pinocchio::srdf::removeCollisionPairs(model, geom_model, srdf_filename);

  // Build the data associated to the geom_model
  // GeometryData 存储：各几何体在世界坐标系中的位姿 + 碰撞检测结果
  GeometryData geom_data(geom_model);

  // ---- 加载 SRDF 中的命名参考配置 ----
  // loadReferenceConfigurations：解析 SRDF 中 <group_state> 标签
  // 将命名配置存入 model.referenceConfigurations 映射
  pinocchio::srdf::loadReferenceConfigurations(
    model,
    srdf_filename);  // the reference configuration stored in the SRDF file is called half_sitting

  // half_sitting：标准半蹲站立配置（膝关节略弯，无碰撞的标准姿态）
  const Model::ConfigVectorType & q = model.referenceConfigurations["half_sitting"];

  // ---- 批量碰撞检测（所有碰撞对）----
  // computeCollisions 内部：forwardKinematics + updateGeometryPlacements + FCL 查询
  // 默认不提前终止（stopAtFirstCollision = false）
  computeCollisions(model, data, geom_model, geom_data, q);

  // Print the status of all the collision pairs
  for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
  {
    const CollisionPair & cp          = geom_model.collisionPairs[k];
    const hpp::fcl::CollisionResult & cr = geom_data.collisionResults[k];

    std::cout << "collision pair: " << cp.first << " , " << cp.second << " - collision: ";
    std::cout << (cr.isCollision() ? "yes" : "no") << std::endl;
  }

  // ---- 提前终止模式：发现第一个碰撞立即返回 ----
  // computeCollisions(..., true) → stopAtFirstCollision=true
  // 适合安全检查（只需知道"是否有碰撞"，不需要知道"哪些碰撞"）
  computeCollisions(model, data, geom_model, geom_data, q, true);

  // ---- 单对碰撞检测 ----
  // 先更新所有几何体位姿（forwardKinematics + 几何传播）
  const PairIndex pair_id = 2;   // 检测第三个碰撞对（0-indexed）
  const Model::ConfigVectorType q_neutral = neutral(model);
  // updateGeometryPlacements：等价于 forwardKinematics + 逐几何体更新
  // 在手动控制碰撞检测流程时使用（与 computeCollisions 的区别）
  updateGeometryPlacements(
    model, data, geom_model, geom_data,
    q_neutral);  // performs a forward kinematics over the whole kinematics model + update the
                 // placement of all the geometries contained inside geom_model

  // 只检测第 pair_id 个碰撞对（比 computeCollisions 更高效，当只关心特定对时）
  computeCollision(geom_model, geom_data, pair_id);

  return 0;
}
