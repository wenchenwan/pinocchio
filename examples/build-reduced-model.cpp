// ============================================================
// 降阶模型（Reduced Model）：锁定指定关节，构建自由度更少的模型
// 应用场景：
//   - 固定手臂，只研究腿部动力学（减少计算量）
//   - 固定上半身，研究下肢平衡控制
//   - 将复杂机器人拆分为子系统分别优化
//
// 两种构建方式：
//   1. 指定"要锁定的关节"列表（直接指定锁定哪些）
//   2. 指定"要保留的关节"列表（间接：锁定所有其他关节）
// ============================================================

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/model.hpp"  // buildReducedModel

#include <iostream>
#include <algorithm>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

// 辅助函数：检查元素是否在 vector 中（C++11 无 contains）
template<typename T>
bool is_in_vector(const std::vector<T> & vector, const T & elt)
{
  return vector.end() != std::find(vector.begin(), vector.end(), elt);
}

int main(int argc, char ** argv)
{
  using namespace pinocchio;

  const std::string urdf_filename =
    (argc <= 1) ? PINOCCHIO_MODEL_DIR
                    + std::string("/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf")
                : argv[1];

  // Load the urdf model（固定基 UR5，6 DOF）
  Model model;
  pinocchio::urdf::buildModel(urdf_filename, model);

  // ====================================================
  // 方式一：指定"要锁定的关节"列表
  // ====================================================
  std::vector<std::string> list_of_joints_to_lock_by_name;
  list_of_joints_to_lock_by_name.push_back("elbow_joint");
  list_of_joints_to_lock_by_name.push_back("wrist_3_joint"); // 顺序无关（内部会排序）
  list_of_joints_to_lock_by_name.push_back("wrist_2_joint");
  list_of_joints_to_lock_by_name.push_back("blabla"); // 不存在的关节名（会被跳过）

  // 将关节名转为 ID（跳过不存在的名称）
  std::vector<JointIndex> list_of_joints_to_lock_by_id;
  for (std::vector<std::string>::const_iterator it = list_of_joints_to_lock_by_name.begin();
       it != list_of_joints_to_lock_by_name.end(); ++it)
  {
    const std::string & joint_name = *it;
    if (model.existJointName(joint_name))  // 验证关节名是否存在
      list_of_joints_to_lock_by_id.push_back(model.getJointId(joint_name));
    else
      std::cout << "joint: " << joint_name << " does not belong to the model" << std::endl;
  }

  // q_rand：锁定关节将被固定在该配置对应的角度
  // 被锁定关节的运动学参数被"烘焙"到模型几何结构中（变为固定偏置）
  Eigen::VectorXd q_rand    = randomConfiguration(model);
  Eigen::VectorXd q_neutral = neutral(model);
  PINOCCHIO_UNUSED_VARIABLE(q_neutral);  // 仅用于演示 neutral() 用法

  std::cout << "\n\nFIRST CASE: BUILD A REDUCED MODEL FROM A LIST OF JOINT TO LOCK" << std::endl;

  // buildReducedModel(model, joints_to_lock, reference_q)
  // 返回降阶模型（njoints 减少 = 锁定关节数，nq/nv 相应减小）
  // reference_q：被锁定关节固定在此配置
  Model reduced_model = buildReducedModel(model, list_of_joints_to_lock_by_id, q_rand);

  // Print the list of joints in the original model
  std::cout << "List of joints in the original model:" << std::endl;
  for (JointIndex joint_id = 1; joint_id < model.joints.size(); ++joint_id)
    std::cout << "\t- " << model.names[joint_id] << std::endl;

  // Print the list of joints in the reduced model
  std::cout << "List of joints in the reduced model:" << std::endl;
  for (JointIndex joint_id = 1; joint_id < reduced_model.joints.size(); ++joint_id)
    std::cout << "\t- " << reduced_model.names[joint_id] << std::endl;

  // ====================================================
  // 方式二：指定"要保留的关节"列表（等价但更直观）
  // 适用场景：已知要保留哪些关节（如只保留腿部关节）
  // ====================================================
  std::cout << "\n\nSECOND CASE: BUILD A REDUCED MODEL FROM A LIST OF JOINT TO KEEP UNLOCKED"
            << std::endl;

  std::vector<std::string> list_of_joints_to_keep_unlocked_by_name;
  list_of_joints_to_keep_unlocked_by_name.push_back("shoulder_pan_joint");
  list_of_joints_to_keep_unlocked_by_name.push_back("shoulder_lift_joint");
  list_of_joints_to_keep_unlocked_by_name.push_back("wrist_1_joint");

  std::vector<JointIndex> list_of_joints_to_keep_unlocked_by_id;
  for (std::vector<std::string>::const_iterator it =
         list_of_joints_to_keep_unlocked_by_name.begin();
       it != list_of_joints_to_keep_unlocked_by_name.end(); ++it)
  {
    const std::string & joint_name = *it;
    if (model.existJointName(joint_name))
      list_of_joints_to_keep_unlocked_by_id.push_back(model.getJointId(joint_name));
    else
      std::cout << "joint: " << joint_name << " does not belong to the model";
  }

  // 转换：将"保留列表"取补集得到"锁定列表"
  // 遍历所有关节（从 1 开始，0=宇宙关节），不在保留列表中的即锁定
  list_of_joints_to_lock_by_id.clear();
  for (JointIndex joint_id = 1; joint_id < model.joints.size(); ++joint_id)
  {
    const std::string joint_name = model.names[joint_id];
    if (is_in_vector(list_of_joints_to_keep_unlocked_by_name, joint_name))
      continue;   // 保留的关节，跳过
    else
      list_of_joints_to_lock_by_id.push_back(joint_id);  // 其余的锁定
  }

  // 构建第二个降阶模型（只保留 shoulder_pan, shoulder_lift, wrist_1）
  Model reduced_model2 = buildReducedModel(model, list_of_joints_to_lock_by_id, q_rand);

  // Print the list of joints in the second reduced model
  std::cout << "List of joints in the second reduced model:" << std::endl;
  for (JointIndex joint_id = 1; joint_id < reduced_model2.joints.size(); ++joint_id)
    std::cout << "\t- " << reduced_model2.names[joint_id] << std::endl;
}
