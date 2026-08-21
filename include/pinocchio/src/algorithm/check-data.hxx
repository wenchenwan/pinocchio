//
// Copyright (c) 2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/check-data.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/check-data.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  // ============================================================
  // checkData：校验一个 Data 是否与 Model 结构匹配（是 model.check(data) 的落地实现）。
  // 逐条断言 Data 的所有缓存数组/矩阵维度、以及派生索引表都对得上 Model 的
  // njoints / nq / nv / nvExtended / nframes。全部通过返回 true，任一不对立即返回 false。
  // 主要用途：算法运行前确认"这个 data 确实是这个 model 用 createData() 造出来的"，
  //           抓住"传了旧的/别的模型的 Data"这类隐蔽 bug（对应不变量：Model 改了必须重建 Data）。
  // 外层 DIAGNOSTIC_PUSH/IGNORED_DEPRECATED/POP：局部压制访问已弃用字段时的告警。
  // ============================================================
  PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
  PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool checkData(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef typename Model::JointModel JointModel;

// 局部宏：条件不成立就立刻 return false（"断言即返回"的短路写法，全过才到末尾 return true）。
// 用局部 #define + 末尾 #undef，省去一长串重复的 if(!x) return false;，且不外泄。
#define CHECK_DATA(a)                                                                              \
  if (!(a))                                                                                        \
    return false;

    // TODO JMinvJt,sDUiJt are never explicitly initialized.
    // TODO impulse_c
    // They are not check neither

    // —— 组①：逐关节数组，长度必须 == njoints（每个关节一格） ——
    CHECK_DATA((int)data.joints.size() == model.njoints);
    CHECK_DATA((int)data.a.size() == model.njoints);
    CHECK_DATA((int)data.a_gf.size() == model.njoints);
    CHECK_DATA((int)data.v.size() == model.njoints);
    CHECK_DATA((int)data.f.size() == model.njoints);
    CHECK_DATA((int)data.oMi.size() == model.njoints);
    CHECK_DATA((int)data.liMi.size() == model.njoints);
    CHECK_DATA((int)data.Ycrb.size() == model.njoints);
    CHECK_DATA((int)data.Yaba.size() == model.njoints);
    CHECK_DATA((int)data.Fcrb.size() == model.njoints);
    for (const typename Data::Matrix6x & F : data.Fcrb)
    {
      CHECK_DATA(F.cols() == model.nv);
    }
    CHECK_DATA((int)data.iMf.size() == model.njoints);
    CHECK_DATA((int)data.iMf.size() == model.njoints);
    CHECK_DATA((int)data.com.size() == model.njoints);
    CHECK_DATA((int)data.vcom.size() == model.njoints);
    CHECK_DATA((int)data.acom.size() == model.njoints);
    CHECK_DATA((int)data.mass.size() == model.njoints);

    // —— 组②：切空间维度的向量/矩阵。注意 J 的列数按 nvExtended（雅可比按扩展维度建），其余按 nv ——
    CHECK_DATA(data.tau.size() == model.nv);
    CHECK_DATA(data.nle.size() == model.nv);
    CHECK_DATA(data.ddq.size() == model.nv);
    CHECK_DATA(data.u.size() == model.nv);
    CHECK_DATA(data.M.rows() == model.nv);
    CHECK_DATA(data.M.cols() == model.nv);
    CHECK_DATA(data.Ag.cols() == model.nv);
    CHECK_DATA(data.U.cols() == model.nv);
    CHECK_DATA(data.U.rows() == model.nv);
    CHECK_DATA(data.D.size() == model.nv);
    CHECK_DATA(data.tmp.size() >= model.nv);
    CHECK_DATA(data.J.cols() == model.nvExtended);
    CHECK_DATA(data.Jcom.cols() == model.nv);
    CHECK_DATA(data.torque_residual.size() == model.nv);
    CHECK_DATA(data.dq_after.size() == model.nv);
    // CHECK_DATA( data.impulse_c.size()== model.nv );

    CHECK_DATA(data.kinematic_hessians.dimension(0) == 6);
    CHECK_DATA(data.kinematic_hessians.dimension(1) == model.nv);
    CHECK_DATA(data.kinematic_hessians.dimension(2) == model.nv);

    // —— 组③：坐标系数 —— oMf 每个 Frame 一格
    CHECK_DATA((int)data.oMf.size() == model.nframes);

    // —— 组④：fromRow 索引表（稀疏 Cholesky 按自由度重排的树），长度按 nvExtended ——
    CHECK_DATA((int)data.lastChild.size() == model.njoints);
    CHECK_DATA((int)data.nvSubtree.size() == model.njoints);
    CHECK_DATA((int)data.parents_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.mimic_parents_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.non_mimic_parents_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.idx_vExtended_to_idx_v_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.nvSubtree_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.start_idx_v_fromRow.size() == model.nvExtended);
    CHECK_DATA((int)data.end_idx_v_fromRow.size() == model.nvExtended);

    // —— 组⑤：Model 并行索引表 ↔ 关节对象自报，必须一致（最能抓 bug 的一组） ——
    // model.nqs/idx_qs/... 是"同一信息在 Model 数组里的副本"，jmodel.nq()/idx_q() 是关节自存的，
    // 二者若不同步（如手工改了 Model 却没改关节）会让下游静默算错，这里逐关节核对。
    for (JointIndex joint_id = 1; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
      const typename Model::JointModel & jmodel = model.joints[joint_id];

      CHECK_DATA(model.nqs[joint_id] == jmodel.nq());
      CHECK_DATA(model.idx_qs[joint_id] == jmodel.idx_q());
      CHECK_DATA(model.nvs[joint_id] == jmodel.nv());
      CHECK_DATA(model.idx_vs[joint_id] == jmodel.idx_v());
      CHECK_DATA(model.nvExtendeds[joint_id] == jmodel.nvExtended());
      CHECK_DATA(model.idx_vExtendeds[joint_id] == jmodel.idx_vExtended());
    }
    // —— 组⑥：树/子树一致性 ——
    for (JointIndex j = 1; int(j) < model.njoints; ++j)
    {
      // c = 子树里 id 最大的后代；Pinocchio 保证"子树 = 连续区间 [j, c]"
      const JointIndex c = model.subtrees[j].back();
      CHECK_DATA((int)c < model.njoints);

      // 累加 j 及其所有后代 [j+1, c] 的 nv，必须 == data.nvSubtree[j]；
      // 同时验证这些后代确实在子树内（parents[d] >= j）。
      int nv = model.joints[j].nv();
      for (JointIndex d = j + 1; d <= c; ++d) // explore all descendant
      {
        CHECK_DATA(model.parents[d] >= j);

        nv += model.joints[d].nv();
      }

      CHECK_DATA(nv == data.nvSubtree[j]);

      // c 之后的关节不能是 j 的后代（父在 [j,c] 之外），保证子树是连续区间
      for (JointIndex d = c + 1; (int)d < model.njoints; ++d)
        CHECK_DATA((model.parents[d] < j) || (model.parents[d] > c));

      CHECK_DATA(
        data.nvSubtree[j] == data.nvSubtree_fromRow[(size_t)model.joints[j].idx_vExtended()]);

      // fromRow 父指针一致性：按关节的 idx_vExtended 行，校验 parents_fromRow 指向父关节的
      // 最后一个扩展 v 行；根(row==0)时父指针为 -1；并按父关节是否 mimic 核对 mimic/non_mimic 版本。
      int row = model.joints[j].idx_vExtended();
      const JointModel & jparent = model.joints[model.parents[j]];
      if (row == 0)
      {
        CHECK_DATA(data.parents_fromRow[(size_t)row] == -1);
        CHECK_DATA(data.mimic_parents_fromRow[(size_t)row] == -1);
        CHECK_DATA(data.non_mimic_parents_fromRow[(size_t)row] == -1);
      }
      else
      {
        CHECK_DATA(
          jparent.idx_vExtended() + jparent.nvExtended() - 1 == data.parents_fromRow[(size_t)row]);
        if (boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(&jparent))
        {
          CHECK_DATA(data.parents_fromRow[(size_t)row] == data.mimic_parents_fromRow[(size_t)row]);
        }
        else
        {
          CHECK_DATA(
            data.parents_fromRow[(size_t)row] == data.non_mimic_parents_fromRow[(size_t)row]);
        }
      }
    }

    // —— 组⑦：mimic 记账校验 —— data.mimic_subtree_joint[k] 应 == 该 mimicking 关节子树里
    //          第一个 nv≠0 的后代（确保 mimic 传动的记账正确）。
    if (model.mimicking_joints.size() != 0)
    {
      for (size_t k = 0; k < model.mimicking_joints.size(); k++)
      {
        // Check the mimic_subtree_joint
        const auto & mimicking_sub = model.subtrees[model.mimicking_joints[k]];
        size_t j = 1;
        JointIndex id_subtree = 0;
        for (; j < mimicking_sub.size(); j++)
        {
          if (model.nvs[mimicking_sub[j]] != 0)
          {
            id_subtree = mimicking_sub[j];
            break;
          }
        }
        CHECK_DATA(id_subtree == data.mimic_subtree_joint[k]);
      }
    }

#undef CHECK_DATA
    return true;
  }
  PINOCCHIO_COMPILER_DIAGNOSTIC_POP

} // namespace pinocchio
