# 可复现流程

## 目标与范围
- 目标：在现有数据与脚本基础上，根据审稿意见修订文章，并开展新的关联分析。
- 当前计划：CBCL total raw 与 SC 的相关；SC 与 S-A axis 的关联；harmonize/ComBat 方法按评论意见更新。
- 约束：不在此文档中自行决策新的 harmonize/ComBat 方案；仅记录已存在流程与待更新项。

## 数据与路径
- 新的数据输入与结果输出分别放在 `data/` 与 `outputs/`。
- 允许直接使用绝对路径以兼容既有数据与集群环境。
- 若使用绝对路径，请在此文档或会话记录中注明数据版本、来源与路径。
- `wd/` 作为历史结果与中间产物目录，保持不动。
- 本仓库不使用 `configs/` 或 `src/` 结构，避免将现有脚本迁移导致路径失效。

## 分析流程（对应现有脚本）
1) 数据整理与 SC 强度提取（`development_script/1st_dataclean`）
   - 合并数据与一致性阈值筛选。
   - 提取大尺度 SC 矩阵并生成汇总数据。
   - 进行站点/批次 harmonize（ComBat，待按评论更新）。
2) 发育模型与导数分析（`development_script/2nd_fitdevelopmentalmodel`）
   - GAM/GAMM 建模，生成拟合曲线与导数。
   - 计算与 S-A 轴的相关与可视化。
3) S-A 轴可视化（`development_script/3rd_plotConnectionalAxis`）
4) 年龄分段与变化率分析（`development_script/4th_changerate_SAcorr`）
5) 认知与心理病理相关分析（`development_script/5th_cognition`、`development_script/6th_pfactor`）

## 待补充说明
- 根据 `docs/research/Comments.pdf` 与 `docs/research/Manurscript_20251112.pdf` 更新 harmonize/ComBat 的描述与使用场景。
- 明确 CBCL total raw 与 SC 相关分析的输入表结构与协变量设置（仅记录，不做实现）。
