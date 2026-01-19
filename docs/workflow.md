# 可复现流程

## 目标与范围
- 目标：在现有数据与脚本基础上，根据审稿意见修订文章，并开展新的关联分析。
- 当前计划：CBCL total raw 与 SC 的相关（含 S-A axis 相关）；harmonize/ComBat 方法按评论意见更新。
- 约束：不在此文档中自行决策新的 harmonize/ComBat 方案；仅记录已存在流程与待更新项。

## 数据与路径
- 新的数据输入与结果输出分别放在 `data/` 与 `outputs/`。
- 允许直接使用绝对路径以兼容既有数据与集群环境。
- 若使用绝对路径，请在此文档或会话记录中注明数据版本、来源与路径。
- 所有新生成的数据与结果必须保存在当前项目路径下。
- `wd/` 作为历史结果与中间产物目录，保持不动。
- 本仓库不使用 `configs/` 或 `src/` 结构，避免将现有脚本迁移导致路径失效。
- `combat_gam/` 用于集中存放 ComBat-GAM 与纵向 ComBat 相关代码及其依赖。
- ComBat 运行日志统一写入 `outputs/logs/combat_gam/`，输出结果保存到 `outputs/results/combat_gam/`。

## ComBat-GAM 运行约定
- 小型测试可直接运行 `combat_gam/scripts/*.sh`；正式任务必须使用 `combat_gam/sbatch/*.sbatch` 提交到 `q_fat_c`。
- ComBat-GAM 使用项目专用环境：`/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/envs/scdevelopment`。
- HCP-D/Chinese 的 GAM 由 `rpy2` 调用 `mgcv::smoothCon` 生成 `s(Age, k=3, bs="tp")` 基函数矩阵，并作为协变量输入 `neuroHarmonize`（不再启用自动平滑）。
- 运行相关依赖应安装在 `scdevelopment` 环境；若包依赖冲突需建立独立虚拟环境，并在此处补充说明。
- ABCD 的 Nonlinear-ComBat-GAM 支持并行：`nlongcombat` 使用 `mclapply`，核数由 `SLURM_CPUS_PER_TASK` 控制；未设置时默认单核。
- 结构连接 R² 方差分解图（Raw vs ComBat）由 `combat_gam/scripts/plot_abcd_variance_decomposition.R` 生成，输出在 `outputs/figures/combat_gam/`，包含：
  - `abcd_variance_decomp_base`（age+sex+meanFD）
  - `abcd_variance_decomp_cognition`（含 cognition）
  - `abcd_variance_decomp_pfactor`（含 p-factor）
- CBCL total raw 的 Raw vs ComBat 方差分解图由 `development_script/1st_dataclean/S4th_plot_abcd_variance_decomp_cbcl.R` 生成，输出为 `abcd_variance_decomp_cbcl_totalraw`；ComBat 输出默认写入 `outputs/results/combat_cbcl/`。
- CBCL total raw 的 siteID R² 明细与汇总由 `development_script/1st_dataclean/S5th_export_cbcl_site_r2_summary.R` 导出，输出为 `outputs/results/combat_cbcl/cbcl_site_r2_raw_vs_combat.csv` 与 `cbcl_site_r2_summary.csv`。
- HCPD 与 Chinese Cohort 的 Raw vs ComBat 方差分解图由 `combat_gam/scripts/plot_hcpd_chinese_variance_decomposition.R` 生成，输出为：
  - `hcpd_variance_decomp`
  - `chinese_variance_decomp`
  - 注：Chinese Cohort 中 `study`（siteID）与年龄分布强烈不均衡时，ComBat-GAM 只能在保留年龄效应的前提下校正批次差异，部分 residual site 效应可能仍存在。

## 分析流程（对应现有脚本）
1) 数据整理与 SC 强度提取（`development_script/1st_dataclean`）
   - 合并数据与一致性阈值筛选。
   - 提取大尺度 SC 矩阵并生成汇总数据。
   - 进行站点/批次 harmonize（ComBat，待按评论更新）；CBCL total raw 可用 `S3rd_combat_controlsite_ABCD_CBCLtotalraw.R` 的测试参数先小样本跑通。
2) 发育模型与导数分析（`development_script/2nd_fitdevelopmentalmodel`）
   - GAM/GAMM 建模，生成拟合曲线与导数。
   - 计算与 S-A 轴的相关与可视化。
3) S-A 轴可视化（`development_script/3rd_plotConnectionalAxis`）
4) 年龄分段与变化率分析（`development_script/4th_changerate_SAcorr`）
5) 认知与心理病理相关分析（`development_script/5th_cognition`、`development_script/6th_pfactor`）
   - CBCL total raw 关联分析脚本见 `development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD.R`。
   - 小样本验证脚本见 `development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD_smalltest.R`。

## Conda 环境运行（CBCL 关联）
- 统一使用本地 conda 环境：`/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/envs/scdevelopment_r41`（R 4.1.3）。
- 小样本测试作业脚本：`sbatch/run_cbcl_assoc_smalltest.sbatch`（脚本内已 `conda activate`）。
- 全量作业脚本：`sbatch/run_cbcl_assoc_full.sbatch`（脚本内已 `conda activate`）。
- 运行前确认 `outputs/logs/` 存在，脚本内会自动创建。

## 待补充说明
- 根据 `docs/research/Comments.pdf` 与 `docs/research/Manurscript_20251112.pdf` 更新 harmonize/ComBat 的描述与使用场景。
- 明确 CBCL total raw 与 SC 相关分析的输入表结构与协变量设置（与 p-factor 保持一致，ComBat 使用 `gamfunction/combat.R`）。

## 常见报错与处理
- `lme4` 载入失败（GLIBC 版本不匹配，指向 `GPFS/.../R/packages/lme4`）：sbatch 环境优先加载用户库导致；在脚本内清空 `R_LIBS_USER`/`R_LIBS` 并显式设置 `.libPaths()` 到 conda 环境库。
- `gratia`/`psych` 载入失败（GLIBC 版本不匹配）：优先确保作业激活 conda 环境并避免用户库路径污染。
- `tidyverse` 报 `readr/forcats/lubridate` 缺失：在 `scdevelopment` 环境安装 `r-tidyverse`（或补装 `r-readr`、`r-forcats`）。
- CBCL ComBat 输出中 `age` 列为嵌套 data.frame：在 CBCL 关联脚本中从 `demopath/DemodfScreenFinal.csv` 按 `scanID` 回填 `age`。
- `gratia` 需与 `ggplot2`/`mgcv` 版本匹配：R 4.1.x 下使用 `gratia_0.8.1`（CRAN Archive）以兼容 `mgcv 1.8`。
- `dplyr`/`rlang` 报 `undefined symbol: R_existsVarInFrame`：通常是 `R_LIBS_USER` 指向旧用户库（如 `/GPFS/.../R/packages`）导致 ABI 不一致；运行前显式设置 `R_LIBS_USER`/`R_LIBS` 到当前 conda 环境库路径，并用 `.libPaths()` 置顶该路径。

## 可复用经验（本次会话）
- **环境隔离**：为复现旧版脚本，建议新建独立 conda 环境（如 `scdevelopment_r41`，R 4.1.3），避免在旧环境上补装导致依赖漂移；必要时用 CRAN Archive 固定包版本（例如 `gratia_0.8.1`）。
- **sbatch 入参约定**：若脚本要求输出路径（如 Chinese ComBat），sbatch 模板应显式传参并将日志写到 `outputs/logs/`，避免默认 `slurm-*.out` 分散且误判为“无日志”。  
- **方差分解图解读**：`facet_grid(..., scales=\"free_y\")` 会放大 ComBat 面板的视觉高度；判断校正效果应结合数值计算（如导出 `siteID` 的 mean/max），避免仅凭堆叠柱形高度下结论。
- **数据/结果不入库**：`demopath/`、`outputs/`、`wd/`、`containers/` 等数据目录应在 `.gitignore` 忽略；若已被追踪，需要 `git rm -r --cached` 从索引移除（不删除磁盘文件）。
