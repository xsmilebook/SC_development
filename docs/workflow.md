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

## 容器运行（Singularity/Apptainer）
- 背景：部分计算节点无法加载在登录节点安装/更新过的 conda R 包（GLIBC 版本不匹配），且计算节点常无法访问 CRAN，导致运行时反复报错（`GLIBC_2.xx not found` / `cannot open URL .../PACKAGES`）。
- 解决：使用项目内 Singularity 容器（R 4.1.3，固定 mgcv/nlme/ecostats/gratia 等关键版本），避免节点差异导致的依赖不一致。
- 构建（集群）：`sbatch sbatch/build_scdevelopment_r41_container.sbatch`，生成镜像 `outputs/containers/scdevelopment_r41.sif`。
- 运行（集群）：`sbatch sbatch/run_hcpd_devmodel_combatgam_CV75_container.sbatch`（与非容器版参数一致，支持 `N_EDGES`/`SKIP_POSTERIOR`）。
- 定义文件：`containers/scdevelopment_r41.def`（可按需扩展依赖，但应尽量保持版本稳定）。
  - 常见构建报错：`could not use fakeroot: no mapping entry found in /etc/subuid for <user>`：说明集群未为该用户配置 fakeroot/subuid。当前构建脚本不使用 `--fakeroot`；若仍失败需联系管理员开启 setuid build 或提供可用的 fakeroot 配置。
  - 常见构建报错：`You must be the root user ... use --remote or --fakeroot`：说明本地 build 需要 root/setuid；当前构建脚本使用 `singularity build --remote`。
  - 常见构建报错：`Error: object ‘attr’ is not exported by 'namespace:xfun'`（`knitr` lazy-load 失败）：通常是 `xfun` API 随 CRAN 更新发生漂移。处理：在容器定义中固定 `xfun==0.52` 与 `knitr==1.43`，并确保安装时不自动拉取 `Suggests` 依赖导致的版本覆盖。
  - 网络：计算节点联网需代理；构建脚本已内置 `http(s)_proxy` 等环境变量导出。

## ComBat-GAM 运行约定
- 小型测试可直接运行 `combat_gam/scripts/*.sh`；正式任务必须使用 `combat_gam/sbatch/*.sbatch` 提交到 `q_fat_c`。
- ComBat-GAM 使用项目专用环境：`/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/envs/scdevelopment`。
- HCP-D/Chinese 推荐使用 neuroHarmonize 原生 `smooth_terms`（GAM）实现：`combat_gam/scripts/run_combat_gam_neuroharmonize_native.py`（支持 `smooth_fx` 固定平滑或 k-fold 选择）。
- 运行相关依赖应安装在 `scdevelopment` 环境；若包依赖冲突需建立独立虚拟环境，并在此处补充说明。
- Reviewer2（Q5）补充：ABCD baseline 的 `age+sex+meanFD` 纵向 ComBat（不保护 cognition）可用 `combat_gam/sbatch/abcd_combat_gam_baseline_age_sex_meanfd.sbatch` 提交，输出 `*combatgam_age_sex_meanfd_baseline.rds`。
- ABCD 的 Nonlinear-ComBat-GAM 支持并行：`nlongcombat` 使用 `mclapply`，核数由 `SLURM_CPUS_PER_TASK` 控制；未设置时默认单核。
- 结构连接 R² 方差分解图（Raw vs ComBat）由 `combat_gam/scripts/plot_abcd_variance_decomposition.R` 生成，输出在 `outputs/figures/combat_gam/`，包含：
  - `abcd_variance_decomp_base`（age+sex+meanFD）
  - `abcd_variance_decomp_cognition`（含 cognition）
  - `abcd_variance_decomp_pfactor`（含 p-factor）
- 方差分解的变量解释量采用 **序列（sequential）R²**：按 `age → sex → mean_fd → (cognition/pfactor/cbcl 等) → site` 顺序逐步加入变量，每一步的解释量为 `R²_k - R²_{k-1}`（不再计算所有子集 Shapley）。
- 集群绘图可直接提交：
  - `combat_gam/sbatch/plot_abcd_variance_decomposition.sbatch`
  - `combat_gam/sbatch/plot_hcpd_chinese_variance_decomposition.sbatch`

## HCP-D/Chinese ComBat-GAM（neuroHarmonize 原生 smooth_terms）
- HCP-D 提交：`sbatch combat_gam/sbatch/hcpd_combat_gam_native.sbatch`
- Chinese 提交：`sbatch combat_gam/sbatch/chinese_combat_gam_native.sbatch`
- 输出默认写入：
  - `outputs/results/combat_gam/hcpd/*combatgam_native.rds`
  - `outputs/results/combat_gam/chinese/*combatgam_native.rds`
- 注：statsmodels GAM 的 B-spline 要求 `smooth_df >= smooth_degree + 1`；当前脚本默认使用 `smooth_degree=3` 且 `smooth_df=4`。
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
   - HCP-D（基于 ComBat-GAM 输出）的可复现运行入口：
     - 输入默认：`outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds`
     - sbatch：`sbatch/run_hcpd_devmodel_combatgam_CV75.sbatch`
     - 产物目录：
       - intermediates：`outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/`
       - results：`outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/`
       - figures：`outputs/figures/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/`
3) S-A 轴可视化（`development_script/3rd_plotConnectionalAxis`）
4) 年龄分段与变化率分析（`development_script/4th_changerate_SAcorr`）
5) 认知与心理病理相关分析（`development_script/5th_cognition`、`development_script/6th_pfactor`）
   - CBCL total raw 关联分析脚本见 `development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD.R`。
   - 小样本验证脚本见 `development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD_smalltest.R`。
   - CBCL 关联分析绘图输出（t-value matrix、S-A rank 散点、分位数组轨迹）写入 `outputs/figures/cbcl_totprob/`。

## Conda 环境运行（CBCL 关联）
- 统一使用本地 conda 环境：`/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/envs/scdevelopment_r41`（R 4.1.3）。
- 小样本测试作业脚本：`sbatch/run_cbcl_assoc_smalltest.sbatch`（脚本内已 `conda activate`）。
- 全量作业脚本：`sbatch/run_cbcl_assoc_full.sbatch`（脚本内已 `conda activate`）。
- 运行前确认 `outputs/logs/` 存在，脚本内会自动创建。
- 两个 sbatch 脚本均会显式设置 `R_LIBS_USER`/`R_PROFILE_USER`，避免作业环境误加载用户库（`/GPFS/.../R/packages`）导致 ABI 不一致报错。

## 待补充说明
- 根据 `docs/research/Comments.pdf` 与 `docs/research/Manurscript_20251112.pdf` 更新 harmonize/ComBat 的描述与使用场景。
- 明确 CBCL total raw 与 SC 相关分析的输入表结构与协变量设置（与 p-factor 保持一致，ComBat 使用 `gamfunction/combat.R`）。

## 常见报错与处理
- `lme4` 载入失败（GLIBC 版本不匹配，指向 `GPFS/.../R/packages/lme4`）：sbatch 环境优先加载用户库导致；在脚本内清空 `R_LIBS_USER`/`R_LIBS` 并显式设置 `.libPaths()` 到 conda 环境库。
- `gratia`/`psych` 载入失败（GLIBC 版本不匹配）：优先确保作业激活 conda 环境并避免用户库路径污染。
- `cli.so`/`farver.so`/`Rcpp.so` 报 `GLIBC_2.xx not found`（HCP-D 发育模型、以及依赖 tidyverse 的脚本中较常见）：通常是 **登录节点与计算节点的 GLIBC 版本不一致**，导致在登录节点安装/更新过的 R 包（即使位于 conda 环境目录）在计算节点无法加载。
  - 处理：将 `R_LIBS_USER` 指向项目内的 R 库目录（例如 `outputs/r_libs/scdevelopment_r41/`），并在 sbatch 中对缺失/载入失败的包执行 `install.packages(..., lib=R_LIBS_USER)`，使其在计算节点上编译安装以匹配该节点的 GLIBC；同时保留 `R_LIBS_SITE=${CONDA_PREFIX}/lib/R/library` 作为回退路径。
  - 现有实现：`sbatch/run_hcpd_devmodel_combatgam_CV75.sbatch` 在运行前会自动创建并使用 `outputs/r_libs/scdevelopment_r41/`，并在“屏蔽系统/用户库”条件下检查与补装依赖。
- `install.packages()` 报 `cannot open URL .../PACKAGES` 或误报 “package is not available for this version of R”（HCP-D 发育模型中曾频繁出现）：通常是 **计算节点无法访问 CRAN**（网络策略/HTTPS 受限），导致无法拉取索引文件。
  - 处理：在登录节点预先构建离线 CRAN 源码仓库（`outputs/r_cran_repo/src/contrib/` + `PACKAGES.gz`），sbatch 在计算节点只从本地 repo 安装依赖（`options(repos=c(CRAN="file://.../outputs/r_cran_repo"))`），无需外网。
  - 现有实现：`sbatch/run_hcpd_devmodel_combatgam_CV75.sbatch` 会检查 `outputs/r_cran_repo/src/contrib/PACKAGES.gz` 是否存在；若存在则用离线 repo 安装缺失/载入失败的包到 `outputs/r_libs/scdevelopment_r41/`。
- `gratia` 报 `there is no package called 'mvnfast'`：这是 `gratia` 的依赖缺失；若仅运行 CBCL 关联分析，已移除 `gamfunction/gamminteraction.R` 对 `gratia` 的依赖以避免该类问题。
- `tidyverse` 报 `readr/forcats/lubridate` 缺失：在 `scdevelopment` 环境安装 `r-tidyverse`（或补装 `r-readr`、`r-forcats`）。
- CBCL ComBat 输出中 `age` 列为嵌套 data.frame：在 CBCL 关联脚本中从 `demopath/DemodfScreenFinal.csv` 按 `scanID` 回填 `age`。
- `gratia` 需与 `ggplot2`/`mgcv` 版本匹配：R 4.1.x 下使用 `gratia_0.8.1`（CRAN Archive）以兼容 `mgcv 1.8`。
- `dplyr`/`rlang` 报 `undefined symbol: R_existsVarInFrame`：通常是 `R_LIBS_USER` 指向旧用户库（如 `/GPFS/.../R/packages`）导致 ABI 不一致；运行前显式设置 `R_LIBS_USER`/`R_LIBS` 到当前 conda 环境库路径，并用 `.libPaths()` 置顶该路径。
  - 注：`sbatch/run_cbcl_assoc_{smalltest,full}.sbatch` 已内置该隔离设置。

## 可复用经验（本次会话）
- **环境隔离**：为复现旧版脚本，建议新建独立 conda 环境（如 `scdevelopment_r41`，R 4.1.3），避免在旧环境上补装导致依赖漂移；必要时用 CRAN Archive 固定包版本（例如 `gratia_0.8.1`）。
- **sbatch 入参约定**：若脚本要求输出路径（如 Chinese ComBat），sbatch 模板应显式传参并将日志写到 `outputs/logs/`，避免默认 `slurm-*.out` 分散且误判为“无日志”。  
- **R ABI 报错排查**：出现 `rlang.so undefined symbol`、`GLIBC` 不匹配等问题时，优先排查作业是否误加载用户库（如 `/GPFS/.../R/packages`）；在 sbatch 内用 `Rscript --vanilla` 并显式设置 `R_LIBS_USER` 到 conda 库可显著降低此类问题。
- **并行与线程**：脚本并行尽量从 `SLURM_CPUS_PER_TASK` 读取核数（如 `mclapply`/`nlongcombat`）；同时限制 `OPENBLAS_NUM_THREADS/OMP_NUM_THREADS/MKL_NUM_THREADS=1` 避免线程过订阅导致性能下降或报错。
- **neuroHarmonize 平滑项**：若通过 `smoothCon(... )$X` 手动构造基函数并作为线性协变量输入，`fx` 的影响通常不会体现在 `X` 上；要让固定/惩罚平滑生效，应使用 neuroHarmonize 原生 `smooth_terms`（`smooth_fx`）。
- **statsmodels spline 约束**：使用 `harmonizationLearn(..., smooth_terms=...)` 时，`smooth_df` 需满足 `smooth_df >= smooth_degree + 1`（例如三次样条 `degree=3` 时 `df>=4`）。
- **方差分解图耗时**：绘图阶段会对每条边重复拟合多次模型（Raw/ComBat × 多个 covariate 集合），通常远慢于一次性 ComBat；序列（sequential）R² 相比全子集方法更可控，且可通过并行显著缩短总耗时。
- **方差分解图解读**：`facet_grid(..., scales=\"free_y\")` 会放大 ComBat 面板的视觉高度；判断校正效果应结合数值计算（如导出 `siteID` 的 mean/max），避免仅凭堆叠柱形高度下结论。
- **数据/结果不入库**：`demopath/`、`outputs/`、`wd/`、`containers/` 等数据目录应在 `.gitignore` 忽略；若已被追踪，需要 `git rm -r --cached` 从索引移除（不删除磁盘文件）。
