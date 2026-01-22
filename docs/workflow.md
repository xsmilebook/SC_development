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
- 构建（集群）：`sbatch sbatch/build_scdevelopment_r41_container.sbatch`，默认生成新镜像 `outputs/containers/scdevelopment_r41_<tag>.sif`（避免覆盖正在使用的旧镜像）；可用 `SIF_TAG` 自定义 `<tag>`。
- 运行（集群）：`sbatch sbatch/run_hcpd_devmodel_combatgam_CV75_container.sbatch`（与非容器版参数一致，支持 `N_EDGES`/`SKIP_POSTERIOR`）；如需指定新镜像，传入 `SIF_PATH=/.../scdevelopment_r41_<tag>.sif`。
  - 补图策略：脚本会在检测到 S3 的 decile 图或 S4 的相关性散点图缺失时，自动对对应步骤启用 `--force=1` 以回填图片（可用 `FORCE_S3=0/1`、`FORCE_S4=0/1` 覆盖）。
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
- ABCD 纵向 Nonlinear-ComBat-GAM（新增变体）：
  - CBCL total problems：`sbatch combat_gam/sbatch/abcd_combat_gam_cbcl.sbatch`，输出 `*combatgam_cbcl.rds`（协变量列：`cbcl_scr_syn_totprob_r`；不做 baseline-only）。
  - NIH Toolbox total cognition（age-corrected，baseline-only）：`sbatch combat_gam/sbatch/abcd_combat_gam_comp_agecorrected_baseline.sbatch`，输出 `*combatgam_comp_agecorrected_baseline.rds`（协变量列：`nihtbx_totalcomp_agecorrected`；仅保留 baseline 与 cognition 方案一致）。
  - 注：ABCD 的 `input_rds` 常仅包含 SC 与基础协变量；若缺少上述表型列，脚本会按 `scanID` 从 `demopath/DemodfScreenFinal.csv` 自动回填后再运行（回填失败会给出明确报错）。
- ABCD 的 Nonlinear-ComBat-GAM 支持并行：`nlongcombat` 使用 `mclapply`，核数由 `SLURM_CPUS_PER_TASK` 控制；未设置时默认单核。

### 非容器作业的 GLIBC 报错处理
- 典型报错：`/lib64/libc.so.6: version 'GLIBC_2.32' not found (required by .../cli.so)`（常由 `dplyr/cli` 等包触发）。
- 根因：在登录节点编译/安装的 R 包会链接登录节点的 GLIBC；当计算节点 GLIBC 更旧时，作业在 `dyn.load()` 阶段直接失败。
- 推荐做法（不使用容器时）：
  - 优先使用已在计算节点验证可用的 conda 环境（本项目默认 `scdevelopment`）。
  - 若必须使用 `scdevelopment_r41` 一类新环境，请在计算节点上重新安装触发问题的 R 包（从源码编译），或在同一 GLIBC 版本的节点上完成安装后再提交作业。
  - `combat_gam/sbatch/plot_*_variance_decomposition.sbatch` 已默认使用 `CONDA_ENV=scdevelopment`（可通过 `CONDA_ENV=scdevelopment_r41` 覆盖）；当出现 GLIBC 报错时请不要在计算节点继续使用会触发报错的环境。
- 结构连接 R² 方差分解图（Raw vs ComBat）由 `combat_gam/scripts/plot_abcd_variance_decomposition.R` 生成，输出在 `outputs/figures/combat_gam/`，包含：
  - `abcd_variance_decomp_base`（age+sex+meanFD）
  - `abcd_variance_decomp_cognition`（含 cognition）
  - `abcd_variance_decomp_pfactor`（含 p-factor）
  - `abcd_variance_decomp_cbcl_totprob`（含 CBCL total problems；读取 `cbcl_scr_syn_totprob_r`，必要时从 `demopath/DemodfScreenFinal.csv` 按 `scanID` 回填）
  - `abcd_variance_decomp_totalcomp_agecorrected`（含 NIH Toolbox total cognition age-corrected；读取 `nihtbx_totalcomp_agecorrected`，必要时从 `demopath/DemodfScreenFinal.csv` 按 `scanID` 回填；Raw 侧按 baseline-only 过滤）
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

## 发育模型脚本的“存在即跳过”约定（HCP-D）
- 为避免重复计算（尤其是 S1 拟合与 S2 posterior derivative），`development_script/2nd_fitdevelopmentalmodel` 的 HCP-D Rscript 版本默认采用 **存在即跳过**：
  - S1：当 `gamresults*.rds` 与 `gammodel*.rds` 对应输出已存在时跳过重新拟合；其他中间文件（如 `plotdatasum.df_*`、`SCdata.diw_*`、`*_scale_TRUE.rds`）同理按文件存在判断跳过。
  - S2：当 `derivative.df*.rds`（以及 posterior 的 `derivative.posterior*.rds`）存在时跳过。
  - S3：当 `plotdatasum_scale_TRUE_SA12.rds` 与关键图（`devcurve_Rsq_fit.ratio.tiff`、`devcurve_meanderv2_fit.Z.tiff`）存在时跳过。
  - S4：当 `SCrank_correlation_summary.csv` 存在时跳过。
- 强制重跑：为任一步脚本增加 `--force=1`（例如：`Rscript .../S3rd_visualizationfitSCcurves_SA12sumSCinvnode_HCPD.R --force=1`）。

## HCP-D 发育模型输出 QC（conda 环境）
- 目的：在 **conda R 环境** 下检查已生成的 `gamresults/gammodel`（raw+scaled）是否一致、是否为 `mgcv::gam` 对象、以及 `mgcv::predict.gam()` 是否可用（含 `se.fit` 路径）。
- 脚本：`development_script/2nd_fitdevelopmentalmodel/QC_check_gam_outputs_HCPD.R`
- 集群提交：`sbatch sbatch/qc_hcpd_devmodel_outputs_condenv.sbatch`（日志：`outputs/logs/qc_hcpd_gamout_%j.log`）
- 输出：`outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/qc/`
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
   - ABCD p-factor（GENERAL；Nonlinear-ComBat-GAM 输出 `*combatgam_pfactor.rds`）S1 复现入口（连续效应 + age-by-pfactor smooth interaction）：
     - sbatch（容器版，72 核）：`sbatch sbatch/run_abcd_pfactor_effect_continuous_container.sbatch`
     - 结果：`outputs/results/6th_pfactor/abcd/pfactor/`
   - ABCD fluid cognition（uncorrected；Nonlinear-ComBat-GAM 输出 `*combatgam_cognition.rds`）可复现入口（原始设定：控制 `age(smooth)+sex+mean_fd`）：
     - sbatch（容器版，72 核）：`sbatch sbatch/run_abcd_cognition_fluid_uncorrected_container.sbatch`
     - 结果：`outputs/results/5th_cognition/abcd/cognition/`
     - 图像：`outputs/figures/5th_cognition/abcd/cognition/`
   - ABCD total cognition（age-corrected，baseline-only；Nonlinear-ComBat-GAM 变体 `*combatgam_comp_agecorrected_baseline.rds`）可复现入口：
     - sbatch（容器版，72 核）：`sbatch sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch`
     - 结果：`outputs/results/5th_cognition/abcd/comp_agecorrected/`
     - 图像（tiff+pdf）：`outputs/figures/5th_cognition/abcd/comp_agecorrected/`
     - 注：为避免 `pandoc` 依赖，复现入口使用 `Rscript`（不走 `rmarkdown::render`）。
     - 协变量设定（更新）：在 `run_abcd_cognition_comp_agecorrected_{S1,S2}.R` 中，新增的 SC–cognition 关联不再调整 `age`（不含 `s(age, ...)`）且不包含 `sex`；当前仅控制 `mean_fd`。
     - 输出命名（避免覆盖）：脚本默认使用 `COG_ASSOC_TAG=meanfd_only` 作为文件名后缀，所有 `*.rds` 与图片均写成 `..._comp_agecorrected_<COG_ASSOC_TAG>.*`；如需并行保留多个变体，可在提交时覆盖，例如 `COG_ASSOC_TAG=v2_meanfd_only sbatch sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch`。
     - 欧氏距离控制项默认读取：`wd/interdataFolder_ABCD/average_EuclideanDistance_12.csv`（可用 `ABCD_EUCLID_CSV` 覆盖）。
     - 并行：脚本使用 `mclapply`（fork）并默认最多使用 60 个 worker（sbatch 仍可申请 72 CPU）；若遇到 `Cannot fork` 会按 60→50→40→30→20→… 自动降档直到可运行。

## CBCL 关联运行
- 默认使用容器镜像：`outputs/containers/scdevelopment_r41.sif`（可用 `SIF_PATH=/.../scdevelopment_r41_<tag>.sif` 指向新构建镜像）。
- CBCL 关联分析默认输入（ABCD Nonlinear-ComBat-GAM 输出）：`outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cbcl.rds`。
- 全量作业脚本（容器版）：`sbatch/run_cbcl_assoc_full.sbatch`（使用 `outputs/containers/scdevelopment_r41.sif`；可用 `SIF_PATH` 指向新构建镜像）。
- 运行前确认 `outputs/logs/` 存在，脚本内会自动创建。
- 若需非容器运行，建议参考 `docs/workflow.md` 中关于 R ABI 隔离的说明，避免作业环境误加载用户库（`/GPFS/.../R/packages`）导致 ABI 不一致报错。

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
- HCP-D 发育模型容器作业报 `Error in match.names(clabs, names(xi)) : names do not match previous names`（`do.call -> rbind`）：通常是某些边级模型拟合失败导致返回对象列名不一致，直接 `rbind` 合并会崩溃。
  - 处理：对边级拟合使用 `tryCatch` 并跳过失败边；合并时用 `dplyr::bind_rows()`；将失败边标签写入 `failed_edges_*.txt` 便于追踪与复跑。
- HCP-D 发育模型容器作业在保存图片时报 `Error in Ops.data.frame(guide_loc, panel_loc) : '==' only defined for equally-sized data frames`（`ggsave -> add_guides`）：
  - 常见触发：容器内 `ggplot2`/`patchwork` 版本组合不稳定（特别是 `ggplot2 4.x`）。
  - 处理：使用已固定到稳定版本的容器重新构建镜像（当前定义中 `ggplot2==3.5.2`、`patchwork==1.3.0`），然后用 `SIF_PATH=/.../scdevelopment_r41_<tag>.sif sbatch sbatch/run_hcpd_devmodel_combatgam_CV75_container.sbatch` 指向新镜像执行。
- ABCD cognition（`gamfunction/gamcog.R`）在并行/集群环境下报 `object 'nbinom2' not found`（来自 `ecostats::anovaPB`）：
  - 原因：`anovaPB` 内部并行或对象序列化在某些节点上不稳定。
  - 处理：在 `gamfunction/gamcog.R` 中强制 `anovaPB(..., ncpus=1)` 并对失败回退为 `NA`（后续按 `p=1` 处理），避免全边失败。
- HCP-D 发育模型容器作业在 S3（可视化）报 `Error in plotdatasum[[i]][, -14] : incorrect number of dimensions` 且伴随 `all scheduled cores encountered errors`：常见于 **上游 S1 跳过部分边导致 `gammodelsum`/`gamresults` 不再满 78 条**，但 S3 仍按固定 `elementnum` 遍历并用固定列号删列。
  - 处理：S3 改为按可用边数迭代（`min(length(gammodelsum), nrow(gamresults))`），对 `plotdata_generate()` 增加 `tryCatch`；删列按响应列名（与 `parcel` 同名）删除，避免依赖固定列位置。
- `plotdatasum_scale_TRUE_SA12.rds` 中出现 `try-error` 且报 `lm object does not have a proper 'qr' component`：通常是 `plotdata_generate()` 内部 `predict(..., se.fit=TRUE)` 在少数模型上失败导致。
  - 处理：`gamfunction/plotdata_generate.R` 对 `se.fit=TRUE` 增加容错，失败时回退到 `se.fit=FALSE`（仅保留 `fit`，CI 列为 NA），确保后续流程可继续运行。
- `ggsave()` 报 `The package "svglite" is required to save as SVG`：当前环境/容器未安装 `svglite`。
  - 处理：HCP-D 的 S3/S4 脚本默认不再输出 `.svg`，改为输出无需额外依赖的 `.pdf`（同时保留 `.tiff`）。
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
- **先隔离再诊断**：计算节点/登录节点差异（GLIBC、CRAN 访问）会把“包缺失/不可用”伪装为各种运行时错误；优先确保作业屏蔽 `R_LIBS_USER/R_LIBS` 并固定 `.libPaths()`，再判断是真缺包还是 ABI/网络问题。
- **容器优先补齐系统依赖**：`nloptr/lme4/gamm4/ragg/textshaping` 等失败常源于系统库缺失（如 `cmake/libwebp-dev/harfbuzz/freetype`），应优先写入 `containers/*.def`，比反复补装 R 包更稳定可复现。
- **版本固定优于“装最新版”**：R 4.1.x 下 CRAN 当前包常提高 R 版本门槛或发生 API 漂移（如 `gamm4`、`xfun/knitr`）；对关键链条使用 CRAN Archive 固定版本可显著降低构建/运行波动。
- **避免嵌套并行**：外层 `mclapply` 已并行时，`ecostats::anovaPB()` 等内部并行应显式收敛到单核（如 `ncpus=1`），避免序列化/worker 环境不一致导致的随机失败。
- **边级结果合并要容错**：边级拟合/预测一旦出现失败对象，`do.call(rbind, ...)` 容易整体崩溃；推荐 `tryCatch` 捕获失败并输出失败边清单，同时用 `dplyr::bind_rows()` 合并。
- **按标签对齐而非按序号**：下游步骤（S3/S4）不要假设永远有 78 条边且顺序固定；对 `SCrank`、欧氏距离、meanSC 等应按 `parcel` 映射对齐以兼容“部分缺失/跳过”的情况。
- **预测显式使用 mgcv 方法**：绘图预测应显式调用 `mgcv::predict.gam()`，并确保 `newdata` 的 factor level 与拟合数据一致；`se.fit=TRUE` 路径不稳定时提供回退（仅输出 `fit`、CI 为 NA）。
- **输出格式减少额外依赖**：`ggsave(.svg)` 依赖 `svglite`，在 HPC/容器环境中常缺；默认采用 `tiff + pdf` 更稳，降低系统依赖与编译链风险。
- **默认“存在即跳过”，保留强制重跑**：按步骤检测输出文件存在即可从失败点续跑并节省计算资源；同时提供 `--force=1` 便于单步重跑。
- **把排错固化为 QC**：用 conda 环境直接读取 RDS 做一致性与 `predict.gam` 可用性检查，可快速区分“模型问题”与“下游计算/绘图问题”，并输出可追踪的 QC 报告供复跑与比较。
- **neuroHarmonize 平滑项**：若通过 `smoothCon(... )$X` 手动构造基函数并作为线性协变量输入，`fx` 的影响通常不会体现在 `X` 上；要让固定/惩罚平滑生效，应使用 neuroHarmonize 原生 `smooth_terms`（`smooth_fx`）。
- **statsmodels spline 约束**：使用 `harmonizationLearn(..., smooth_terms=...)` 时，`smooth_df` 需满足 `smooth_df >= smooth_degree + 1`（例如三次样条 `degree=3` 时 `df>=4`）。
- **方差分解图耗时**：绘图阶段会对每条边重复拟合多次模型（Raw/ComBat × 多个 covariate 集合），通常远慢于一次性 ComBat；序列（sequential）R² 相比全子集方法更可控，且可通过并行显著缩短总耗时。
- **方差分解图解读**：`facet_grid(..., scales=\"free_y\")` 会放大 ComBat 面板的视觉高度；判断校正效果应结合数值计算（如导出 `siteID` 的 mean/max），避免仅凭堆叠柱形高度下结论。
- **数据/结果不入库**：`demopath/`、`outputs/`、`wd/`、`containers/` 等数据目录应在 `.gitignore` 忽略；若已被追踪，需要 `git rm -r --cached` 从索引移除（不删除磁盘文件）。
