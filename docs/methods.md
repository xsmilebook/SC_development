# 方法学定义（基于现有脚本与研究文档）

## 数据与样本
- 主要数据集：ABCD 与 HCP-D（详见 `README.md` 中的数据来源说明）。
- 本文档基于现有脚本与 `docs/research/` 中的研究文档整理；需在后续逐条对照审稿意见更新细节。

## 结构连接（SC）强度提取与一致性阈值
- 细粒度 SC 网络基于 Schaefer-400 分区，连接权重为按节点体积归一化的纤维束计数。
- 排除边缘系统间连接，最终得到 376×376 网络中的 70,786 条边。
- 计算每条边的变异系数（CV），使用一致性阈值筛选伪连接；主分析采用 P75th CV 阈值（参考 Baum et al., 2020）。
- 见 `development_script/1st_dataclean/S1st_mergedata_*.R`。

## 大尺度 SC 矩阵与分辨率
- 将细粒度连接汇总为大尺度 SC 矩阵（例如 12×12），每条连接为经节点体积归一化后的汇总强度。
- 主分析分辨率为 12（78 条连接）；并在 7 与 17 分辨率上做敏感性分析。
- 见 `development_script/1st_dataclean/S2nd_mergedata_SA_ds_sumSC_*.R`。

## Harmonize/ComBat
- 现有脚本使用 ComBat 对多站点 SC 强度进行校正。
- 分别针对发育模型、认知模型与 p-factor 模型进行校正，协变量包含年龄、性别、头动等（认知与 p-factor 模型各自加入对应表型）。
- 见 `development_script/1st_dataclean/S3rd_combat_controlsite_*.R` 与 `gamfunction/combat.R`。
- 具体 harmonize/ComBat 方案需以 `docs/research/Comments.pdf` 与 `docs/research/Manurscript_20251112.pdf` 的修改意见为准，本文件不作自行决策。

## 发育模型（GAM/GAMM）
- HCP-D 使用 GAM，ABCD 使用 GAMM；每条边单独建模。
- 因变量为边的 SC 强度，年龄作为平滑项；协变量包括性别与头动。
- 平滑函数结点数 k=3，采用 REML 估计平滑参数。
- 为降低初始强度差异影响，将每条边的 SC 强度除以最小年龄处的拟合值（HCP-D：8 岁；ABCD：8.9 岁）得到比例化指标。
- 见 `development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_*.R`。

## 导数与后验导数
- 在年龄范围内等间隔采样 1,000 个年龄点，计算一阶导数并保存结果。
- 从每条边的平滑函数后验中抽样 1,000 次以估计导数的稳定性。
- 见 `development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_*.R`。

## S-A 连接轴与相关分析
- S-A 连接轴用于描述大尺度连接的梯度排序，并用于连接分组与可视化。
- 对发育参数与 S-A 连接轴排名进行 Spearman 相关分析；在相关前回归掉欧氏距离的影响。
- 见 `development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_*.R` 与 `development_script/2nd_fitdevelopmentalmodel/computeEuclidistance.R`。
- S-A 轴的可视化与分区文件见 `development_script/3rd_plotConnectionalAxis/`。

## 需要补充的来源信息
- 来自 `docs/research/Comments.pdf`、`docs/research/Manurscript_20251112.pdf` 与 `docs/research/SupplementaryMaterials20251112.pdf` 的方法学细节尚未逐条摘录，需在后续补充。
