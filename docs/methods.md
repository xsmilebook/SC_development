# 方法学定义（基于现有脚本与研究文档）

## 数据与样本
- 主要数据集：HCP-D（发现集）、ABCD（独立验证与行为关联）、中国队列（跨文化泛化）。
- 本文档基于现有脚本与 `docs/research/` 中的研究文档整理；需在后续逐条对照审稿意见更新细节。

### HCP-D
- 数据来源：NDA Lifespan HCP-D release 2.0。
- 初始样本 652 名，排除：dMRI 不完整、解剖异常、年龄小于 8 岁、头动过大（FD > 均值 + 3SD）。
- 最终样本：590 名（273 男，8.1–21.9 岁）。

### ABCD
- 数据来源：ABCD release 5.1（影像来自 fast-track portal 2022-06）。
- 仅使用 SIEMENS 扫描以减少厂商偏差：基线 5,803 次、2 年随访 4,547 次（含 dMRI、场图、T1WI）。
- 主要排除：不符合官方影像纳入标准、dMRI 缺失/失败、语言与健康相关排除、早产/低体重、头动过大等。
- 最终纳入：7,104 次扫描（基线 3,949 次，8.9–11.0 岁；随访 3,155 次，10.6–13.8 岁）。

### 中国队列（EFNY / devCCNP / SAND）
- 初始 723 次扫描（含 dMRI 与 T1WI）。
- 主要排除：严重感官/神经问题、早产/低体重、年龄或性别缺失、处理失败、年龄小于 6 岁、头动过大。
- 最终纳入：609 次扫描（6.1–22.0 岁），其中 EFNY 152 次、devCCNP 312 次（59 名多次测量）、SAND 145 次。

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

## 行为与表型评估
- 认知：HCP-D 与 ABCD 使用 NIH Toolbox Cognition Battery，五个任务的未年龄校正标准分（均值 100，SD 15）平均后再标准化得到 fluid composite。
- 精神病理：ABCD 使用 CBCL 量表，基于既有双因子模型估计 p-factor，并在全样本上进行跨时间点的度量不变性约束。

## 影像采集与预处理
- 采集：HCP-D、ABCD、EFNY、SAND 使用 3T SIEMENS；devCCNP 使用 GE Discovery MR750 3T。
- T1WI：HCP-D 使用最小化处理数据；ABCD 与 EFNY 使用 SIEMENS 标准化 T1WI；中国队列其余为原始 T1WI。
- 预处理：QSIPrep 0.16.0，用于 T1WI 与 dMRI；T1WI 进行强度非均匀校正、去除非脑组织、MNI 标准化。
- FreeSurfer 7.1.1：表面与组织分割；HCP-D 使用 HCP pipeline 工作流，其余使用 recon-all。
- dMRI：多步骤预处理（MP-PCA 去噪、Gibbs 去环、B1 校正、eddy 纠正、与 T1WI 配准等）；devCCNP 与 SAND 无场图，未做相位编码畸变校正。

## 结构连接重建与权重计算
- 纤维追踪：MRtrix3 iFOD2，ACT + HSVS 约束；每次扫描生成 1,000 万条流线，长度 30–250 mm。
- 方向模型：多壳使用 MSMT-CSD；单壳使用 SS3T-CSD。
- 权重：SIFT2 计算每条流线权重。
- 系统划分：Schaefer-400 沿 S-A 轴排序，去除边缘系统，剩余 376 区域并分为 12 个系统（每系统 31–32 区），并在 7 与 17 系统上做敏感性分析。
- 连边定义：流线端点 2 mm 半径搜索近邻灰质节点；CV>75th 的连边被剔除。
- SC 强度：连边流线数 × SIFT2 权重 / 连接系统体积均值。

## 敏感性分析（补充材料）
- 替代网络规模：7 系统与 17 系统（S-A 轴划分）。
- 替代分区：Yeo-7 与 Yeo-17（基于 Schaefer-400），去除边缘系统后重建大尺度 SC。
- TractSeg：基于 72 条主束重建连接，以减少全局追踪噪声影响。
- 距离控制：回归掉系统间欧氏距离后再评估与 S-A 轴的关联。
- 额外协变量：SES 与 ICV。
- 轴定义替代：S-A 连接轴由“系统轴秩之积”替代“系统轴秩平方和”。

## 需要补充的来源信息
- 来自 `docs/research/Comments.pdf`、`docs/research/Manurscript_20251112.pdf` 与 `docs/research/SupplementaryMaterials20251112.pdf` 的方法学细节尚未逐条摘录，需在后续补充。
