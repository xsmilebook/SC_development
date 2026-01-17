# ComBat 实现计划（按审稿意见）

## 目标与范围
- 目标：按审稿意见更新批次校正流程（ComBat-GAM / Nonlinear-ComBat-GAM），不改变既有分析结构。
- 范围：仅规划与记录实现步骤、协变量配置与输出约定；不在本文档中做方法选择或代码实现。

## 方法选择与对应数据集
- HCP-D 与 Chinese Cohort：使用 ComBat-GAM（neuroHarmonize/`neuroCombat.R`）。
- ABCD：使用 Nonlinear-ComBat-GAM（`nonlinearlongcombat.R` + `neuroCombat.R`）。

## 基础模型与参数约束
- 批次变量：采集站点（site）。
- 协变量模型：年龄以平滑项建模；其余协变量为线性项。
- 平滑设定：`bs=tp`、`fx=TRUE`、`k=3`。
- GAMM 设定与 `gamfunction/gammsmooth.R` 保持一致（`gamm4` + `REML=TRUE` + `random=~(1|subID)`）。
- 公式依据：按 `docs/research/Manurscript_20251112.pdf` “Correction for multi-site batch effects” 部分的公式与描述执行（年龄平滑 + 性别 + 平均头动；认知与 p-factor 分别加入）。

## ABCD 的三套数据与协变量
- 方案 A（基础）：`age + sex + meanFD`。
- 方案 B（认知）：`age + sex + cognition + meanFD`。
- 方案 C（精神病理）：`age + sex + p-factor + meanFD`。

> 说明：`cognition` 为 NIH Toolbox fluid composite，`p-factor` 为 CBCL 双因子模型的一般因子；具体变量名与取值需与现有数据表一致。

## 输入、输出与路径约定
- 新数据：放入 `data/`，结果输出到 `outputs/`，历史 `wd/` 不改动。
- 允许使用绝对路径，但需在会话记录中注明数据版本与路径。
- 本次使用既有数据，仅调整方法流程与批次校正策略。

## 数据与路径清单（待补充）
- HCP-D：`/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_HCPD`（来源：`development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R`）
- ABCD：`/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD`（来源：`development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD.R`）
- Chinese Cohort：`/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort`
- 若以上路径在集群上缺失或迁移，需在此处补充新的绝对路径（含 ibmgpfs/GPFS 前缀）。

### 具体输入文件（用于 ComBat）
- HCP-D SC 数据：`/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds`
- HCP-D 协变量表：`/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment/demopath/HCPD_demo_behav.csv`
- ABCD SC 数据：`/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds`
- ABCD 协变量表：`/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment/demopath/DemodfScreenFinal.csv`
- Chinese Cohort SC 数据：`/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/<待确认文件名>`
- Chinese Cohort 协变量表：`/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/<待确认文件名>`

## 实现步骤（规划）
1) **梳理输入表与变量名**
   - 明确每个数据集的样本 ID、站点变量、年龄、性别、头动、SES、认知、p-factor 的字段名与缺失策略。
   - 记录在 `docs/workflow.md` 与会话记录中。
2) **HCP-D / Chinese：ComBat-GAM 计划**
   - 基于 `neuroCombat.R` 生成 ComBat-GAM 方案与设计矩阵。
   - 设计矩阵含：站点 batch、年龄平滑项（`bs=tp`、`fx=TRUE`、`k=3`）、性别、平均头动；认知与 p-factor 不并入该路径。
   - 输出：校正后的 SC 连接强度矩阵与对应日志。
3) **ABCD：Nonlinear-ComBat-GAM 计划**
   - 使用 `nonlinearlongcombat.R` 进行非线性与纵向结构的批次校正。
   - ABCD 合并基线与随访后统一校正（纵向方法一致）。
   - 按三套协变量分别运行，保持相同的年龄平滑设定（`bs=tp`、`fx=TRUE`、`k=3`）。
   - 调用 GAMM 的实现位置与 `gamfunction/gammsmooth.R` 保持一致。
   - 输出三套结果，各自标记对应协变量方案。
4) **一致性检查**
   - 确认校正后数据维度与原始 SC 连接数一致（例如 12×12 系统对应 78 条连接）。
   - 记录每个方案的样本量与排除情况。

## 交付物（文档与产物）
- `outputs/`：按数据集与协变量方案分层保存校正结果。
- `docs/workflow.md`：记录批次校正的输入字段、参数设定与输出路径。
- `PROGRESS.md` 与 `docs/sessions/`：记录执行与版本。

## 待确认问题
- 无。
