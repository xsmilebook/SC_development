# AGENTS.md

在 `sc_connectome_trajectories` 中进行 AI 协作的操作规则。
本仓库以 `ARCHITECTURE.md` 为稳定结构的唯一来源，除非明确要求，不得提出结构性变更。

## 基础流程

1) 在规划或编辑前先阅读 `ARCHITECTURE.md`。
2) 仅依据 `PLAN.md` 或用户明确指令执行更改。
3) 若更改影响使用方式或结构，更新 `README.md` 与/或 `docs/workflow.md`。
4) 每完成一组更改后更新 `PROGRESS.md`。
5) 每次会话记录到 `docs/sessions/`（按日期命名）。
6) 若修改了代码或文档，需用 git 提交（自行撰写提交信息）。
7) 如使用 `create-plan` 技能，需将计划写入根目录 `PLAN.md`。

## 范围与约束

- 默认范围为仅文档，除非用户明确要求工程更改。
- 不修改 `data/` 与 `outputs/` 下的运行时产物。
- 未经明确要求不得更改目录名称或结构。
- 保持最小化改动；避免与需求无关的重构。

## 工程约定（仅在被要求时）

- 集群 GPU 工作遵循 `docs/cluster_gpu_usage.md`（Slurm + Singularity + 不在节点上修改环境）。
- 以 `python -m scripts.<entry>` 作为执行入口。
- 避免临时 `sys.path` 修改；如不可避免，仅限入口脚本并说明原因。
- 所有路径必须来自 `configs/`。

## 文档约定

- 使用严谨、科学的语言。
- `docs/` 内内容用中文撰写，文件名保持英文。
- 根目录文件名保持英文。
- sc_connectome_trajectories 相关说明集中在 `docs/workflow.md` 与 `configs/paths.yaml`（数据集部分）。
