# 集群使用说明

## Slurm 提交约定
- 仅在需要时使用 `sbatch` 提交任务，队列使用 `q_fat_c`。
- 集群不限制运行时间，内存由系统动态分配，**禁止**在 `sbatch` 脚本中设置 `--time` 或 `--mem`。
- 单节点可申请 72 核（通过 `--cpus-per-task` 设置）。
- 日志输出路径保持统一，便于检索。

## 并行执行约定
- ABCD 的 Nonlinear-ComBat-GAM 在 R 内部使用 `mclapply` 并行，核数读取 `SLURM_CPUS_PER_TASK`。
