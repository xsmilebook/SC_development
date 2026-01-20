# 集群使用说明

## Slurm 提交约定
- 仅在需要时使用 `sbatch` 提交任务，队列使用 `q_fat_c`。
- 集群不限制运行时间，内存由系统动态分配，**禁止**在 `sbatch` 脚本中设置 `--time` 或 `--mem`。
- 单节点可申请 72 核（通过 `--cpus-per-task` 设置）。
- 日志输出路径保持统一，便于检索。

## 并行执行约定
- ABCD 的 Nonlinear-ComBat-GAM 在 R 内部使用 `mclapply` 并行，核数读取 `SLURM_CPUS_PER_TASK`。
- 并行任务需限制 BLAS 线程（如 `OPENBLAS_NUM_THREADS=1`、`OMP_NUM_THREADS=1`），避免线程过载与创建失败。

## 容器（Singularity/3.7.0）
- 使用前：`module load singularity/3.7.0`
- 常见构建报错：`FATAL: You must be the root user ... you can use --remote or --fakeroot to build from a Singularity recipe file`
  - 解释：当前集群不允许普通用户直接本地构建（需要 root / setuid），且若未为用户配置 `/etc/subuid`，`--fakeroot` 也不可用。
  - 处理：使用 `--remote` 构建（由远端 builder 完成 root 权限构建），或联系管理员配置 fakeroot/subuid。
- 网络：计算节点联网通常需要代理；容器构建脚本应在 sbatch 内设置 `http_proxy/https_proxy/...`（若不需要可忽略）。
