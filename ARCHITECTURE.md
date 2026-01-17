# Project Structure

```
SCDevelopment/                  # 项目根目录
├── README.md                   # 项目简介与入口
├── AGENTS.md                   # AI 协作与修改约束
├── PROGRESS.md                 # 进度记录
├── ARCHITECTURE.md             # 结构说明（本文）
│
├── docs/                       # 方法、流程与决策记录（索引见 docs/README.md）
│   ├── README.md               # 文档集合索引
│   ├── workflow.md             # 可复现流程
│   ├── methods.md              # 方法学细节
│   ├── research/               # 论文、补充材料、审稿意见等来源
│   └── sessions/               # 会话记录（yy/mm/dd.md）
│
├── development_script/         # 论文主分析与敏感性分析脚本
├── gamfunction/                # 分析函数库（R）
├── dMRIprocessing/             # 预处理与重建脚本
├── demopath/                   # 人口学与行为数据表
├── combat_gam/                 # ComBat-GAM/纵向 ComBat 代码与外部依赖
├── data/                       # 新数据输入（不纳入版本控制）
├── outputs/                    # 新结果输出（不纳入版本控制）
└── wd/                         # 历史中间结果、模型输出与可视化数据（保持不动）
```

## Notes

- 本仓库不使用 `configs/`、`src/`、`scripts/` 等通用结构，按现有分析脚本组织。
- 为兼容历史数据与共享集群，允许使用绝对路径；路径与数据版本需在文档或会话记录中说明。
- 结构变更需明确指令；默认仅更新说明文档以匹配现有布局。
