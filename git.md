# 空间插值任务：9km → 1km 多变量插值与评估

本项目旨在实现从 9km 分辨率气象场到 1km 分辨率的高精度空间插值，并对多种插值算法进行系统性比较与验证。插值变量包括：**温度、降水、地表太阳辐照度、10米风速**。

## 📌 任务概述

- **插值目标**：将 9km 网格数据插值至 1km 分辨率  
- **插值变量**：2m 温度、总降水、地表太阳辐照、10m 风速  
- **插值算法实现与比较**：
  - 双线性插值（Bilinear Interpolation）
  - 反距离加权（Inverse Distance Weighting, IDW）
  - 克里金插值（Kriging）【自选算法】
- **验证策略**：
  - **自验证**：随机点交叉验证 + 区域误差分析（新疆、四川、山西、江苏）
  - **外部验证**：与 ERA5/ERA5-Land 再分析数据对比
- **可视化**：插值结果、误差分布、ERA5 对比图

---

## 由于时间及内存原因，部分任务仅选择新疆地区，ERA5数据因下载队列问题只有温度、uv风。

## 📂 项目结构

```bash
├── data/                     # 原始 9km 数据与 ERA5 数据（示例或链接）
├── scripts/
│   ├── interpolation.py      # 插值算法实现（三种方法）
│   ├── validation.py         # 自验证与 ERA5 验证
│   ├── visualization.py      # 绘图脚本（对比图、误差图等）
│   └── utils.py              # 工具函数（读取、重采样、指标计算等）
├── results/
│   ├── interpolation/
│   ├── validation/
│   └── figures/              # 所有生成的图像（含 README 中展示图）
├── environment.yml           # 推荐的 Conda 环境配置
└── README.md