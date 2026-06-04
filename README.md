# 浮游生物群落生态调查分析

基于 2025 年 5 月三次野外调查的水质与浮游生物数据的定量生态学分析。

## 目录结构

```
bio/
├── README.md                 # 项目说明
├── requirements.txt          # Python 依赖
├── .gitignore
│
├── data/                     # 数据文件
│   ├── water.txt             # 清洗后的水质数据（TSV）
│   └── data.txt              # 清洗后的生物样本数据（TSV）
│
├── scripts/                  # 分析脚本
│   ├── clean.py              # 数据清洗脚本（文档+验证）
│   ├── analysis.py           # 基础多样性指数计算
│   ├── advanced_analysis.py  # PCA/NMDS/聚类/功能群/指示物种等
│   ├── mantel_simper.py      # Mantel检验 + SIMPER分析
│   └── gen_report.py         # 综合Word报告生成
│
└── output/                   # 分析输出
    ├── figures/              # 7张PNG分析图
    ├── results/              # 文本结果文件
    └── report/               # 正式Word报告
```

## 数据概况

- **采样点**：黎照湖、香雪海、梦川、后山水池、菜根潭（5个）
- **调查批次**：2025-05-16、2025-05-23、2025-05-31（3次）
- **总样本数**：15个（5样点 × 3次）
- **记录物种**：~50种，分属13个功能类群
- **水质参数**：TDS、盐度、pH、电导率、温度

## 分析方法

| 方法 | 说明 | 脚本 |
|------|------|------|
| Shannon-Wiener H' | 多样性指数 | analysis.py |
| Simpson D | 多样性指数 | analysis.py |
| Pielou J' | 均匀度指数 | analysis.py |
| Margalef d | 丰富度指数 | analysis.py |
| PCA | 水质变量主成分分析 | advanced_analysis.py |
| NMDS | 非度量多维尺度排序 | advanced_analysis.py |
| 层次聚类 | Ward聚类（Bray-Curtis距离） | advanced_analysis.py |
| 功能群分析 | 按门/类群汇总 | advanced_analysis.py |
| Beta多样性 | Bray-Curtis相异矩阵 | advanced_analysis.py |
| 指示物种 | IndVal分析 | advanced_analysis.py |
| K-优势度曲线 | 累积优势度分布 | advanced_analysis.py |
| Mantel Test | 群落-环境距离相关性 | mantel_simper.py |
| SIMPER | 物种差异贡献分解 | mantel_simper.py |

## 快速开始

### 环境要求

- Python 3.10+
- pip 包管理器

### 安装依赖

```bash
pip install -r requirements.txt
```

### 运行分析

所有脚本从项目根目录运行：

```bash
# 1. 基础多样性分析
python scripts/analysis.py

# 2. 高级生态分析（生成7张图）
python scripts/advanced_analysis.py

# 3. Mantel检验 + SIMPER分析
python scripts/mantel_simper.py

# 4. 生成综合Word报告
python scripts/gen_report.py
```

### 输出文件

- `output/results/analysis_result.txt` — 多样性指数表 + 物种多度排名
- `output/results/mantel_simper_result.txt` — Mantel检验 + SIMPER详细结果
- `output/figures/*.png` — 7张分析图表
- `output/report/浮游生物群落生态调查报告.docx` — 综合报告

## 主要发现

### 环境梯度
- PCA 前两个主成分累计解释 **83.8%** 的水质方差
- PC1（61.6%）：离子浓度轴（TDS + 盐度 + 电导率）
- PC2（22.2%）：水化学-温度轴（pH + 温度）

### 群落特征
- 绿藻门占绝对优势（**39.7%**）—— 盘星藻、小球藻为优势种
- 梦川物种最丰富（39种），后山水池个体最多（283）
- NMDS Stress = 0.077，群落结构拟合良好

### 群落驱动力
- **时间距离**是群落差异的最强解释因素（Mantel r = 0.402，p < 0.001）
- 空间异质性效应不显著（p = 0.670）
- TDS 和电导率对群落有次要但显著的影响（p < 0.05）

### 关键物种
- **盘星藻**驱动时间维度的群落演替（SIMPER 贡献率 14.4%）
- **剑水蚤**驱动后山水池的空间异质性（SIMPER 贡献率 12-15%）
- **轮虫**是梦川的完美指示种（IndVal = 1.000）

## 数据说明

### 模糊计数映射
野外镜检中部分物种密度过高无法精确计数，使用模糊描述：

| 原始表述 | 映射值 | 依据 |
|----------|--------|------|
| 若干 | 3 | 视野内少量（约1-5个体） |
| 很多 | 8 | 视野内大量（难以逐个数） |
| 大量 | 15 | 视野内密集（含"若干（大量）"） |

> ⚠️ 此为近似处理，正式引用时应说明映射规则，建议辅以敏感性分析。

### 数据清洗
- 日期标准化：`第一次5.16` → `2025-5.16`
- 物种-数量拆分：正则匹配 `物种名+数字`
- 录入错误修正：`中心刚硅藻` → `中心纲硅藻`
- 异常标注清理：去除 `（疑似）` 前缀
- 单位剥离：水质表头统一

## 引用
如需引用本分析结果，请说明数据来源和分析方法，并注明模糊计数映射方案及其不确定性。
