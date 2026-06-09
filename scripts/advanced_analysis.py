# -*- coding: utf-8 -*-
"""
高级生态学定量分析
====================
基于 5 采样点 × 3 次调查 的水质+生物数据

分析方法：
  1. PCA — 水质变量降维 + 与多样性指数的关联
  2. 度量多维排序（Sammon映射）— 基于物种组成 Bray-Curtis 距离
  3. 聚类分析 — 样点物种组成的层次聚类
  4. 相关性分析 — 多样性指数 vs 环境因子矩阵
  5. 功能群分析 — 按门/类群汇总
  6. Beta 多样性 — Bray-Curtis 相异矩阵
  7. 指示物种分析 — 各采样点特征种
"""
import math
import warnings
from collections import defaultdict
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 配置中文字体（Windows 系统）
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'KaiTi', 'FangSong']
matplotlib.rcParams['axes.unicode_minus'] = False  # 修复负号显示

warnings.filterwarnings('ignore')

# ============================================================
# 0. 数据读取
# ============================================================
import os
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE = os.path.join(PROJECT_ROOT, "data", "data.txt")
WATER_FILE = os.path.join(PROJECT_ROOT, "data", "water.txt")
OUT_DIR = os.path.join(PROJECT_ROOT, "output", "figures")

# 模糊计数映射
COUNT_MAP = {"若干": 27, "很多": 8, "大量": 15}

# --- 读取生物数据 ---
bio_rows = []
with open(DATA_FILE, "r", encoding="utf-8") as f:
    lines = f.readlines()
for line in lines[1:]:
    parts = line.strip().split("\t")
    if len(parts) == 4:
        date, loc, sp, cnt = parts
        n = int(cnt) if cnt.isdigit() else COUNT_MAP.get(cnt, 1)
        bio_rows.append((date, loc, sp, n))

# --- 读取水质数据 ---
water_rows = []
with open(WATER_FILE, "r", encoding="utf-8") as f:
    lines = f.readlines()
for line in lines[1:]:
    parts = line.strip().split("\t")
    if len(parts) == 7:
        water_rows.append((parts[0], parts[1],
                           float(parts[2]), float(parts[3]),
                           float(parts[4]), float(parts[5]),
                           float(parts[6])))

# 样本标签 (15个)
SAMPLES = sorted(set((r[0], r[1]) for r in bio_rows), key=lambda x: (x[0], x[1]))
SAMPLE_LABELS = [f"{d}\n{l}" for d, l in SAMPLES]
SAMPLE_IDS = [f"{d}_{l}" for d, l in SAMPLES]
N = len(SAMPLES)  # 15

# 构建物种矩阵: (N samples) × (M species)
all_species = sorted(set(r[2] for r in bio_rows))
M = len(all_species)
sp_idx = {sp: i for i, sp in enumerate(all_species)}

species_mat = np.zeros((N, M))
for date, loc, sp, n in bio_rows:
    sample_key = (date, loc)
    si = SAMPLE_IDS.index(f"{date}_{loc}")
    sj = sp_idx[sp]
    species_mat[si, sj] += n

# 构建水质矩阵
env_vars = ['TDS', '盐度', 'pH', '电导率', '温度']
env_mat = np.zeros((N, 5))
for date, loc, tds, sal, ph, cond, temp in water_rows:
    key = f"{date}_{loc}"
    if key in SAMPLE_IDS:
        si = SAMPLE_IDS.index(key)
        env_mat[si] = [tds, sal, ph, cond, temp]

# 多样性指数
def calc_indices(ab_vec):
    nz = ab_vec[ab_vec > 0]
    S = len(nz)
    N_total = nz.sum()
    H = 0.0
    if N_total > 0:
        H = -sum((x/N_total) * math.log(x/N_total) for x in nz)
    D = 0.0
    if N_total > 0:
        D = 1.0 - sum((x/N_total)**2 for x in nz)
    J = H / math.log(S) if S > 1 else 1.0
    Mg = (S-1) / math.log(N_total) if N_total > 1 else 0.0
    return S, int(N_total), H, D, J, Mg

div_indices = np.array([calc_indices(species_mat[i]) for i in range(N)])
div_labels = ['丰富度S', '总多度N', "Shannon H'", 'Simpson D', "Pielou J'", 'Margalef d']

# ============================================================
# 1. PCA 主成分分析（水质变量）
# ============================================================
print("=" * 60)
print("1. PCA 主成分分析")
print("=" * 60)

env_std = (env_mat - env_mat.mean(axis=0)) / env_mat.std(axis=0, ddof=1)
cov = np.cov(env_std.T)
eigvals, eigvecs = eigh(cov)
# 从大到小排序
order = np.argsort(-eigvals)
eigvals = eigvals[order]
eigvecs = eigvecs[:, order]

# 各主成分方差解释
var_exp = eigvals / eigvals.sum() * 100
cum_var = np.cumsum(var_exp)
print(f"  特征值: {eigvals}")
print(f"  方差解释率(%): {', '.join(f'{v:.1f}' for v in var_exp)}")
print(f"  累计方差解释率(%): {', '.join(f'{v:.1f}' for v in cum_var)}")

# PC 载荷
pc_scores = env_std @ eigvecs

print("\n  主成分载荷 (PC1-PC3):")
for i, var in enumerate(env_vars):
    loadings = ", ".join(f"PC{j+1}={eigvecs[i, j]:+.3f}" for j in range(min(3, len(eigvecs))))
    print(f"    {var:6s}: {loadings}")

# PC1-PC2 与多样性指数的相关性
print("\n  PC1 与多样性指数的 Pearson 相关:")
for j, label in enumerate(div_labels):
    r, p = stats.pearsonr(pc_scores[:, 0], div_indices[:, j])
    sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    print(f"    {label:12s}: r={r:+.3f}, p={p:.3f} {sig}")

# ---- PCA 图 ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 图1: PCA 样本得分 biplot
ax = axes[0]
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
loc_list = sorted(set(l[1] for l in SAMPLES))
loc_colors = {l: colors[i % len(colors)] for i, l in enumerate(loc_list)}
markers = {'2026-5.16': 'o', '2026-5.23': 's', '2026-5.31': '^'}

for i, (date, loc) in enumerate(SAMPLES):
    ax.scatter(pc_scores[i, 0], pc_scores[i, 1],
               c=loc_colors[loc], marker=markers[date],
               s=120, edgecolors='k', linewidths=0.5, zorder=5)
    ax.annotate(f"{loc[:2]}", (pc_scores[i, 0], pc_scores[i, 1]),
                textcoords="offset points", xytext=(5, 5), fontsize=8)

# 画载荷向量
for j, var in enumerate(env_vars):
    ax.arrow(0, 0, eigvecs[j, 0] * 3.5, eigvecs[j, 1] * 3.5,
             head_width=0.08, head_length=0.12, fc='grey', ec='grey', alpha=0.7)
    ax.text(eigvecs[j, 0] * 3.8, eigvecs[j, 1] * 3.8, var, fontsize=9, color='grey')

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)', fontsize=11)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)', fontsize=11)
ax.set_title('PCA: 水质变量排序', fontsize=13)
ax.axhline(0, color='grey', lw=0.5)
ax.axvline(0, color='grey', lw=0.5)

# 图例
from matplotlib.lines import Line2D
l_leg = [Line2D([0], [0], marker='o', color='w', markerfacecolor=loc_colors[l],
                markersize=8, label=l) for l in loc_list]
d_leg = [Line2D([0], [0], marker=m, color='w', markerfacecolor='grey',
                markersize=8, label=d) for d, m in markers.items()]
leg1 = ax.legend(handles=l_leg, title='采样点', loc='upper left', fontsize=8)
ax.add_artist(leg1)
ax.legend(handles=d_leg, title='日期', loc='lower left', fontsize=8)

# 图2: PC1 vs Shannon H'
ax = axes[1]
for i, (date, loc) in enumerate(SAMPLES):
    ax.scatter(pc_scores[i, 0], div_indices[i, 2],
               c=loc_colors[loc], marker=markers[date],
               s=120, edgecolors='k', linewidths=0.5)
    ax.annotate(loc[:2], (pc_scores[i, 0], div_indices[i, 2]),
                textcoords="offset points", xytext=(5, 5), fontsize=8)

r_val, p_val = stats.pearsonr(pc_scores[:, 0], div_indices[:, 2])
ax.set_xlabel('PC1 (水质第一主成分)', fontsize=11)
ax.set_ylabel("Shannon H'", fontsize=11)
ax.set_title(f"PC1 vs Shannon 多样性 (r={r_val:+.3f}, p={p_val:.3f})", fontsize=13)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/pca_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: pca_analysis.png")

# ============================================================
# 2. 度量多维排序分析（Sammon 映射）
# ============================================================
print("\n" + "=" * 60)
print("2. 度量多维排序（Sammon 映射）")
print("=" * 60)

def bray_curtis(x, y):
    """Bray-Curtis 相异度"""
    return np.sum(np.abs(x - y)) / np.sum(x + y) if np.sum(x + y) > 0 else 0

# Bray-Curtis 距离矩阵
bc_dist = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        bc_dist[i, j] = bray_curtis(species_mat[i], species_mat[j])

print(f"  Bray-Curtis 相异度矩阵 ({N}×{N}):")
print(f"    均值={bc_dist[np.triu_indices(N,1)].mean():.3f}")
print(f"    范围=[{bc_dist[np.triu_indices(N,1)].min():.3f}, {bc_dist[np.triu_indices(N,1)].max():.3f}]")

# Sammon 映射（加权度量MDS，PCoA初始化 + 加权梯度下降）
def simple_nmds(dist_mat, dims=2, max_iter=200, eps=1e-6):
    """Sammon 映射: PCoA 初始化 + 加权梯度下降 (w=1/d)"""
    n = dist_mat.shape[0]
    # 初始化: 用 PCoA (即经典MDS)
    d2 = dist_mat ** 2
    j_mat = np.eye(n) - np.ones((n, n)) / n
    b_mat = -0.5 * j_mat @ d2 @ j_mat
    eigvals, eigvecs = eigh(b_mat)
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order][:dims]
    eigvecs = eigvecs[:, order][:, :dims]
    # 取正特征值
    pos_mask = eigvals > 0
    init = np.zeros((n, dims))
    for d in range(dims):
        if d < len(eigvals) and eigvals[d] > 0:
            init[:, d] = eigvecs[:, d] * np.sqrt(eigvals[d])

    # 梯度下降优化 Stress
    config = init.copy()
    lr = 0.01
    prev_stress = float('inf')
    for it in range(max_iter):
        # 当前配置的距离
        cfg_dist = np.sqrt(np.sum((config[:, None, :] - config[None, :, :]) ** 2, axis=2))
        cfg_dist += 1e-10

        # 梯度
        grad = np.zeros_like(config)
        stress = 0
        total = 0
        for i in range(n):
            for j in range(i+1, n):
                if dist_mat[i, j] > 0:
                    diff = cfg_dist[i, j] - dist_mat[i, j]
                    weight = 1.0 / dist_mat[i, j]  # 权重: 距离远的点对惩罚更轻
                    grad_ij = diff * weight * (config[i] - config[j]) / cfg_dist[i, j]
                    grad[i] += grad_ij
                    grad[j] -= grad_ij
                    stress += diff**2 * weight
                    total += dist_mat[i, j]

        stress = stress / total if total > 0 else 0

        if abs(prev_stress - stress) < eps:
            break
        prev_stress = stress
        config -= lr * grad

    return config, stress

nmds_coords, nmds_stress = simple_nmds(bc_dist)
print(f"  度量多维排序 Stress = {nmds_stress:.4f}")

# Sammon 排序图
fig, ax = plt.subplots(figsize=(8, 7))
for i, (date, loc) in enumerate(SAMPLES):
    ax.scatter(nmds_coords[i, 0], nmds_coords[i, 1],
               c=loc_colors[loc], marker=markers[date],
               s=150, edgecolors='k', linewidths=0.5, zorder=5)
    ax.annotate(f"{date[-4:]}\n{loc}", (nmds_coords[i, 0], nmds_coords[i, 1]),
                textcoords="offset points", xytext=(6, 6), fontsize=7)

# 连接同一采样点的三次调查
for loc in loc_list:
    pts = np.array([[nmds_coords[i, 0], nmds_coords[i, 1]]
                    for i, (d, l) in enumerate(SAMPLES) if l == loc])
    if len(pts) == 3:
        ax.plot(pts[:, 0], pts[:, 1], '--', color=loc_colors[loc], alpha=0.4, lw=1)

ax.set_xlabel('MDS1', fontsize=11)
ax.set_ylabel('MDS2', fontsize=11)
ax.set_title(f'度量多维排序 (Sammon, Stress={nmds_stress:.4f})', fontsize=13)
ax.legend(handles=l_leg, title='采样点', loc='best', fontsize=8)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/nmds_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: nmds_analysis.png")

# ============================================================
# 3. 聚类分析
# ============================================================
print("\n" + "=" * 60)
print("3. 层次聚类分析")
print("=" * 60)

condensed = squareform(bc_dist)
Z = linkage(condensed, method='ward')

fig, ax = plt.subplots(figsize=(10, 5))
dn = dendrogram(Z, labels=SAMPLE_LABELS, ax=ax,
                leaf_font_size=9, color_threshold=0.5 * max(Z[:, 2]))
ax.set_title('Ward 层次聚类 (基于 Bray-Curtis 相异度)', fontsize=13)
ax.set_ylabel('距离', fontsize=11)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/cluster_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: cluster_analysis.png")

# 聚类结果
clusters = fcluster(Z, t=4, criterion='maxclust')
print(f"  4 类聚类结果:")
for i, (date, loc) in enumerate(SAMPLES):
    print(f"    {date} {loc:6s} → 类 {clusters[i]}")

# ============================================================
# 4. 相关性分析
# ============================================================
print("\n" + "=" * 60)
print("4. 多样性指数 vs 环境因子 相关性矩阵")
print("=" * 60)

all_indices = np.hstack([div_indices, env_mat])
all_labels = div_labels + env_vars

corr_mat = np.zeros((len(all_labels), len(all_labels)))
pval_mat = np.zeros((len(all_labels), len(all_labels)))

for i in range(len(all_labels)):
    for j in range(len(all_labels)):
        r, p = stats.pearsonr(all_indices[:, i], all_indices[:, j])
        corr_mat[i, j] = r
        pval_mat[i, j] = p

# 打印关键相关
print("\n  多样性指数 × 环境因子 (Pearson r):")
print(f"  {'':12s}", end="")
for ev in env_vars:
    print(f"  {ev:>6s}", end="")
print()
for di in range(len(div_labels)):
    print(f"  {div_labels[di]:12s}", end="")
    for ei in range(5):
        r = corr_mat[di, 6 + ei]
        p = pval_mat[di, 6 + ei]
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else " ")
        print(f"  {r:+5.2f}{sig}", end="")
    print()

# 相关性热力图
fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(corr_mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
ax.set_xticks(range(len(all_labels)))
ax.set_yticks(range(len(all_labels)))
ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=8)
ax.set_yticklabels(all_labels, fontsize=8)

for i in range(len(all_labels)):
    for j in range(len(all_labels)):
        color = 'white' if abs(corr_mat[i, j]) > 0.6 else 'black'
        ax.text(j, i, f'{corr_mat[i,j]:.2f}', ha='center', va='center',
                fontsize=7, color=color)

plt.colorbar(im, ax=ax, shrink=0.8)
ax.set_title('Pearson 相关矩阵 (多样性指数 + 环境因子)', fontsize=13)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/correlation_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: correlation_heatmap.png")

# ============================================================
# 5. 功能群分析
# ============================================================
print("\n" + "=" * 60)
print("5. 功能群（门/类群）组成分析")
print("=" * 60)

PHYLUM_MAP = {
    '盘星藻': '绿藻门', '栅藻': '绿藻门', '小球藻': '绿藻门', '空球藻': '绿藻门',
    '实球藻': '绿藻门', '团藻': '绿藻门', '衣藻': '绿藻门', '新月藻': '绿藻门',
    '鼓藻': '绿藻门', '角星鼓藻': '绿藻门', '美丽鼓藻': '绿藻门', '三角鼓藻': '绿藻门',
    '四胞藻': '绿藻门', '四角藻': '绿藻门', '网球藻': '绿藻门', '双星藻': '绿藻门',
    '集星藻': '绿藻门', '十字藻': '绿藻门', '微芒藻': '绿藻门', '空星藻': '绿藻门',
    '空星藻属': '绿藻门', '盘藻': '绿藻门', '星形冠盘藻': '绿藻门',
    '月牙藻': '绿藻门', '月牙藻集群': '绿藻门', '卵囊藻': '绿藻门',
    '多形裸藻': '绿藻门',  # 裸藻在有些分类系统中归绿藻门
    '水绵': '绿藻门',
    '羽文纲硅藻': '硅藻门', '中心纲硅藻': '硅藻门', '中心刚硅藻': '硅藻门',
    '针杆藻': '硅藻门', '桥弯藻': '硅藻门', '小环藻': '硅藻门',
    '小环藻（中心纲硅藻': '硅藻门', '星杆藻': '硅藻门', '楔形藻': '硅藻门',
    '等鞭金藻科': '金藻门', '金藻门锥囊藻': '金藻门', '锥囊藻': '金藻门',
    '鱼鳞藻': '金藻门', '黄丝藻': '黄藻门',
    '甲藻': '甲藻门', '角甲藻': '甲藻门', '甲藻残片': '甲藻门', '角藻': '甲藻门',
    '裸藻': '裸藻门', '扁裸藻': '裸藻门', '囊裸藻': '裸藻门',
    '梭形裸藻': '裸藻门', '梭状裸藻': '裸藻门', '袋鞭藻': '裸藻门',
    '鳞孔藻': '裸藻门',
    '颤藻': '蓝藻门', '念珠藻': '蓝藻门', '微囊藻': '蓝藻门',
    '色球藻': '蓝藻门', '平裂藻': '蓝藻门', '隐球藻': '蓝藻门',
    '卵形隐藻': '隐藻门', '隐藻': '隐藻门',
    '眼虫': '裸藻门',  # 眼虫 = Euglena → 裸藻门
    '轮虫': '轮虫动物', '线虫': '线虫动物', '水螨': '节肢动物',
    '剑水蚤': '节肢动物', '剑水蚤残骸': '节肢动物',
    '水蚤': '节肢动物', '裸腹蚤': '节肢动物', '鞭虫卵': '线虫动物',
    '草履虫': '原生动物', '表壳虫': '原生动物', '太阳虫': '原生动物',
    '肾形虫': '原生动物', '纤毛虫': '原生动物', '变形虫': '原生动物',
    '尾滴虫': '原生动物', '喇叭虫': '原生动物',
    '硅藻': '硅藻门',  # 泛称
}

def get_phylum(species):
    if species in PHYLUM_MAP:
        return PHYLUM_MAP[species]
    # 模糊匹配
    for key, val in PHYLUM_MAP.items():
        if key in species or species in key:
            return val
    return '其他'

# 按样本汇总功能群
phylum_list = sorted(set(list(PHYLUM_MAP.values()) + ['其他']))
phylums_by_sample = np.zeros((N, len(phylum_list)))
for i, (date, loc) in enumerate(SAMPLES):
    for sp, n in zip(all_species, species_mat[i]):
        if n > 0:
            ph = get_phylum(sp)
            pi = phylum_list.index(ph)
            phylums_by_sample[i, pi] += n

# 功能群堆叠柱状图
fig, ax = plt.subplots(figsize=(14, 6))
bottom = np.zeros(N)
phylum_colors = plt.cm.tab20(np.linspace(0, 1, len(phylum_list)))
for pi, ph in enumerate(phylum_list):
    vals = phylums_by_sample[:, pi]
    ax.bar(range(N), vals, bottom=bottom, label=ph, color=phylum_colors[pi], edgecolor='white', lw=0.3)
    bottom += vals

ax.set_xticks(range(N))
ax.set_xticklabels(SAMPLE_LABELS, fontsize=7)
ax.set_ylabel('映射多度', fontsize=11)
ax.set_title('各样本功能群组成（按门/类群）', fontsize=13)
ax.legend(loc='upper left', fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/functional_groups.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: functional_groups.png")

# 各门占比
total_by_phylum = phylums_by_sample.sum(axis=0)
print("\n  各功能群总多度:")
for pi, ph in enumerate(phylum_list):
    if total_by_phylum[pi] > 0:
        pct = total_by_phylum[pi] / total_by_phylum.sum() * 100
        print(f"    {ph:10s}: {total_by_phylum[pi]:6.0f} ({pct:4.1f}%)")

# ============================================================
# 6. Beta 多样性热力图
# ============================================================
print("\n" + "=" * 60)
print("6. Beta 多样性 (Bray-Curtis 相异矩阵)")
print("=" * 60)

fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(bc_dist, cmap='YlOrRd', aspect='equal')
ax.set_xticks(range(N))
ax.set_yticks(range(N))
ax.set_xticklabels(SAMPLE_LABELS, rotation=45, ha='right', fontsize=7)
ax.set_yticklabels(SAMPLE_LABELS, fontsize=7)
plt.colorbar(im, ax=ax, shrink=0.8, label='Bray-Curtis 相异度')

for i in range(N):
    for j in range(N):
        ax.text(j, i, f'{bc_dist[i,j]:.2f}', ha='center', va='center',
                fontsize=5, color='white' if bc_dist[i,j] > 0.5 else 'black')

ax.set_title('Beta 多样性: Bray-Curtis 相异矩阵', fontsize=13)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/beta_diversity.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: beta_diversity.png")

# 组内 vs 组间
same_site = []
diff_site = []
for i in range(N):
    for j in range(i+1, N):
        if SAMPLES[i][1] == SAMPLES[j][1]:
            same_site.append(bc_dist[i, j])
        else:
            diff_site.append(bc_dist[i, j])
print(f"\n  同一采样点内相异度: {np.mean(same_site):.3f} ± {np.std(same_site):.3f} (n={len(same_site)})")
print(f"  不同采样点间相异度: {np.mean(diff_site):.3f} ± {np.std(diff_site):.3f} (n={len(diff_site)})")
t, p = stats.ttest_ind(same_site, diff_site)
print(f"  t-test: t={t:.3f}, p={p:.4f} {'**' if p < 0.01 else ('*' if p < 0.05 else 'ns')}")

# ============================================================
# 7. 指示物种分析 (IndVal 简化版)
# ============================================================
print("\n" + "=" * 60)
print("7. 指示物种分析 (各采样点特征种)")
print("=" * 60)

for loc in loc_list:
    loc_indices = [i for i, (d, l) in enumerate(SAMPLES) if l == loc]
    other_indices = [i for i, (d, l) in enumerate(SAMPLES) if l != loc]

    indicators = []
    for sj, sp in enumerate(all_species):
        in_loc = species_mat[loc_indices, sj].sum()
        in_other = species_mat[other_indices, sj].sum()
        total = in_loc + in_other
        if total < 5:
            continue
        # 特异性 A = 在该样点中的多度 / 该物种总多度
        A = in_loc / total
        # 保真度 B = 出现该物种的该样点样本数 / 该样点样本总数
        presence_in = sum(1 for i in loc_indices if species_mat[i, sj] > 0)
        B = presence_in / len(loc_indices)
        # IndVal = A × B
        indval = A * B
        if indval > 0.4:
            indicators.append((sp, indval, in_loc, presence_in))

    indicators.sort(key=lambda x: -x[1])
    print(f"\n  {loc} (Top 5 指示物种):")
    for sp, iv, ab, pr in indicators[:5]:
        print(f"    {sp:16s}  IndVal={iv:.3f}  多度={ab:.0f}  出现{pr}/3次")

# ============================================================
# 8. K-优势度曲线
# ============================================================
print("\n" + "=" * 60)
print("8. K-优势度曲线")
print("=" * 60)

fig, ax = plt.subplots(figsize=(10, 6))
for i, (date, loc) in enumerate(SAMPLES):
    nz = species_mat[i][species_mat[i] > 0]
    if len(nz) == 0:
        continue
    nz_sorted = np.sort(nz)[::-1]
    cum_pct = np.cumsum(nz_sorted) / nz_sorted.sum() * 100
    rank_pct = np.arange(1, len(cum_pct)+1) / len(nz_sorted) * 100
    ax.plot(rank_pct, cum_pct, color=loc_colors[loc],
            linestyle=['-', '--', '-.'][['2026-5.16', '2026-5.23', '2026-5.31'].index(date)],
            alpha=0.7, label=f"{date[-4:]} {loc}")

ax.set_xlabel('物种累积百分比 (%)', fontsize=11)
ax.set_ylabel('累积优势度 (%)', fontsize=11)
ax.set_title('K-优势度曲线', fontsize=13)
ax.legend(fontsize=6, ncol=3, loc='lower right')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/k_dominance.png", dpi=150, bbox_inches='tight')
plt.close()
print("  → 图保存: k_dominance.png")

# ============================================================
# 汇总
# ============================================================
print("\n" + "=" * 60)
print("所有分析完成！")
print(f"共生成 7 张图:")
print("  1. pca_analysis.png        — PCA 主成分分析")
print("  2. nmds_analysis.png       — 度量多维排序（Sammon 映射）")
print("  3. cluster_analysis.png    — 层次聚类树状图")
print("  4. correlation_heatmap.png — 相关性热力图")
print("  5. functional_groups.png   — 功能群组成")
print("  6. beta_diversity.png      — Beta 多样性热力图")
print("  7. k_dominance.png         — K-优势度曲线")
print("=" * 60)
