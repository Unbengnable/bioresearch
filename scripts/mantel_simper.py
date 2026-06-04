# -*- coding: utf-8 -*-
"""
Mantel Test & SIMPER 分析
=========================
1. Mantel Test — 检验群落相异矩阵与环境距离矩阵的相关性
2. SIMPER — 识别对采样点间差异贡献最大的物种
"""
import math
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr

# ============================================================
# 0. 数据加载（复用 advanced_analysis 的逻辑）
# ============================================================
import os
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE = os.path.join(PROJECT_ROOT, "data", "data.txt")
WATER_FILE = os.path.join(PROJECT_ROOT, "data", "water.txt")
COUNT_MAP = {"若干": 27, "很多": 8, "大量": 15}

# --- 生物数据 ---
bio_rows = []
with open(DATA_FILE, "r", encoding="utf-8") as f:
    lines = f.readlines()
for line in lines[1:]:
    parts = line.strip().split("\t")
    if len(parts) == 4:
        d, l, sp, cnt = parts
        n = int(cnt) if cnt.isdigit() else COUNT_MAP.get(cnt, 1)
        bio_rows.append((d, l, sp, n))

# --- 水质数据 ---
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

SAMPLES = sorted(set((r[0], r[1]) for r in bio_rows), key=lambda x: (x[0], x[1]))
SAMPLE_IDS = [f"{d}_{l}" for d, l in SAMPLES]
N = len(SAMPLES)

all_species = sorted(set(r[2] for r in bio_rows))
M = len(all_species)
sp_idx = {sp: i for i, sp in enumerate(all_species)}

species_mat = np.zeros((N, M))
for date, loc, sp, n in bio_rows:
    si = SAMPLE_IDS.index(f"{date}_{loc}")
    species_mat[si, sp_idx[sp]] += n

env_mat = np.zeros((N, 5))
env_vars = ['TDS', '盐度', 'pH', '电导率', '温度']
for date, loc, tds, sal, ph, cond, temp in water_rows:
    key = f"{date}_{loc}"
    if key in SAMPLE_IDS:
        si = SAMPLE_IDS.index(key)
        env_mat[si] = [tds, sal, ph, cond, temp]

# 按采样点分组（聚合三次调查）
LOCATIONS = sorted(set(l for _, l in SAMPLES))
loc_groups = {l: [] for l in LOCATIONS}
for i, (d, l) in enumerate(SAMPLES):
    loc_groups[l].append(i)

# ============================================================
# 1. Bray-Curtis 群落相异矩阵
# ============================================================
def bray_curtis(x, y):
    s = np.sum(np.abs(x - y))
    t = np.sum(x + y)
    return s / t if t > 0 else 0.0

bc_mat = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        bc_mat[i, j] = bray_curtis(species_mat[i], species_mat[j])
bc_condensed = squareform(bc_mat)

# ============================================================
# 2. 环境距离矩阵（Euclidean，标准化后）
# ============================================================
env_std = (env_mat - env_mat.mean(axis=0)) / env_mat.std(axis=0, ddof=1)
env_dist_mat = squareform(pdist(env_std, metric='euclidean'))
env_dist_condensed = squareform(env_dist_mat)  # 转为 condensed 形式供 mantel_test

print("=" * 70)
print("一、Mantel Test — 群落相异度 vs 环境距离")
print("=" * 70)

def mantel_test(dist1, dist2, n_perm=9999, seed=42):
    """
    Mantel 检验：两个距离矩阵之间的相关性
    dist1, dist2: 压缩形式的下三角矩阵 (n*(n-1)/2,)
    n_perm: 置换次数
    返回: (r, p值, 置换分布)
    """
    # 观测到的 Pearson 相关系数
    r_obs, _ = pearsonr(dist1, dist2)

    # 置换检验
    rng = np.random.default_rng(seed)
    count = 0
    r_perm = np.zeros(n_perm)

    for k in range(n_perm):
        # 随机置换 dist2 的标签
        perm_idx = rng.permutation(len(dist2))
        rk, _ = pearsonr(dist1, dist2[perm_idx])
        r_perm[k] = rk
        if rk >= r_obs:
            count += 1

    p_value = (count + 1) / (n_perm + 1)
    return r_obs, p_value, r_perm

r_mantel, p_mantel, r_perm = mantel_test(bc_condensed, env_dist_condensed, n_perm=9999)

print(f"\n  全模型 Mantel Test:")
print(f"    r = {r_mantel:+.4f}")
print(f"    p = {p_mantel:.4f} ({'***' if p_mantel<0.001 else ('**' if p_mantel<0.01 else ('*' if p_mantel<0.05 else 'ns'))})")
print(f"    置换次数: 9999")
print(f"    置换分布: mean={r_perm.mean():.4f}, std={r_perm.std():.4f}")
print(f"    95%CI = [{np.percentile(r_perm, 2.5):.4f}, {np.percentile(r_perm, 97.5):.4f}]")

# ---- 3.1 偏 Mantel Test（分别控制每个环境变量） ----
print(f"\n  环境变量逐一检验 (Partial Mantel by single variable):")

env_dist_single = {}
for k in range(5):
    single = pdist(env_std[:, k:k+1], metric='euclidean')  # 已为 condensed 形式
    env_dist_single[env_vars[k]] = single

for var in env_vars:
    d_env = env_dist_single[var]
    r, p, _ = mantel_test(bc_condensed, d_env, n_perm=9999, seed=12)
    sig = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else 'ns'))
    print(f"    {var:6s}: r={r:+.4f}, p={p:.4f} {sig}")

# ---- 3.2 时间 Mantel（仅时间距离）----
print(f"\n  时间距离 Mantel Test:")
time_dist = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        # 日期差（按批次）
        di = ['2026-5.16', '2026-5.23', '2026-5.31'].index(SAMPLES[i][0])
        dj = ['2026-5.16', '2026-5.23', '2026-5.31'].index(SAMPLES[j][0])
        time_dist[i, j] = abs(di - dj)
time_condensed = squareform(time_dist)

r_time, p_time, _ = mantel_test(bc_condensed, time_condensed, n_perm=9999, seed=24)
sig = '***' if p_time<0.001 else ('**' if p_time<0.01 else ('*' if p_time<0.05 else 'ns'))
print(f"    r={r_time:+.4f}, p={p_time:.4f} {sig}")

# ---- 3.3 空间距离 Mantel ----
print(f"\n  空间距离 Mantel Test（若有 GPS 坐标可替换为真实距离，当前用二值指示）:")
# 没有真实坐标，用"是否同采样点"二值矩阵
spatial_dist = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        spatial_dist[i, j] = 0.0 if SAMPLES[i][1] == SAMPLES[j][1] else 1.0
spatial_condensed = squareform(spatial_dist)

r_space, p_space, _ = mantel_test(bc_condensed, spatial_condensed, n_perm=9999, seed=36)
sig = '***' if p_space<0.001 else ('**' if p_space<0.01 else ('*' if p_space<0.05 else 'ns'))
print(f"    r={r_space:+.4f}, p={p_space:.4f} {sig}")

# ============================================================
# 4. SIMPER 分析
# ============================================================
print("\n" + "=" * 70)
print("二、SIMPER — 采样点间差异贡献分析")
print("=" * 70)

def simper_between_groups(mat_a, mat_b, species_names, max_species=10):
    """
    SIMPER: 计算两组样本间 Bray-Curtis 差异中各物种的贡献百分比

    参数:
      mat_a: 组A的物种矩阵 (n_a × M)
      mat_b: 组B的物种矩阵 (n_b × M)
      species_names: 物种名列表
      max_species: 返回贡献最大的前N个物种

    返回: [(物种名, 平均贡献, 贡献百分比, 累积百分比), ...]
    """
    M = mat_a.shape[1]
    group_a_mean = mat_a.mean(axis=0)
    group_b_mean = mat_b.mean(axis=0)

    # 两组间的平均 Bray-Curtis 差异
    n_a, n_b = mat_a.shape[0], mat_b.shape[0]
    total_bc = 0.0
    pairwise = []
    for i in range(n_a):
        for j in range(n_b):
            bc = bray_curtis(mat_a[i], mat_b[j])
            pairwise.append(bc)
            total_bc += bc
    avg_bc = total_bc / (n_a * n_b)

    # 每个物种的贡献
    species_contrib = []
    for j in range(M):
        # 物种 j 对平均 Bray-Curtis 差异的贡献
        avg_abund_diff = abs(group_a_mean[j] - group_b_mean[j])
        species_contrib.append((j, avg_abund_diff))

    # 按贡献从大到小排序
    species_contrib.sort(key=lambda x: -x[1])

    total_diff = sum(sc[1] for sc in species_contrib)

    results = []
    cum_pct = 0
    for k, (j, contrib) in enumerate(species_contrib):
        pct = contrib / total_diff * 100 if total_diff > 0 else 0
        cum_pct += pct
        results.append((species_names[j], contrib, pct, cum_pct,
                        group_a_mean[j], group_b_mean[j]))
        if k >= max_species - 1 and cum_pct > 70:
            break

    return avg_bc, results

# ---- 对各采样点对做 SIMPER ----
print("\n  采样点两两比较（各三次调查聚合）:\n")

# 聚合每个采样点三次调查的平均
loc_mats = {}
for loc in LOCATIONS:
    indices = loc_groups[loc]
    loc_mats[loc] = species_mat[indices]

for ia in range(len(LOCATIONS)):
    for ib in range(ia + 1, len(LOCATIONS)):
        la, lb = LOCATIONS[ia], LOCATIONS[ib]
        avg_bc, contribs = simper_between_groups(
            loc_mats[la], loc_mats[lb], all_species, max_species=12
        )
        print(f"  {la:6s} vs {lb:6s}  ──  平均 Bray-Curtis = {avg_bc:.3f}")
        print(f"  {'─' * 50}")
        print(f"  {'物种':16s}  {'贡献值':>8s}  {'贡献%':>7s}  {'累积%':>7s}  {'平均多度('+la+')':>12s}  {'平均多度('+lb+')':>12s}")
        print(f"  {'─' * 50}")
        for sp, contrib, pct, cum, a_mean, b_mean in contribs[:10]:
            print(f"  {sp:16s}  {contrib:8.3f}  {pct:6.1f}%  {cum:6.1f}%  {a_mean:12.1f}  {b_mean:12.1f}")
        print()

# ---- 按日期分组 SIMPER ----
print("\n" + "─" * 70)
print("  时间分组 SIMPER（三次调查两两比较，聚合所有采样点）:\n")

dates_sorted = sorted(set(d for d, _ in SAMPLES))
for ia in range(len(dates_sorted)):
    for ib in range(ia + 1, len(dates_sorted)):
        da, db = dates_sorted[ia], dates_sorted[ib]
        indices_a = [i for i, (d, l) in enumerate(SAMPLES) if d == da]
        indices_b = [i for i, (d, l) in enumerate(SAMPLES) if d == db]
        mat_a = species_mat[indices_a]
        mat_b = species_mat[indices_b]

        avg_bc, contribs = simper_between_groups(mat_a, mat_b, all_species, max_species=10)
        print(f"  {da} vs {db}  ──  平均 Bray-Curtis = {avg_bc:.3f}")
        print(f"  {'─' * 50}")
        print(f"  {'物种':16s}  {'贡献值':>8s}  {'贡献%':>7s}  {'累积%':>7s}")
        print(f"  {'─' * 50}")
        for sp, contrib, pct, cum, _, _ in contribs[:8]:
            print(f"  {sp:16s}  {contrib:8.3f}  {pct:6.1f}%  {cum:6.1f}%")
        print()

# ============================================================
# 5. 保存完整结果到文件
# ============================================================
out_file = os.path.join(PROJECT_ROOT, "output", "results", "mantel_simper_result.txt")
with open(out_file, "w", encoding="utf-8") as f:
    f.write("=" * 70 + "\n")
    f.write("Mantel Test & SIMPER 分析结果\n")
    f.write("=" * 70 + "\n\n")

    # Mantel
    f.write("一、Mantel Test\n\n")
    f.write(f"全模型: r={r_mantel:.4f}, p={p_mantel:.4f}\n")
    f.write(f"置换分布 (n=9999): mean={r_perm.mean():.4f}, std={r_perm.std():.4f}\n\n")

    f.write("单变量:\n")
    for var in env_vars:
        d_env = env_dist_single[var]
        r, p, _ = mantel_test(bc_condensed, d_env, n_perm=9999, seed=12)
        f.write(f"  {var}: r={r:+.4f}, p={p:.4f}\n")

    f.write(f"\n时间距离: r={r_time:+.4f}, p={p_time:.4f}\n")
    f.write(f"空间距离: r={r_space:+.4f}, p={p_space:.4f}\n")

    # SIMPER
    f.write("\n\n二、SIMPER 分析\n\n")
    for ia in range(len(LOCATIONS)):
        for ib in range(ia + 1, len(LOCATIONS)):
            la, lb = LOCATIONS[ia], LOCATIONS[ib]
            avg_bc, contribs = simper_between_groups(
                loc_mats[la], loc_mats[lb], all_species, max_species=15
            )
            f.write(f"{la} vs {lb} (Bray-Curtis={avg_bc:.3f})\n")
            f.write(f"{'物种':20s} {'贡献%':>7s} {'累积%':>7s}\n")
            for sp, _, pct, cum, _, _ in contribs:
                f.write(f"{sp:20s} {pct:6.1f}% {cum:6.1f}%\n")
            f.write("\n")

    # 时间分组
    f.write("\n时间分组:\n\n")
    for ia in range(len(dates_sorted)):
        for ib in range(ia + 1, len(dates_sorted)):
            da, db = dates_sorted[ia], dates_sorted[ib]
            indices_a = [i for i, (d, l) in enumerate(SAMPLES) if d == da]
            indices_b = [i for i, (d, l) in enumerate(SAMPLES) if d == db]
            avg_bc, contribs = simper_between_groups(
                species_mat[indices_a], species_mat[indices_b], all_species, max_species=15
            )
            f.write(f"{da} vs {db} (Bray-Curtis={avg_bc:.3f})\n")
            f.write(f"{'物种':20s} {'贡献%':>7s} {'累积%':>7s}\n")
            for sp, _, pct, cum, _, _ in contribs:
                f.write(f"{sp:20s} {pct:6.1f}% {cum:6.1f}%\n")
            f.write("\n")

print(f"\n完整结果已保存至: {out_file}")
print("=" * 70)
