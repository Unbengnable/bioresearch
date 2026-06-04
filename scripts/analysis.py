# -*- coding: utf-8 -*-
"""
生物多样性分析：多度、丰富度、多样性指数
基于清洗后的 data.txt
"""
import math
from collections import defaultdict

# ============================================================
# 1. 读取数据
# ============================================================
import os
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_file = os.path.join(PROJECT_ROOT, "data", "data.txt")
rows = []
with open(data_file, "r", encoding="utf-8") as f:
    lines = f.readlines()

header = lines[0].strip().split("\t")
for line in lines[1:]:
    parts = line.strip().split("\t")
    if len(parts) == 4:
        rows.append(parts)  # [date, location, species, count]

# ============================================================
# 2. 模糊计数 → 数值映射
# ============================================================
COUNT_MAP = {
    "若干": 27,
    "很多": 8,
    "大量": 15,
}

def to_numeric(count_str):
    """将数量转为数值，非数字按映射表转换"""
    if count_str.isdigit():
        return int(count_str)
    return COUNT_MAP.get(count_str, 1)  # 未知类型按 1 计

# 按 (日期, 采样点) 分组
groups = defaultdict(list)
for date, loc, species, count in rows:
    n = to_numeric(count)
    groups[(date, loc)].append((species, count, n))

# ============================================================
# 3. 计算各指数
# ============================================================
def shannon(abundances):
    """Shannon-Wiener 多样性指数 H' = -Σ(pi * ln(pi))"""
    total = sum(abundances)
    if total == 0:
        return 0.0
    h = 0.0
    for n in abundances:
        if n > 0:
            p = n / total
            h -= p * math.log(p)
    return h

def simpson(abundances):
    """Simpson 多样性指数 D = 1 - Σ(pi²)"""
    total = sum(abundances)
    if total == 0:
        return 0.0
    d = sum((n / total) ** 2 for n in abundances)
    return 1.0 - d

def simpson_inverse(abundances):
    """Simpson 倒数指数 1/D = 1 / Σ(pi²)"""
    total = sum(abundances)
    if total == 0:
        return 0.0
    d = sum((n / total) ** 2 for n in abundances)
    if d == 0:
        return float("inf")
    return 1.0 / d

def pielou(shannon_h, richness):
    """Pielou 均匀度指数 J' = H' / ln(S)"""
    if richness <= 1:
        return 1.0
    return shannon_h / math.log(richness)

def margalef(richness, total_abundance):
    """Margalef 丰富度指数 d = (S - 1) / ln(N)"""
    if total_abundance <= 1:
        return 0.0
    return (richness - 1) / math.log(total_abundance)

# ============================================================
# 4. 逐样点逐日期计算
# ============================================================
print("=" * 90)
print("  采样日期      采样点      丰富度S  总多度N  Shannon H'  Simpson D  Pielou J'  Margalef")
print("=" * 90)

# 排序
sorted_keys = sorted(groups.keys(), key=lambda x: (x[0], x[1]))

summary = []
for (date, loc) in sorted_keys:
    entries = groups[(date, loc)]
    richness = len(set(e[0] for e in entries))  # 独特物种数
    total_n = sum(e[2] for e in entries)         # 总个体数（映射后）
    abundances = [e[2] for e in entries]

    h = shannon(abundances)
    d = simpson(abundances)
    j = pielou(h, richness)
    mg = margalef(richness, total_n)

    summary.append((date, loc, richness, total_n, h, d, j, mg))
    print(f"  {date:12s} {loc:6s}  {richness:6d}  {total_n:8d}  {h:10.4f}  {d:9.4f}  {j:9.4f}  {mg:8.4f}")

# ============================================================
# 5. 按采样点汇总（三次调查合并）
# ============================================================
print()
print("=" * 90)
print("  【按采样点汇总】")
print("=" * 90)

loc_aggregate = defaultdict(lambda: {"species": set(), "total_n": 0, "samples": 0})
for (date, loc), entries in groups.items():
    agg = loc_aggregate[loc]
    agg["samples"] += 1
    for species, _, n in entries:
        agg["species"].add(species)
        agg["total_n"] += n

for loc in ["黎照湖", "香雪海", "梦川", "后山水池", "菜根谭"]:
    agg = loc_aggregate[loc]
    s = len(agg["species"])
    n = agg["total_n"]
    print(f"  {loc:6s}  总物种数={s:3d}  总个体数(映射)={n:5d}  平均每次物种数={s/3:.1f}")

# ============================================================
# 6. 按日期汇总（五个采样点合并）
# ============================================================
print()
print("=" * 90)
print("  【按调查日期汇总】")
print("=" * 90)

date_aggregate = defaultdict(lambda: {"species": set(), "total_n": 0})
for (date, loc), entries in groups.items():
    agg = date_aggregate[date]
    for species, _, n in entries:
        agg["species"].add(species)
        agg["total_n"] += n

for date in sorted(date_aggregate.keys()):
    agg = date_aggregate[date]
    print(f"  {date}  总物种数={len(agg['species']):3d}  总个体数(映射)={agg['total_n']:5d}")

# ============================================================
# 7. 多度排名（全部调查累计）
# ============================================================
print()
print("=" * 90)
print("  【物种多度排名 Top 20（按映射后总数量）】")
print("=" * 90)

species_total = defaultdict(int)
for (date, loc), entries in groups.items():
    for species, _, n in entries:
        species_total[species] += n

ranked = sorted(species_total.items(), key=lambda x: -x[1])
for i, (sp, n) in enumerate(ranked[:20], 1):
    print(f"  {i:2d}. {sp:16s}  {n:5d}")

# ============================================================
# 8. 保存结果到文件
# ============================================================
out_file = os.path.join(PROJECT_ROOT, "output", "results", "analysis_result.txt")
with open(out_file, "w", encoding="utf-8") as f:
    f.write("=" * 90 + "\n")
    f.write("生物多样性分析结果\n")
    f.write("注：模糊计数映射规则：若干→27, 很多→8, 大量→15\n")
    f.write("=" * 90 + "\n\n")

    # 逐样点表
    f.write("采样日期\t采样点\t丰富度S\t总多度N\tShannon H'\tSimpson D\tPielou J'\tMargalef d\n")
    for row in summary:
        f.write("\t".join(str(x) for x in row) + "\n")

    f.write("\n物种多度排名\n")
    f.write("物种\t映射总数量\n")
    for sp, n in ranked:
        f.write(f"{sp}\t{n}\n")

print(f"\n结果已保存至: {out_file}")
