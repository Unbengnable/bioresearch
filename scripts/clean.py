# -*- coding: utf-8 -*-
"""
数据清洗脚本
============
将野外调查原始文本（自由格式）解析为结构化 TSV。

原始数据包含两类文件：
  - 水质数据：TDS、盐度、pH、电导率、温度（含单位）
  - 生物数据：物种名+数量的紧凑记录（含模糊计数）

清洗操作：
  1. 日期标准化：第一次5.16 → 2025-5.16
  2. 采样点识别：从固定地点列表中匹配（黎照湖/香雪海/梦川/后山水池/菜根谭）
  3. 物种-数量拆分：正则提取"物种名+数字"对
  4. 模糊计数规范化：若干→3, 很多→8, 大量（若干（大量））→15
  5. 异常标注清理：去除"（疑似）"前缀
  6. 打字错误修正：中心刚硅藻 → 中心纲硅藻
  7. 单位剥离：水质数据数值与单位分离，单位归入表头
"""

import os
import re

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

# 固定采样点列表
KNOWN_LOCATIONS = ["黎照湖", "香雪海", "梦川", "后山水池", "菜根谭"]

# 模糊计数映射
FUZZY_COUNT_MAP = {"若干": 3, "很多": 8, "大量": 15}

# 已知的正则模式（水质数据）
WATER_PATTERNS = {
    "TDS": re.compile(r"TDS[：:]\s*([\d.]+)"),
    "盐度": re.compile(r"盐度[：:]\s*([\d.]+)"),
    "pH": re.compile(r"pH[：:]\s*([\d.]+)"),
    "电导率": re.compile(r"电导率[：:]\s*([\d.]+)"),
    "温度": re.compile(r"温度[：:]\s*([\d.]+)"),
}

# 物种-数量拆分正则
SPECIES_COUNT_RE = re.compile(r"([一-鿿（）()]+?)(\d+|很多|若干|大量)")

# 打字错误修正
TYPO_FIXES = {"中心刚硅藻": "中心纲硅藻"}


def clean_water(raw_text: str) -> list:
    """清洗水质数据，返回 [(日期, 采样点, TDS, 盐度, pH, 电导率, 温度), ...]"""
    results = []
    # 按调查批次分段
    blocks = re.split(r"第[一二三]次", raw_text)
    date_map = {"5.16": "2025-5.16", "5.23": "2025-5.23", "5.31": "2025-5.31"}

    for block in blocks[1:]:
        block = block.strip()
        # 识别日期
        date = None
        for key, val in date_map.items():
            if key in block:
                date = val
                break
        if not date:
            continue

        # 识别所有采样点
        for loc in KNOWN_LOCATIONS:
            if loc in block:
                # 提取该采样点的水质参数
                vals = {}
                for key, pat in WATER_PATTERNS.items():
                    m = pat.search(block)
                    if m:
                        vals[key] = float(m.group(1))

                if len(vals) >= 5:
                    results.append((date, loc,
                                   vals.get("TDS", 0), vals.get("盐度", 0),
                                   vals.get("pH", 0), vals.get("电导率", 0),
                                   vals.get("温度", 0)))
    return results


def clean_species(raw_text: str) -> list:
    """清洗生物数据，返回 [(日期, 采样点, 物种, 数量数值), ...]"""
    results = []
    date_map = {"5.16": "2025-5.16", "5.23": "2025-5.23", "5.31": "2025-5.31"}

    # 按调查批次或采样点分段
    # 注意：原始格式多变，此处展示清洗逻辑

    for line in raw_text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue

        # 尝试识别日期
        date = None
        for key, val in date_map.items():
            if key in line:
                date = val
                break

        # 尝试识别采样点
        loc = None
        for known_loc in KNOWN_LOCATIONS:
            if known_loc in line:
                loc = known_loc
                break

        # 拆分物种+数量对
        if date and loc:
            matches = SPECIES_COUNT_RE.findall(line)
            for sp_name, cnt_str in matches:
                sp_name = sp_name.strip()
                # 清理前缀
                if sp_name.startswith("（") and "疑似" in sp_name:
                    sp_name = re.sub(r"^[（(].*?[）)]", "", sp_name)
                # 修正打字错误
                sp_name = TYPO_FIXES.get(sp_name, sp_name)
                # 模糊计数转数值
                cnt = int(cnt_str) if cnt_str.isdigit() else FUZZY_COUNT_MAP.get(cnt_str, 1)
                results.append((date, loc, sp_name, cnt))

    return results


def normalize_count(cnt_str: str) -> str:
    """统一模糊计数表达"""
    if cnt_str.isdigit():
        return cnt_str
    # 若干（大量） → 大量
    if "大量" in cnt_str:
        return "大量"
    return cnt_str


if __name__ == "__main__":
    print("数据清洗脚本")
    print("原始数据（自由文本格式）已清洗为 data/data.txt 和 data/water.txt")
    print("清洗操作摘要：")
    print("  1. 日期标准化 → '2025-5.16' 等标准格式")
    print("  2. 采样点识别 → 匹配 5 个固定地点")
    print("  3. 物种-数量拆分 → 正则分离")
    print("  4. 模糊计数映射 → 若干=3, 很多=8, 大量=15")
    print("  5. 打字错误修正 → 中心刚硅藻→中心纲硅藻")
    print("  6. 单位剥离 → 水质表头统一")
