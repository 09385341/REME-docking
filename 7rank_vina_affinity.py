#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''成功解析 vina.log 中的 affinity，并结合 pocket_info.txt 中的口袋信息，生成蛋白质对接亲和力排名。'''
import os
import sys
import csv
import re
from typing import List, Dict, Optional

def parse_vina_log(log_path: str) -> List[float]:
    """
    解析 vina.log，返回所有模式（mode）的 affinity 列表（float）。
    解析策略：
    - 表格行形如：<mode> <affinity> <rmsd l.b.> <rmsd u.b.>
    - 识别以数字开头且至少 4 列的行，第二列作为 affinity。
    """
    affinities: List[float] = []
    try:
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                tokens = s.split()
                # 形如: 1 -6.9 0.000 0.000
                if len(tokens) >= 4 and tokens[0].isdigit():
                    try:
                        aff = float(tokens[1])
                        affinities.append(aff)
                    except ValueError:
                        continue
    except FileNotFoundError:
        pass
    return affinities

def parse_pocket_info(info_path: str) -> Dict[str, Optional[str]]:
    """
    解析 pocket_info.txt 中常见的口袋信息。
    优先提取常见键：center_x/y/z, size_x/y/z（不区分大小写，支持 :, = 分隔）
    同时保留原始文本到 'raw' 字段。未找到则为 None。
    """
    result = {
        "center_x": None,
        "center_y": None,
        "center_z": None,
        "size_x": None,
        "size_y": None,
        "size_z": None,
        "raw": None,
    }
    try:
        with open(info_path, "r", encoding="utf-8", errors="ignore") as f:
            raw = f.read()
            result["raw"] = raw
    except FileNotFoundError:
        return result

    def extract_float(label_variants: List[str]) -> Optional[float]:
        # 在整段文本中查找 key [:=] value（float）
        for lv in label_variants:
            pattern = rf"{lv}\s*[:=]\s*([\-+]?\d+(?:\.\d+)?)"
            m = re.search(pattern, result["raw"], flags=re.IGNORECASE)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    pass
        return None

    # 常见别名/写法
    result["center_x"] = extract_float(["center_x", "center x", "box_center_x", "centerx"])
    result["center_y"] = extract_float(["center_y", "center y", "box_center_y", "centery"])
    result["center_z"] = extract_float(["center_z", "center z", "box_center_z", "centerz"])

    result["size_x"] = extract_float(["size_x", "size x", "box_size_x", "sizex"])
    result["size_y"] = extract_float(["size_y", "size y", "box_size_y", "sizey"])
    result["size_z"] = extract_float(["size_z", "size z", "box_size_z", "sizez"])

    return result

def collect_best_affinities_with_pockets(results_dir: str) -> List[Dict]:
    """
    遍历 results_dir 下的蛋白质子文件夹；对每个蛋白质：
    - 遍历 pocket_* 子文件夹；
    - 解析 vina.log，找出该口袋绝对值最大的 affinity；
    - 在所有口袋中选取绝对值最大的口袋作为该蛋白质的最佳对接结果；
    返回列表，每项包含：
    {
        "protein": 子文件夹名,
        "pocket": 口袋子文件夹名,
        "best_affinity": float,  # 原始值（通常为负）
        "abs_affinity": float,
        "pocket_info": {center_x, center_y, center_z, size_x, size_y, size_z, raw},
        "pocket_path": 口袋目录的相对路径（相对于 results_dir）
    }
    """
    ranking = []
    if not os.path.isdir(results_dir):
        print(f"错误：未找到目录 {results_dir}")
        return ranking

    for protein_entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
        if not protein_entry.is_dir():
            continue
        protein_name = protein_entry.name

        best_record_for_protein = None  # 存储该蛋白质的最佳（abs 最大）
        # 遍历 pocket_* 子目录
        try:
            pocket_entries = sorted(os.scandir(protein_entry.path), key=lambda e: e.name)
        except PermissionError:
            print(f"警告：无法访问目录 {protein_entry.path}")
            continue

        has_any_pocket = False
        for pocket_entry in pocket_entries:
            if not pocket_entry.is_dir():
                continue
            pocket_dir_name = pocket_entry.name
            # 仅考虑以 'pocket' 开头的目录（常见为 pocket_1, pocket_2, ...）
            if not pocket_dir_name.lower().startswith("pocket"):
                continue
            has_any_pocket = True

            log_path = os.path.join(pocket_entry.path, "vina.log")
            affinities = parse_vina_log(log_path)
            if not affinities:
                # 该口袋没有可用 affinity，跳过
                # 仅提示一次，避免刷屏
                print(f"提示：{protein_name}/{pocket_dir_name} 未解析到 affinity（缺失或格式不匹配）：{log_path}")
                continue

            # 该口袋的最佳 affinity（按绝对值）
            pocket_best_affinity = max(affinities, key=lambda x: abs(x))
            pocket_abs = abs(pocket_best_affinity)

            # 解析口袋信息
            info_path = os.path.join(pocket_entry.path, "pocket_info.txt")
            pocket_info = parse_pocket_info(info_path)

            record = {
                "protein": protein_name,
                "pocket": pocket_dir_name,
                "best_affinity": pocket_best_affinity,
                "abs_affinity": pocket_abs,
                "pocket_info": pocket_info,
                # 存储相对路径，便于查看
                "pocket_path": os.path.relpath(pocket_entry.path, results_dir),
            }

            if (best_record_for_protein is None) or (pocket_abs > best_record_for_protein["abs_affinity"]):
                best_record_for_protein = record

        if not has_any_pocket:
            print(f"警告：{protein_name} 未发现 pocket_* 子文件夹。")
        elif best_record_for_protein is None:
            print(f"警告：{protein_name} 的 pocket_* 中均未找到有效的 vina.log。")
        else:
            ranking.append(best_record_for_protein)

    return ranking

def print_ranking(ranking: List[Dict]) -> None:
    """
    打印排名到标准输出，包含蛋白质、最佳口袋名及 affinity。
    """
    if not ranking:
        print("未找到任何可用的对接结果。")
        return
    print("{:<6} {:<30} {:<14} {:>14} {:>14}".format("Rank", "Protein", "Pocket", "BestAffinity", "AbsValue"))
    print("-" * 90)
    for i, item in enumerate(ranking, start=1):
        print("{:<6} {:<30} {:<14} {:>14.3f} {:>14.3f}".format(
            i, item["protein"], item["pocket"], item["best_affinity"], item["abs_affinity"]
        ))

def save_csv(ranking: List[Dict], out_path: str) -> None:
    """
    将排名保存为 CSV 文件，包含口袋的常见参数与原始信息。
    """
    if not ranking:
        return
    fieldnames = [
        "Rank",
        "Protein",
        "Pocket",
        "BestAffinity(kcal/mol)",
        "AbsValue",
        "PocketPath",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "pocket_info_raw",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, item in enumerate(ranking, start=1):
            pi = item.get("pocket_info", {}) or {}
            row = {
                "Rank": i,
                "Protein": item["protein"],
                "Pocket": item["pocket"],
                "BestAffinity(kcal/mol)": f"{item['best_affinity']:.3f}",
                "AbsValue": f"{item['abs_affinity']:.3f}",
                "PocketPath": item.get("pocket_path", ""),
                "center_x": "" if pi.get("center_x") is None else f"{pi['center_x']:.3f}",
                "center_y": "" if pi.get("center_y") is None else f"{pi['center_y']:.3f}",
                "center_z": "" if pi.get("center_z") is None else f"{pi['center_z']:.3f}",
                "size_x": "" if pi.get("size_x") is None else f"{pi['size_x']:.3f}",
                "size_y": "" if pi.get("size_y") is None else f"{pi['size_y']:.3f}",
                "size_z": "" if pi.get("size_z") is None else f"{pi['size_z']:.3f}",
                "pocket_info_raw": pi.get("raw") or "",
            }
            writer.writerow(row)
    print(f"已保存结果到: {out_path}")

def main():
    # 默认读取脚本同级目录下的 results/
    if len(sys.argv) >= 2:
        results_dir = sys.argv[1]
    else:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        results_dir = os.path.join(base_dir, "results")

    ranking = collect_best_affinities_with_pockets(results_dir)
    # 按绝对值从大到小排序
    ranking.sort(key=lambda x: x["abs_affinity"], reverse=True)

    print_ranking(ranking)

    out_csv = os.path.join(os.getcwd(), "enzyme_affinity_ranking_with_pockets.csv")
    save_csv(ranking, out_csv)

if __name__ == "__main__":
    main()

'''PS D:\Project\pythonProject\pythonProject\reactionsearch> & D:/conda_envs/reme-try/python.exe d:/Project/pythonProject/pythonProject/reactionsearch/rank_vina_affinity.py
Rank   Protein                        Pocket           BestAffinity       AbsValue
------------------------------------------------------------------------------------------
1      Q8W404_alphafold               pocket_1               -9.500          9.500
2      A0A1L7NRS0_alphafold           pocket_1               -8.800          8.800
3      Q8W405_alphafold               pocket_1               -8.800          8.800
4      B9A1Q4_alphafold               pocket_1               -8.700          8.700
5      Q93YP7_alphafold               pocket_2               -8.500          8.500
6      A0A168S819_alphafold           pocket_1               -8.300          8.300
7      F9UM18_7EBO                    pocket_3               -8.300          8.300
8      Q0QFI7_alphafold               pocket_1               -8.300          8.300
9      A0A168S850_alphafold           pocket_1               -8.100          8.100
10     I1NJA2_alphafold               pocket_2               -8.100          8.100
11     P76469_2VWT                    pocket_3               -8.000          8.000
12     Q8GHB2_2XLY                    pocket_1               -8.000          8.000
13     Q9RF52_alphafold               pocket_1               -7.800          7.800
14     P76469_2VWS                    pocket_4               -7.700          7.700
15     Q8GHB2_2XM5                    pocket_1               -7.700          7.700
16     Q8GHB2_2XLQ                    pocket_1               -7.600          7.600
17     Q8GHB2_2XM7                    pocket_1               -7.600          7.600
18     Q6Q874_alphafold               pocket_1               -7.500          7.500
19     A0A1B1WA14_alphafold           pocket_1               -7.100          7.100
20     A0A1B1WA13_alphafold           pocket_1               -7.000          7.000
21     Q53B70_alphafold               pocket_2               -6.400          6.400
22     Q93XE6_alphafold               pocket_2               -5.900          5.900
已保存结果到: D:\Project\pythonProject\pythonProject\reactionsearch\enzyme_affinity_ranking_with_pockets.csv'''