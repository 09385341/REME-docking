#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
成功运行后会在脚本同目录生成两个 CSV 文件：centroids_all_pdb.csv 和 centroids_all_alphafold.csv，包含所有处理的 PDB 文件的口袋质心信息。
batch_detect_pockets_dual.py

功能：
- 自动读取脚本同一目录下的 protein_structures/pdb 与 protein_structures/alphafold 文件夹中所有 .pdb 文件（非递归）。
- 调用 pyKVFinder，计算每个输入 PDB 的口袋质心。
- 仅输出两个汇总 CSV 到脚本同目录：
    - centroids_all_pdb.csv
    - centroids_all_alphafold.csv
  CSV 格式为：
    source_pdb,cavity_name,label,centroid_x,centroid_y,centroid_z
- 不再生成伪原子 PDB，不再生成每个 PDB 的单独 CSV。

可选参数（保持与原脚本类似的 KVFinder 参数配置）：
    --step           网格间距（Å）
    --probe-in       Probe in（Å）
    --probe-out      Probe out（Å）
    --volume-cutoff  体积阈值（Å^3）
    --no-depth       不计算深度
    --hydropathy     计算疏水性（可能较慢）
    --surface-only   优先使用 surface 网格（若可用，否则退回 cavities 网格）

用法：
    python batch_detect_pockets_dual.py
"""

from __future__ import annotations
import argparse
import pathlib
import sys
import traceback
import numpy as np
from typing import List, Tuple, Dict, Optional

try:
    import pyKVFinder
except Exception as e:
    print("错误：无法导入 pyKVFinder，请先安装或确保在 PYTHONPATH 中。")
    print(e)
    sys.exit(1)


def voxel_index_to_coord(vertices: np.ndarray, step: float, ix: int, iy: int, iz: int) -> np.ndarray:
    P1, P2, P3, P4 = np.asarray(vertices, dtype=float)
    v2 = P2 - P1
    v3 = P3 - P1
    v4 = P4 - P1

    def safe_unit(v):
        norm = np.linalg.norm(v)
        if norm < 1e-8:
            return v
        return v / norm

    vx = safe_unit(v2)
    vy = safe_unit(v3)
    vz = safe_unit(v4)

    coord = P1 + (ix + 0.5) * step * vx + (iy + 0.5) * step * vy + (iz + 0.5) * step * vz
    return coord


def get_cavity_name_from_label(label: int) -> str:
    index = int(label) - 2
    first = chr(65 + (index // 26) % 26)
    second = chr(65 + (index % 26))
    return f"K{first}{second}"


def find_pdb_files_in_dir(input_dir: pathlib.Path) -> List[pathlib.Path]:
    # 非递归查找 .pdb 文件（大小写不敏感）
    if not input_dir.exists() or not input_dir.is_dir():
        return []
    return [p for p in sorted(input_dir.iterdir()) if p.is_file() and p.suffix.lower() == ".pdb"]


def compute_centroids(results, surface_only: bool = True) -> List[Tuple[str, int, float, float, float]]:
    """
    从 pyKVFinder 结果计算所有口袋的质心列表。
    返回：[(cavity_name, label, x, y, z), ...]
    """
    vertices = getattr(results, "_vertices", None)
    step = getattr(results, "_step", None)
    if vertices is None or step is None:
        raise RuntimeError("结果对象缺少 _vertices 或 _step，无法从栅格坐标映射到三维坐标。")

    vertices = np.asarray(vertices, dtype=float)
    step = float(step)

    cavities_grid = np.asarray(results.cavities)
    surface_grid = getattr(results, "surface", None)

    grid = np.asarray(surface_grid) if (surface_only and surface_grid is not None) else cavities_grid
    mask = grid >= 2

    if not np.any(mask):
        return []

    ix_list, iy_list, iz_list = np.where(mask)
    labels = grid[ix_list, iy_list, iz_list].astype(int)

    cavity_coords: Dict[int, List[np.ndarray]] = {}
    for ix, iy, iz, lab in zip(ix_list, iy_list, iz_list, labels):
        coord = voxel_index_to_coord(vertices, step, int(ix), int(iy), int(iz))
        cavity_coords.setdefault(int(lab), []).append(coord)

    centroids: List[Tuple[str, int, float, float, float]] = []
    for lab in sorted(cavity_coords.keys()):
        pts = np.asarray(cavity_coords[lab])
        centroid = pts.mean(axis=0)
        cavity_name = get_cavity_name_from_label(lab)
        centroids.append((cavity_name, lab, float(centroid[0]), float(centroid[1]), float(centroid[2])))

    return centroids


def run_kvfinder_on_file(pdb_path: pathlib.Path, args) -> Optional[List[Tuple[str, int, float, float, float]]]:
    """
    调用 pyKVFinder 处理单个 PDB，返回质心列表；失败时返回 None。
    """
    try:
        results = pyKVFinder.run_workflow(
            input=str(pdb_path),
            ligand=None,
            vdw=None,
            box=None,
            step=args.step,
            probe_in=args.probe_in,
            probe_out=args.probe_out,
            removal_distance=2.4,
            volume_cutoff=args.volume_cutoff,
            ligand_cutoff=5.0,
            include_depth=args.depth,
            include_hydropathy=args.hydropathy,
            hydrophobicity_scale="KyteDoolittle",
            surface="SES",
            ignore_backbone=False,
            nthreads=None,
            verbose=False,
        )
    except TypeError:
        try:
            results = pyKVFinder.run_workflow(str(pdb_path))
        except Exception:
            print("run_workflow 调用失败（降级后仍失败），打印 traceback：")
            traceback.print_exc()
            return None
    except Exception:
        print("调用 run_workflow 时发生异常，跳过此 PDB，打印 traceback：")
        traceback.print_exc()
        return None

    if results is None:
        print(f"pyKVFinder 返回 None 或未检测到口袋：{pdb_path.name}")
        return []

    try:
        return compute_centroids(results, surface_only=args.surface_only)
    except Exception:
        print(f"为 {pdb_path.name} 计算质心失败，打印 traceback：")
        traceback.print_exc()
        return None


def process_dir(input_dir: pathlib.Path, out_csv_path: pathlib.Path, args) -> None:
    """
    处理指定目录下的所有 .pdb 文件，将质心汇总到 out_csv_path。
    """
    pdb_files = find_pdb_files_in_dir(input_dir)

    with out_csv_path.open("w", encoding="utf-8") as fh:
        fh.write("source_pdb,cavity_name,label,centroid_x,centroid_y,centroid_z\n")

    if not pdb_files:
        print(f"在 {input_dir} 中未找到任何 .pdb 文件。已写入空的 CSV header 到 {out_csv_path.name}")
        return

    print(f"发现 {len(pdb_files)} 个 PDB 文件于 {input_dir}，开始处理...")

    appended = 0
    for pdb_path in pdb_files:
        basename = pdb_path.stem
        print("\n" + "=" * 60)
        print(f"处理: {pdb_path.name}")

        centroids = run_kvfinder_on_file(pdb_path, args)
        if centroids is None:
            continue

        if len(centroids) == 0:
            print(f"{basename} 未检测到口袋（或结果为空）。")
            continue

        try:
            with out_csv_path.open("a", encoding="utf-8") as fh:
                for cname, lab, x, y, z in centroids:
                    fh.write(f"{basename},{cname},{lab},{x:.3f},{y:.3f},{z:.3f}\n")
                    appended += 1
        except Exception:
            print("追加到汇总 CSV 时发生异常：")
            traceback.print_exc()

    print(f"\n目录 {input_dir.name} 处理完成。共写入 {appended} 条质心记录到：{out_csv_path.resolve()}")


def main():
    parser = argparse.ArgumentParser(
        description="自动读取脚本同目录下 protein_structures/pdb 与 protein_structures/alphafold 的 PDB，使用 pyKVFinder 输出两个汇总 CSV。"
    )
    parser.add_argument("--step", type=float, default=0.6, help="网格间距 Å (默认: 0.6)")
    parser.add_argument("--probe-in", type=float, default=1.4, help="Probe in (Å)")
    parser.add_argument("--probe-out", type=float, default=4.0, help="Probe out (Å)")
    parser.add_argument("--volume-cutoff", type=float, default=50.0, help="体积阈值 Å^3")
    parser.add_argument("--no-depth", dest="depth", action="store_false", help="不计算深度")
    parser.add_argument("--hydropathy", action="store_true", help="计算疏水性（可能慢）")
    parser.add_argument("--surface-only", dest="surface_only", action="store_true", help="优先使用 surface 网格（若可用）")
    args = parser.parse_args()

    base_dir = pathlib.Path(__file__).resolve().parent
    protein_structures_dir = base_dir / "protein_structures"
    pdb_dir = protein_structures_dir / "pdb"
    alphafold_dir = protein_structures_dir / "alphafold"

    out_csv_pdb = base_dir / "centroids_all_pdb.csv"
    out_csv_alphafold = base_dir / "centroids_all_alphafold.csv"

    print("脚本目录：", base_dir)
    print("目标目录：", protein_structures_dir)
    print("pdb 子目录：", pdb_dir)
    print("alphafold 子目录：", alphafold_dir)
    print("输出 CSV：", out_csv_pdb.name, "与", out_csv_alphafold.name)

    any_dir_exists = False

    if pdb_dir.exists() and pdb_dir.is_dir():
        any_dir_exists = True
        process_dir(pdb_dir, out_csv_pdb, args)
    else:
        print(f"警告：未找到目录 {pdb_dir}，将写入空的 CSV header 至 {out_csv_pdb.name}")
        with out_csv_pdb.open("w", encoding="utf-8") as fh:
            fh.write("source_pdb,cavity_name,label,centroid_x,centroid_y,centroid_z\n")

    if alphafold_dir.exists() and alphafold_dir.is_dir():
        any_dir_exists = True
        process_dir(alphafold_dir, out_csv_alphafold, args)
    else:
        print(f"警告：未找到目录 {alphafold_dir}，将写入空的 CSV header 至 {out_csv_alphafold.name}")
        with out_csv_alphafold.open("w", encoding="utf-8") as fh:
            fh.write("source_pdb,cavity_name,label,centroid_x,centroid_y,centroid_z\n")

    if not any_dir_exists:
        print("错误：未找到 protein_structures/pdb 或 protein_structures/alphafold 目录。请检查目录结构。")
        sys.exit(2)

    print("\n全部处理完成。")
    print("centroids_all_pdb.csv 路径：", out_csv_pdb.resolve())
    print("centroids_all_alphafold.csv 路径：", out_csv_alphafold.resolve())
    print("结束。")


if __name__ == "__main__":
    main()

'''PS C:\Users\Jiaya> & D:/conda_envs/reme-try/python.exe "d:/Project/pythonProject/pythonProject/reactionsearch - 副本/3batch_detect_pockets_dual.py"
脚本目录： D:\Project\pythonProject\pythonProject\reactionsearch - 副本
目标目录： D:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures
pdb 子目录： D:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb
alphafold 子目录： D:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold
输出 CSV： centroids_all_pdb.csv 与 centroids_all_alphafold.csv
发现 38 个 PDB 文件于 D:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb，开始处理...

============================================================
处理: A0A0H4SN47_6YC8.pdb
Warning: No cavities detected, returning None!
pyKVFinder 返回 None 或未检测到口袋：A0A0H4SN47_6YC8.pdb
A0A0H4SN47_6YC8 未检测到口袋（或结果为空）。

============================================================
处理: C5H429_5DXX.pdb

============================================================
处理: C5H429_5DXY.pdb

============================================================
处理: C5H429_5DY2.pdb

============================================================
处理: C5H429_5DY3.pdb
Warning: No cavities detected, returning None!
pyKVFinder 返回 None 或未检测到口袋：C5H429_5DY3.pdb
C5H429_5DY3 未检测到口袋（或结果为空）。

============================================================
处理: F9UM18_7EBO.pdb

============================================================
处理: P04964_1DAP.pdb

============================================================
处理: P04964_1F06.pdb

============================================================
处理: P04964_2DAP.pdb

============================================================
处理: P04964_3DAP.pdb

============================================================
处理: P04964_5LOA.pdb

============================================================
处理: P04964_5LOC.pdb

============================================================
处理: P51110_2IOD.pdb

============================================================
处理: P51110_2NNL.pdb

============================================================
处理: P51110_3BXX.pdb

============================================================
处理: P51110_3C1T.pdb

============================================================
处理: P76469_2VWS.pdb

============================================================
处理: P76469_2VWT.pdb

============================================================
处理: Q65CJ7_3BA1.pdb

============================================================
处理: Q65CJ7_3BAZ.pdb
Warning: No cavities detected, returning None!
pyKVFinder 返回 None 或未检测到口袋：Q65CJ7_3BAZ.pdb
Q65CJ7_3BAZ 未检测到口袋（或结果为空）。

============================================================
处理: Q8GHB2_2XLQ.pdb

============================================================
处理: Q8GHB2_2XLY.pdb

============================================================
处理: Q8GHB2_2XM5.pdb

============================================================
处理: Q8GHB2_2XM7.pdb

============================================================
处理: Q93TU6_1O8U.pdb

============================================================
处理: Q9BXD5_6ARH.pdb

============================================================
处理: Q9R9V9_5CPL.pdb

============================================================
处理: Q9R9V9_5CPM.pdb

============================================================
处理: Q9R9V9_5CPN.pdb

============================================================
处理: Q9R9V9_5CPO.pdb
Warning: No cavities detected, returning None!
pyKVFinder 返回 None 或未检测到口袋：Q9R9V9_5CPO.pdb
Q9R9V9_5CPO 未检测到口袋（或结果为空）。

============================================================
处理: Q9R9V9_5N6Q.pdb

============================================================
处理: Q9R9V9_8A8I.pdb

============================================================
处理: Q9R9V9_8AU8.pdb

============================================================
处理: Q9R9V9_8AU9.pdb

============================================================
处理: Q9R9V9_8AUF.pdb

============================================================
处理: Q9R9V9_8AUG.pdb

============================================================
处理: Q9R9V9_8AUH.pdb

============================================================
处理: Q9R9V9_8AUI.pdb

目录 pdb 处理完成。共写入 121 条质心记录到：D:\Project\pythonProject\pythonProject\reactionsearch - 副本\centroids_all_pdb.csv
发现 29 个 PDB 文件于 D:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold，开始处理...

============================================================
处理: A0A0A0PXZ7_alphafold.pdb

============================================================
处理: A0A0C6ENR6_alphafold.pdb

============================================================
处理: A0A168S819_alphafold.pdb

============================================================
处理: A0A168S850_alphafold.pdb

============================================================
处理: A0A1B1WA13_alphafold.pdb

============================================================
处理: A0A1B1WA14_alphafold.pdb

============================================================
处理: A0A1L7NRS0_alphafold.pdb

============================================================
处理: B9A1Q4_alphafold.pdb

============================================================
处理: F7YTT4_alphafold.pdb

============================================================
处理: I1NJA2_alphafold.pdb

============================================================
处理: M4DU26_alphafold.pdb

============================================================
处理: P40952_alphafold.pdb

============================================================
处理: P44430_alphafold.pdb

============================================================
处理: P51102_alphafold.pdb

============================================================
处理: P51103_alphafold.pdb

============================================================
处理: P51104_alphafold.pdb

============================================================
处理: P51105_alphafold.pdb

============================================================
处理: P51107_alphafold.pdb

============================================================
处理: Q0QFI7_alphafold.pdb

============================================================
处理: Q53B70_alphafold.pdb

============================================================
处理: Q6Q874_alphafold.pdb

============================================================
处理: Q6UAQ7_alphafold.pdb

============================================================
处理: Q8W404_alphafold.pdb

============================================================
处理: Q8W405_alphafold.pdb

============================================================
处理: Q8YNV6_alphafold.pdb

============================================================
处理: Q93XE6_alphafold.pdb

============================================================
处理: Q93YP7_alphafold.pdb

============================================================
处理: Q9AJC6_alphafold.pdb

============================================================
处理: Q9RF52_alphafold.pdb

目录 alphafold 处理完成。共写入 73 条质心记录到：D:\Project\pythonProject\pythonProject\reactionsearch - 副本\centroids_all_alphafold.csv

全部处理完成。
centroids_all_pdb.csv 路径： D:\Project\pythonProject\pythonProject\reactionsearch - 副本\centroids_all_pdb.csv
centroids_all_alphafold.csv 路径： D:\Project\pythonProject\pythonProject\reactionsearch - 副本\centroids_all_alphafold.csv
结束。'''