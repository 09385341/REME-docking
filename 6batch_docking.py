'''成功运行后，结果文件夹结构如下：
results/
    A0A1B1WA13/
        pocket_1/
            vina-output.pdbqt
            vina.log
            pocket_info.txt
        pocket_2/
            vina-output.pdbqt
            vina.log
            pocket_info.txt
'''
import subprocess
import sys
from pathlib import Path


def read_config(config_path: Path) -> dict:
    cfg = {}
    if not config_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_path}")
    with config_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                v = v.strip()
                if k:
                    v = v.strip().strip('"').strip("'")
                    cfg[k] = v
    return cfg


def parse_pockets(pockets_path: Path):
    """
    解析 all_pockets.txt，每一行例如：
    A0A1B1WA13_alphafold\t--center_x 15.9 --center_y 5.4 --center_z -11.8 --size_x 20.0 --size_y 20.0 --size_z 20.0
    返回列表：[(protein_name, box_dict, original_tokens), ...]
    """
    pockets = []
    if not pockets_path.exists():
        raise FileNotFoundError(f"口袋文件不存在: {pockets_path}")
    with pockets_path.open("r", encoding="utf-8") as f:
        for idx, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) < 2:
                print(f"警告: 第{idx}行格式异常，已跳过: {line}")
                continue
            protein_name, args_str = parts[0], parts[1]
            tokens = args_str.split()
            box = {}
            i = 0
            while i < len(tokens):
                token = tokens[i]
                if token.startswith("--"):
                    key = token[2:]
                    if i + 1 < len(tokens):
                        value = tokens[i + 1]
                        box[key] = value
                        i += 2
                    else:
                        print(f"警告: 第{idx}行参数 {token} 缺少值")
                        i += 1
                else:
                    i += 1

            required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
            missing = [r for r in required if r not in box]
            if missing:
                print(f"警告: 第{idx}行缺少必要参数: {missing}，已跳过")
                continue
            pockets.append((protein_name, box, tokens))
    return pockets


def report_path(label: str, p: Path, must_exist: bool = False, hint: str = "") -> bool:
    p = Path(p)
    exists = p.exists()
    status = "找到" if exists else "未找到"
    extra = f" | {hint}" if hint else ""
    print(f"[路径检查] {label}: {status} -> {p}{extra}")
    if must_exist and not exists:
        raise FileNotFoundError(f"{label} 不存在: {p}")
    return exists


def resolve_path(base_dir: Path, raw_value: str | None) -> Path | None:
    if raw_value is None:
        return None
    s = str(raw_value).strip().strip('"').strip("'")
    if not s:
        return None
    p = Path(s)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def run_cmd(cmd, cwd=None):
    """
    执行命令并输出。
    兼容 Windows 下外部程序输出非 UTF-8 的情况（比如 vina.exe）。
    """
    print(f"\n执行命令: {' '.join(cmd)}")

    # text=True 时 Python 会按 encoding 解码 stdout。
    # vina.exe 在 Windows 上经常输出本地代码页（GBK/CP936），会导致 utf-8 解码报错。
    # 这里优先 utf-8，失败时自动回退到系统本地编码（mbcs），并且替换无法解码的字符。
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="strict",
        )
        out = result.stdout
    except UnicodeDecodeError:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="mbcs",   # Windows 本地代码页（通常是 cp936）
            errors="replace",
        )
        out = result.stdout

    print(out)
    if result.returncode != 0:
        raise RuntimeError(f"命令执行失败，返回码: {result.returncode}")


def main():
    base_dir = Path(__file__).parent.resolve()
    print("[工作目录] base_dir =", base_dir)

    config_path = base_dir / "config.txt"
    pockets_path = base_dir / "all_pockets.txt"
    proteins_root = base_dir / "protein_structures_done"
    alphafold_dir = proteins_root / "alphafold"
    pdb_dir = proteins_root / "pdb"
    tmp_root = base_dir / "tmp"
    tmp_common = base_dir / "tmp_common"
    results_root = base_dir / "results"

    report_path("config.txt", config_path, must_exist=True)
    report_path("all_pockets.txt", pockets_path, must_exist=True)
    report_path("protein_structures_done/", proteins_root, must_exist=False,
                hint="若不存在，将导致后续受体文件全部无法找到")
    report_path("alphafold 目录", alphafold_dir, must_exist=False)
    report_path("pdb 目录", pdb_dir, must_exist=False)

    cfg = read_config(config_path)

    amdock_dir_cfg = cfg.get("AMDock_dir")
    if amdock_dir_cfg:
        amdock_dir_path = resolve_path(base_dir, amdock_dir_cfg)
        print(f"[配置] AMDock_dir(来自config) = {amdock_dir_cfg}")
    else:
        amdock_dir_path = (base_dir / "AMDock-win-master").resolve()
        print(f"[配置] 未在 config.txt 指定 AMDock_dir，默认使用: {amdock_dir_path}")

    print(f"[路径解析] 最终使用的 AMDock_dir = {amdock_dir_path}")
    report_path("AMDock_dir", amdock_dir_path, must_exist=True)

    python_path = amdock_dir_path / "python.exe"
    prepare_receptor_py = amdock_dir_path / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
    prepare_ligand_py = amdock_dir_path / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"

    vina_cfg = cfg.get("vina")
    if vina_cfg:
        vina_path = resolve_path(base_dir, vina_cfg)
        print(f"[配置] vina(来自config) = {vina_cfg}")
    else:
        vina_path = amdock_dir_path / "Lib" / "site-packages" / "AMDock" / "programs" / "vina.exe"
        print(f"[配置] 未指定 vina，默认使用 AMDock 自带: {vina_path}")

    print(f"[路径解析] 最终使用的 vina = {vina_path}")

    exhaustiveness = cfg.get("exhaustiveness", "20")
    energy_range = cfg.get("energy_range", "100")
    num_modes = cfg.get("num_modes", "10")

    ligand_file = cfg.get("ligand_pdb", "output_with_charges.pdb")
    ligand_path = base_dir / ligand_file

    report_path("AMDock python.exe", python_path, must_exist=True)
    report_path("prepare_receptor4.py", prepare_receptor_py, must_exist=True)
    report_path("prepare_ligand4.py", prepare_ligand_py, must_exist=True)
    report_path("vina.exe", vina_path, must_exist=True)
    report_path("ligand pdb 文件", ligand_path, must_exist=True, hint=f"来自 config: ligand_pdb={ligand_file}")

    tmp_root.mkdir(exist_ok=True)
    tmp_common.mkdir(exist_ok=True)
    results_root.mkdir(exist_ok=True)
    report_path("tmp/", tmp_root, must_exist=True)
    report_path("tmp_common/", tmp_common, must_exist=True)
    report_path("results/", results_root, must_exist=True)

    ligand_pdbqt = tmp_common / "ligand.pdbqt"
    if report_path("共用配体 pdbqt", ligand_pdbqt, must_exist=False):
        print("\n已存在共用配体 pdbqt，跳过生成:", ligand_pdbqt)
    else:
        print("\n第1步，准备配体的 pdbqt 文件（一次性）...")
        run_cmd([str(python_path), str(prepare_ligand_py), "-l", str(ligand_path), "-o", str(ligand_pdbqt)])
        report_path("共用配体 pdbqt(生成后)", ligand_pdbqt, must_exist=True)

    pockets = parse_pockets(pockets_path)
    print(f"[口袋解析] 共解析到 {len(pockets)} 条口袋记录")
    if not pockets:
        print("未解析到任何口袋信息，退出。")
        return

    pocket_index_map = {}
    receptor_ready = {}

    for protein_name, box, _tokens in pockets:
        idx = pocket_index_map.get(protein_name, 0) + 1
        pocket_index_map[protein_name] = idx
        print(f"\n==== 开始对接: {protein_name} | 口袋 #{idx} ====")

        if "alphafold" in protein_name.lower():
            receptor_src = alphafold_dir / f"{protein_name}.pdb"
            src_hint = "来源目录=alphafold"
        else:
            receptor_src = pdb_dir / f"{protein_name}.pdb"
            src_hint = "来源目录=pdb"

        if not report_path("受体 PDB 文件", receptor_src, must_exist=False, hint=src_hint):
            print(f"警告: 找不到蛋白质 PDB 文件，跳过该口袋: {protein_name} pocket#{idx}")
            continue

        protein_tmp = tmp_root / protein_name
        protein_tmp.mkdir(exist_ok=True)
        protein_out_base = results_root / protein_name
        protein_out_base.mkdir(exist_ok=True)
        report_path("protein_tmp 目录", protein_tmp, must_exist=True)
        report_path("protein_results 目录", protein_out_base, must_exist=True)

        pocket_out = protein_out_base / f"pocket_{idx}"
        pocket_out.mkdir(exist_ok=True)
        report_path("pocket_out 目录", pocket_out, must_exist=True)

        receptor_pdbqt = protein_tmp / "receptor.pdbqt"
        report_path("receptor.pdbqt(预期位置)", receptor_pdbqt, must_exist=False)

        if not receptor_ready.get(protein_name) or not receptor_pdbqt.exists():
            print("第2步，生成蛋白的 pdbqt 文件...")
            run_cmd([
                str(python_path), str(prepare_receptor_py),
                "-r", str(receptor_src),
                "-o", str(receptor_pdbqt),
                "-A", "checkhydrogens",
                "-U", "nphs"
            ])
            receptor_ready[protein_name] = True
            report_path("receptor.pdbqt(生成后)", receptor_pdbqt, must_exist=True)
        else:
            print("已存在受体 pdbqt，复用:", receptor_pdbqt)

        print("第3步，执行 AutoDock Vina ...")
        vina_out = pocket_out / "vina-output.pdbqt"
        vina_log = pocket_out / "vina.log"

        pocket_info = pocket_out / "pocket_info.txt"
        with pocket_info.open("w", encoding="utf-8") as f:
            f.write(f"protein={protein_name}\n")
            f.write(f"pocket_index={idx}\n")
            f.write(f"center_x={box['center_x']}\n")
            f.write(f"center_y={box['center_y']}\n")
            f.write(f"center_z={box['center_z']}\n")
            f.write(f"size_x={box['size_x']}\n")
            f.write(f"size_y={box['size_y']}\n")
            f.write(f"size_z={box['size_z']}\n")
            f.write(f"exhaustiveness={exhaustiveness}\n")
            f.write(f"energy_range={energy_range}\n")
            f.write(f"num_modes={num_modes}\n")
        report_path("pocket_info.txt(写入后)", pocket_info, must_exist=True)

        cmd = [
            str(vina_path),
            "--exhaustiveness", str(exhaustiveness),
            "--energy_range", str(energy_range),
            "--num_modes", str(num_modes),
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(box["center_x"]),
            "--center_y", str(box["center_y"]),
            "--center_z", str(box["center_z"]),
            "--size_x", str(box["size_x"]),
            "--size_y", str(box["size_y"]),
            "--size_z", str(box["size_z"]),
            "--out", str(vina_out),
            "--log", str(vina_log)
        ]

        report_path("vina 输出文件(将生成)", vina_out, must_exist=False)
        report_path("vina 日志(将生成)", vina_log, must_exist=False)

        run_cmd(cmd)

        report_path("vina 输出文件(生成后)", vina_out, must_exist=False)
        report_path("vina 日志(生成后)", vina_log, must_exist=False)

        print(f"完成: {protein_name} 口袋#{idx} -> {vina_out}")

    print("\n全部对接完成！结果位于:", results_root)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("发生错误:", e)
        sys.exit(1)