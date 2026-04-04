import os
import glob
from pymol import cmd

def preprocess(input_path, output_path):
    print(f"正在处理: {input_path} 输出到 {output_path}")
    cmd.reinitialize()
    cmd.load(input_path, "prot")
    cmd.remove("resn HOH")
    cmd.remove("hetatm and not resn HOH")
    cmd.h_add("prot")
    cmd.save(output_path, "prot")
    cmd.delete("prot")

def batch_preprocess(input_dir, output_dir):
    print(f"检查文件夹: {input_dir}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    print(f"找到 {len(pdb_files)} 个PDB文件")
    for pdb_file in pdb_files:
        filename = os.path.basename(pdb_file)
        output_path = os.path.join(output_dir, filename)
        preprocess(pdb_file, output_path)

if __name__ == "__main__":
    print("启动批量预处理脚本……")
    # 获取当前脚本所在目录，确保无论在哪运行脚本都能定位到数据
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

    src_root = os.path.join(ROOT_DIR, "protein_structures")
    dst_root = os.path.join(ROOT_DIR, "protein_structures_done")
    for subdir in ["alphafold", "pdb"]:
        input_dir = os.path.join(src_root, subdir)
        output_dir = os.path.join(dst_root, subdir)
        batch_preprocess(input_dir, output_dir)