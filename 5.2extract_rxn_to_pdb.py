from rdkit import Chem
from rdkit.Chem import AllChem
import os

# 定义分子的 SMILES 字符串(输入反应物的SMILES)
smiles = 'CCOC(=O)C1=NO[C@]2(CCC(=O)CCO2)[C@@H]1OCc1ccccc1'
mol = Chem.MolFromSmiles(smiles)

# 添加氢原子
mol = Chem.AddHs(mol)

# 进行能量最小化
AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
AllChem.MMFFOptimizeMolecule(mol)

# 计算 Gasteiger 部分电荷
AllChem.ComputeGasteigerCharges(mol)

# 打印每个原子的部分电荷
for atom in mol.GetAtoms():
    print(f"Atom: {atom.GetSymbol()}, Gasteiger Charge: {atom.GetProp('_GasteigerCharge')}")

# 保存带有部分电荷的分子为 PDB 文件（写到脚本所在目录）
script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, 'output_with_charges.pdb')


print("write to:", out_path)
pdb_text = Chem.MolToPDBBlock(mol)
with open(out_path, "w", encoding="utf-8", newline="\n") as f:
    f.write(pdb_text)
print("PDB written to:", out_path)