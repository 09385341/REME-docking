from rdkit import Chem
from rdkit.Chem import rdChemReactions

# 反应 SMILES，反应物和产物用 '>>' 分隔
rxn_smiles = 'Nc1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)n1.Nc1nc2c(c(=O)[nH]1)N[C@H]1C3=C(S[MoH](=O)(=O)(O)S3)[C@@H](COP(=O)(O)O)O[C@H]1N2>>Nc1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@H]4Nc5nc(N)[nH]c(=O)c5N[C@H]4C4=C3S[MoH](=O)(=O)(O)S4)[C@@H](O)[C@H]2O)c(=O)n1.O=P(O)(O)OP(=O)(O)O'
# 创建反应对象
rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)

# 写入 RXN 文件
with open("Rhea31335.rxn", "w") as f:
    f.write(rdChemReactions.ReactionToRxnBlock(rxn))