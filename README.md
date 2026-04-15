# 基于分子对接的候选酶活性与选择性预测流程

本项目针对目标反应，通过 REME 平台筛选候选酶，获取蛋白结构，进行分子对接并评估结合亲和力，为酶的实验选择提供计算依据。

## 研究背景

1. **候选酶初步筛选**  
   利用 REME 平台，基于多种反应相似性评价方法，筛选反应机理上合理的酶 EC 编号，构建候选酶集。

2. **蛋白结构获取与预处理**  
   从 PDB 或 AlphaFold 获取候选酶结构，进行去水、去配体、加氢、能量最小化等预处理。

3. **分子对接预测活性与选择性**  
   将反应物分子对接至酶活性口袋，通过结合亲和力打分排序，辅助酶的选择。

## 工作流程与脚本说明

所有脚本均放在项目根目录下运行，按顺序执行：

| 脚本 | 功能 |
| :--- | :--- |
| `0_edgedriver_auto_fetch.py` | 自动获取 EdgeDriver，为爬虫做准备 |
| `0_SMILES_to_rxn.py` | 将反应的SMILES转为rxn文件，方便上传至REME |
| `1_inputreaction.py` | 读取 RXN 文件，自动提交至 REME 网站，爬取获得 UniProt ID |
| `2_download_structures.py` | 根据 UniProt ID 自动下载 PDB 结构（有晶体结构从 RCSB PDB 获取，无则下载 AlphaFold 预测结构） |
| `3_batch_detect_pockets_dual.py` | 调用 pyKVFinder 批量计算所有 PDB 的口袋质心，输出 `centroids_all_pdb.csv` 和 `centroids_all_alphafold.csv` |
| `4_scripts_merge_pockets_Version2.py` | 合并两个 CSV 文件，生成 `all_pockets.txt`（保留来源蛋白名称） |
| `5.1_pymol_preprocess_proteins.py` | 对接前使用 PyMOL 预处理蛋白结构 |
| `5.2_extract_rxn_to_pdb.py` | 手动输入 SMILES，将反应物小分子预处理为对接格式 |
| `6_batch_docking.py` | 批量执行 AutoDock Vina 分子对接 |
| `7_rank_vina_affinity.py` | 解析 vina.log 中的亲和力，结合口袋信息生成排名 |
| `config.txt` | 对接信息记录，做对接不可或缺 |

## 注意

1. **如果成功执行了一个任务后，做另一个反应任务**  
   必须确保tmp、tmp_common文件夹都已经移除，如果存在则可能使用上一个任务生成的小分子配体、蛋白质的pdbqt文件做对接。
   必须确保result文件夹都已经移除，如果存在则对接后评分会引入之前任务的对接结果与分数。
   
2.**如果不知道反应的SMILES怎么写**  
   手动打开REME网站（https://reme.biodesign.ac.cn/），点击SMILES Generator，在画布上自己画出反应，然后保存就可以直接得到rxn文件，就不用再执行`0_SMILES_to_rxn.py`。

## 安装与环境配置

本项目使用 conda 环境管理依赖。请确保已安装 [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 或 [Anaconda](https://www.anaconda.com/)。

```bash
# 克隆仓库
git clone https://github.com/09385341/REME-docking.git
cd your-repo-name

# 创建环境
conda env create -f environment.yml

# 激活环境
conda activate reme-try
