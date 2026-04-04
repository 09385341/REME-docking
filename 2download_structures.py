"""
蛋白质结构自动下载工具
自动检测并下载蛋白质的实验结构(PDB)或AlphaFold预测结构

v4 修正（按你的最新要求）：
- 只过滤以下 6 个“确定不是 UniProt ID”的条目：
  ===、RDKIT、RXNFP、DRFP、(COSINE)、(EUCLIDEAN)
- 不做任何“正则判断/长度判断”的额外过滤，避免误伤真实 UniProt ID
- 同时对 UniProt ID 去重（保序）
- 其余逻辑不变：优先 UniProt idmapping 批量 -> 抓不到再备用；PDB优先，失败再 AlphaFold
"""

import requests
import json
import os
import time
import sys
from typing import List, Dict, Optional, Tuple
from urllib.parse import quote
import argparse
import re

# 新增：导入Bio模块用于ExPASy查询
try:
    from Bio import ExPASy, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False
    print("⚠️  Bio模块未安装，将跳过ExPASy查询功能")


# -------------------------
# 仅过滤指定的 6 个 token + 去重（保序）
# -------------------------
EXCLUDE_TOKENS = {"===", "RDKIT", "RXNFP", "DRFP", "(COSINE)", "(EUCLIDEAN)"}


def dedupe_preserve_order(items: List[str]) -> List[str]:
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def sanitize_uniprot_ids(raw_ids: List[str]) -> List[str]:
    """
    只做两件事：
    1) upper/strip
    2) 过滤 EXCLUDE_TOKENS 这 6 个
    3) 去重（保序）

    不做任何额外的格式校验，避免误过滤真实ID
    """
    cleaned = []
    for x in raw_ids:
        t = (x or "").strip().upper()
        if not t:
            continue
        if t in EXCLUDE_TOKENS:
            continue
        cleaned.append(t)
    return dedupe_preserve_order(cleaned)


class ProteinDownloader:
    def __init__(self, output_dir: str = "./protein_structures"):
        """
        初始化下载器
        """
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if not os.path.isabs(output_dir):
            self.output_dir = os.path.join(script_dir, output_dir)
        else:
            self.output_dir = output_dir

        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })

        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "pdb"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "alphafold"), exist_ok=True)

        self.stats = {
            'total': 0,
            'pdb_found': 0,
            'alphafold_found': 0,
            'not_found': 0,
            'errors': 0
        }

        self.results = []

        # v4：批量 UniProt->PDB 映射缓存
        self._bulk_pdb_map: Dict[str, List[str]] = {}

    # -------------------------
    # v4：批量 UniProt idmapping（一次提交全部）
    # -------------------------
    def _uniprot_idmapping_run(self, ids: List[str], from_db: str, to_db: str) -> Optional[str]:
        try:
            url = "https://rest.uniprot.org/idmapping/run"
            data = {
                "from": from_db,
                "to": to_db,
                "ids": ",".join(ids)
            }
            resp = self.session.post(url, data=data, timeout=30)
            if not resp.ok:
                print(f"  ❌ UniProt idmapping/run 失败: HTTP {resp.status_code}, 响应: {resp.text[:200]}")
                return None
            return resp.json().get("jobId")
        except Exception as e:
            print(f"  ❌ UniProt idmapping/run 异常: {e}")
            return None

    def _uniprot_idmapping_poll_results(self, job_id: str, timeout_sec: int = 180, poll_interval: float = 2.0) -> Optional[Dict]:
        try:
            details_url = f"https://rest.uniprot.org/idmapping/details/{job_id}"
            results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"

            start = time.time()
            details = None
            while time.time() - start < timeout_sec:
                d = self.session.get(details_url, timeout=30)
                if d.ok:
                    details = d.json()
                    job_status = details.get("jobStatus") or details.get("status")
                    if job_status and str(job_status).upper() in ("FINISHED", "COMPLETE", "COMPLETED"):
                        break
                    if details.get("results") is not None or details.get("failedIds") is not None:
                        break
                time.sleep(poll_interval)

            r = self.session.get(results_url, timeout=60)
            if r.ok:
                return r.json()

            return details
        except Exception as e:
            print(f"  ❌ UniProt idmapping 轮询异常: {e}")
            return None

    def bulk_prepare_pdb_mapping(self, uniprot_ids: List[str]):
        """
        v4：批量预查询 UniProt->PDB 映射，并缓存
        """
        self._bulk_pdb_map = {}
        if not uniprot_ids:
            return

        print("\n🚀 预处理：使用UniProt ID Mapping（批量）查询 PDB 映射...")

        # 确保不包含那 6 个 token
        ids = sanitize_uniprot_ids(uniprot_ids)
        if not ids:
            print("  ⚠️  没有可用于批量映射的ID（可能全是被排除的token）")
            return

        job_id = self._uniprot_idmapping_run(ids, from_db="UniProtKB_AC-ID", to_db="PDB")
        if not job_id:
            print("  ⚠️  批量映射任务提交失败，将在单个ID处理时走备用方案")
            return

        print(f"  ✅ 已提交映射任务 jobId={job_id}，等待结果...")

        data = self._uniprot_idmapping_poll_results(job_id)
        if not data:
            print("  ⚠️  批量映射未获得结果，将在单个ID处理时走备用方案")
            return

        results = data.get("results") or []
        for item in results:
            fr = (item.get("from") or "").strip().upper()
            to = (item.get("to") or "").strip().upper()
            if not fr or not to:
                continue
            self._bulk_pdb_map.setdefault(fr, []).append(to)

        for k, v in list(self._bulk_pdb_map.items()):
            self._bulk_pdb_map[k] = dedupe_preserve_order([x.strip().upper() for x in v if x.strip()])

        print(f"  📌 批量映射完成：获得 {len(self._bulk_pdb_map)} 个UniProt ID 的PDB映射")

    # -------------------------
    # 下载逻辑（保持原逻辑）
    # -------------------------
    def download_all_pdb_structures(self, uniprot_id: str, pdb_ids: List[str]) -> List[str]:
        downloaded_files = []

        for pdb_id in pdb_ids:
            download_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            filename = os.path.join(self.output_dir, "pdb", f"{uniprot_id}_{pdb_id}.pdb")

            print(f"  📥 下载PDB结构: {pdb_id}")
            if self.download_file(download_url, filename):
                print(f"  ✅ 下载成功: {filename}")
                downloaded_files.append(filename)
            else:
                print(f"  ❌ 下载失败: {pdb_id}")

        return downloaded_files

    def download_file(self, url: str, filename: str) -> bool:
        try:
            response = self.session.get(url, timeout=60, stream=True)
            response.raise_for_status()

            with open(filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            return True

        except Exception as e:
            print(f"  ❌ 下载失败: {e}")
            return False

    # -------------------------
    # PDB 查询：优先批量映射缓存 -> 备用RCSB -> 备用ExPASy
    # -------------------------
    def check_pdb_structure(self, uniprot_id: str) -> Tuple[bool, Optional[List[str]], Optional[Dict]]:
        uniprot_id = uniprot_id.strip().upper()

        # 如果是那 6 个 token，直接不查（保险）
        if uniprot_id in EXCLUDE_TOKENS:
            return False, None, None

        # 1) 优先批量映射缓存
        if uniprot_id in self._bulk_pdb_map:
            pdb_ids = self._bulk_pdb_map.get(uniprot_id) or []
            if pdb_ids:
                print(f"    ✅（批量映射缓存）找到PDB IDs: {pdb_ids}")

                details = None
                try:
                    detail_response = self.session.get(
                        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_ids[0]}",
                        timeout=30
                    )
                    if detail_response.ok:
                        details = detail_response.json()
                except Exception as e:
                    print(f"  ⚠️  获取PDB详细信息失败: {e}")

                return True, pdb_ids, details

        # 2) 备用：RCSB Search API（保留你的原逻辑）
        pdb_found, pdb_ids, pdb_details = self._check_pdb_rcsb(uniprot_id)
        if pdb_found:
            return pdb_found, pdb_ids, pdb_details

        # 3) 备用：ExPASy（保留你的原逻辑）
        return self._check_pdb_expasy(uniprot_id)

    def _check_pdb_rcsb(self, uniprot_id: str) -> Tuple[bool, Optional[List[str]], Optional[Dict]]:
        try:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "pager": {"start": 0, "rows": 10}
                }
            }

            print(f"    查询参数: {uniprot_id}")
            response = self.session.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=query,
                timeout=30,
                headers={'Content-Type': 'application/json'}
            )

            print(f"    PDB API响应状态: {response.status_code}")
            if not response.ok:
                print(f"  ❌ PDB查询失败: HTTP {response.status_code}, 响应: {response.text[:200]}")
                return False, None, None

            data = response.json()
            print(f"    PDB查询结果数量: {len(data.get('result_set', []))}")

            if data.get('result_set') and len(data['result_set']) > 0:
                pdb_ids = [entry['identifier'] for entry in data['result_set']]
                print(f"    找到PDB IDs: {pdb_ids}")

                details = None
                if pdb_ids:
                    try:
                        detail_response = self.session.get(
                            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_ids[0]}",
                            timeout=30
                        )
                        if detail_response.ok:
                            details = detail_response.json()
                    except Exception as e:
                        print(f"  ⚠️  获取PDB详细信息失败: {e}")

                return True, pdb_ids, details

            return False, None, None

        except Exception as e:
            print(f"  ❌ PDB查询异常: {e}")
            return False, None, None

    def _check_pdb_expasy(self, uniprot_id: str) -> Tuple[bool, Optional[List[str]], Optional[Dict]]:
        if not BIO_AVAILABLE:
            print("  ❌ Bio模块不可用，跳过ExPASy查询")
            return False, None, None

        try:
            print(f"  🧬 使用ExPASy查询UniProt记录: {uniprot_id}")

            handle = ExPASy.get_sprot_raw(uniprot_id)
            record = SeqIO.read(handle, 'swiss')

            pdb_ids = []
            for ref in record.dbxrefs:
                if 'PDB' in ref:
                    pdb_id = ref.split(':')[1]
                    pdb_ids.append(pdb_id)

            if pdb_ids:
                print(f"    ExPASy方法找到PDB IDs: {pdb_ids}")
                return True, pdb_ids, None
            else:
                print("    ExPASy方法未找到PDB引用")
                return False, None, None

        except Exception as e:
            print(f"  ❌ ExPASy查询异常: {e}")
            return False, None, None

    # -------------------------
    # AlphaFold（保持原逻辑）
    # -------------------------
    def check_alphafold_structure(self, uniprot_id: str) -> Tuple[bool, Optional[str], Optional[Dict]]:
        try:
            response = self.session.get(
                f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}",
                timeout=30
            )

            if response.ok:
                data = response.json()
                if data and len(data) > 0:
                    entry = data[0]
                    pdb_url = entry.get('pdbUrl')
                    return True, pdb_url, entry

            return False, None, None

        except Exception as e:
            print(f"  ❌ AlphaFold查询异常: {e}")
            return False, None, None

    # -------------------------
    # 单个ID处理（保持原逻辑，但排除那 6 个 token）
    # -------------------------
    def process_uniprot_id(self, uniprot_id: str) -> Optional[Dict]:
        uniprot_id = (uniprot_id or "").strip().upper()

        # 只跳过你指定的 6 个
        if uniprot_id in EXCLUDE_TOKENS:
            print(f"\n⏭️  跳过非UniProt条目: {uniprot_id}")
            return None

        print(f"\n🔍 处理 {uniprot_id}...")

        result = {
            'uniprot_id': uniprot_id,
            'status': 'not_found',
            'source': None,
            'filename': None,
            'pdb_ids': None,
            'details': {}
        }

        print("  📊 检查PDB（优先UniProt idmapping批量缓存）...")
        pdb_found, pdb_ids, pdb_details = self.check_pdb_structure(uniprot_id)

        if pdb_found and pdb_ids:
            print(f"  ✅ 找到PDB结构: {', '.join(pdb_ids)}")

            print(f"  📥 开始下载 {len(pdb_ids)} 个PDB结构...")
            downloaded_files = self.download_all_pdb_structures(uniprot_id, pdb_ids)

            if downloaded_files:
                print(f"  ✅ 成功下载 {len(downloaded_files)} 个PDB结构")
                result.update({
                    'status': 'found_pdb',
                    'source': 'PDB',
                    'filename': downloaded_files,
                    'pdb_ids': pdb_ids,
                    'download_count': len(downloaded_files),
                    'details': {
                        'title': pdb_details.get('struct', {}).get('title') if pdb_details else None,
                        'resolution': pdb_details.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0] if pdb_details else None,
                        'method': pdb_details.get('exptl', [{}])[0].get('method') if pdb_details else None
                    }
                })
                self.stats['pdb_found'] += 1
                return result
            else:
                print("  ❌ 所有PDB结构下载失败，尝试AlphaFold...")

        print("  🧬 检查AlphaFold...")
        af_found, af_url, af_details = self.check_alphafold_structure(uniprot_id)

        if af_found and af_url:
            print("  ✅ 找到AlphaFold结构")

            filename = os.path.join(self.output_dir, "alphafold", f"{uniprot_id}_alphafold.pdb")

            print("  📥 下载AlphaFold结构")
            if self.download_file(af_url, filename):
                print(f"  ✅ 下载成功: {filename}")
                result.update({
                    'status': 'found_alphafold',
                    'source': 'AlphaFold',
                    'filename': filename,
                    'details': {
                        'confidence': af_details.get('confidenceScore') if af_details else None,
                        'model_version': af_details.get('modelCreatedDate') if af_details else None
                    }
                })
                self.stats['alphafold_found'] += 1
                return result
            else:
                print("  ❌ AlphaFold下载失败")
                result['status'] = 'error'
                self.stats['errors'] += 1
                return result

        print("  ❌ 未找到任何结构")
        self.stats['not_found'] += 1
        return result

    # -------------------------
    # 批量处理：只过滤 6 个 token + 去重，然后批量映射，然后逐个处理
    # -------------------------
    def process_uniprot_list(self, uniprot_ids: List[str], delay: float = 1.0) -> List[Dict]:
        uniprot_ids = sanitize_uniprot_ids(uniprot_ids)

        print(f"🚀 开始处理 {len(uniprot_ids)} 个蛋白质...")
        print(f"📁 输出目录: {os.path.abspath(self.output_dir)}")

        self.stats['total'] = len(uniprot_ids)

        self.bulk_prepare_pdb_mapping(uniprot_ids)

        for i, uniprot_id in enumerate(uniprot_ids, 1):
            print(f"\n{'='*50}")
            print(f"进度: {i}/{len(uniprot_ids)}")

            result = self.process_uniprot_id(uniprot_id)
            if result is not None:
                self.results.append(result)

            if i < len(uniprot_ids):
                time.sleep(delay)

        return self.results

    # -------------------------
    # 报告/摘要（保持原逻辑）
    # -------------------------
    def save_report(self, filename: str = None) -> str:
        if filename is None:
            filename = os.path.join(self.output_dir, "download_report.txt")

        with open(filename, 'w', encoding='utf-8') as f:
            f.write("蛋白质结构下载报告\n")
            f.write("=" * 50 + "\n\n")

            f.write("统计信息:\n")
            f.write(f"总计: {self.stats['total']}\n")
            f.write(f"PDB结构: {self.stats['pdb_found']}\n")
            f.write(f"AlphaFold结构: {self.stats['alphafold_found']}\n")
            f.write(f"未找到: {self.stats['not_found']}\n")
            f.write(f"错误: {self.stats['errors']}\n\n")

            f.write("详细结果:\n")
            f.write("-" * 50 + "\n")

            for result in self.results:
                f.write(f"\nUniProt ID: {result['uniprot_id']}\n")
                f.write(f"状态: {result['status']}\n")
                f.write(f"数据源: {result['source'] or 'N/A'}\n")

                if result['filename']:
                    if isinstance(result['filename'], list):
                        f.write(f"下载文件数量: {len(result['filename'])}\n")
                        for i, fn in enumerate(result['filename'], 1):
                            f.write(f"文件{i}: {fn}\n")
                    else:
                        f.write(f"文件: {result['filename']}\n")
                else:
                    f.write("文件: N/A\n")

                if result['pdb_ids']:
                    f.write(f"PDB IDs: {', '.join(result['pdb_ids'])}\n")

                if result.get('download_count'):
                    f.write(f"成功下载数量: {result['download_count']}\n")

                if result['details']:
                    for key, value in result['details'].items():
                        if value is not None:
                            f.write(f"{key}: {value}\n")

                f.write("-" * 30 + "\n")

        return filename

    def print_summary(self):
        print(f"\n{'='*60}")
        print("📊 处理完成！统计摘要:")
        print(f"{'='*60}")
        print(f"总计: {self.stats['total']}")
        print(f"✅ PDB结构: {self.stats['pdb_found']}")
        print(f"🧬 AlphaFold结构: {self.stats['alphafold_found']}")
        print(f"❌ 未找到: {self.stats['not_found']}")
        print(f"⚠️  错误: {self.stats['errors']}")
        print(f"\n📁 文件保存在: {os.path.abspath(self.output_dir)}")

        success_rate = ((self.stats['pdb_found'] + self.stats['alphafold_found']) / self.stats['total'] * 100) if self.stats['total'] > 0 else 0
        print(f"🎯 成功率: {success_rate:.1f}%")


def read_uniprot_ids_from_file(filename: str) -> List[str]:
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read()

        uniprot_ids = re.split(r'[,;\s\n]+', content)
        uniprot_ids = [uid.strip().upper() for uid in uniprot_ids if uid.strip()]

        # 只过滤你指定的 6 个 + 去重（保序）
        return sanitize_uniprot_ids(uniprot_ids)

    except Exception as e:
        print(f"❌ 读取文件失败: {e}")
        return []


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input_file = os.path.join(script_dir, 'uniprot_ids.txt')

    print(f"🔍 脚本目录: {script_dir}")
    print(f"📁 查找输入文件: {default_input_file}")

    if os.path.exists(default_input_file):
        print(f"✅ 找到输入文件: {default_input_file}")
        uniprot_ids = read_uniprot_ids_from_file(default_input_file)
    else:
        print(f"❌ 未找到输入文件: {default_input_file}")
        print("请在脚本目录下创建 uniprot_ids.txt 文件，每行一个UniProt ID")
        return

    if not uniprot_ids:
        print("❌ 输入文件中没有有效的UniProt ID（可能全是被排除的token）")
        return

    print(f"📋 将处理 {len(uniprot_ids)} 个UniProt ID（已去重，并排除指定token）:")
    for uid in uniprot_ids[:10]:
        print(f"  - {uid}")
    if len(uniprot_ids) > 10:
        print(f"  ... 还有 {len(uniprot_ids) - 10} 个")

    downloader = ProteinDownloader("protein_structures")

    try:
        downloader.process_uniprot_list(uniprot_ids, 1.0)

        report_file = downloader.save_report()
        print(f"\n📄 详细报告已保存: {report_file}")

        downloader.print_summary()

    except KeyboardInterrupt:
        print("\n❌ 用户中断操作")
        downloader.print_summary()
    except Exception as e:
        print(f"\n❌ 程序异常: {e}")
        downloader.print_summary()


if __name__ == "__main__":
    main()

'''PS C:\Users\Jiaya> & D:/conda_envs/reme-try/python.exe "d:/Project/pythonProject/pythonProject/reactionsearch - 副本/download_structures.py"
🔍 脚本目录: d:\Project\pythonProject\pythonProject\reactionsearch - 副本
📁 查找输入文件: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\uniprot_ids.txt
✅ 找到输入文件: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\uniprot_ids.txt
📋 将处理 39 个UniProt ID（已去重，并排除指定token）:
  - F9UM18
  - A0A1B1WA13
  - A0A1B1WA14
  - Q93YP7
  - Q8W405
  - Q8W404
  - Q0QFI7
  - A0A168S819
  - A0A168S850
  - Q6Q874
  ... 还有 29 个
🚀 开始处理 39 个蛋白质...
📁 输出目录: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures

🚀 预处理：使用UniProt ID Mapping（批量）查询 PDB 映射...
  ✅ 已提交映射任务 jobId=ddkfApZTTR，等待结果...
  📌 批量映射完成：获得 7 个UniProt ID 的PDB映射

==================================================
进度: 1/39

🔍 处理 F9UM18...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['7EBO']
  ✅ 找到PDB结构: 7EBO
  📥 开始下载 1 个PDB结构...
  📥 下载PDB结构: 7EBO
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\F9UM18_7EBO.pdb
  ✅ 成功下载 1 个PDB结构

==================================================
进度: 2/39

🔍 处理 A0A1B1WA13...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A1B1WA13
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A1B1WA13
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A1B1WA13_alphafold.pdb

==================================================
进度: 3/39

🔍 处理 A0A1B1WA14...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A1B1WA14
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A1B1WA14
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A1B1WA14_alphafold.pdb

==================================================
进度: 4/39

🔍 处理 Q93YP7...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q93YP7
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q93YP7
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q93YP7_alphafold.pdb

==================================================
进度: 5/39

🔍 处理 Q8W405...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q8W405
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q8W405
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q8W405_alphafold.pdb

==================================================
进度: 6/39

🔍 处理 Q8W404...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q8W404
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q8W404
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q8W404_alphafold.pdb

==================================================
进度: 7/39

🔍 处理 Q0QFI7...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q0QFI7
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q0QFI7
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q0QFI7_alphafold.pdb

==================================================
进度: 8/39

🔍 处理 A0A168S819...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A168S819
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A168S819
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A168S819_alphafold.pdb

==================================================
进度: 9/39

🔍 处理 A0A168S850...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A168S850
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A168S850
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A168S850_alphafold.pdb

==================================================
进度: 10/39

🔍 处理 Q6Q874...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q6Q874
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q6Q874
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q6Q874_alphafold.pdb

==================================================
进度: 11/39

🔍 处理 P76469...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['2VWS', '2VWT']
  ✅ 找到PDB结构: 2VWS, 2VWT
  📥 开始下载 2 个PDB结构...
  📥 下载PDB结构: 2VWS
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P76469_2VWS.pdb
  📥 下载PDB结构: 2VWT
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P76469_2VWT.pdb
  ✅ 成功下载 2 个PDB结构

==================================================
进度: 12/39

🔍 处理 B9A1Q4...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: B9A1Q4
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: B9A1Q4
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\B9A1Q4_alphafold.pdb

==================================================
进度: 13/39

🔍 处理 I1NJA2...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: I1NJA2
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: I1NJA2
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\I1NJA2_alphafold.pdb

==================================================
进度: 14/39

🔍 处理 A0A1L7NRS0...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A1L7NRS0
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A1L7NRS0
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A1L7NRS0_alphafold.pdb

==================================================
进度: 15/39

🔍 处理 Q53B70...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q53B70
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q53B70
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q53B70_alphafold.pdb

==================================================
进度: 16/39

🔍 处理 Q93XE6...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q93XE6
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q93XE6
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q93XE6_alphafold.pdb

==================================================
进度: 17/39

🔍 处理 Q8GHB2...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['2XLQ', '2XLY', '2XM5', '2XM7']
  ✅ 找到PDB结构: 2XLQ, 2XLY, 2XM5, 2XM7
  📥 开始下载 4 个PDB结构...
  📥 下载PDB结构: 2XLQ
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q8GHB2_2XLQ.pdb
  📥 下载PDB结构: 2XLY
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q8GHB2_2XLY.pdb
  📥 下载PDB结构: 2XM5
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q8GHB2_2XM5.pdb
  📥 下载PDB结构: 2XM7
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q8GHB2_2XM7.pdb
  ✅ 成功下载 4 个PDB结构

==================================================
进度: 18/39

🔍 处理 Q9AJC6...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q9AJC6
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q9AJC6
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q9AJC6_alphafold.pdb

==================================================
进度: 19/39

🔍 处理 C5H429...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['5DXX', '5DXY', '5DY2', '5DY3']
  ✅ 找到PDB结构: 5DXX, 5DXY, 5DY2, 5DY3
  📥 开始下载 4 个PDB结构...
  📥 下载PDB结构: 5DXX
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\C5H429_5DXX.pdb
  📥 下载PDB结构: 5DXY
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\C5H429_5DXY.pdb
  📥 下载PDB结构: 5DY2
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\C5H429_5DY2.pdb
  📥 下载PDB结构: 5DY3
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\C5H429_5DY3.pdb
  ✅ 成功下载 4 个PDB结构

==================================================
进度: 20/39

🔍 处理 P40952...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P40952
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P40952
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P40952_alphafold.pdb

==================================================
进度: 21/39

🔍 处理 Q9R9V9...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['5CPL', '5CPM', '5CPN', '5CPO', '5N6Q', '8A8I', '8AU8', '8AU9', '8AUF', '8AUG', '8AUH', '8AUI']
  ✅ 找到PDB结构: 5CPL, 5CPM, 5CPN, 5CPO, 5N6Q, 8A8I, 8AU8, 8AU9, 8AUF, 8AUG, 8AUH, 8AUI
  📥 开始下载 12 个PDB结构...
  📥 下载PDB结构: 5CPL
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_5CPL.pdb
  📥 下载PDB结构: 5CPM
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_5CPM.pdb
  📥 下载PDB结构: 5CPN
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_5CPN.pdb
  📥 下载PDB结构: 5CPO
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_5CPO.pdb
  📥 下载PDB结构: 5N6Q
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_5N6Q.pdb
  📥 下载PDB结构: 8A8I
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8A8I.pdb
  📥 下载PDB结构: 8AU8
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AU8.pdb
  📥 下载PDB结构: 8AU9
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AU9.pdb
  📥 下载PDB结构: 8AUF
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AUF.pdb
  📥 下载PDB结构: 8AUG
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AUG.pdb
  📥 下载PDB结构: 8AUH
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AUH.pdb
  📥 下载PDB结构: 8AUI
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9R9V9_8AUI.pdb
  ✅ 成功下载 12 个PDB结构

==================================================
进度: 22/39

🔍 处理 A0A0H4SN47...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['6YC8']
  ✅ 找到PDB结构: 6YC8
  📥 开始下载 1 个PDB结构...
  📥 下载PDB结构: 6YC8
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\A0A0H4SN47_6YC8.pdb
  ✅ 成功下载 1 个PDB结构

==================================================
进度: 23/39

🔍 处理 Q93TU6...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    ✅（批量映射缓存）找到PDB IDs: ['1O8U']
  ✅ 找到PDB结构: 1O8U
  📥 开始下载 1 个PDB结构...
  📥 下载PDB结构: 1O8U
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q93TU6_1O8U.pdb
  ✅ 成功下载 1 个PDB结构

==================================================
进度: 24/39

🔍 处理 Q8YNV6...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q8YNV6
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q8YNV6
  ❌ ExPASy查询异常: No records found in handle
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q8YNV6_alphafold.pdb

==================================================
进度: 25/39

🔍 处理 Q9BXD5...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q9BXD5
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q9BXD5
    ExPASy方法找到PDB IDs: ['6ARH', '6ARH']
  ✅ 找到PDB结构: 6ARH, 6ARH
  📥 开始下载 2 个PDB结构...
  📥 下载PDB结构: 6ARH
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9BXD5_6ARH.pdb
  📥 下载PDB结构: 6ARH
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q9BXD5_6ARH.pdb
  ✅ 成功下载 2 个PDB结构

==================================================
进度: 26/39

🔍 处理 P04964...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P04964
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P04964
    ExPASy方法找到PDB IDs: ['1DAP', '1F06', '2DAP', '3DAP', '5LOA', '5LOC', '1DAP', '1F06', '2DAP', '3DAP', '5LOA', '5LOC']
  ✅ 找到PDB结构: 1DAP, 1F06, 2DAP, 3DAP, 5LOA, 5LOC, 1DAP, 1F06, 2DAP, 3DAP, 5LOA, 5LOC
  📥 开始下载 12 个PDB结构...
  📥 下载PDB结构: 1DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_1DAP.pdb
  📥 下载PDB结构: 1F06
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_1F06.pdb
  📥 下载PDB结构: 2DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_2DAP.pdb
  📥 下载PDB结构: 3DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_3DAP.pdb
  📥 下载PDB结构: 5LOA
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_5LOA.pdb
  📥 下载PDB结构: 5LOC
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_5LOC.pdb
  📥 下载PDB结构: 1DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_1DAP.pdb
  📥 下载PDB结构: 1F06
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_1F06.pdb
  📥 下载PDB结构: 2DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_2DAP.pdb
  📥 下载PDB结构: 3DAP
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_3DAP.pdb
  📥 下载PDB结构: 5LOA
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_5LOA.pdb
  📥 下载PDB结构: 5LOC
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P04964_5LOC.pdb
  ✅ 成功下载 12 个PDB结构

==================================================
进度: 27/39

🔍 处理 F7YTT4...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: F7YTT4
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: F7YTT4
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\F7YTT4_alphafold.pdb

==================================================
进度: 28/39

🔍 处理 Q65CJ7...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q65CJ7
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q65CJ7
    ExPASy方法找到PDB IDs: ['3BA1', '3BAZ', '3BA1', '3BAZ']
  ✅ 找到PDB结构: 3BA1, 3BAZ, 3BA1, 3BAZ
  📥 开始下载 4 个PDB结构...
  📥 下载PDB结构: 3BA1
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q65CJ7_3BA1.pdb
  📥 下载PDB结构: 3BAZ
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q65CJ7_3BAZ.pdb
  📥 下载PDB结构: 3BA1
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q65CJ7_3BA1.pdb
  📥 下载PDB结构: 3BAZ
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\Q65CJ7_3BAZ.pdb
  ✅ 成功下载 4 个PDB结构

==================================================
进度: 29/39

🔍 处理 P51102...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51102
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51102
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P51102_alphafold.pdb

==================================================
进度: 30/39

🔍 处理 P51103...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51103
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51103
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P51103_alphafold.pdb

==================================================
进度: 31/39

🔍 处理 P51104...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51104
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51104
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P51104_alphafold.pdb

==================================================
进度: 32/39

🔍 处理 P51110...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51110
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51110
    ExPASy方法找到PDB IDs: ['2IOD', '2NNL', '3BXX', '3C1T', '2IOD', '2NNL', '3BXX', '3C1T']
  ✅ 找到PDB结构: 2IOD, 2NNL, 3BXX, 3C1T, 2IOD, 2NNL, 3BXX, 3C1T
  📥 开始下载 8 个PDB结构...
  📥 下载PDB结构: 2IOD
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_2IOD.pdb
  📥 下载PDB结构: 2NNL
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_2NNL.pdb
  📥 下载PDB结构: 3BXX
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_3BXX.pdb
  📥 下载PDB结构: 3C1T
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_3C1T.pdb
  📥 下载PDB结构: 2IOD
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_2IOD.pdb
  📥 下载PDB结构: 2NNL
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_2NNL.pdb
  📥 下载PDB结构: 3BXX
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_3BXX.pdb
  📥 下载PDB结构: 3C1T
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\pdb\P51110_3C1T.pdb
  ✅ 成功下载 8 个PDB结构

==================================================
进度: 33/39

🔍 处理 P51105...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51105
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51105
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P51105_alphafold.pdb

==================================================
进度: 34/39

🔍 处理 P51107...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P51107
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P51107
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P51107_alphafold.pdb

==================================================
进度: 35/39

🔍 处理 Q6UAQ7...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: Q6UAQ7
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: Q6UAQ7
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\Q6UAQ7_alphafold.pdb

==================================================
进度: 36/39

🔍 处理 A0A0A0PXZ7...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A0A0PXZ7
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A0A0PXZ7
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A0A0PXZ7_alphafold.pdb

==================================================
进度: 37/39

🔍 处理 A0A0C6ENR6...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: A0A0C6ENR6
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: A0A0C6ENR6
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\A0A0C6ENR6_alphafold.pdb

==================================================
进度: 38/39

🔍 处理 M4DU26...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: M4DU26
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: M4DU26
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\M4DU26_alphafold.pdb

==================================================
进度: 39/39

🔍 处理 P44430...
  📊 检查PDB（优先UniProt idmapping批量缓存）...
    查询参数: P44430
    PDB API响应状态: 400
  ❌ PDB查询失败: HTTP 400, 响应: {
  "status" : 400,
  "message" : "JSON schema validation failed for query: {\"query\":{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_polymer_entity_container_identif
  🧬 使用ExPASy查询UniProt记录: P44430
    ExPASy方法未找到PDB引用
  🧬 检查AlphaFold...
  ✅ 找到AlphaFold结构
  📥 下载AlphaFold结构
  ✅ 下载成功: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\alphafold\P44430_alphafold.pdb

📄 详细报告已保存: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures\download_report.txt

============================================================
📊 处理完成！统计摘要:
============================================================
总计: 39
✅ PDB结构: 11
🧬 AlphaFold结构: 28
❌ 未找到: 0
⚠️  错误: 0

📁 文件保存在: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\protein_structures
🎯 成功率: 100.0%'''