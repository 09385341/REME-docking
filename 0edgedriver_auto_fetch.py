#成功获取 Edge 版本并下载对应的 EdgeDriver，支持 Windows 32/64 位和 ARM64 架构。
import os
import zipfile
import platform
import requests
import winreg

def get_edge_version_windows() -> str | None:
    """
    从注册表读取 Edge 版本（优先当前用户，其次本机）
    返回如 '144.0.3719.115'，失败返回 None
    """
    keys = [
        (winreg.HKEY_CURRENT_USER, r"Software\\Microsoft\\Edge\\BLBeacon"),
        (winreg.HKEY_LOCAL_MACHINE, r"Software\\Microsoft\\Edge\\BLBeacon"),
        (winreg.HKEY_LOCAL_MACHINE, r"Software\\WOW6432Node\\Microsoft\\Edge\\BLBeacon"),
    ]
    for hive, subkey in keys:
        try:
            with winreg.OpenKey(hive, subkey) as k:
                val, _ = winreg.QueryValueEx(k, "version")
                if isinstance(val, str) and val.strip():
                    return val.strip()
        except OSError:
            continue
    return None

def get_windows_arch() -> str:
    """
    返回 'win64' / 'win32' / 'arm64'
    """
    m = platform.machine().lower()
    arch, _ = platform.architecture()
    if 'arm' in m:
        return 'arm64'
    if arch == '32bit':
        return 'win32'
    return 'win64'

def build_candidate_urls(version: str, arch: str) -> list[str]:
    """
    构造可能可用的官方下载链接（按优先级尝试）
    1) msedgewebdriverstorage.blob.core.windows.net
    2) msedgedriver.azureedge.net
    """
    filename = {
        'win64': 'edgedriver_win64.zip',
        'win32': 'edgedriver_win32.zip',
        'arm64': 'edgedriver_arm64.zip',
    }[arch]
    return [
        f"https://msedgewebdriverstorage.blob.core.windows.net/edgewebdriver/{version}/{filename}",
        f"https://msedgedriver.azureedge.net/{version}/{filename}",
    ]

def download_file(url: str, dst: str, timeout: int = 30) -> bool:
    """
    下载文件到 dst，成功返回 True，否则 False（不抛异常，便于依次尝试多个 URL）
    """
    try:
        with requests.get(url, timeout=timeout, stream=True) as r:
            if r.status_code != 200:
                return False
            total = int(r.headers.get("Content-Length", "0")) or None
            with open(dst, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 256):
                    if chunk:
                        f.write(chunk)
        return True
    except Exception:
        return False

def download_edgedriver_to_script_dir() -> str:
    """
    将 EdgeDriver 下载并解压到脚本所在目录，返回 msedgedriver.exe 的绝对路径。
    若失败抛出异常。
    """
    version = get_edge_version_windows()
    if not version:
        raise RuntimeError("无法自动检测 Edge 版本。请在 Edge 地址栏输入 edge://version 查看版本，并手动下载。")

    arch = get_windows_arch()
    urls = build_candidate_urls(version, arch)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(script_dir, exist_ok=True)

    zip_name = f"msedgedriver_{version}_{arch}.zip"
    zip_path = os.path.join(script_dir, zip_name)

    # 逐个 URL 尝试下载
    ok = False
    for url in urls:
        print(f"尝试下载 EdgeDriver：{url}")
        ok = download_file(url, zip_path, timeout=40)
        if ok:
            print(f"已下载：{zip_path}")
            break
        print("下载失败，尝试下一个地址...")

    if not ok:
        raise RuntimeError("所有下载地址均失败。可能是 DNS/网络限制或版本不存在。")

    # 解压到脚本目录
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(script_dir)
    print(f"已解压到脚本目录：{script_dir}")

    # 定位 msedgedriver.exe
    driver_path = os.path.join(script_dir, "msedgedriver.exe")
    if not os.path.exists(driver_path):
        # 某些包把可执行文件放在子目录，递归查找
        for root, _, files in os.walk(script_dir):
            if "msedgedriver.exe" in files:
                driver_path = os.path.join(root, "msedgedriver.exe")
                break

    if not os.path.exists(driver_path):
        raise RuntimeError("解压后未找到 msedgedriver.exe，请检查压缩包内容。")

    print(f"EdgeDriver 就绪：{driver_path}")
    return os.path.abspath(driver_path)

if __name__ == "__main__":
    try:
        p = download_edgedriver_to_script_dir()
        print("成功:", p)
    except Exception as e:
        print("失败：", e)
'''PS C:\Users\Jiaya> & D:/conda_envs/reme-try/python.exe d:/Project/pythonProject/pythonProject/reactionsearch/edgedriver_auto_fetch_Version2.py
尝试下载 EdgeDriver：https://msedgewebdriverstorage.blob.core.windows.net/edgewebdriver/144.0.3719.115/edgedriver_win64.zip
已下载：d:\Project\pythonProject\pythonProject\reactionsearch\msedgedriver_144.0.3719.115_win64.zip
已解压到脚本目录：d:\Project\pythonProject\pythonProject\reactionsearch
EdgeDriver 就绪：d:\Project\pythonProject\pythonProject\reactionsearch\msedgedriver.exe
成功: d:\Project\pythonProject\pythonProject\reactionsearch\msedgedriver.exe'''