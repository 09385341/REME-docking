import os
import shutil
import time
import re
import traceback
from datetime import datetime

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.action_chains import ActionChains

# ------------------ 配置项 ------------------
# 点击 Execute 后，在出现选项卡前后都适当延时，确保结果稳定
RESULTS_LOAD_DELAY = 8.0               # 选项卡出现后的基础等待（秒）
BEFORE_TAB_CLICK_EXTRA_DELAY = 6.0     # 准备点击选项卡之前再多等（秒）
POST_TAB_CLICK_WAIT = 5.0              # 点击选项卡之后的等待（秒），再开始抓取

# 打开每个 protein 详情窗口后，为了等内容加载更充分
DETAIL_WINDOW_INITIAL_WAIT = 3.5       # 切到详情窗口后的初次等待（秒）
AFTER_OPEN_ALL_LINKS_WAIT = 3.0        # 当前页点击所有 protein 链接后，等待窗口全部打开（秒）
UNIPROT_EXTRA_WAIT_IF_EMPTY = 3.0      # 如果初次没找到 Uniprot 链接，再额外等待后重试（秒）

WINDOW_TITLE_BLACKLIST = [
    "brenda enzyme database",
    "creative commons",
    "deed - attribution 4.0 international - creative commons",
    "网络违法犯罪信息举报网站",
    "中央网络安全和信息化委员会办公室",
    "违法和不良信息举报中心",
    "cyberpolice",
]
LINK_HREF_BLACKLIST = [
    "brenda-enzymes.org",
    "brenda.bio",
    "creativecommons.org",
    "12377.cn",
    "cac.gov.cn",
    "cyberpolice.cn",
    "cyberpolice",
]

# ------------------ 基础工具与���动 ------------------
def snapshot(driver, prefix="error"):
    """在当前目录保存截图，便于定位问题"""
    try:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"{prefix}_{ts}.png")
        if driver:
            driver.save_screenshot(path)
        print(f"[DEBUG] 已保存截图：{path}")
    except Exception as e:
        print(f"[DEBUG] 保存截图失败：{e}")

def create_edge_driver(options: Options) -> webdriver.Edge:
    """优先使用本地 EdgeDriver，找不到则让 Selenium Manager 处理"""
    candidates = []
    for env_name in ("MSEDGEDRIVER", "EDGE_DRIVER_PATH"):
        p = os.environ.get(env_name)
        if p:
            candidates.append(p)
    for name in ("msedgedriver", "msedgedriver.exe"):
        p = shutil.which(name)
        if p:
            candidates.append(p)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_driver = os.path.join(script_dir, "msedgedriver.exe")
    if os.path.exists(local_driver):
        candidates.append(local_driver)
    normalized = []
    for path in candidates:
        if os.path.isdir(path):
            normalized.append(os.path.join(path, "msedgedriver.exe"))
        else:
            normalized.append(path)
    for path in normalized:
        if path and os.path.exists(path):
            print(f"使用本地 EdgeDriver: {path}")
            try:
                return webdriver.Edge(service=Service(path), options=options)
            except Exception as e:
                print(f"使用本地 EdgeDriver 失败：{e}")
    print("本地未找到驱动，尝试使用 Selenium Manager ...")
    return webdriver.Edge(options=options)

def find_rxn_file():
    """查找脚本所在目录中的 .rxn 文件"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for file in os.listdir(script_dir):
        if file.endswith('.rxn'):
            return os.path.join(script_dir, file)
    return None

def get_total_pages(driver):
    """尽力获取分页总页数；失败时回落到动态翻页模式"""
    try:
        page_items = driver.find_elements(
            By.XPATH,
            "//li[contains(@class, 'ant-pagination-item') and not(contains(@class, 'ant-pagination-item-active'))]//a"
        )
        if page_items:
            nums = []
            for item in page_items:
                t = item.text.strip()
                if t.isdigit():
                    nums.append(int(t))
            if nums:
                total = max(nums)
                print(f"检测到总页数：{total}")
                return total
        active_page = driver.find_element(By.XPATH, "//li[contains(@class, 'ant-pagination-item-active')]")
        if active_page:
            print("只检测到当前页，将使用动态翻页模式")
            return None
    except Exception as e:
        print(f"获取总页数时发生错误：{e}")
    print("无法检测总页数，将使用动态翻页模式")
    return None

# ------------------ 复选框操作（参考 Version4） ------------------
def _get_card_element(wait, card_title_text):
    title_el = wait.until(EC.presence_of_element_located(
        (By.XPATH, f"//div[normalize-space()='{card_title_text}']")
    ))
    card_el = title_el.find_element(By.XPATH, "./ancestor::div[contains(@class,'ant-card')]")
    return card_el

def _get_checkbox_wrapper_in_card(card_el, label_text):
    return card_el.find_element(
        By.XPATH,
        f".//label[contains(@class,'ant-checkbox-wrapper')][.//span[normalize-space()='{label_text}']]"
    )

def _checkbox_is_checked(wrapper_el):
    """更鲁棒的选中判断：input + wrapper类 + 内部span类 + aria-checked"""
    try:
        input_el = wrapper_el.find_element(By.CSS_SELECTOR, "input.ant-checkbox-input[type='checkbox']")
        if input_el.is_selected():
            return True
        aria = (input_el.get_attribute("aria-checked") or "").lower()
        if aria == "true":
            return True
    except Exception:
        pass
    cls = (wrapper_el.get_attribute("class") or "")
    if "ant-checkbox-wrapper-checked" in cls:
        return True
    try:
        spans = wrapper_el.find_elements(By.CSS_SELECTOR, "span.ant-checkbox")
        for s in spans:
            if "ant-checkbox-checked" in (s.get_attribute("class") or ""):
                return True
    except Exception:
        pass
    return False

def _click_wrapper(driver, wrapper_el):
    driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", wrapper_el)
    time.sleep(0.15)
    try:
        wrapper_el.click()
    except Exception:
        driver.execute_script("arguments[0].click();", wrapper_el)
    time.sleep(0.15)

def _ensure_checkbox_state(driver, wrapper_el, desired_checked=True, max_tries=3):
    tries = 0
    while tries < max_tries and _checkbox_is_checked(wrapper_el) != desired_checked:
        _click_wrapper(driver, wrapper_el)
        tries += 1
    return _checkbox_is_checked(wrapper_el) == desired_checked

def clear_all_method_checkboxes(driver, wait):
    """将 RDKit/RXNFP/DRFP 三个卡片的相关复选框全部取消"""
    try:
        rdkit_card = _get_card_element(wait, "RDKit")
        for label in ["Select"]:
            try:
                w = _get_checkbox_wrapper_in_card(rdkit_card, label)
                _ensure_checkbox_state(driver, w, desired_checked=False)
            except Exception:
                pass
        rxnfp_card = _get_card_element(wait, "RXNFP")
        for label in ["Cosine", "Euclidean distance"]:
            try:
                w = _get_checkbox_wrapper_in_card(rxnfp_card, label)
                _ensure_checkbox_state(driver, w, desired_checked=False)
            except Exception:
                pass
        drfp_card = _get_card_element(wait, "DRFP")
        for label in ["Cosine", "Euclidean distance"]:
            try:
                w = _get_checkbox_wrapper_in_card(drfp_card, label)
                _ensure_checkbox_state(driver, w, desired_checked=False)
            except Exception:
                pass
        print("已清空所有方法复选框。")
    except Exception as e:
        print(f"清空复选框时出错：{e}")

def select_method_once(driver, wait, card_title, labels_to_check):
    """只勾选指定卡片的标签，并校验状态"""
    clear_all_method_checkboxes(driver, wait)
    try:
        card = _get_card_element(wait, card_title)
        for label in labels_to_check:
            w = _get_checkbox_wrapper_in_card(card, label)
            ok = _ensure_checkbox_state(driver, w, desired_checked=True)
            if ok:
                print(f"勾选成功：{card_title} - {label}")
            else:
                print(f"[WARN] 勾选可能未成功：{card_title} - {label}（将继续尝试执行）")
    except Exception as e:
        print(f"[WARN] 设置 {card_title} 复选框异常：{e}")
        snapshot(driver, "select_method_error")
    time.sleep(0.2)

# ------------------ Excute 按钮与选项卡（强化 RXNFP euclidean 点击） ------------------
def wait_execute_and_click(driver, timeout=25):
    """等待 Excute 按钮可点击且不禁用，然后点击"""
    try:
        btn = WebDriverWait(driver, timeout).until(EC.presence_of_element_located(
            (By.XPATH, "//button[span[normalize-space()='Excute']]")
        ))
        disabled = (btn.get_attribute("disabled") is not None) or ("ant-btn-disabled" in (btn.get_attribute("class") or ""))
        if disabled:
            WebDriverWait(driver, timeout).until(
                lambda d: (btn.get_attribute("disabled") is None) and ("ant-btn-disabled" not in (btn.get_attribute("class") or ""))
            )
        driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", btn)
        time.sleep(0.2)
        try:
            btn.click()
        except Exception:
            driver.execute_script("arguments[0].click();", btn)
        print("已点击 Excute 按钮。")
        return True
    except Exception as e:
        print(f"[WARN] 点击 Excute 失败：{e}")
        snapshot(driver, "execute_click_error")
        return False

def collect_tab_texts(driver):
    """返回当前页所有 role=tab 的可见文本（列表）"""
    tabs = driver.find_elements(By.XPATH, "//div[@role='tab']")
    texts = [t.text.strip() for t in tabs if t.text.strip()]
    print(f"[DEBUG] 检测到结果选项卡：{texts}")
    return texts, tabs

def click_results_tab(driver, algo, metric=None, timeout=35, retries=2):
    """
    鲁棒点击结果选项卡，优先按精确/同义匹配：
    例如 'RDKit'、'RXNFP (cosine)'、'RXNFP (euclidean)'、'RXNFP (Euclidean distance)'、'DRFP (euclidean)' 等。
    """
    WebDriverWait(driver, timeout).until(EC.presence_of_element_located((By.XPATH, "//div[@role='tab']")))
    # 在选项卡出现后进一步等待，确保列表加载充分
    time.sleep(RESULTS_LOAD_DELAY + BEFORE_TAB_CLICK_EXTRA_DELAY)

    metric_synonyms = {
        "euclidean": ["euclidean", "euclidean distance"],
        "cosine": ["cosine"]
    }
    target_algo = (algo or "").strip().lower()
    target_metric = (metric or "").strip().lower()
    metric_candidates = metric_synonyms.get(target_metric, [target_metric] if target_metric else [])

    for attempt in range(retries + 1):
        texts, tabs = collect_tab_texts(driver)
        norm_texts = [t.strip().lower() for t in texts]

        # 1) 精确匹配 "algo (metric)" 文本
        exact_targets = []
        if target_algo and target_metric:
            for m in metric_candidates:
                exact_targets.append(f"{target_algo} ({m})")
        elif target_algo:
            exact_targets.append(target_algo)

        for i, nt in enumerate(norm_texts):
            if nt in exact_targets:
                tab_el = tabs[i]
                print(f"[DEBUG] 精确匹配到选项卡：{texts[i]}")
                driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", tab_el)
                time.sleep(0.4)
                try:
                    tab_el.click()
                except Exception:
                    driver.execute_script("arguments[0].click();", tab_el)
                return texts[i]

        # 2) 同义词包含匹配：algo + (任一 metric 变体)
        if target_algo and metric_candidates:
            for i, nt in enumerate(norm_texts):
                if (target_algo in nt) and any(m in nt for m in metric_candidates):
                    tab_el = tabs[i]
                    print(f"[DEBUG] 同义匹配到选项卡：{texts[i]}")
                    driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", tab_el)
                    time.sleep(0.4)
                    try:
                        tab_el.click()
                    except Exception:
                        driver.execute_script("arguments[0].click();", tab_el)
                    return texts[i]

        # 3) 回退评分匹配（algo 权重更高，metric 次之）
        def score(text):
            t = text.strip().lower()
            s = 0
            if target_algo and target_algo in t: s += 10
            if target_metric and any(m in t for m in metric_candidates): s += 5
            return s
        scored = sorted([(score(txt), i, txt) for i, txt in enumerate(texts)], reverse=True)
        if scored and scored[0][0] > 0:
            best_idx = scored[0][1]
            best_tab = tabs[best_idx]
            print(f"[DEBUG] 评分匹配点击选项卡：{texts[best_idx]}（algo={target_algo}, metric={target_metric}）")
            driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", best_tab)
            time.sleep(0.4)
            try:
                best_tab.click()
            except Exception:
                driver.execute_script("arguments[0].click();", best_tab)
            return texts[best_idx]

        # 等待后重试
        if attempt < retries:
            print(f"[WARN] 未找到匹配的选项卡（尝试 {attempt+1}/{retries}），稍后重试...")
            time.sleep(2.0)

    print(f"[ERROR] 未找到匹配的选项卡（algo={algo}, metric={metric}）。")
    snapshot(driver, f"tab_not_found_{algo}_{metric}".replace(" ", "_"))
    return None

def go_back_to_selection(driver):
    """点击返回/Back 或浏览器后退，然后等待选择页出现"""
    clicked = False
    selectors = [
        "//button[.//span[normalize-space()='返回']]", "//button[.//span[normalize-space()='Back']]",
        "//button[normalize-space()='返回']", "//button[normalize-space()='Back']",
        "//a[.//span[normalize-space()='返回']]", "//a[.//span[normalize-space()='Back']]",
    ]
    for xp in selectors:
        try:
            btn = WebDriverWait(driver, 5).until(EC.element_to_be_clickable((By.XPATH, xp)))
            driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", btn)
            time.sleep(0.2)
            btn.click()
            clicked = True
            break
        except Exception:
            continue
    if not clicked:
        print("未找到返回按钮，使用浏览器后退。")
        driver.back()
    WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//div[normalize-space()='RDKit']")))
    print("已返回到选择页面。")

# ------------------ 黑名单工具 ------------------
def link_is_blacklisted(link):
    """根据 href/text 过滤掉明显的不相关链接，避免点击"""
    try:
        href = (link.get_attribute("href") or "").lower()
        text = (link.text or "").lower()
        for pat in LINK_HREF_BLACKLIST:
            if pat in href:
                return True
        text_pats = ["brenda", "creative commons", "举报", "网络安全", "中央网络安全", "不良信息", "cyberpolice"]
        return any(pat in text for pat in text_pats)
    except Exception:
        return False

def title_is_blacklisted(title):
    """窗口标题黑名单匹配（包含判断）"""
    t = (title or "").lower()
    for pat in WINDOW_TITLE_BLACKLIST:
        if pat.lower() in t:
            return True
    return False

# ------------------ 结果页抓取（不点击 Uniprot 链接；不去重） ------------------
def process_results_for_method(driver, algo, metric, categorized_ids):
    """
    点击对应选项卡后：
    - 过滤并点击当前页面的 protein 链接（避开黑名单域）
    - 在详情窗口中查找 Uniprot 链接，记录其 text 作为 Uniprot ID（不点击链接，且不去重）
    - 避免处理黑名单标题的窗口
    - 翻页直到最后
    """
    desired_key = f"{algo}" + (f" ({metric})" if metric else "")
    print(f"\n开始处理选项卡：{desired_key}")

    clicked_text = click_results_tab(driver, algo, metric)
    if not clicked_text:
        print(f"[WARN] 未能点击到 {desired_key}，跳过该轮抓取。")
        return

    # 点击选项卡后再等一段时间，确保结果稳定
    time.sleep(POST_TAB_CLICK_WAIT)

    current_page = 1
    total_pages = get_total_pages(driver)
    print(f"[{desired_key}] {'总页数 '+str(total_pages) if total_pages else '动态翻页模式'}")

    categorized_ids.setdefault(desired_key, [])

    while True:
        try:
            print(f"[{desired_key}] 正在处理第 {current_page} 页")

            main_window = driver.current_window_handle
            protein_links = driver.find_elements(
                By.XPATH,
                "//a[contains(@class, 'ant-typography') and contains(@class, 'css-zbbqwj')]"
            )
            print(f"[{desired_key}] 本页 protein 链接数量（原始）：{len(protein_links)}")

            # 预过滤不相关链接（避免点击打开无用窗口）
            filtered_links = []
            for link in protein_links:
                if link_is_blacklisted(link):
                    try:
                        print(f"[{desired_key}] 跳过黑名单链接：text='{link.text}' href='{link.get_attribute('href')}'")
                    except Exception:
                        print(f"[{desired_key}] 跳过黑名单链接（无法读取文本或href）")
                    continue
                filtered_links.append(link)

            print(f"[{desired_key}] 本页 protein 链接数量（过滤后）：{len(filtered_links)}")

            # 点击过滤后的 protein 链接以打开详情窗口
            for i, link in enumerate(filtered_links):
                try:
                    driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", link)
                    time.sleep(0.12)
                    link.click()
                    time.sleep(0.7)
                except Exception as e:
                    print(f"[{desired_key}] 点击第 {i+1} 个 protein 链接失败：{e}")

            # 等待所有窗口打开完成后再处理
            time.sleep(AFTER_OPEN_ALL_LINKS_WAIT)

            # 处理已打开的详情窗口
            all_windows = driver.window_handles
            protein_windows = [w for w in all_windows if w != main_window]
            print(f"[{desired_key}] 打开详情窗口数：{len(protein_windows)}")

            for i, window in enumerate(protein_windows):
                try:
                    driver.switch_to.window(window)
                    # 切到详情窗口后，先等待加载充分
                    time.sleep(DETAIL_WINDOW_INITIAL_WAIT)

                    title = driver.title
                    print(f"[{desired_key}] 当前详情窗口标题：{title}")

                    # 黑名单标题窗口直接关闭
                    if title_is_blacklisted(title):
                        print(f"[{desired_key}] 跳过黑名单窗口：{title}")
                        driver.close()
                        continue

                    # 查找 Uniprot 链接（只记录，不点击）
                    uniprot_links = driver.find_elements(
                        By.XPATH,
                        "//a[contains(@href, 'https://www.uniprot.org/uniprotkb/') and @target='_blank']"
                    )
                    # 如果初次为0，额外等待后重试一次
                    if len(uniprot_links) == 0:
                        time.sleep(UNIPROT_EXTRA_WAIT_IF_EMPTY)
                        uniprot_links = driver.find_elements(
                            By.XPATH,
                            "//a[contains(@href, 'https://www.uniprot.org/uniprotkb/') and @target='_blank']"
                        )

                    print(f"[{desired_key}] 在当前详情窗口找到 {len(uniprot_links)} 个 Uniprot 链接")

                    for j, uniprot_link in enumerate(uniprot_links):
                        try:
                            uniprot_id = uniprot_link.text.strip()
                            if uniprot_id:
                                # 不去重：发现就记录
                                categorized_ids[desired_key].append(uniprot_id)
                                print(f"[{desired_key}] 记录 Uniprot ID: {uniprot_id}")
                        except Exception as e:
                            print(f"[{desired_key}] 记录第 {j+1} 个 Uniprot 链接失败：{e}")

                    driver.close()
                except Exception as e:
                    print(f"[{desired_key}] 处理第 {i+1} 个详情窗口失败：{e}")

            driver.switch_to.window(main_window)
            print(f"[{desired_key}] 第 {current_page} 页处理完成")

            # 翻页逻辑
            if total_pages and current_page >= total_pages:
                print(f"[{desired_key}] 已处理完所有 {total_pages} 页")
                break
            try:
                next_btn = driver.find_element(
                    By.XPATH,
                    "//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"
                )
                parent_li = next_btn.find_element(By.XPATH, "..")
                if "ant-pagination-disabled" in parent_li.get_attribute("class"):
                    print(f"[{desired_key}] 下一页按钮已禁用，最后一页")
                    break
                driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", next_btn)
                next_btn.click()
                current_page += 1
                time.sleep(1.4)
            except Exception as e:
                print(f"[{desired_key}] 没有可用的下一页按钮：{e}")
                break

        except Exception as e:
            print(f"[{desired_key}] 处理页面 {current_page} 失败：{e}")
            snapshot(driver, f"process_error_{desired_key.replace(' ', '_')}")
            break

# ------------------ 主流程 ------------------
def main():
    options = Options()
    options.add_argument("--start-maximized")
    driver = create_edge_driver(options)

    categorized_ids = {}  # 分类保存 Uniprot ID（允许重复）

    try:
        driver.set_page_load_timeout(60)
        driver.get("https://reme.biodesign.ac.cn/")
        wait = WebDriverWait(driver, 40)

        # Step 1: SMILES Generator
        wait.until(EC.element_to_be_clickable((By.XPATH, "//button[span[text()='SMILES Generator']]"))).click()

        # Step 2: Open Structure
        wait.until(EC.element_to_be_clickable((By.XPATH, "//button[@data-testid='open-file-button']"))).click()

        # Step 3: 上传 .rxn
        rxn_file_path = find_rxn_file()
        if not rxn_file_path:
            print("当前目录中未找到 .rxn 文件")
            return
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "input[type='file']"))).send_keys(rxn_file_path)

        # Step 4: Add to Canvas
        wait.until(EC.element_to_be_clickable((By.XPATH, "//input[@value='Add to Canvas']"))).click()
        time.sleep(2)

        # Step 5: 点击画布中心（模拟用户操作）
        window_size = driver.get_window_size()
        center_x = window_size['width'] // 2
        center_y = window_size['height'] // 2
        actions = ActionChains(driver)
        steps = 20
        step_x = center_x / steps
        step_y = center_y / steps
        for _ in range(steps):
            actions.move_by_offset(step_x, step_y).perform()
            time.sleep(2 / steps)
        actions.click().perform()
        time.sleep(2)

        # Step 6: OK
        wait.until(EC.element_to_be_clickable((By.XPATH, "//button[span[text()='OK']]"))).click()

        # Step 7: Search
        wait.until(EC.element_to_be_clickable((By.XPATH, "//button[span[text()='Search']]"))).click()

        # Step 8: Yes
        wait.until(EC.element_to_be_clickable((By.XPATH, "//button[span[text()='Yes']]"))).click()

        # 逐次运行五个方法（复选框 + 选项卡；抓取）
        runs = [
            {"card": "RDKit", "labels": ["Select"], "algo": "RDKit", "metric": None},
            {"card": "RXNFP", "labels": ["Cosine"], "algo": "RXNFP", "metric": "cosine"},
            {"card": "RXNFP", "labels": ["Euclidean distance"], "algo": "RXNFP", "metric": "euclidean"},
            {"card": "DRFP",  "labels": ["Cosine"], "algo": "DRFP",  "metric": "cosine"},
            {"card": "DRFP",  "labels": ["Euclidean distance"], "algo": "DRFP",  "metric": "euclidean"},
        ]

        for run in runs:
            print("\n==============================")
            print(f"开始运行：{run['card']} - {', '.join(run['labels'])}")

            # Step 9: 勾选当次方法（鲁棒校验）
            select_method_once(driver, wait, run['card'], run['labels'])

            # Step 10: 点击 Excute（等待按钮可用）
            if not wait_execute_and_click(driver, timeout=25):
                print("[WARN] 本轮未能点击 Excute，尝试继续（可能页面未响应）")

            # 等待至少出现一个结果选项卡，并额外等待页面加载充分
            WebDriverWait(driver, 35).until(EC.presence_of_element_located((By.XPATH, "//div[@role='tab']")))
            time.sleep(RESULTS_LOAD_DELAY + BEFORE_TAB_CLICK_EXTRA_DELAY)

            # Step 11/12: 点击对应选项卡并抓取 Uniprot ID（不再点击链接；允许重复）
            process_results_for_method(driver, run["algo"], run["metric"], categorized_ids)

            # 返回到选择页
            go_back_to_selection(driver)
            time.sleep(0.8)

        # 输出分类结果
        script_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(script_dir, "uniprot_ids_by_method.txt")
        total_count = sum(len(v) for v in categorized_ids.values())

        with open(out_path, "w", encoding="utf-8") as f:
            for k, ids in categorized_ids.items():
                f.write(f"=== {k} ===\n")
                if ids:
                    for uid in ids:
                        f.write(uid + "\n")
                else:
                    f.write("(no results)\n")
                f.write("\n")
        print(f"\n成功生成分类txt文件：{out_path}")
        print(f"总共保存了 {total_count} 个 Uniprot ID（按算法-度量分类，包含重复）")

        print("操作完成！按五个方法分别抓取并分类输出。")

    except Exception as e:
        print("发生错误:", e)
        print("详细堆栈：")
        traceback.print_exc()
        snapshot(driver, "fatal_error")
    finally:
        print("脚本运行完毕，浏览器保持打开状态。")
        input("按下回车以退出脚本，但浏览器仍将保持打开状态。")

if __name__ == "__main__":
    main()

'''PS C:\Users\Jiaya> C:/Users/Jiaya/anaconda3/Scripts/activate
PS C:\Users\Jiaya> conda activate reme-try
PS C:\Users\Jiaya> & D:/conda_envs/reme-try/python.exe "d:/Project/pythonProject/pythonProject/reactionsearch - 副本/1inputreaction_Version4.py"
使用本地 EdgeDriver: d:\Project\pythonProject\pythonProject\reactionsearch - 副本\msedgedriver.exe

DevTools listening on ws://127.0.0.1:57419/devtools/browser/579949d5-0d43-4500-b413-11fbf74faeca
[7120:38248:0214/003535.711:ERROR:chrome\browser\importer\edge_china_browsers\edge_qqbrowser_importer_utils_win.cc:163] QQBrowser user data path not found.
[7120:25108:0214/003536.680:ERROR:chrome\browser\task_manager\providers\fallback_task_provider.cc:126] Every renderer should have at least one task provided by a primary task provider. If a "Renderer" fallback task is shown, it is a bug. If you have repro steps, please file a new bug and tag it as a dependency of crbug.com/739782.

==============================
开始运行：RDKit - Select
已清空所有方法复选框。
勾选成功：RDKit - Select
已点击 Excute 按钮。

开始处理选项卡：RDKit
[DEBUG] 检测到结果选项卡：['Introduction', 'RDKit']
[DEBUG] 精确匹配到选项卡：RDKit
检测到总页数：2
[RDKit] 总页数 2
[RDKit] 正在处理第 1 页
[RDKit] 本页 protein 链接数量（原始）：15
[RDKit] 跳过黑名单链接：text='BRENDA' href='https://www.brenda-enzymes.org/'
[RDKit] 跳过黑名单链接：text='CC BY 4.0' href='https://creativecommons.org/licenses/by/4.0/'
[RDKit] 跳过黑名单链接：text='Reporting Center for Illegal and Harmful Information' href='https://www.12377.cn/'
[RDKit] 跳过黑名单链接：text='Cyberspace Administration of China' href='http://www.cac.gov.cn/'
[RDKit] 跳过黑名单链接：text='Internet Crime Info Reporting Center' href='https://cyberpolice.mps.gov.cn/wfjb'
[RDKit] 本页 protein 链接数量（过滤后）：10
[RDKit] 打开详情窗口数：10
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 0 个 Uniprot 链接
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 3 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: F9UM18
[RDKit] 记录 Uniprot ID: A0A1B1WA13
[RDKit] 记录 Uniprot ID: A0A1B1WA14
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 6 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: Q93YP7
[RDKit] 记录 Uniprot ID: Q8W405
[RDKit] 记录 Uniprot ID: Q8W404
[RDKit] 记录 Uniprot ID: Q0QFI7
[RDKit] 记录 Uniprot ID: A0A168S819
[RDKit] 记录 Uniprot ID: A0A168S850
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 1 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: Q6Q874
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 1 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: P76469
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 2 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: A0A1B1WA13
[RDKit] 记录 Uniprot ID: A0A1B1WA14
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 3 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: B9A1Q4
[RDKit] 记录 Uniprot ID: I1NJA2
[RDKit] 记录 Uniprot ID: A0A1L7NRS0
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 2 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: Q53B70
[RDKit] 记录 Uniprot ID: Q93XE6
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 1 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: Q8GHB2
[RDKit] 当前详情窗口标题：REME
[RDKit] 在当前详情窗口找到 2 个 Uniprot 链接
[RDKit] 记录 Uniprot ID: Q53B70
[RDKit] 记录 Uniprot ID: Q93XE6
[RDKit] 第 1 页处理完成
[RDKit] 没有可用的下一页按钮：Message: no such element: Unable to locate element: {"method":"xpath","selector":"//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"}
  (Session info: MicrosoftEdge=144.0.3719.115); For documentation on this error, please visit: https://www.selenium.dev/documentation/webdriver/troubleshooting/errors#nosuchelementexception
Stacktrace:
Symbols not available. Dumping unresolved backtrace:
        0x7ff7d5285045
        0x7ff7d52850a4
        0x7ff7d55d28c4
        0x7ff7d5026a9a
        0x7ff7d5026d35
        0x7ff7d5064b07
        0x7ff7d501d6b7
        0x7ff7d50628f7
        0x7ff7d501cbcc
        0x7ff7d501be26
        0x7ff7d501c9f3
        0x7ff7d50e0214
        0x7ff7d50dc685
        0x7ff7d50ece73
        0x7ff7d52a0378
        0x7ff7d52a8c46
        0x7ff7d528c724
        0x7ff7d528c869
        0x7ff7d527a702
        0x7ffce0b3e8d7
        0x7ffce0f6c40c

未找到返回按钮，使用浏览器后退。
已返回到选择页面。

==============================
开始运行：RXNFP - Cosine
已清空所有方法复选框。
勾选成功：RXNFP - Cosine
已点击 Excute 按钮。

开始处理选项卡：RXNFP (cosine)
[DEBUG] 检测到结果选项卡：['Introduction', 'RXNFP (cosine)']
[DEBUG] 精确匹配到选项卡：RXNFP (cosine)
检测到总页数：2
[RXNFP (cosine)] 总页数 2
[RXNFP (cosine)] 正在处理第 1 页
[RXNFP (cosine)] 本页 protein 链接数量（原始）：15
[RXNFP (cosine)] 跳过黑名单链接：text='BRENDA' href='https://www.brenda-enzymes.org/'
[RXNFP (cosine)] 跳过黑名单链接：text='CC BY 4.0' href='https://creativecommons.org/licenses/by/4.0/'
[RXNFP (cosine)] 跳过黑名单链接：text='Reporting Center for Illegal and Harmful Information' href='https://www.12377.cn/'
[RXNFP (cosine)] 跳过黑名单链接：text='Cyberspace Administration of China' href='http://www.cac.gov.cn/'
[RXNFP (cosine)] 跳过黑名单链接：text='Internet Crime Info Reporting Center' href='https://cyberpolice.mps.gov.cn/wfjb'
[RXNFP (cosine)] 本页 protein 链接数量（过滤后）：10
[RXNFP (cosine)] 打开详情窗口数：10
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q9AJC6
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 3 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: C5H429
[RXNFP (cosine)] 记录 Uniprot ID: P40952
[RXNFP (cosine)] 记录 Uniprot ID: Q9R9V9
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: A0A0H4SN47
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 2 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q93TU6
[RXNFP (cosine)] 记录 Uniprot ID: Q8YNV6
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q9AJC6
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q9BXD5
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 2 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: P04964
[RXNFP (cosine)] 记录 Uniprot ID: F7YTT4
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q65CJ7
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 10 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: P51102
[RXNFP (cosine)] 记录 Uniprot ID: P51103
[RXNFP (cosine)] 记录 Uniprot ID: P51104
[RXNFP (cosine)] 记录 Uniprot ID: P51110
[RXNFP (cosine)] 记录 Uniprot ID: P51105
[RXNFP (cosine)] 记录 Uniprot ID: P51107
[RXNFP (cosine)] 记录 Uniprot ID: Q6UAQ7
[RXNFP (cosine)] 记录 Uniprot ID: A0A0A0PXZ7
[RXNFP (cosine)] 记录 Uniprot ID: A0A0C6ENR6
[RXNFP (cosine)] 记录 Uniprot ID: M4DU26
[RXNFP (cosine)] 当前详情窗口标题：REME
[RXNFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (cosine)] 记录 Uniprot ID: Q65CJ7
[RXNFP (cosine)] 第 1 页处理完成
[RXNFP (cosine)] 没有可用的下一页按钮：Message: no such element: Unable to locate element: {"method":"xpath","selector":"//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"}
  (Session info: MicrosoftEdge=144.0.3719.115); For documentation on this error, please visit: https://www.selenium.dev/documentation/webdriver/troubleshooting/errors#nosuchelementexception
Stacktrace:
Symbols not available. Dumping unresolved backtrace:
        0x7ff7d5285045
        0x7ff7d52850a4
        0x7ff7d55d28c4
        0x7ff7d5026a9a
        0x7ff7d5026d35
        0x7ff7d5064b07
        0x7ff7d501d6b7
        0x7ff7d50628f7
        0x7ff7d501cbcc
        0x7ff7d501be26
        0x7ff7d501c9f3
        0x7ff7d50e0214
        0x7ff7d50dc685
        0x7ff7d50ece73
        0x7ff7d52a0378
        0x7ff7d52a8c46
        0x7ff7d528c724
        0x7ff7d528c869
        0x7ff7d527a702
        0x7ffce0b3e8d7
        0x7ffce0f6c40c

未找到返回按钮，使用浏览器后退。
已返回到选择页面。

==============================
开始运行：RXNFP - Euclidean distance
已清空所有方法复选框。
勾选成功：RXNFP - Euclidean distance
已点击 Excute 按钮。

开始处理选项卡：RXNFP (euclidean)
[DEBUG] 检测到结果选项卡：['Introduction', 'RXNFP (euclidean)']
[DEBUG] 精确匹配到选项卡：RXNFP (euclidean)
检测到总页数：2
[RXNFP (euclidean)] 总页数 2
[RXNFP (euclidean)] 正在处理第 1 页
[RXNFP (euclidean)] 本页 protein 链接数量（原始）：15
[RXNFP (euclidean)] 跳过黑名单链接：text='BRENDA' href='https://www.brenda-enzymes.org/'
[RXNFP (euclidean)] 跳过黑名单链接：text='CC BY 4.0' href='https://creativecommons.org/licenses/by/4.0/'
[RXNFP (euclidean)] 跳过黑名单链接：text='Reporting Center for Illegal and Harmful Information' href='https://www.12377.cn/'
[RXNFP (euclidean)] 跳过黑名单链接：text='Cyberspace Administration of China' href='http://www.cac.gov.cn/'
[RXNFP (euclidean)] 跳过黑名单链接：text='Internet Crime Info Reporting Center' href='https://cyberpolice.mps.gov.cn/wfjb'
[RXNFP (euclidean)] 本页 protein 链接数量（过滤后）：10
[RXNFP (euclidean)] 打开详情窗口数：10
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q65CJ7
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 2 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q93TU6
[RXNFP (euclidean)] 记录 Uniprot ID: Q8YNV6
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q9AJC6
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 3 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: C5H429
[RXNFP (euclidean)] 记录 Uniprot ID: P40952
[RXNFP (euclidean)] 记录 Uniprot ID: Q9R9V9
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q9BXD5
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q9AJC6
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 10 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: P51102
[RXNFP (euclidean)] 记录 Uniprot ID: P51103
[RXNFP (euclidean)] 记录 Uniprot ID: P51104
[RXNFP (euclidean)] 记录 Uniprot ID: P51110
[RXNFP (euclidean)] 记录 Uniprot ID: P51105
[RXNFP (euclidean)] 记录 Uniprot ID: P51107
[RXNFP (euclidean)] 记录 Uniprot ID: Q6UAQ7
[RXNFP (euclidean)] 记录 Uniprot ID: A0A0A0PXZ7
[RXNFP (euclidean)] 记录 Uniprot ID: A0A0C6ENR6
[RXNFP (euclidean)] 记录 Uniprot ID: M4DU26
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: P44430
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 2 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: P04964
[RXNFP (euclidean)] 记录 Uniprot ID: F7YTT4
[RXNFP (euclidean)] 当前详情窗口标题：REME
[RXNFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[RXNFP (euclidean)] 记录 Uniprot ID: Q65CJ7
[RXNFP (euclidean)] 第 1 页处理完成
[RXNFP (euclidean)] 没有可用的下一页按钮：Message: no such element: Unable to locate element: {"method":"xpath","selector":"//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"}
  (Session info: MicrosoftEdge=144.0.3719.115); For documentation on this error, please visit: https://www.selenium.dev/documentation/webdriver/troubleshooting/errors#nosuchelementexception
Stacktrace:
Symbols not available. Dumping unresolved backtrace:
        0x7ff7d5285045
        0x7ff7d52850a4
        0x7ff7d55d28c4
        0x7ff7d5026a9a
        0x7ff7d5026d35
        0x7ff7d5064b07
        0x7ff7d501d6b7
        0x7ff7d50628f7
        0x7ff7d501cbcc
        0x7ff7d501be26
        0x7ff7d501c9f3
        0x7ff7d50e0214
        0x7ff7d50dc685
        0x7ff7d50ece73
        0x7ff7d52a0378
        0x7ff7d52a8c46
        0x7ff7d528c724
        0x7ff7d528c869
        0x7ff7d527a702
        0x7ffce0b3e8d7
        0x7ffce0f6c40c

未找到返回按钮，使用浏览器后退。
已返回到选择页面。

==============================
开始运行：DRFP - Cosine
已清空所有方法复选框。
勾选成功：DRFP - Cosine
已点击 Excute 按钮。

开始处理选项卡：DRFP (cosine)
[DEBUG] 检测到结果选项卡：['Introduction', 'RXNFP (cosine)']
[DEBUG] 评分匹配点击选项卡：RXNFP (cosine)（algo=drfp, metric=cosine）
检测到总页数：2
[DRFP (cosine)] 总页数 2
[DRFP (cosine)] 正在处理第 1 页
[DRFP (cosine)] 本页 protein 链接数量（原始）：15
[DRFP (cosine)] 跳过黑名单链接：text='BRENDA' href='https://www.brenda-enzymes.org/'
[DRFP (cosine)] 跳过黑名单链接：text='CC BY 4.0' href='https://creativecommons.org/licenses/by/4.0/'
[DRFP (cosine)] 跳过黑名单链接：text='Reporting Center for Illegal and Harmful Information' href='https://www.12377.cn/'
[DRFP (cosine)] 跳过黑名单链接：text='Cyberspace Administration of China' href='http://www.cac.gov.cn/'
[DRFP (cosine)] 跳过黑名单链接：text='Internet Crime Info Reporting Center' href='https://cyberpolice.mps.gov.cn/wfjb'
[DRFP (cosine)] 本页 protein 链接数量（过滤后）：10
[DRFP (cosine)] 打开详情窗口数：10
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 10 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: P51102
[DRFP (cosine)] 记录 Uniprot ID: P51103
[DRFP (cosine)] 记录 Uniprot ID: P51104
[DRFP (cosine)] 记录 Uniprot ID: P51110
[DRFP (cosine)] 记录 Uniprot ID: P51105
[DRFP (cosine)] 记录 Uniprot ID: P51107
[DRFP (cosine)] 记录 Uniprot ID: Q6UAQ7
[DRFP (cosine)] 记录 Uniprot ID: A0A0A0PXZ7
[DRFP (cosine)] 记录 Uniprot ID: A0A0C6ENR6
[DRFP (cosine)] 记录 Uniprot ID: M4DU26
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 3 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: C5H429
[DRFP (cosine)] 记录 Uniprot ID: P40952
[DRFP (cosine)] 记录 Uniprot ID: Q9R9V9
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q9AJC6
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q65CJ7
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 2 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: P04964
[DRFP (cosine)] 记录 Uniprot ID: F7YTT4
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q65CJ7
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: A0A0H4SN47
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q9BXD5
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q9AJC6
[DRFP (cosine)] 当前详情窗口标题：REME
[DRFP (cosine)] 在当前详情窗口找到 2 个 Uniprot 链接
[DRFP (cosine)] 记录 Uniprot ID: Q93TU6
[DRFP (cosine)] 记录 Uniprot ID: Q8YNV6
[DRFP (cosine)] 第 1 页处理完成
[DRFP (cosine)] 没有可用的下一页按钮：Message: no such element: Unable to locate element: {"method":"xpath","selector":"//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"}
  (Session info: MicrosoftEdge=144.0.3719.115); For documentation on this error, please visit: https://www.selenium.dev/documentation/webdriver/troubleshooting/errors#nosuchelementexception
Stacktrace:
Symbols not available. Dumping unresolved backtrace:
        0x7ff7d5285045
        0x7ff7d52850a4
        0x7ff7d55d28c4
        0x7ff7d5026a9a
        0x7ff7d5026d35
        0x7ff7d5064b07
        0x7ff7d501d6b7
        0x7ff7d50628f7
        0x7ff7d501cbcc
        0x7ff7d501be26
        0x7ff7d501c9f3
        0x7ff7d50e0214
        0x7ff7d50dc685
        0x7ff7d50ece73
        0x7ff7d52a0378
        0x7ff7d52a8c46
        0x7ff7d528c724
        0x7ff7d528c869
        0x7ff7d527a702
        0x7ffce0b3e8d7
        0x7ffce0f6c40c

未找到返回按钮，使用浏览器后退。
已返回到选择页面。

==============================
开始运行：DRFP - Euclidean distance
已清空所有方法复选框。
勾选成功：DRFP - Euclidean distance
已点击 Excute 按钮。

开始处理选项卡：DRFP (euclidean)
[DEBUG] 检测到结果选项卡：['Introduction', 'RXNFP (euclidean)']
[DEBUG] 评分匹配点击选项卡：RXNFP (euclidean)（algo=drfp, metric=euclidean）
检测到总页数：2
[DRFP (euclidean)] 总页数 2
[DRFP (euclidean)] 正在处理第 1 页
[DRFP (euclidean)] 本页 protein 链接数量（原始）：15
[DRFP (euclidean)] 跳过黑名单链接：text='BRENDA' href='https://www.brenda-enzymes.org/'
[DRFP (euclidean)] 跳过黑名单链接：text='CC BY 4.0' href='https://creativecommons.org/licenses/by/4.0/'
[DRFP (euclidean)] 跳过黑名单链接：text='Reporting Center for Illegal and Harmful Information' href='https://www.12377.cn/'
[DRFP (euclidean)] 跳过黑名单链接：text='Cyberspace Administration of China' href='http://www.cac.gov.cn/'
[DRFP (euclidean)] 跳过黑名单链接：text='Internet Crime Info Reporting Center' href='https://cyberpolice.mps.gov.cn/wfjb'
[DRFP (euclidean)] 本页 protein 链接数量（过滤后）：10
[DRFP (euclidean)] 打开详情窗口数：10
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q9AJC6
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q9AJC6
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q9BXD5
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q65CJ7
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 2 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: P04964
[DRFP (euclidean)] 记录 Uniprot ID: F7YTT4
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 10 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: P51102
[DRFP (euclidean)] 记录 Uniprot ID: P51103
[DRFP (euclidean)] 记录 Uniprot ID: P51104
[DRFP (euclidean)] 记录 Uniprot ID: P51110
[DRFP (euclidean)] 记录 Uniprot ID: P51105
[DRFP (euclidean)] 记录 Uniprot ID: P51107
[DRFP (euclidean)] 记录 Uniprot ID: Q6UAQ7
[DRFP (euclidean)] 记录 Uniprot ID: A0A0A0PXZ7
[DRFP (euclidean)] 记录 Uniprot ID: A0A0C6ENR6
[DRFP (euclidean)] 记录 Uniprot ID: M4DU26
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 2 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q93TU6
[DRFP (euclidean)] 记录 Uniprot ID: Q8YNV6
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 3 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: C5H429
[DRFP (euclidean)] 记录 Uniprot ID: P40952
[DRFP (euclidean)] 记录 Uniprot ID: Q9R9V9
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: P44430
[DRFP (euclidean)] 当前详情窗口标题：REME
[DRFP (euclidean)] 在当前详情窗口找到 1 个 Uniprot 链接
[DRFP (euclidean)] 记录 Uniprot ID: Q65CJ7
[DRFP (euclidean)] 第 1 页处理完成
[DRFP (euclidean)] 没有可用的下一页按钮：Message: no such element: Unable to locate element: {"method":"xpath","selector":"//a[contains(@class, 'ant-pagination-item-link') and (contains(@aria-label,'Next Page') or contains(@aria-label,'Next'))]"}
  (Session info: MicrosoftEdge=144.0.3719.115); For documentation on this error, please visit: https://www.selenium.dev/documentation/webdriver/troubleshooting/errors#nosuchelementexception
Stacktrace:
Symbols not available. Dumping unresolved backtrace:
        0x7ff7d5285045
        0x7ff7d52850a4
        0x7ff7d55d28c4
        0x7ff7d5026a9a
        0x7ff7d5026d35
        0x7ff7d5064b07
        0x7ff7d501d6b7
        0x7ff7d50628f7
        0x7ff7d501cbcc
        0x7ff7d501be26
        0x7ff7d501c9f3
        0x7ff7d50e0214
        0x7ff7d50dc685
        0x7ff7d50ece73
        0x7ff7d52a0378
        0x7ff7d52a8c46
        0x7ff7d528c724
        0x7ff7d528c869
        0x7ff7d527a702
        0x7ffce0b3e8d7
        0x7ffce0f6c40c

未找到返回按钮，使用浏览器后退。
已返回到选择页面。

成功生成分类txt文件：d:\Project\pythonProject\pythonProject\reactionsearch - 副本\uniprot_ids_by_method.txt
总共保存了 113 个 Uniprot ID（按算法-度量分类，包含重复）
操作完成！按五个方法分别抓取并分类输出。
脚本运行完毕，浏览器保持打开状态。
按下回车以退出脚本，但浏览器仍将保持打开状态。'''