#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
图 3 - 社区网络可视化（统一 UI/UX 规范）
- 背景：#F8F9FA（极浅灰）
- 字体：Arial（无衬线）
- 社区配色：Tableau 10；若社区数 > 10，则剩余部分使用 Viridis 离散采样
- 关键蛋白（Top hubs）：Orange/Gold #F39C12
- Tooltip（悬停提示）必须包含：蛋白 Symbol、全称、功能描述、连接度数（Degree）
- 缩放：支持鼠标滚轮缩放
- 动态网络：使用 Pyvis，并统一 physics（物理引擎）参数

输入文件（STRING v12.0，小鼠 10090）：
- 10090.protein.links.v12.0.txt.gz
- 10090.protein.info.v12.0.txt.gz
（可选）- 10090.protein.aliases.v12.0.txt.gz（本图不必需）
"""

import os
import pandas as pd
import numpy as np
import networkx as nx

# -----------------------------
# 0) 路径与统一 UI 参数
# -----------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

LINKS_GZ = os.path.join(DATA_DIR, "10090.protein.links.v12.0.txt.gz")
INFO_GZ  = os.path.join(DATA_DIR, "10090.protein.info.v12.0.txt.gz")

SCORE_CUTOFF = 900          # 分数阈值（可改为 400/700/900）
MAX_NODES_TO_PLOT = 1000    # 绘图节点上限（建议 1500~3000）
TOP_HUBS = 25               # 关键蛋白数量（按度最高 Top N）-> 橙色强调
TOP_LABELS = 25             # 显示标签（label）的节点数（只给少数点打字，避免糊）
RANDOM_SEED = 42

OUT_HTML = os.path.join(os.path.dirname(__file__), "..", "figures", f"fig3_community_network_th{SCORE_CUTOFF}.html")
OUT_CSV  = os.path.join(os.path.dirname(__file__), "..", "outputs", f"community_assignments_th{SCORE_CUTOFF}.csv")
OUT_GEXF = os.path.join(os.path.dirname(__file__), "..", "outputs", f"fig3_network_th{SCORE_CUTOFF}.gexf")

# 创建输出目录
os.makedirs(os.path.dirname(OUT_HTML), exist_ok=True)
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
os.makedirs(os.path.dirname(OUT_GEXF), exist_ok=True)

# ---- 颜色规范（严格）----
BG_COLOR = "#F8F9FA"        # 背景色
TEXT_COLOR = "#2C3E50"      # 科学蓝/深色文字
EDGE_COLOR = "rgba(44,62,80,0.22)"  # 边的淡色（低对比，避免抢视觉）
HUB_COLOR = "#F39C12"       # 关键蛋白强调色（橙/金）

# ---- 统一 physics 参数（全组复制同一份，不要个人随意改）----
PHYSICS_OPTIONS = {
    "enabled": True,
    "solver": "forceAtlas2Based",
    "forceAtlas2Based": {
        "gravitationalConstant": -60,
        "centralGravity": 0.012,
        "springLength": 120,
        "springConstant": 0.08,
        "damping": 0.45,
        "avoidOverlap": 0.7
    },
    "maxVelocity": 40,
    "minVelocity": 0.1,
    "timestep": 0.5,
    "stabilization": {
        "enabled": True,
        "iterations": 1200,
        "updateInterval": 25
    }
}

# -----------------------------
# 1) 读取 protein.info：构建 id -> (symbol / full name / description) 映射
# -----------------------------
def load_info(info_path: str):
    info = pd.read_csv(info_path, sep="\t", compression="gzip")

    # --- 列名鲁棒处理 ---
    # STRING 常见列：protein_external_id, preferred_name, annotation（描述）
    id_col = "protein_external_id" if "protein_external_id" in info.columns else info.columns[0]
    symbol_col = "preferred_name" if "preferred_name" in info.columns else (info.columns[1] if len(info.columns) > 1 else id_col)

    # annotation/description 列兜底识别
    desc_candidates = [c for c in info.columns if ("annot" in c.lower() or "desc" in c.lower())]
    desc_col = desc_candidates[0] if desc_candidates else None

    # 构建映射：protein_id -> symbol
    id2symbol = dict(zip(info[id_col].astype(str), info[symbol_col].astype(str)))

    # 构建映射：protein_id -> 描述（annotation）
    if desc_col is not None:
        id2desc = dict(zip(info[id_col].astype(str), info[desc_col].astype(str)))
    else:
        id2desc = {}

    # "全称"在 STRING info 中通常没有独立列，annotation 往往更像"名称+功能"的综合描述
    # 这里做一个实用约定：
    # - full_name：优先用 annotation（若有）；否则用 symbol
    # - function_desc：annotation（若有），否则留空
    id2full = {}
    for pid in id2symbol:
        ann = id2desc.get(pid, "")
        if isinstance(ann, str) and ann.strip():
            id2full[pid] = ann.strip()
        else:
            id2full[pid] = id2symbol.get(pid, pid)

    return id2symbol, id2full, id2desc


# -----------------------------
# 2) 读取 links 并建图（按阈值过滤边）
# -----------------------------
def build_graph(links_path: str, score_cutoff: int) -> nx.Graph:
    df = pd.read_csv(
        links_path,
        sep=r"\s+",
        compression="gzip",
        engine="python"
    )
    # 标准列名一般为：protein1 protein2 combined_score
    df.columns = ["protein1", "protein2", "combined_score"]
    df = df[df["combined_score"] >= score_cutoff].copy()

    G = nx.Graph()
    # 将 combined_score（0-1000）缩放为 weight（0-1）
    for p1, p2, s in df.itertuples(index=False):
        if p1 != p2:
            G.add_edge(str(p1), str(p2), score=int(s), weight=float(s)/1000.0)
    return G


# -----------------------------
# 3) 选择用于绘图的子图（避免太大/太乱）
# -----------------------------
def choose_plot_subgraph(G: nx.Graph, max_nodes: int) -> nx.Graph:
    if G.number_of_nodes() == 0:
        return G

    # 取最大连通子图（Largest Connected Component）
    lcc = max(nx.connected_components(G), key=len)
    H = G.subgraph(lcc).copy()

    if H.number_of_nodes() <= max_nodes:
        return H

    # 若仍过大：按度排序，取度最高的 top_nodes 构建诱导子图（Induced Subgraph）
    deg = dict(H.degree())
    top_nodes = sorted(deg, key=deg.get, reverse=True)[:max_nodes]
    H2 = H.subgraph(top_nodes).copy()
    return H2


# -----------------------------
# 4) 社区检测（Louvain）
# -----------------------------
def louvain_partition(H: nx.Graph, seed: int = 42) -> dict:
    import community as community_louvain  # python-louvain
    part = community_louvain.best_partition(H, weight="weight", random_state=seed)
    return part


# -----------------------------
# 5) 社区配色（Tableau10 + Viridis 超出部分）- 严格统一
# -----------------------------
def build_comm_colors(comm_ids):
    # 使用 plotly 内置调色板（跨平台一致）
    from plotly.colors import qualitative, sample_colorscale

    tableau10 = qualitative.T10  # 10 种颜色
    comm_ids_sorted = sorted(comm_ids)

    mapping = {}
    if len(comm_ids_sorted) <= 10:
        # 社区数 <= 10：直接用 Tableau10
        for i, cid in enumerate(comm_ids_sorted):
            mapping[cid] = tableau10[i]
    else:
        # 社区数 > 10：前 10 个用 Tableau10，剩余用 Viridis 离散采样补齐
        for i, cid in enumerate(comm_ids_sorted[:10]):
            mapping[cid] = tableau10[i]
        remain = comm_ids_sorted[10:]
        viridis_colors = sample_colorscale("Viridis", [i/(max(1, len(remain)-1)) for i in range(len(remain))])
        for cid, col in zip(remain, viridis_colors):
            mapping[cid] = col

    return mapping


# -----------------------------
# 6) 用 Pyvis 生成互动网络图（统一 UI/UX + Tooltip 字段齐全）
# -----------------------------
def export_pyvis(H: nx.Graph, part: dict, id2symbol: dict, id2full: dict, id2desc: dict,
                 out_html: str, out_csv: str, out_gexf: str):

    from pyvis.network import Network
    import json

    # 计算每个节点的度（Degree）：用于 tooltip 与 hub 识别
    deg = dict(H.degree())
    top_hubs = set(sorted(deg, key=deg.get, reverse=True)[:TOP_HUBS])
    top_labels = set(sorted(deg, key=deg.get, reverse=True)[:TOP_LABELS])

    # 获取社区编号集合，并为每个社区分配颜色
    comm_ids = set(part.values())
    comm_color = build_comm_colors(comm_ids)

    # 导出社区分配结果 CSV（用于后续分析/表格/复现）
    assign = pd.DataFrame({
        "protein_id": list(H.nodes()),
        "symbol": [id2symbol.get(n, n) for n in H.nodes()],
        "community": [part.get(n, -1) for n in H.nodes()],
        "degree": [deg.get(n, 0) for n in H.nodes()],
        "annotation": [id2desc.get(n, "") for n in H.nodes()],
    })
    assign.to_csv(out_csv, index=False, encoding="utf-8-sig")

    # 导出 GEXF（可导入 Cytoscape/Gephi 做论文级静态排版）
    for n in H.nodes():
        H.nodes[n]["symbol"] = id2symbol.get(n, n)
        H.nodes[n]["community"] = int(part.get(n, -1))
        H.nodes[n]["degree"] = int(deg.get(n, 0))
        H.nodes[n]["annotation"] = id2desc.get(n, "")
    nx.write_gexf(H, out_gexf)

    # 创建 Pyvis 网络对象（浅色背景、统一字体）
    net = Network(
        height="850px",
        width="100%",
        bgcolor=BG_COLOR,
        font_color=TEXT_COLOR,
        directed=False,
        notebook=False
    )

    # 设置交互与视觉选项：缩放、悬停提示、节点/边样式、物理引擎参数
    options = {
        "interaction": {
            "hover": True,
            "tooltipDelay": 120,
            "zoomView": True,
            "dragView": True,
            "navigationButtons": False
        },
        "nodes": {
            "shape": "dot",
            "borderWidth": 1,
            "font": {
                "face": "Arial",
                "size": 14,
                "color": TEXT_COLOR,
                "bold": False
            }
        },
        "edges": {
            "color": EDGE_COLOR,
            "smooth": {"enabled": True, "type": "dynamic"},
            "width": 0.7
        },
        "physics": PHYSICS_OPTIONS
    }
    net.set_options(json.dumps(options))

    # 添加节点：社区上色 + hub 强调 + tooltip 字段齐全
    for n in H.nodes():
        symbol = id2symbol.get(n, n)
        full_name = id2full.get(n, symbol)
        desc = id2desc.get(n, "")
        degree = deg.get(n, 0)
        comm = part.get(n, -1)

        # hub（关键蛋白）强调：橙色 + 更大点
        if n in top_hubs:
            color = HUB_COLOR
            size = 18
        else:
            color = comm_color.get(comm, "#3498DB")  # 若意外缺失则回退到蓝色
            size = 10

        # 只给少量节点显示 label，避免整体糊成一片
        label = symbol if n in top_labels else ""

        # Tooltip：必须包含 Symbol、全称、功能描述、Degree（再附带社区编号）
        title = (
            f"<div style='font-family:Arial;color:{TEXT_COLOR};line-height:1.35;'>"
            f"<b>Symbol</b>: {symbol}<br>"
            f"<b>Full name</b>: {full_name}<br>"
            f"<b>Function</b>: {desc if desc else full_name}<br>"
            f"<b>Degree</b>: {degree}<br>"
            f"<b>Community</b>: {comm}"
            f"</div>"
        )

        net.add_node(
            n,
            label=label,
            title=title,
            color=color,
            size=size
        )

    # 添加边：悬停显示 combined_score
    for u, v, data in H.edges(data=True):
        score = data.get("score", None)
        etitle = (
            f"<div style='font-family:Arial;color:{TEXT_COLOR};'>"
            f"<b>combined_score</b>: {score}"
            f"</div>"
        ) if score is not None else ""
        net.add_edge(u, v, title=etitle, value=data.get("weight", 0.5))

    # 设置页面标题区（标题加粗，比正文大 2pt；并说明高亮规则）
    net.heading = (
        f"<h2 style='font-family:Arial;color:{TEXT_COLOR};"
        f"font-weight:700;margin:20px 0 10px 0;'>"
        f"Figure 3. Community Network (STRING 10090) | cutoff={SCORE_CUTOFF}"
        f"</h2>"
        f"<div style='font-family:Arial;color:{TEXT_COLOR};margin-bottom:12px;'>"
        f"<span style='font-weight:700;'>Highlight</span>: Top {TOP_HUBS} hubs in "
        f"<span style='color:{HUB_COLOR};font-weight:700;'>Orange/Gold</span>."
        f"</div>"
    )

    # 输出 HTML
    net.save_graph(out_html)

    # 后处理：给页面增加约 10% 的 padding，并统一背景色
    try:
        with open(out_html, "r", encoding="utf-8") as f:
            html = f.read()

        # 在 <body> 后注入容器 padding（只替换一次）
        if "<body>" in html:
            html = html.replace(
                "<body>",
                f"<body style='margin:0;background:{BG_COLOR};'>"
                f"<div style='padding:5% 5% 5% 5%;'>",
                1
            )
            html = html.replace("</body>", "</div></body>", 1)

        with open(out_html, "w", encoding="utf-8") as f:
            f.write(html)
    except Exception:
        pass


def main():
    # 检查数据文件是否存在
    if not os.path.exists(LINKS_GZ):
        print(f"错误: 找不到文件 {LINKS_GZ}")
        print("请确保已下载 STRING 数据文件到正确的数据目录中")
        return

    if not os.path.exists(INFO_GZ):
        print(f"错误: 找不到文件 {INFO_GZ}")
        print("请确保已下载 STRING 数据文件到正确的数据目录中")
        return

    # 1) 读取映射表（id -> symbol/full/desc）
    id2symbol, id2full, id2desc = load_info(INFO_GZ)

    # 2) 建图（按 combined_score 阈值过滤）
    print(f"[INFO] Reading links and building graph (cutoff={SCORE_CUTOFF}) ...")
    G = build_graph(LINKS_GZ, SCORE_CUTOFF)
    print(f"[INFO] Raw graph: nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")

    # 3) 选择用于绘图的子图（最大连通子图；若过大则取 top-degree 诱导子图）
    H = choose_plot_subgraph(G, MAX_NODES_TO_PLOT)
    print(f"[INFO] Plot graph: nodes={H.number_of_nodes()}, edges={H.number_of_edges()}")

    # 4) 社区检测（Louvain）
    print("[INFO] Running Louvain community detection ...")
    part = louvain_partition(H, RANDOM_SEED)
    n_comm = len(set(part.values()))
    print(f"[INFO] Communities found: {n_comm}")

    # 5) 导出：HTML 互动图 + CSV 社区表 + GEXF 网络文件
    print("[INFO] Exporting Pyvis HTML + CSV + GEXF ...")
    export_pyvis(H, part, id2symbol, id2full, id2desc, OUT_HTML, OUT_CSV, OUT_GEXF)

    print("\n[DONE]")
    print("HTML :", OUT_HTML)
    print("CSV  :", OUT_CSV)
    print("GEXF :", OUT_GEXF)


if __name__ == "__main__":
    main()