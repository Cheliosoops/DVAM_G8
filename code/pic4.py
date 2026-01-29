#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
方向四：相互作用证据分布（STRING 多通道子得分）
输出：
1) fig4_evidence_share_th{cutoff}.html
   - 各 evidence 通道"非零边占比"柱状图（交互）
2) fig4_evidence_radar_keyproteins_th{cutoff}.html
   - Top hubs 关键蛋白在各证据通道的平均得分雷达图（交互）
3) evidence_summary_th{cutoff}.csv
   - 证据通道统计汇总表

UI/UX 规范（全组统一）：
- 背景：#F8F9FA
- 字体：Arial（sans-serif）
- 主色：#2C3E50 / #3498DB
- 强调色（关键蛋白）：#F39C12
- 聚类/多类别配色：Tableau10 或 Viridis（此脚本中：雷达图多条曲线用 Tableau10，并把第一名 hub 强制金色）
"""

import os
import numpy as np
import pandas as pd

import plotly.graph_objects as go
from plotly.colors import qualitative

# -----------------------------
# 0) 路径与统一 UI 参数（注意：你要求的 Windows 路径）
# -----------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

DETAILED_GZ = os.path.join(DATA_DIR, "10090.protein.links.detailed.v12.0.txt.gz")
INFO_GZ     = os.path.join(DATA_DIR, "10090.protein.info.v12.0.txt.gz")

# 若你已经跑过方向三，会生成这个文件；脚本会优先用它来选 top hubs
COMMUNITY_ASSIGN_CSV = os.path.join(DATA_DIR, "community_assignments_th700.csv")  # 如果 cutoff 改了，下面会自动替换

SCORE_CUTOFF = 700          # 400/700/900 可调整（建议与你方向三一致）
TOP_HUBS = 20               # "关键蛋白候选"数量（通常用于表格/强调）
RADAR_TOPN = 8              # 雷达图展示的关键蛋白数量（建议 5~10，太多会乱）

# 输出
OUT_SUMMARY_CSV = os.path.join(os.path.dirname(__file__), "..", "outputs", f"evidence_summary_th{SCORE_CUTOFF}.csv")
OUT_BAR_HTML    = os.path.join(os.path.dirname(__file__), "..", "figures", f"fig4_evidence_share_th{SCORE_CUTOFF}.html")
OUT_RADAR_HTML  = os.path.join(os.path.dirname(__file__), "..", "figures", f"fig4_evidence_radar_keyproteins_th{SCORE_CUTOFF}.html")

# 创建输出目录
os.makedirs(os.path.dirname(OUT_SUMMARY_CSV), exist_ok=True)
os.makedirs(os.path.dirname(OUT_BAR_HTML), exist_ok=True)
os.makedirs(os.path.dirname(OUT_RADAR_HTML), exist_ok=True)

# ---- UI 颜色（严格统一）----
BG_COLOR   = "#F8F9FA"
TEXT_COLOR = "#2C3E50"
MAIN_BLUE  = "#3498DB"
HUB_COLOR  = "#F39C12"
GRAY_NS    = "#B0B0B0"   # 语义色：非显著可用（本方向通常不涉及 p 值，但预留）

# ---- 字体统一 ----
FONT_FAMILY = "Arial"

# -----------------------------
# 1) 读取 protein.info：构建 id -> symbol/annotation
# -----------------------------
def load_info(info_path: str):
    info = pd.read_csv(info_path, sep="\t", compression="gzip")

    # STRING 常见列名：protein_external_id, preferred_name, annotation
    id_col = "protein_external_id" if "protein_external_id" in info.columns else info.columns[0]
    sym_col = "preferred_name" if "preferred_name" in info.columns else (info.columns[1] if len(info.columns) > 1 else id_col)

    desc_candidates = [c for c in info.columns if ("annot" in c.lower() or "desc" in c.lower())]
    desc_col = desc_candidates[0] if desc_candidates else None

    id2symbol = dict(zip(info[id_col].astype(str), info[sym_col].astype(str)))
    id2desc = dict(zip(info[id_col].astype(str), info[desc_col].astype(str))) if desc_col else {}

    return id2symbol, id2desc


# -----------------------------
# 2) 读取 links.detailed 并过滤阈值
# -----------------------------
def load_detailed_edges(detailed_path: str, cutoff: int) -> pd.DataFrame:
    """
    STRING links.detailed 一般包含：
    protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score

    这里做鲁棒处理：如果列名不同，会尝试自动识别/重命名。
    """
    df = pd.read_csv(detailed_path, sep=r"\s+", compression="gzip", engine="python")

    # 若文件没有 header（极少见），需要手动命名；这里做兜底
    if "protein1" not in df.columns or "protein2" not in df.columns:
        # 常见长度为 10 列（2 + 7 + 1）
        if df.shape[1] == 10:
            df.columns = [
                "protein1", "protein2",
                "neighborhood", "fusion", "cooccurence", "coexpression",
                "experimental", "database", "textmining",
                "combined_score"
            ]
        else:
            raise ValueError(f"无法识别列名/列数：当前 df.shape={df.shape}，请检查文件是否为 STRING links.detailed。")

    # 确保 combined_score 存在
    if "combined_score" not in df.columns:
        # 有时会叫 combined_score 或 combined
        cands = [c for c in df.columns if "combined" in c.lower()]
        if not cands:
            raise ValueError("未找到 combined_score 列，请检查 links.detailed 文件格式。")
        df = df.rename(columns={cands[0]: "combined_score"})

    # 按阈值过滤边
    df = df[df["combined_score"].astype(int) >= cutoff].copy()

    # 统一 protein id 为 str
    df["protein1"] = df["protein1"].astype(str)
    df["protein2"] = df["protein2"].astype(str)

    return df


# -----------------------------
# 3) 证据通道统计：非零占比 + 分布统计
# -----------------------------
def summarize_evidence(df: pd.DataFrame) -> pd.DataFrame:
    # 识别证据通道列（排除 protein1/protein2/combined_score）
    ignore = {"protein1", "protein2", "combined_score"}
    evidence_cols = [c for c in df.columns if c not in ignore]

    records = []
    total_edges = len(df)

    for col in evidence_cols:
        s = df[col].astype(float)
        nonzero = (s > 0).sum()
        nonzero_ratio = nonzero / total_edges if total_edges > 0 else 0.0

        s_nz = s[s > 0]
        rec = {
            "evidence_channel": col,
            "total_edges": total_edges,
            "nonzero_edges": int(nonzero),
            "nonzero_ratio": float(nonzero_ratio),
            "mean_nonzero_score": float(s_nz.mean()) if len(s_nz) else 0.0,
            "median_nonzero_score": float(s_nz.median()) if len(s_nz) else 0.0,
            "max_score": float(s.max()) if len(s) else 0.0
        }
        records.append(rec)

    out = pd.DataFrame(records).sort_values("nonzero_ratio", ascending=False).reset_index(drop=True)
    return out


# -----------------------------
# 4) 关键蛋白选择（优先用方向三输出；否则从边集计算度）
# -----------------------------
def get_top_hubs(df_edges: pd.DataFrame, assign_csv_path: str, cutoff: int, topk: int) -> list:
    """
    返回 protein_id 列表（按度从高到低）
    """
    # 如果有方向三的 community_assignments 文件：直接用里面的 degree 排序
    candidate_path = assign_csv_path.replace("th700", f"th{cutoff}")
    if os.path.exists(candidate_path):
        ass = pd.read_csv(candidate_path)
        if "protein_id" in ass.columns and "degree" in ass.columns:
            ass = ass.sort_values("degree", ascending=False)
            return ass["protein_id"].astype(str).head(topk).tolist()

    # 否则：从 df_edges 统计度（无向图：端点出现次数）
    deg = pd.concat([df_edges["protein1"], df_edges["protein2"]]).value_counts()
    return deg.head(topk).index.astype(str).tolist()


# -----------------------------
# 5) 为关键蛋白计算"各证据通道平均得分"（用于雷达图）
# -----------------------------
def compute_protein_evidence_profile(df_edges: pd.DataFrame, protein_ids: list) -> pd.DataFrame:
    """
    对每个 protein，统计其 incident edges 的各 evidence score 均值（0~1000）
    """
    ignore = {"protein1", "protein2", "combined_score"}
    evidence_cols = [c for c in df_edges.columns if c not in ignore]

    # 把边展开到端点视角：每条边对两个端点各贡献一条记录
    df1 = df_edges[["protein1"] + evidence_cols].rename(columns={"protein1": "protein"})
    df2 = df_edges[["protein2"] + evidence_cols].rename(columns={"protein2": "protein"})
    df_end = pd.concat([df1, df2], ignore_index=True)

    # 只保留关键蛋白
    df_end = df_end[df_end["protein"].isin(protein_ids)].copy()

    # degree：每个蛋白出现次数（即 incident edges 数）
    deg = df_end.groupby("protein").size().rename("degree").reset_index()

    # 各证据通道均值
    mean_scores = df_end.groupby("protein")[evidence_cols].mean().reset_index()

    out = mean_scores.merge(deg, on="protein", how="left")
    # 按 degree 排序
    out = out.sort_values("degree", ascending=False).reset_index(drop=True)
    return out


# -----------------------------
# 6) 绘图（Plotly）：证据占比柱状图
# -----------------------------
def plot_evidence_share(summary_df: pd.DataFrame, out_html: str, cutoff: int):
    # X：通道，Y：非零占比
    x = summary_df["evidence_channel"].tolist()
    y = (summary_df["nonzero_ratio"] * 100).tolist()  # 转成百分比

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=x,
            y=y,
            marker_color=MAIN_BLUE,
            hovertemplate=(
                "<b>Evidence</b>: %{x}<br>"
                "<b>Non-zero edges</b>: %{customdata[0]}<br>"
                "<b>Total edges</b>: %{customdata[1]}<br>"
                "<b>Ratio</b>: %{y:.2f}%<extra></extra>"
            ),
            customdata=np.stack([summary_df["nonzero_edges"], summary_df["total_edges"]], axis=1)
        )
    )

    fig.update_layout(
        title=dict(
            text=f"<b>Figure 4A. Evidence Channel Coverage (cutoff={cutoff})</b>",
            font=dict(family=FONT_FAMILY, size=18, color=TEXT_COLOR)  # 标题比轴标签大 2pt（下面轴标签 16）
        ),
        font=dict(family=FONT_FAMILY, size=16, color=TEXT_COLOR),
        paper_bgcolor=BG_COLOR,
        plot_bgcolor=BG_COLOR,
        margin=dict(l=80, r=80, t=90, b=80),  # ~10% padding 的效果
        xaxis=dict(title="Evidence channel（证据通道）"),
        yaxis=dict(title="Non-zero edge ratio（%）", rangemode="tozero"),
    )

    fig.write_html(out_html, include_plotlyjs="cdn")
    print("[OK] 写出证据占比柱状图：", out_html)


# -----------------------------
# 7) 绘图（Plotly）：关键蛋白雷达图
# -----------------------------
def plot_radar(prof_df: pd.DataFrame, id2symbol: dict, id2desc: dict, out_html: str, cutoff: int, topn: int):
    """
    prof_df：包含 protein + evidence_cols + degree
    """
    # 证据通道列
    ignore = {"protein", "degree"}
    evidence_cols = [c for c in prof_df.columns if c not in ignore]

    # 雷达图的角度标签
    theta = evidence_cols + [evidence_cols[0]]  # 闭合

    # 配色：Tableau10（严格统一），并把第一名 hub 强制为金色
    tableau10 = qualitative.T10
    fig = go.Figure()

    # 只展示 topn 条曲线
    show = prof_df.head(topn).copy()

    for i, row in show.iterrows():
        pid = row["protein"]
        symbol = id2symbol.get(pid, pid)
        desc = id2desc.get(pid, "")
        degree = int(row["degree"])

        r = [float(row[c]) for c in evidence_cols]
        r = r + [r[0]]  # 闭合

        # 颜色规则：第一条（金色强调），其余用 Tableau10 循环
        if i == 0:
            color = HUB_COLOR
        else:
            color = tableau10[(i - 1) % len(tableau10)]

        # hover 信息：Symbol / 描述 / Degree / 各通道值
        hover = (
            f"<b>Symbol</b>: {symbol}<br>"
            f"<b>Description</b>: {desc}<br>"
            f"<b>Degree</b>: {degree}<br>"
            f"<extra></extra>"
        )

        fig.add_trace(
            go.Scatterpolar(
                r=r,
                theta=theta,
                mode="lines+markers",
                name=f"{symbol} (deg={degree})",
                line=dict(color=color, width=2),
                marker=dict(size=5, color=color),
                hovertemplate=hover
            )
        )

    fig.update_layout(
        title=dict(
            text=f"<b>Figure 4B. Key Proteins Evidence Profile (cutoff={cutoff})</b>",
            font=dict(family=FONT_FAMILY, size=18, color=TEXT_COLOR)
        ),
        font=dict(family=FONT_FAMILY, size=16, color=TEXT_COLOR),
        paper_bgcolor=BG_COLOR,
        plot_bgcolor=BG_COLOR,
        margin=dict(l=90, r=90, t=90, b=80),
        polar=dict(
            bgcolor=BG_COLOR,
            radialaxis=dict(
                visible=True,
                range=[0, 1000],   # STRING 子得分通常 0~1000
                tickfont=dict(family=FONT_FAMILY, size=14, color=TEXT_COLOR),
                gridcolor="rgba(44,62,80,0.15)"
            ),
            angularaxis=dict(
                tickfont=dict(family=FONT_FAMILY, size=14, color=TEXT_COLOR),
                gridcolor="rgba(44,62,80,0.15)"
            )
        ),
        legend=dict(
            bgcolor="rgba(248,249,250,0.6)",
            bordercolor="rgba(44,62,80,0.15)",
            borderwidth=1
        )
    )

    fig.write_html(out_html, include_plotlyjs="cdn")
    print("[OK] 写出关键蛋白雷达图：", out_html)


def main():
    # 基本检查
    if not os.path.exists(DETAILED_GZ):
        print(f"错误: 找不到文件 {DETAILED_GZ}")
        print("请确保已下载 STRING 数据文件到正确的数据目录中")
        return
        
    if not os.path.exists(INFO_GZ):
        print(f"错误: 找不到文件 {INFO_GZ}")
        print("请确保已下载 STRING 数据文件到正确的数据目录中")
        return

    # 1) 读取 info（用于 tooltip/名称）
    id2symbol, id2desc = load_info(INFO_GZ)

    # 2) 读取 detailed 边并按 cutoff 过滤
    print(f"[INFO] 读取 links.detailed 并按 cutoff={SCORE_CUTOFF} 过滤 ...")
    df = load_detailed_edges(DETAILED_GZ, SCORE_CUTOFF)
    print(f"[INFO] 过滤后边数：{len(df):,}")

    # 3) 证据通道统计汇总
    summary = summarize_evidence(df)
    summary.to_csv(OUT_SUMMARY_CSV, index=False, encoding="utf-8-sig")
    print("[OK] 写出证据统计表：", OUT_SUMMARY_CSV)

    # 4) 选 top hubs（优先用方向三输出；否则从 df 算）
    hubs = get_top_hubs(df, COMMUNITY_ASSIGN_CSV, SCORE_CUTOFF, TOP_HUBS)
    print(f"[INFO] Top hubs（前 {TOP_HUBS}）：", hubs[:10], "..." if len(hubs) > 10 else "")

    # 5) 计算关键蛋白证据画像（用于雷达图）
    prof = compute_protein_evidence_profile(df, hubs)
    # 6) 画证据占比柱状图
    plot_evidence_share(summary, OUT_BAR_HTML, SCORE_CUTOFF)
    # 7) 画关键蛋白雷达图
    plot_radar(prof, id2symbol, id2desc, OUT_RADAR_HTML, SCORE_CUTOFF, RADAR_TOPN)

    print("\n[DONE]")
    print("Evidence summary CSV:", OUT_SUMMARY_CSV)
    print("Bar chart HTML      :", OUT_BAR_HTML)
    print("Radar chart HTML    :", OUT_RADAR_HTML)


if __name__ == "__main__":
    main()