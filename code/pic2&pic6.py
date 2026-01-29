import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ==========================================
# 1. UI/UX 全局视觉规范配置
# ==========================================
COLOR_MAIN_DARK = "#2C3E50"  # Scientific Blue/Dark
COLOR_MAIN_LIGHT = "#3498DB"  # Scientific Teal/Blue
COLOR_HIGHLIGHT = "#F39C12"  # Orange/Gold
COLOR_BG = "#FFFFFF"  # 纯白背景
FONT_FAMILY = "Arial, Roboto, sans-serif"


def apply_layout_style(fig, title_text, x_title, y_title):
    """应用统一的图表布局样式"""
    fig.update_layout(
        title={
            'text': f"<b>{title_text}</b>",
            'y': 0.95, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top',
            'font': {'size': 18}
        },
        paper_bgcolor=COLOR_BG,
        plot_bgcolor=COLOR_BG,
        font=dict(family=FONT_FAMILY, color=COLOR_MAIN_DARK),
        margin=dict(l=80, r=80, t=80, b=80),
        hovermode="closest"
    )
    # 【修改点1】全局去除网格线 (showgrid=False)
    fig.update_xaxes(
        title_text=x_title,
        showline=True, linewidth=1, linecolor='lightgrey', mirror=True,
        showgrid=False  # 不显示网格线
    )
    fig.update_yaxes(
        title_text=y_title,
        showline=True, linewidth=1, linecolor='lightgrey', mirror=True,
        showgrid=False  # 不显示网格线
        # gridcolor='#F0F0F0' # 去掉网格颜色设置
    )


# ==========================================
# 2. 数据读取与清洗 (指定绝对路径)
# ==========================================
print("--- 阶段 1: 数据加载 ---")

# 定义文件路径 (使用 r"" 防止转义字符报错)
path_info = r"C:\Users\lxr\Desktop\10090.protein.info.v12.0.txt.gz"
path_links = r"C:\Users\lxr\Desktop\10090.protein.links.v12.0.txt.gz"

try:
    # A. 读取 ID 映射文件 (Info)
    print(f"正在读取映射表: {path_info} ...")
    df_info = pd.read_csv(
        path_info,
        sep='\t',
        compression='gzip',
        usecols=['#string_protein_id', 'preferred_name']
    )
    # 构建字典
    id_to_name = dict(zip(df_info['#string_protein_id'], df_info['preferred_name']))

    # B. 读取网络连边文件 (Links)
    print(f"正在读取连边数据: {path_links} ...")
    df_links = pd.read_csv(
        path_links,
        sep=' ',
        compression='gzip'
    )

    # C. 数据融合与清洗
    print("正在清洗与映射数据...")
    # 过滤 (保留 score >= 400)
    df_clean = df_links[df_links['combined_score'] >= 400].copy()

    # 将 ID 替换为基因名
    df_clean['node1'] = df_clean['protein1'].map(id_to_name).fillna(df_clean['protein1'])
    df_clean['node2'] = df_clean['protein2'].map(id_to_name).fillna(df_clean['protein2'])

    # 提取最终 DataFrame
    df = df_clean[['node1', 'node2', 'combined_score']]
    print(f" 数据准备完毕! 包含 {len(df)} 条连边")

except FileNotFoundError as e:
    print(f"\n 错误: 找不到文件。请检查路径是否正确：\n{e.filename}")
    exit()
except Exception as e:
    print(f"\n 读取错误: {e}")
    print("如果文件其实已经解压了(不是.gz)，请去掉代码里的 compression='gzip' 参数。")
    exit()

# ==========================================
# 3. 分析方向 2: 关键节点识别 (Lollipop Chart)
# ==========================================
print("\n--- 阶段 2: 生成关键节点图 (Top 20 Hubs) ---")

G = nx.from_pandas_edgelist(df, 'node1', 'node2')
degree_dict = nx.degree_centrality(G)

# 获取 Top 20
top20_df = pd.DataFrame(list(degree_dict.items()), columns=['Protein', 'Degree']) \
    .sort_values('Degree', ascending=True).tail(20)

fig_hub = go.Figure()

# 绘制棒 (Lines)
for i, row in top20_df.iterrows():
    is_top3 = i in top20_df.index[-3:]
    line_color = COLOR_HIGHLIGHT if is_top3 else COLOR_MAIN_LIGHT

    fig_hub.add_shape(type='line',
                      x0=0, y0=row['Protein'],
                      x1=row['Degree'], y1=row['Protein'],
                      line=dict(color=line_color, width=3)
                      )

# 绘制糖 (Markers)
# 普通 Top 20
fig_hub.add_trace(go.Scatter(
    x=top20_df['Degree'][:-3], y=top20_df['Protein'][:-3],
    mode='markers', name='Key Proteins',
    marker=dict(color=COLOR_MAIN_LIGHT, size=12),
    hovertemplate="<b>%{y}</b><br>Degree: %{x:.4f}<extra></extra>"
))

# 核心 Top 3
fig_hub.add_trace(go.Scatter(
    x=top20_df['Degree'][-3:], y=top20_df['Protein'][-3:],
    mode='markers', name='Top 3 Hubs',
    marker=dict(color=COLOR_HIGHLIGHT, size=16, line=dict(color='black', width=1)),
    hovertemplate="<b>%{y}</b><br>Degree: %{x:.4f}<br>Rank: Top 3<extra></extra>"
))

apply_layout_style(fig_hub, "Top 20 Hub Proteins (Degree Centrality)", "Degree Centrality Score", "")

# 【修改点3】调整横坐标间距
fig_hub.update_xaxes(
    rangemode="tozero",
    tickmode='linear',  # 强制使用线性刻度
    dtick=0.02,  # 设置刻度间隔为 0.02 (根据数据范围调整)
)
fig_hub.show()

# ==========================================
# 4. 分析方向 6: 阈值敏感性分析 (Dual-Axis Chart)
# ==========================================
print("\n--- 阶段 3: 生成阈值敏感性分析图 ---")

# 定义阈值范围
thresholds = range(400, 951, 50)
node_counts = []
edge_counts = []

print("正在计算不同阈值下的网络拓扑...")
for thresh in thresholds:
    sub_df = df[df['combined_score'] >= thresh]
    nodes = set(sub_df['node1']).union(set(sub_df['node2']))

    node_counts.append(len(nodes))
    edge_counts.append(len(sub_df))

# 绘图 - 双轴
fig_sens = make_subplots(specs=[[{"secondary_y": True}]])

# 左轴: 节点数
fig_sens.add_trace(
    go.Scatter(
        x=list(thresholds), y=node_counts,
        name="Nodes (节点数)",
        mode='lines+markers',
        line=dict(color=COLOR_MAIN_LIGHT, width=3),
        marker=dict(size=8)
    ), secondary_y=False
)

# 右轴: 边数
fig_sens.add_trace(
    go.Scatter(
        x=list(thresholds), y=edge_counts,
        name="Edges (连边数)",
        mode='lines+markers',
        line=dict(color=COLOR_MAIN_DARK, width=3, dash='dot'),
        marker=dict(symbol='diamond', size=8, color=COLOR_MAIN_DARK)
    ), secondary_y=True
)

apply_layout_style(fig_sens, "Network Sensitivity Analysis", "Confidence Score Threshold", "Number of Nodes")
# 右轴本来就不显示网格，这里再次确认
fig_sens.update_yaxes(title_text="Number of Edges", secondary_y=True, showgrid=False)

# 【修改点2】去掉了添加虚线参考线的代码
# for t in [400, 700, 900]:
#     fig_sens.add_vline(x=t, line_width=1, line_dash="dash", line_color="gray", opacity=0.5)

fig_sens.show()

print("\n 所有图表已生成完毕！")