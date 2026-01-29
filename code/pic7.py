import os
import pandas as pd
import holoviews as hv
from holoviews import opts
import networkx as nx
import numpy as np

# --- 1. 初始化引擎 ---
# 必须先执行这一步，否则无法生成交互图表
hv.extension('bokeh')
print("✅ 绘图引擎初始化成功。")

# --- 2. 参数配置 ---
KEYWORD = 'ribosomal'  # 搜索关键词
TOP_N = 30  # 弦图节点数（建议20-40，太多会乱）
SCORE_MIN = 400  # 相互作用置信度阈值

# --- 3. 加载并处理数据 ---
print("正在读取本地数据...")

# 更新数据路径
data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
df_info = pd.read_csv(os.path.join(data_dir, '10090.protein.info.v12.0.txt.gz'), sep='\t')
df_links = pd.read_csv(os.path.join(data_dir, '10090.protein.links.v12.0.txt.gz'), sep=' ')

# 获取 ID 列名（适配您的文件：#string_protein_id）
id_col = df_info.columns[0]

# 筛选包含关键词的蛋白 (注意您的文件列名是 'annotation')
print(f"正在根据关键词 '{KEYWORD}' 筛选核心蛋白...")
subset_info = df_info[df_info['annotation'].str.contains(KEYWORD, case=False, na=False)].head(TOP_N)

if subset_info.empty:
    print("❌ 未匹配到任何蛋白，请检查关键词或文件内容。")
else:
    # 建立映射和 ID 集合
    id_map = dict(zip(subset_info[id_col], subset_info['preferred_name']))
    target_ids = set(subset_info[id_col])

    # 提取这些蛋白之间的连边（诱导子图）
    sub_links = df_links[
        (df_links['protein1'].isin(target_ids)) &
        (df_links['protein2'].isin(target_ids)) &
        (df_links['combined_score'] >= SCORE_MIN)
        ].copy()

    # ID 转为易读名称
    sub_links['source'] = sub_links['protein1'].map(id_map)
    sub_links['target'] = sub_links['protein2'].map(id_map)

    # --- 4. 构建并美化弦图 ---
    print("正在构建网络并渲染弦图...")

    # 定义节点数据集（确保所有筛选出的蛋白都在圆周上）
    nodes = hv.Dataset(subset_info['preferred_name'].unique(), 'index')
    # 定义边数据集
    edges = sub_links[['source', 'target', 'combined_score']]

    # 创建弦图对象
    chord = hv.Chord((edges, nodes))

    # 应用视觉规范
    # 注意：为了解决之前的报错，删除了 bg_fill，改用 hooks 或直接不设背景（默认白底）
    chord.opts(
        opts.Chord(
            width=800,
            height=800,
            title=f"Interaction Subnetwork: {KEYWORD.capitalize()} Proteins",

            # 节点样式
            node_color='index',
            cmap='Category20',
            labels='index',
            label_text_font='Arial',

            # 连线（弦）样式
            edge_color='combined_score',
            edge_cmap='viridis',  # 高分连线更亮（黄），低分偏暗（紫）
            edge_alpha=0.7,

            # 交互功能
            tools=['hover'],
            selection_policy='nodes'
        )
    ).opts(
        # 这种方式设置背景色更稳健，如果还报错可删除此行
        opts.Chord(hooks=[lambda plot, element: plot.state.update(background_fill_color="#F8F9FA")])
    )

    # --- 5. 保存结果 ---
    output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"subnetwork_chord_{KEYWORD}.html")
    
    hv.save(chord, output_file)
    print(f"✨ 成功！请在文件夹中打开 [{output_file}] 查看效果。")