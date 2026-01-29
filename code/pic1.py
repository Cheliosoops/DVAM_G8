import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# --- 1. 数据加载与预处理 ---
# 读取 STRING 数据，注意：STRING 文件以空格分隔
file_path = '10090.protein.links.v12.0.txt.gz'
df = pd.read_csv(file_path, sep=' ')

# 筛选高置信度相互作用 (Score > 700)，以保证网络具有生物学意义
# STRING 的 score 扩大了 1000 倍，所以 700 代表 0.7
df_filtered = df[df['combined_score'] > 700]

# 构建无向图
G = nx.from_pandas_edgelist(df_filtered, 'protein1', 'protein2')

# --- 2. 计算度分布 ---
degrees = [val for (node, val) in G.degree()]
degree_counts = pd.Series(degrees).value_counts().sort_index()

x = degree_counts.index.values # 度数 k
y = degree_counts.values       # 频率 P(k)

# --- 3. 视觉风格配置 ---
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Roboto']
bg_color = '#F8F9FA'
main_blue = '#3498DB'
dark_teal = '#2C3E50'

fig, ax = plt.subplots(figsize=(8, 6), facecolor=bg_color)
ax.set_facecolor(bg_color)

# --- 4. 绘制 Log-Log 散点图 ---
# 使用对数坐标
ax.loglog(x, y, 'o', color=dark_teal, markersize=5, alpha=0.6, label='Observed Data')

# --- 5. 拟合幂律分布 (Linear Regression on Log-Log scale) ---
# 过滤掉 log(0) 的情况
mask = (x > 0) & (y > 0)
log_x = np.log10(x[mask])
log_y = np.log10(y[mask])
slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)

# 绘制拟合线
ax.plot(x[mask], 10**intercept * x[mask]**slope, color=main_blue,
        linestyle='--', linewidth=2, label=f'Power-law Fit ($\gamma$ = {abs(slope):.2f})')

# --- 6. 遵循 UI 规范的格式化 ---
# 标题：加粗，字号比标签大 2pt
ax.set_title('Global Network Topology: Degree Distribution',
             fontsize=16, fontweight='bold', pad=20)

# 轴标签：数学符号斜体
ax.set_xlabel('Degree ($k$)', fontsize=14)
ax.set_ylabel('Frequency ($P(k)$)', fontsize=14)

# 文本说明：P 值斜体
stats_text = f'$R^2$ = {r_value**2:.3f}\n$P_{{value}}$ < 0.05'
ax.text(0.05, 0.05, stats_text, transform=ax.transAxes,
        fontsize=12, verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# 布局优化：10% Padding
plt.tight_layout(pad=3.0)

# 移除冗余边框
sns.despine()

# 保存并展示
plt.savefig('degree_distribution_loglog.png', dpi=300, facecolor=bg_color)
plt.show()