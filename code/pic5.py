import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

# --- 1. 加载本地数据 ---
print("正在加载本地数据...")
# 读取信息文件
df_info = pd.read_csv('10090.protein.info.v12.0.txt.gz', sep='\t')
# 读取富集项文件
df_terms = pd.read_csv('10090.protein.enrichment.terms.v12.0.txt.gz', sep='\t')

# 自动获取 ID 列名（通常是第一列 #string_protein_id）
info_id_col = df_info.columns[0]
term_id_col = df_terms.columns[0]

# --- 2. 构造“功能相关”的输入基因集 ---
# 修正点：df_info 对应的描述列名为 'annotation'
print("正在筛选功能相关的测试基因集（以 ribosomal 为关键词）...")
# 使用 annotation 列进行关键词筛选
test_genes = df_info[df_info['annotation'].str.contains('ribosomal', case=False, na=False)]
my_proteins = test_genes[info_id_col].head(100).tolist()

if len(my_proteins) < 5:
    print("关键词筛选出的基因太少，尝试更换为 'Guanine' 或其他关键词。")
else:
    # --- 3. 本地计算富集 ---
    total_proteins_in_bg = df_terms[term_id_col].nunique()
    my_proteins_in_bg = [p for p in my_proteins if p in df_terms[term_id_col].values]
    n_sample = len(my_proteins_in_bg)

    # 筛选 Process (生物过程)
    # 注意：STRING 文件中 category 可能包含 'GO Biological Process'，这里确保匹配
    df_go = df_terms[df_terms['category'].str.contains('Process', na=False)]

    results = []
    # 在 df_terms 中，描述列名确实是 'description'
    for (term_desc, term_id), group in df_go.groupby(['description', 'term']):
        proteins_in_term = set(group[term_id_col])
        k = len([p for p in my_proteins_in_bg if p in proteins_in_term])

        if k >= 3:
            M = total_proteins_in_bg
            n = len(proteins_in_term)
            p_val = hypergeom.sf(k - 1, M, n, n_sample)
            results.append({'Term': term_desc, 'Count': k, 'P-value': p_val, 'Gene_Ratio': k / n})

    # --- 4. 绘图 (严格遵循 UI 规范) ---
    if results:
        res_df = pd.DataFrame(results).sort_values('P-value').head(15)

        # 设置全局字体
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

        # 创建画布，设置背景色 #F8F9FA
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='#F8F9FA')
        ax.set_facecolor('#F8F9FA')

        # 气泡图绘制
        # c: 使用 -log10(P-value) 映射颜色，符合语义（数值越大越显著）
        # cmap: 使用预设的 Viridis 配色方案
        scatter = ax.scatter(x=res_df['Gene_Ratio'],
                             y=res_df['Term'],
                             s=res_df['Count'] * 30,  # 调整气泡大小比例
                             c=-np.log10(res_df['P-value']),
                             cmap='viridis',
                             alpha=0.8,
                             edgecolors='white',
                             linewidth=0.5)

        # 颜色条设置，P 值使用斜体
        cbar = plt.colorbar(scatter)
        cbar.set_label('Significance: $-\log_{10}(P_{italic})$', fontsize=12)

        # 标题与标签：标题加粗并比轴标签大 2pt
        ax.set_title('GO Biological Process Enrichment', fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel('Gene Ratio (Hits / Term Size)', fontsize=14)
        ax.set_ylabel('Functional Categories', fontsize=14)

        # 布局优化：留出边距防止溢出
        plt.tight_layout(pad=3.0)

        plt.savefig('local_enrichment_corrected.png', dpi=300)
        plt.show()
        print(f"成功！已生成气泡图。")
    else:
        print("未发现显著富集项，请检查输入基因集。")