# DVAM_G8

# 蛋白质相互作用网络可视化分析项目

## 📚 项目简介

本项目基于STRING数据库的小鼠蛋白质相互作用网络数据，通过多种可视化技术对蛋白质网络的结构和功能关系进行深入分析。我们实现了社区检测、相互作用证据分布分析、以及核心蛋白子网络可视化等功能，揭示了蛋白质间复杂的生物学关系。

**数据来源**: STRING Database (https://string-db.org/)  
**物种**: 小家鼠 (Mus musculus, taxon ID: 10090)  
**数据版本**: v12.0

**研究内容**:
- 蛋白质网络社区结构分析
- 相互作用证据通道分布
- 核心蛋白识别与可视化
- 核糖体蛋白子网络分析
- 功能富集分析

## 🎯 研究目标

1. **网络拓扑分析**: 识别蛋白质网络中的社区结构，揭示功能模块
2. **关键蛋白识别**: 通过度中心性分析发现网络中的核心蛋白
3. **证据权重解析**: 分析不同证据通道对蛋白质相互作用的贡献
4. **特定功能网络可视化**: 聚焦核糖体等特定功能蛋白群的相互作用模式
5. **交互式可视化**: 提供直观的交互式图表，便于生物学家探索数据

## 📁 项目结构

```
protein_network_visualization/
│
├── data/                                    # 数据文件
│   ├── 10090.protein.links.v12.0.txt.gz    # 蛋白质相互作用链接数据
│   ├── 10090.protein.links.detailed.v12.0.txt.gz # 详细相互作用证据
│   ├── 10090.protein.info.v12.0.txt.gz     # 蛋白质信息（注释、别名等）
│   └── ...                                 # 其他数据文件
│
├── code/                                    # 源代码
│   ├── pic1.py                             # 网络拓扑分析 - 度分布
│   ├── pic2&pic6.py                        # 网络分析和聚类系数
│   ├── pic3.py                             # 社区网络可视化
│   ├── pic4.py                             # 相互作用证据分布分析
│   ├── pic5.py                             # GO功能富集分析
│   └── pic7.py                             # 核糖体蛋白子网络可视化
│
├── Visualization Assets/                    # 生成的可视化图表
│   ├── fig3_community_network_th900.html   # 社区网络图（阈值900）
│   ├── fig4_evidence_radar_keyproteins_th700.html # 证据雷达图
│   ├── fig4_evidence_share_th700.html      # 证据分布柱状图
│   ├── subnetwork_chord_ribosomal.html     # 核糖体蛋白弦图
│   ├── pic1.png                           # 度分布图
│   ├── pic5.png                           # 功能富集分析图
│   └── pic7.png                           # 特定功能子网络图
│
├── README.md                               # 项目说明文档
├── requirements.txt                        # Python依赖包列表
├── index.html                              # 综合可视化展示网页
└── Visual Analysis Report.pdf              # 分析报告
```

## 🔧 环境配置

### 依赖包

完整的依赖列表请参见 requirements.txt 文件。

### 安装步骤

1. **克隆或下载项目**
```bash
git clone <your-repo-url>
cd protein_network_visualization
```

2. **创建虚拟环境（推荐）**
```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
```

3. **安装依赖**
```bash
pip install -r requirements.txt
```

4. **下载STRING数据文件**
```bash
# 下载必需的数据文件到 data/ 目录
# 10090.protein.links.v12.0.txt.gz
# 10090.protein.links.detailed.v12.0.txt.gz
# 10090.protein.info.v12.0.txt.gz
```

## 🚀 使用方法

### 1. 运行所有可视化脚本

按顺序运行代码目录下的脚本以生成所有可视化结果：

```bash
# 网络拓扑分析 - 度分布
python code/pic1.py

# 网络分析和聚类系数
python code/pic2&pic6.py

# 社区网络可视化
python code/pic3.py

# 相互作用证据分析
python code/pic4.py

# GO功能富集分析
python code/pic5.py

# 核糖体蛋白子网络可视化
python code/pic7.py
```

### 2. 查看综合可视化结果

运行所有脚本后，可以通过浏览器打开[index.html](file:///d:%5C%E7%A0%94%E7%A9%B6%E7%94%9F%5C%E6%95%B0%E6%8D%AE%E5%8F%AF%E8%A7%86%E5%8C%96/index.html)查看综合可视化结果。

## 📊 主要发现

### 1. 网络拓扑特征

- **社区结构**: 蛋白质网络呈现出明显的社区结构，反映了生物学上的功能模块
- **枢纽蛋白**: 发现了一批高度连接的核心蛋白，可能在细胞功能中发挥关键作用
- **连接密度**: 不同阈值下的网络表现出不同的拓扑特征

### 2. 相互作用证据分布

- **主导证据类型**: 实验和数据库证据是最主要的相互作用支持
- **共表达证据**: 在某些功能模块中占比较高，表明转录调控的相关性
- **文本挖掘证据**: 补充了大量潜在的相互作用

### 3. 核心蛋白功能

- **高连接性蛋白**: 多为管家基因产物，参与基本细胞过程
- **功能富集**: 核心蛋白在特定生物学通路中显著富集
- **疾病关联**: 部分核心蛋白与人类疾病相关联

### 4. 特定功能网络

- **核糖体网络**: 形成了紧密的功能模块，反映了蛋白质合成的协同性
- **代谢网络**: 展现出高度的组织性和调控性
- **信号传导网络**: 包含多个枢纽节点，体现信号传递的集中性

## 📈 可视化图表说明

### 图1: 网络拓扑分析 - 度分布
- 展示蛋白质网络的度分布特征
- 使用对数坐标绘制，揭示网络的无标度特性

### 图3: 社区网络图
- 展示蛋白质网络的社区结构
- 不同颜色代表不同功能模块
- 黄色节点标识枢纽蛋白（度最高的前25个）
- 支持缩放、拖拽和悬停查看详细信息

### 图4A: 证据通道分布图
- 柱状图展示各证据通道的非零边占比
- 悬停可查看具体数值和统计信息
- 帮助理解相互作用数据的可靠性来源

### 图4B: 关键蛋白证据雷达图
- 雷达图展示关键蛋白在各证据通道的表现
- 金色线条突出显示最重要的枢纽蛋白
- 揭示不同蛋白的相互作用特点

### 图5: GO功能富集分析
- 气泡图展示显著富集的功能类别
- 横轴表示基因比例，纵轴为功能类别
- 气泡大小和颜色表示富集基因数量和统计显著性

### 弦图: 特定功能子网络
- 弦图展示核糖体蛋白等特定功能群体的相互作用
- 弦的宽度反映相互作用强度
- 便于观察局部网络结构

## 💡 生物学意义

1. **功能预测**: 通过网络邻近性预测未知蛋白的功能
2. **药物靶点识别**: 识别网络中的关键节点作为潜在药物靶点
3. **疾病机制探索**: 理解疾病相关蛋白在网络中的位置和作用
4. **进化保守性分析**: 比较不同物种网络的相似性

## 🧪 技术亮点

1. **多层次可视化**: 从全局网络到局部子网络的多层次展示
2. **交互性设计**: 提供丰富的交互功能，便于深入探索
3. **标准化色彩**: 统一的色彩方案，增强图表一致性
4. **性能优化**: 针对大规模网络的可视化优化

## 👥 团队成员

- [徐永祺] - 学号: 2501212905
- [刘祥瑞] - 学号: 2501212928
- [赵张欣悦] - 学号: 2501212909

## 📝 参考文献

1. STRING database: https://string-db.org/
2. Szklarczyk, D., et al. "The STRING database in 2023: protein-protein association networks and beyond." Nucleic Acids Research (2023).
3. Borgatti, S.P., & Everett, M.G. "A Graph-Theoretic Perspective on Centrality." Social Networks (2006).
4. Brandes, U., et al. "On variants of shortest-path betweenness centrality and their generic computation." Social Networks (2008).

## 📄 许可证

本项目仅供学术研究使用

---

**最后更新**: 2026年1月

**数据版本**: STRING v12.0
