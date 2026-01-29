# Data Download Guide / 数据下载指南

本项目分析所使用的蛋白质相互作用（PPI）数据及功能富集数据均来源于 **STRING Database (v12.0)**。目标物种为 **小鼠 (Mus musculus)**，其 Taxonomy ID 为 `10090`。

## 1. 核心数据文件 (Data Files)

你可以通过下表中的链接直接获取原始数据文件：

| 文件名 | 内容描述 | 下载链接 (Direct Link) |
| :--- | :--- | :--- |
| **10090.protein.links.v12.0.txt.gz** | 蛋白质相互作用网络数据，包含各蛋白质对的综合得分（Combined Score）。 | [Download](https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz) |
| **10090.protein.links.detailed.v12.0.txt.gz** | 详细互作分数，包含实验、数据库、文本挖掘等不同维度的独立评分。 | [Download](https://stringdb-downloads.org/download/protein.links.detailed.v12.0/10090.protein.links.detailed.v12.0.txt.gz) |
| **10090.protein.enrichment.terms.v12.0.txt.gz** | 蛋白质对应的功能富集术语列表（涵盖 GO, KEGG, Reactome 等）。 | [Download](https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0/10090.protein.enrichment.terms.v12.0.txt.gz) |

## 2. 快速获取 (Terminal Usage)

如果你正在使用 Linux/macOS 环境或远程服务器，可以运行以下命令在当前目录下自动创建文件夹并下载数据：

```bash
# 创建数据存放目录
mkdir -p data && cd data

# 下载所有必要文件
wget [https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz)
wget [https://stringdb-downloads.org/download/protein.links.detailed.v12.0/10090.protein.links.detailed.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.links.detailed.v12.0/10090.protein.links.detailed.v12.0.txt.gz)
wget [https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0/10090.protein.enrichment.terms.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0/10090.protein.enrichment.terms.v12.0.txt.gz)

# 校验文件大小 (可选)
ls -lh
