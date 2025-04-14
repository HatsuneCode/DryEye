## Part1： 对于差异分析的程度分类及画柱状图表示
DEG = read.table('/mnt/f/Project/GSE252984.stat/DEG.DvsCtrl.txt', sep = '\t', header = T, check.names = F)
DEG = DEG[which(DEG$p_val_adj < .05 & abs(DEG$avg_log2FC) > log2(1.5)),]

# 1. 添加 group 分类（不带换行）
DEG$group <- NA
DEG$group[DEG$upDown == "Up" & DEG$avg_log2FC > 2] <- "Highly Upregulated"
DEG$group[DEG$upDown == "Up" & DEG$avg_log2FC <= 2 & DEG$avg_log2FC > 1] <- "Mildly Upregulated"
DEG$group[DEG$upDown == "Up" & DEG$avg_log2FC <= 1 & DEG$avg_log2FC > log2(1.5)] <- "lowly Upregulated"
DEG$group[DEG$upDown == "Down" & DEG$avg_log2FC < -2] <- "Highly Downregulated"
DEG$group[DEG$upDown == "Down" & DEG$avg_log2FC >= -2 & DEG$avg_log2FC < -1] <- "Mildly Downregulated"
DEG$group[DEG$upDown == "Down" & DEG$avg_log2FC >= -1 & DEG$avg_log2FC < -log2(1.5)] <- "lowly Downregulated"

# 2. 创建带注释版本（含换行描述）用于显示
DEG$group_label <- recode(DEG$group_simple,
                          "Highly Up" = "Highly Up (log2FC > 2)",
                          "Mildly Up" = "Mildly Up (1 < log2FC ≤ 2)",
                          "Lowly Up" = "Lowly Up (log2(1.5) < log2FC ≤ 1)",
                          "Highly Down" = "Highly Down (log2FC < -2)",
                          "Mildly Down" = "Mildly Down (-2 ≤ log2FC < -1)",
                          "Lowly Down" = "Lowly Down (-1 ≤ log2FC < -log2(1.5))"
)

# 3. 设置因子顺序
DEG$group_simple <- factor(DEG$group_simple, levels = c(
  "Highly Up", "Mildly Up", "Lowly Up",
  "Highly Down", "Mildly Down", "Lowly Down"
))

# 4. 设置柱子颜色
deg_colors_labeled <- c(
  "Highly Up (log2FC > 2)" = "#cc0000",    # 深红
  "Mildly Up (1 < log2FC ≤ 2)" = "#ff6666",# 中红
  "Lowly Up (log2(1.5) < log2FC ≤ 1)" = "#ffcccc", # 淡红
  
  "Highly Down (log2FC < -2)" = "#0000cc",     # 深蓝
  "Mildly Down (-2 ≤ log2FC < -1)" = "#6666ff",# 中蓝
  "Lowly Down (-1 ≤ log2FC < -log2(1.5))" = "#ccccff"  # 淡蓝
)

# 5. 重新统计基因数量
count_df <- DEG %>%
  group_by(group_simple) %>%
  summarise(count = n())

# 6. 创建柱状图
p <- ggplot(DEG, aes(x = group_simple, fill = group_label)) +
  geom_bar() +
  geom_text(data = count_df, aes(x = group_simple, y = count, label = count),
            vjust = -0.5, size = 5, inherit.aes = FALSE) +
  scale_fill_manual(values = deg_colors_labeled, name = "DEG Category") +
  theme_minimal() +
  labs(
    title = "Differential Gene Expression Classification",
    x = "Category",
    y = "Gene Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  )

# 7. 输出图片
png("/mnt/f/Project/GSE252984.stat/stat_box.legend_ranges.png", width = 8, height = 6, units = 'in', res = 300)
p
dev.off()

## Part2： 每一类逐一富集分析、画treeplot、并且挑出每个Term中的前五基因
# 1. 确保 gene 是字符型
DEG$gene <- as.character(DEG$gene)

# 2. 按 group 分组提取
high_up <- DEG$gene[DEG$group == "Highly Upregulated"]
mild_up <- DEG$gene[DEG$group == "Mildly Upregulated"]
low_up  <- DEG$gene[DEG$group == "lowly Upregulated"]

high_down <- DEG$gene[DEG$group == "Highly Downregulated"]
mild_down <- DEG$gene[DEG$group == "Mildly Downregulated"]
low_down  <- DEG$gene[DEG$group == "lowly Downregulated"]

# 3. 指定每一个group名字
gene_lists <- list(
  "Highly_Upregulated" = high_up,
  "Mildly_Upregulated" = mild_up,
  "Lowly_Upregulated" = low_up,
  "Highly_Downregulated" = high_down,
  "Mildly_Downregulated" = mild_down,
  "Lowly_Downregulated" = low_down
)

# 4. 创建保存目录（确保路径存在）
dir.create("/mnt/f/Project/GSE252984.stat/pathway_top5_genes", showWarnings = FALSE)
dir.create("/mnt/f/Project/GSE252984.stat/treeplot", showWarnings = FALSE)

# 5. 循环处理每一类（先画treeplot，在挑每一个Term中的前5个基因）
for (name in names(gene_lists)) {
  genes <- gene_lists[[name]]
  
  # 1. GO enrichment
  go <- Enrich(genes, 'GOBP')
  
  # 2. 去掉 regulation 和 somatic recombination 的描述
  go@result %<>% subset(!grepl('regulation|somatic recom', Description))
  
  # 3. 如果结果为空，跳过
  if (nrow(go@result) == 0) {
    cat(paste0("[", name, "] No enriched terms after filtering.\n"))
    next
  }
  
  # 4. 相似性计算
  go.sm = pairwise_termsim(go, method = 'Wang', semData = godata(annoDb = 'org.Mm.eg.db', ont = 'BP'))

  # 5. 画图，添加 title
  p <- treeplot(go.sm, cluster.params = list(label_format_cladelab = 20),
                fontsize = 4, showCategory = 25) +
    ggtitle(paste("GO Treeplot -", name))

  # 6. 保存图片
  ggsave(
    filename = paste0("/mnt/f/Project/GSE252984.stat/treeplot/", name, ".tree.png"),
    plot = p,
    width = 16,
    height = 7
  )

  cat(paste0("[", name, "] Treeplot saved.\n"))
  
  # 7. 提取 treeplot 中使用的 top N term（比如前25个）
  top_terms <- go.sm@result[1: 25, c("Description", "geneID")]
  
  # 8. 拆分 geneID 字符串为基因向量，提取前5个
  top_terms$Top5Genes <- sapply(top_terms$geneID, function(x) {
    genes <- unlist(strsplit(x, "/"))
    paste(genes[1: min(5, length(genes))], collapse = ", ")
  })
  
  # 9. 拆成5列（可选）
  top_genes_df <- data.frame(
    Description = top_terms$Description,
    do.call(rbind, strsplit(top_terms$Top5Genes, ",\\s*"))
  )

  # 10. 保存到一个Excel文件中
  colnames(top_genes_df)[2: 6] <- paste0("Gene", 1: 5)
  write.xlsx(top_genes_df, file = paste0("/mnt/f/Project/GSE252984.stat/pathway_top5_genes/", name, "_top25_terms_top5genes.xlsx"), rowNames = FALSE)
}
