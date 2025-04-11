RNAseq.checkDupRow = function(expr, method = 'mean') {
  ## process symbols
  gene  = sub('^ENS.*?_', '', rownames(expr))
  dgene = unique(gene[duplicated(gene)])
  ## check duplicated symbols
  if (length(dgene)) {
    expr1 = expr[!gene %in% dgene,]
    expr2 = do.call(rbind, lapply(dgene, function(g) {
      e = expr[gene == g, , drop = F]
      if (method == 'mean')
        t(setNames(data.frame(colMeans(e)), g))
    }))
    expr = rbind(expr1, expr2)
    rm(expr1, expr2)
  }
  rm(gene, dgene)
  ## restore symbol
  rownames(expr) = sub('^ENS.*?_', '', rownames(expr))
  expr
}

RNAseq.Normalize = function(expr, log2 = T, method = 'DESeq2') {
  if (method == 'DESeq2') {
    suppressMessages(library(DESeq2))
    dds  = DESeqDataSetFromMatrix(
      countData = expr,
      colData = data.frame(row.names = colnames(expr), samples = factor(colnames(expr))),
      design = ~ samples)
    dds  = estimateSizeFactors(dds)
    expr = counts(dds, normalize = T)
  }
  if (log2) expr = log2(expr + 1)
  expr
}

expr = read.table('/mnt/f/Project/GSE252984.stat/4.expected_count.xls', sep = '\t', header = T, row.names = 1, check.names = F)
expr = RNAseq.checkDupRow(expr)
expr = round(expr)
write.table(cbind(Gene = rownames(expr), expr), '/mnt/f/Project/GSE252984.stat/Raw.counts.txt', sep = '\t', row.names = F)
expr = RNAseq.Normalize(expr)
write.table(cbind(Gene = rownames(expr), expr), '/mnt/f/Project/GSE252984.stat/Normalized.counts.txt', sep = '\t', row.names = F)

# PCA part
PCA = function(expr, ...) {
  pca = data.frame(prcomp(t(expr), ...)$x, check.rows = F)
  pca$sample = rownames(pca)
  pca$sample = factor(pca$sample, pca$sample)
  pca
}

expr.f = expr[order(apply(expr, 1, mad), decreasing = T)[1:3e3],]
pca = PCA(expr.f)
pca$group = sub('.$', '', pca$sample)

ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 4, max.overlaps = 50) +
  labs(title = "PCA of Top 3000 Most Variable Genes",
       x = "PC1", y = "PC2", color = "Group") +
  scale_color_brewer(palette = "Set1") +
  theme_classic(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave("/mnt/f/Project/GSE252984.stat/PCA.Top3000.png", width = 9, height = 7, dpi = 300)

expr.f = expr[order(apply(expr, 1, mad), decreasing = T),]
pca = PCA(expr.f)
pca$group = sub('.$', '', pca$sample)

ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 4, max.overlaps = 50) +
  labs(title = "PCA of All Variable Genes",
       x = "PC1", y = "PC2", color = "Group") +
  scale_color_brewer(palette = "Set1") +
  theme_classic(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave("/mnt/f/Project/GSE252984.stat/PCA.All.png", width = 9, height = 7, dpi = 300)


# 差异分析
RNAseq.DESeq2 = function(expr, pos = NULL, neg = NULL, name = NULL, exp_cut = 10, cut.method = 'inter') {
  ## cut.method: inter, union
  suppressMessages(library(DESeq2))
  ## expr split
  exprP = if (length(pos)) expr[, colnames(expr) %in% pos, drop = F] else 
    expr[, !colnames(expr) %in% neg, drop = F]
  exprN = if (length(neg)) expr[, colnames(expr) %in% neg, drop = F] else 
    expr[, !colnames(expr) %in% pos, drop = F]
  ## expr type
  type = paste(paste(colnames(exprP), collapse = ','), 'vs', paste(colnames(exprN), collapse = ',') )
  if (!length(name)) name = type
  message('DEG: ', type)
  ## condition control ~ treatment
  condition = factor( c(rep('Neg', ncol(exprN)), rep('Pos', ncol(exprP))), c('Neg', 'Pos') )
  ## counts
  expr  = cbind(exprN, exprP)
  expr  = expr[rowSums(expr) > 0, , drop = F]
  ## cut-off
  if (length(exp_cut)) {
    gene = lapply(c('Pos', 'Neg'), function(g) {
      expr.f = expr[, condition == g, drop = F]
      rownames(expr.f)[apply(expr.f, 1, function(i) !any(i < exp_cut) )]
    } )
    if ('inter' %in% cut.method) gene = Reduce(intersect, gene)
    if ('union' %in% cut.method) gene = Reduce(union, gene)
    message('--> valid gene number: ',  length(gene), ' <--')
    expr = expr[gene, , drop = F]
  }
  ## split pos/neg
  exprP = expr[, condition == 'Pos', drop = F]
  exprN = expr[, condition == 'Neg', drop = F]
  ## meta
  meta = data.frame(row.names = colnames(expr), condition)
  ## DESeq2
  dds = DESeqDataSetFromMatrix(countData = expr, colData = meta, design = ~ condition)
  dds = DESeq(dds)
  dds = data.frame(results(dds), check.names = F)
  ## output
  data.frame(p_val = dds$pvalue, avg_log2FC = dds$log2FoldChange, 
             pct.1 = apply(exprP, 1, function(i) sum(i > 0)/ncol(exprP) ),
             pct.2 = apply(exprN, 1, function(i) sum(i > 0)/ncol(exprN) ),
             p_val_adj = dds$padj, gene = rownames(dds), 
             average = rowMeans(expr),  median = apply(expr,  1, median), 
             posAvg  = rowMeans(exprP), posMed = apply(exprP, 1, median),
             negAvg  = rowMeans(exprN), negMed = apply(exprN, 1, median),
             type    = name, 
             upDown  = factor(ifelse(dds$log2FoldChange > 0, 'Up', 'Down'), c('Up', 'Down')), 
             row.names = NULL )
}

raw = read.table('/mnt/f/Project/GSE252984.stat/Raw.counts.txt', sep = '\t', header = T, check.names = F, row.names = 1)
pos = c('BAC1', 'BAC2', 'BAC3', 'BAC4')
neg = c('PBS1', 'PBS2', 'PBS3', 'PBS4')
deg = RNAseq.DESeq2(raw, pos, neg)
write.table(deg, '/mnt/f/Project/GSE252984.stat/DEG.DvsCtrl.txt', sep = '\t', row.names = F)

# heatmap of genes number
pref = '/mnt/f/Project/GSE252984.stat/stat.DEG'
dir.create(pref, F)
para.DEG = list(
  ## p .01
  'FC1.5 Padj.01'  = c(1.5, .01),       'FC1.5 P.01'  = c(1.5, .01), 
  'FC1 Padj.01'    = c(1, .01),         'FC1 P.01'    = c(1, .01), 
  'FC0.58 Padj.01' = c(log2(1.5), .01), 'FC0.58 P.01' = c(log2(1.5), .01),
  'FC0.26 Padj.01' = c(log2(1.2), .01), 'FC0.26 P.01' = c(log2(1.2), .01),
  'FC0 Padj.01'    = c(0, .01),         'FC0 P.01'    = c(0, .01),
  ## p .05
  'FC1.5 Padj.05'  = c(1.5, .05),       'FC1.5 P.05'  = c(1.5, .05),
  'FC1 Padj.05'    = c(1, .05),         'FC1 P.05'    = c(1, .05), 
  'FC0.58 Padj.05' = c(log2(1.5), .05), 'FC0.58 P.05' = c(log2(1.5), .05),
  'FC0.26 Padj.05' = c(log2(1.2), .05), 'FC0.26 P.05' = c(log2(1.2), .05),
  'FC0 Padj.05'    = c(0, .05),         'FC0 P.05'    = c(0, .05)
)

DEG = read.table('/mnt/f/Project/GSE252984.stat/DEG.DvsCtrl.txt', sep = '\t', header = T, check.names = F)
# 计算 DEG 数量
DEGsig = do.call(rbind, lapply(seq(para.DEG), function(i) {
  fc = para.DEG[[i]][1]
  pv = para.DEG[[i]][2]
  nm = names(para.DEG)[i]
  adj = grepl('Padj', nm)
  deg = if (adj) DEG[which(DEG$p_val_adj < pv), ] else DEG[which(DEG$p_val < pv), ]
  up  = deg[deg$avg_log2FC >  fc, ]
  dn  = deg[deg$avg_log2FC < -fc, ]
  data.frame(FC = fc, PV = pv, Up = nrow(up), Down = nrow(dn), ADJ = adj, Name = nm)
}))

# 循环绘图：4种组合
for (adj_val in c(TRUE, FALSE)) {
  for (pv_val in c(0.05, 0.01)) {
    
    pd = DEGsig[DEGsig$ADJ == adj_val & DEGsig$PV == pv_val, c("Name", "Up", "Down")]
    if (nrow(pd) == 0) next  # 如果这一组没有数据，就跳过
    
    print(pd$Name)
    rownames(pd) = sub(' .*', '', pd$Name)
    print(pd$Name)
    pd = t(pd[rev(seq(nrow(pd))),-1])
    pd[2, ] = -pd[2, ]  # Down 设为负数方便绘图显示
    
    # 创建热图
    p = Heatmap(pd,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                name = 'DEG.n',
                col = colorRamp2(c(min(pd), 0, max(pd)), c('blue', 'white', 'red')),
                cell_fun = function(j, i, x, y, w, h, col) {
                  grid.text(abs(pd[i, j]), x, y)
                })
    
    # 构建文件名
    fname = paste0("DEG.adj", tolower(adj_val), ".p", format(pv_val, nsmall = 2), ".png")
    fpath = file.path(pref, fname)
    message("Saving heatmap: ", fpath)
    
    # 保存图片
    png(fpath, width = 8, height = 4.5, units = "in", res = 300)
    draw(p, 
         column_title = paste0("DEG Heatmap (adj=", adj_val, ", p=", pv_val, ")"),
         column_title_gp = gpar(fontsize = 20, fontface = "bold", fontfamily = "serif"))
    dev.off()
  }
}
