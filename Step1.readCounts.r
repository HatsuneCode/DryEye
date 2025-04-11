## Step1: read counts ##
expr = read.table('/mnt/f/Project/GSE252984.stat/4.expected_count.xls', sep = '\t', header = T, row.names = 1, check.names = F)
expr = RNAseq.checkDupRow(expr)
expr = round(expr)
write.table(cbind(Gene = rownames(expr), expr), '/mnt/f/Project/GSE252984.stat/Raw.counts.txt', sep = '\t', row.names = F)
expr = RNAseq.Normalize(expr)
write.table(cbind(Gene = rownames(expr), expr), '/mnt/f/Project/GSE252984.stat/Normalized.counts.txt', sep = '\t', row.names = F)

## PCA
expr.f = expr[order(apply(expr, 1, mad), decreasing = T)[1:3e3], ]
pca = PCA(expr.f)
pca$group = sub('.$', '', pca$sample)
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 4, max.overlaps = 50) +
  labs(title = "PCA of Top3k Variable Genes",
       x = 'PC1', y = 'PC2', color = 'Group') +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    axis.title = element_text(size = 16), axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave('/mnt/f/Project/GSE252984.stat/PCA.Top3k.png', width = 9, height = 7, dpi = 300)

# 差异分析
raw = read.table('/mnt/f/Project/GSE252984.stat/Raw.counts.txt', sep = '\t', header = T, check.names = F, row.names = 1)
pos = c('BAC1', 'BAC2', 'BAC3', 'BAC4')
neg = c('PBS1', 'PBS2', 'PBS3', 'PBS4')
DEG = RNAseq.DESeq2(raw, pos, neg)
write.table(DEG, '/mnt/f/Project/GSE252984.stat/DEG.DvsCtrl.txt', sep = '\t', row.names = F)

# heatmap of genes number
pref = '/mnt/f/Project/GSE252984.stat/stat.DEG'
dir.create(pref, F)
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
for (adj_val in c(TRUE, FALSE)) for (pv_val in c(0.05, 0.01)) {
    pd = DEGsig[DEGsig$ADJ == adj_val & DEGsig$PV == pv_val, c('Name', 'Up', 'Down')]
    if (!nrow(pd)) next  # 如果这一组没有数据，就跳过
    rownames(pd) = sub(' .*', '', pd$Name)
    pd = t(pd[rev(seq(nrow(pd))), -1])
    pd[2, ] = -pd[2, ]  # Down 设为负数方便绘图显示
    # 创建热图
    lim = max(abs(min(pd)), max(pd))
    p = Heatmap(pd,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                name = 'DEG.n',
                col = colorRamp2(c(-lim, 0, lim), c('blue', 'white', 'red')),
                cell_fun = function(j, i, x, y, w, h, col) grid.text(abs(pd[i, j]), x, y) )
    # 构建文件名
    fname = paste0('DEG.adj', tolower(adj_val), '.p', format(pv_val, nsmall = 2), '.png')
    fpath = file.path(pref, fname)
    message('Saving heatmap: ', fpath)
    # 保存图片
    png(fpath, width = 8, height = 4.5, units = 'in', res = 300)
    draw(p, 
         column_title = paste0('DEG Heatmap (adj=', adj_val, ', p=', pv_val, ')'),
         column_title_gp = gpar(fontsize = 20, fontface = 'bold'))
    dev.off()
}
