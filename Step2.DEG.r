# 差异分析
raw = read.table('Raw.counts.txt', sep = '\t', header = T, check.names = F, row.names = 1)
pos = c('D1', 'D2')
neg = c('Ctrl1', 'Ctrl2')
deg = RNAseq.DESeq2(raw, pos, neg)
write.table(deg, '/mnt/f/Project/GSE244787.stat/DEG.DvsCtrl.txt', sep = '\t', row.names = F)

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
