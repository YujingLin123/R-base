library(VennDiagram)
library(RColorBrewer)

#setwd("~/DEGoverlap")

venn.plot <- draw.pairwise.venn(
area1 = 2336, #第一个集大小
area2 = 5130, #第二个集大小
cross.area = 1558, #两个集的交集大小
category = c("NOR", "DML2"), #两个集命名
fill = c("#FB8072","#FDB462"), #集对应的圈填充颜色"blue", "yellow"
lty = "blank", #圆周的线条类型
col = c("white"),
cex = 2, #韦恩图各部分面积标签注释字体大小
cat.cex = 2, #集名字体大小
cat.pos = c(0, 0), #集标签围绕圆的相对位置，0-360°，默认0°在12点钟方向
cat.dist = c(0.03,0.06), #集标签相对于圆位置远近
cat.just = list(c(0, 0), c(0, 0)),
ext.pos = 0, #圈外面积标签位置
ext.dist = -0.05,
ext.length = 0.85, #圈外面积标签连接线长度
ext.line.lwd = 2,
ext.line.lty = "dashed", #圈外面积标签连接线类型
alpha=0.4,
euler.d=T #没有交集，是否分开
);

tiff(filename = "Venn_diagram.tiff");
grid.draw(venn.plot);
dev.off()
