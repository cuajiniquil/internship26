library("ggtree")
library("ape")
library("patchwork")
library("ggplot2")

#=============================================
# TREES & NEWICK FORMAT:
#=============================================
layout(matrix(1:2,2,1,byrow = FALSE))

trextrm <- read.tree(text = "((z:4.1,x:0.5):2, y:2, (w:1,(u:0.9, v:01):1):3);")
plot(trextrm)
edgelabels(trextrm$edge.length, bg = NA, frame = "none", adj= c(0,2))
mtext("((z:4.1,x:2):0.5, y:2, (w:1,(u:0.9, v:01):1):3);", adj = 0)

tr <- read.tree(text = "(((A:5,B:5):2, C:7):4,(D:7,E:7):4);")
plot(tr, show.tip.label = TRUE)
mtext("(((A:5,B:5):2, C:7):4,(D:7,E:7):4);", side = 1, adj = 0)

is.rooted(tr)
is.rooted(trextrm)

rootrex <- root(trextrm, outgroup = c("x","z"), resolve.root = TRUE)
plot(rootrex)
is.rooted(rootrex)

# PLOT HELPER FUNCTION
trplt <- function(tr,txt){
  ggtree(tr, size = 1.5) +   # size = épaisseur des branches
    geom_tiplab(size = 8) +
    geom_nodelab(size = 8, nudge_x = 0.5) +
    annotate("text",
           x = -Inf, y = Inf,
           label = txt,
           size = 8,          # même taille que tes nodelabels
           hjust = -1,         # aligné à droite
           vjust = 2)
}

#=============================================
# POST-ORDER:
#=============================================
tr <- read.tree(text = "(((1:5,2:5):2, 4:7):4,(6:7,7:7):4);")
#plot(tr)

#report: 
tr$node.label <- c("9", "5", "3", "8")
po <- trplt(tr, "B")
ggsave("postordervf.png", po, width = 6, height = 9, dpi = 200)
#nodelabels(c("9","5","3","8"), 6:9, adj = -1, bg = NA, frame = "none", cex = 2) 


#presentation:
ggtree(tr) + geom_text2(aes(label = node))
custom_labels <- data.frame(
  node  = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  labpo = c("1","2","4","6","7","9","5","3","8")  
)

plotpo <- ggtree(tr) %<+% custom_labels +
  geom_text2(aes(label = labpo), hjust = -1, color = "darkgreen", size = 12) +
  geom_point(color = "darkgreen", size = 6)


plotpo

colorpo <- function(nodes, color) {
  list(
  geom_text2(aes(label = labpo, subset = labpo %in% nodes), color = color, hjust = -1, size = 12),
  geom_point2(aes(subset = labpo %in% nodes), color = color, size = 6)
  )
}
  
step1 <- plotpo + colorpo(c(1), "darkred")
s2 <- step1 + colorpo(c(1), "black") + colorpo(c(2), "darkred") 
s3 <- s2 + colorpo(c(2), "black") + colorpo(c(3), "darkred")
s4 <- s3 + colorpo(c(3), "black") + colorpo(c(4), "darkred")
s5 <- s4 + colorpo(c(4), "black") + colorpo(c(5), "darkred")
s6 <- s5 + colorpo(c(5), "black") + colorpo(c(6), "darkred")
s7 <- s6 + colorpo(c(6), "black") + colorpo(c(7), "darkred")
s8 <- s7 + colorpo(c(7), "black") + colorpo(c(8), "darkred")
s9 <- s8 + colorpo(c(8), "black") + colorpo(c(9), "darkred")
s10 <- s9 + colorpo(c(9), "black")
s10

#vibe coded text panel
legend_panel <- ggplot() +
  annotate("text", x = 0, y = 0.8, label = "Postorder traversal",
           size = 5, fontface = "bold", hjust = 0) +
  annotate("point", x = 0, y = 0.5, color = "darkgreen", size = 4) +
  annotate("text", x = 0.1, y = 0.5, label = "unvisited node",
           hjust = 0, size = 4) +
  annotate("point", x = 0, y = 0.3, color = "darkred", size = 4) +
  annotate("text", x = 0.1, y = 0.3, label = "current node",
           hjust = 0, size = 4) +
  annotate("point", x = 0, y = 0.1, color = "black", size = 4) +
  annotate("text", x = 0.1, y = 0.1, label = "visited node",
           hjust = 0, size = 4) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

legend_panel

wrap_plots(c(list(legend_panel), list(
  step1 + ggtitle("A"),
  s2    + ggtitle("B"),
  s3    + ggtitle("C"),
  s4    + ggtitle("D"),
  s5    + ggtitle("E"),
  s6    + ggtitle("F"),
  s7    + ggtitle("G"),
  s8    + ggtitle("H")
)), nrow = 3)

#=============================================
# PRE-ORDER:
#=============================================

#report: 

tr2 <- read.tree(text = "(((4:5,5:5):2, 6:7):4,(8:7,9:7):4);")
#plot(tr2)
#nodelabels(c("1","2","3","7"), 6:9, adj = -1, bg = NA, frame = "none") 

tr2$node.label <- c("1", "2", "3", "7")
pe <- trplt(tr2,"A")
ggsave("preordervf.png", pe, width = 6, height = 9, dpi = 200)

custom_labels <- data.frame(
  node  = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  labpo = c("4","5","6","8","9","1","2","3","7")  
)

#presentation: 

plotpre <- ggtree(tr2) %<+% custom_labels +
  geom_text2(aes(label = labpo), hjust = -1, size = 6) +
  geom_point(size = 3)

plotpre


#=============================================
# LEVEL-ORDER:
#=============================================

#report: 

tr3 <- read.tree(text = "(((8:5,9:5):2, 5:7):4,(6:7,7:7):4);")
#plot(tr3)
#nodelabels(c("1","2","4","3"), 6:9, adj = -1, bg = NA, frame = "none") 

tr3$node.label <- c("1", "2", "4", "3")
lv <- trplt(tr3,"C")
ggsave("levelordervf.png", lv, width = 6, height = 9, dpi = 200)

custom_labels <- data.frame(
  node  = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  labpo = c("8","9","5","6","7","1","2","4","3")  
)

#presentation: 

plotlv <- ggtree(tr3) %<+% custom_labels +
  geom_text2(aes(label = labpo), hjust = -1, size = 6) +
  geom_point(size = 3)

plotlv

#=============================================
# ULTRAMETRICITY:
#=============================================
set.seed(273)
coal <- ggtree(rcoal(5), size = 1.5) +   # size = épaisseur des branches
  geom_tiplab(size = 8)
het <- ggtree(rtree(5), size = 1.5) +   # size = épaisseur des branches
  geom_tiplab(size = 8) 

ggsave("coalescentr.png", coal, width = 6, height = 9, dpi = 200)
ggsave("heterochronictr.png", het, width = 6, height = 9, dpi = 200)

