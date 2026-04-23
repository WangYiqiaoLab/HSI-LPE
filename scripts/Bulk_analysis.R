
## Introduction

There are two part result in the report. The first part is to see if LPE would affect the proliferation and differentiation in neural organoid. The second part is to see how LPE affect the cellular response of primary microglia with the presence of Amyloid-beta peptide.

```{r Load librires}
#| echo: false
#| include: false

library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)

library(gprofiler2)
library(correlationAnalyzeR)
library(clusterProfiler)
#library(xlsx)

library(plotly)
library(DOSE)
library(EnhancedVolcano)
library(biomaRt)
library(htmltools)
library(DT)
library(ReactomePA)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GSVA)
library(UniProt.ws)
#GOBP_list <- getTERM2GENE(GSEA_Type = "GO:BP", Species = "mmusculus") 
#KEGG_list <- getTERM2GENE(GSEA_Type = "KEGG", Species = "mmusculus")


```

## Neurosophere

```{r Load the data organoid}
#| include: false
#| echo: false
#metadata
meta <- read.csv("Bulk/studydesign.csv",sep = ",")

meta <- meta[1:8,]
# data
DATA <-read.delim("../gene_count.xls", sep = "\t", header=TRUE,row.names = "gene_id")
DATA_fpkm <- read.delim("../gene_fpkm.xls", sep = "\t", header=TRUE,row.names = "gene_id")
DATA.org <- read.delim("../gene_count.xls", sep = "\t", header=TRUE, row.names = "gene_id")
DATA.org$gene_id <- rownames(DATA.org)


DATA <- DATA[, meta$Sample]
#rownames(DATA) <- rownames(DATA.org)
```

```{r dds obj organoid, echo=FALSE}
#| include: false
#| echo: false
data <- DESeqDataSetFromMatrix(countData = DATA,
                              colData = meta,
                              design = ~ Group)
```

### PCA

```{r boxplot organoid }
#| column: screen-inset-shaded
#| layout-nrow: 1

par(cex.axis=0.8, mar=c(10, 3, 3, 2))
boxplot( log2(counts(data)+1), las = 2,
  main = "Original Log2 Counts")

#par(cex.axis=0.8, mar=c(10, 3, 3, 2))
#boxplot( assay(rlog(dds)), las = 2,
#  main = "DESeq2 (rlog) Normailized Counts")

vst <- vst(data, blind = F)

par(cex.axis=0.8, mar=c(10, 3, 3, 2))
boxplot( assay(vst), las = 2,
  main = "VST Normailized Counts")

par(cex.axis=0.8, mar=c(10, 3, 3, 2))
boxplot( log2(DATA_fpkm+1), las = 2,
  main = "Log FPKM Counts")

```

```{r PCA organoid, echo=FALSE}
#raw counts
DATA.pca <- prcomp(t(log(DATA+1)), center = T)
#DATA.pca <- prcomp(t(DATA), center = T)
DATA.pca.sum <- as.data.frame(summary(DATA.pca)$importance)
DATA.pca <- as.data.frame(DATA.pca$x)
DATA.pca$Group <-meta$Group
#DATA.pca$Group <- factor(DATA.pca$Group, levels = unique(c("Ctrl", "6_Hours","1_Day", "3_Days")))
DATA.pca$Sample <- rownames(DATA.pca)

pca_raw <-   ggplot(as.data.frame(DATA.pca),
  aes(x = PC1, y = PC2, colour = Group)
  ) +
  geom_point(shape=19, size=4, alpha = 0.7)+
  geom_text(aes(label = Sample), position = position_jitter(width =2, height = 4)) +
  geom_hline(yintercept = 0, colour = "gray65", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray65", linetype = "dashed") +
  ggtitle("PCA on raw count") +
  xlab(paste0("PC1 (",DATA.pca.sum$PC1[2]*100,"%)"))+
  ylab(paste0("PC2 (",DATA.pca.sum$PC2[2]*100,"%)"))+
  #ylab("PC3 2.1%")+
  theme_classic()

#vst counts
DATA.pca <- prcomp(t(assay(vst)), center = T)
#DATA.pca <- prcomp(t(DATA), center = T)
DATA.pca.sum <- as.data.frame(summary(DATA.pca)$importance)
DATA.pca <- as.data.frame(DATA.pca$x)
DATA.pca$Group <-meta$Group
#DATA.pca$Group <- factor(DATA.pca$Group, levels = unique(c("Ctrl", "6_Hours","1_Day", "3_Days")))
DATA.pca$Sample <- rownames(DATA.pca)

pca_vst <-   ggplot(as.data.frame(DATA.pca),
  aes(x = PC1, y = PC2, colour = Group)
  ) +
  geom_point(shape=16, size=4)+
  geom_text(aes(label = Sample), position = position_jitter(width =2, height = 4)) +
  geom_hline(yintercept = 0, colour = "gray65", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray65", linetype = "dashed") +
  ggtitle("PCA on vst normalized count") +
  xlab(paste0("PC1 (",DATA.pca.sum$PC1[2]*100,"%)"))+
  ylab(paste0("PC2 (",DATA.pca.sum$PC2[2]*100,"%)"))+
  #ylab("PC3 2.1%")+
  theme_classic()


```

```{r PCAplot organoid }
#| column: screen-inset-shaded
#| layout-nrow: 1

pca_raw
pdf( "pca_vst.pdf", width = 8, height = 5)
pca_vst
dev.off()
pca_fpkm

rm(pca_raw,pca_vst,pca_fpkm,DATA.pca,DATA.pca.sum)
```

```{r Set the Condition, echo=FALSE}
Con1 = "LPE"


Con2 = "Control"

```

### DGE

```{r make a DESeq object organoid}
#| echo: false
#| include: false

data <- DESeqDataSetFromMatrix(countData = DATA,
                              colData = meta,
                              design = ~ Group )

data <- data[rowSums(counts(data)>=10) >=3,]

data <- DESeq(data)
res <- results(data, contrast=c("Group",Con1,Con2))

#res<-lfcShrink(data, coef =  "Group_LPE_vs_Control")


res$gene_name<-DATA.org[match(sapply(strsplit(rownames(res),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]
res$gene_des<-DATA.org[match(sapply(strsplit(rownames(res),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_description"]
res<-as.data.frame(res)

res$gene_name[is.na(res$gene_name)] <- rownames(res)[is.na(res$gene_name)]


```


```{r get top gene names}

diffexpressed <- ifelse(
    res$log2FoldChange > 1 & res$threshold_OE == TRUE, "#B2182B",
    ifelse(
      res$log2FoldChange < 1 & res$threshold_OE == TRUE, "#2166AC",
        'grey'))
diffexpressed[is.na(diffexpressed)] <- "grey"
names(diffexpressed)[diffexpressed =="#B2182B"] <- "UP"
names(diffexpressed)[diffexpressed =="#2166AC"] <- "DOWN"
names(diffexpressed)[diffexpressed == "grey"] <- "no Sig"


## --- Top 10 UP ---
up <- res[res$diffexpressed == "UP", ]                   # keep only upregulated
top10_up <- up[order(up$padj), ]$gene_name[1:10]           # 10 smallest padj

## --- Top 10 DOWN ---
down <- res[res$diffexpressed == "DOWN", ]                        # keep only downregulated
top10_down <- down[order(down$padj), ]$gene_name[1:10]
```

```{r table deg neural}
deg <- res[res$threshold_OE %in% TRUE,]
deg <- deg[,c(7, 1:10)]
browsable(
  tagList(
    tags$h4("DEGs"),
    datatable(deg)
))
```

```{r}
write.csv(res, "DGE.csv")
```

##### Customized Volcano

```{r Customized Volcano neural}
#| echo: false


vol <- EnhancedVolcano(res, 
                x ="log2FoldChange", 
                y = "padj" , 
                subtitle = paste0( Con1, " vs ", Con2, "   P_adj < 0.05 & abs(LogFC) >= 1"),
                selectLab = c("Wnt4",  "Sorl1"),
                #selectLab =c(c("Sv2a","Nrcam","Nomo1","Grm5","Itgb8","Erbb2","Met","Atp2b1")),
                #selectLab =data_pro$gene_name,
                lab = res$gene_name, 
                pCutoff = 0.05,
                colAlpha = 0.6,
                pointSize = 2.5 ,
                ylab = bquote(~-Log[10] ~ italic(Padj)),
                drawConnectors = T,
                widthConnectors = 0.4,
                max.overlaps = 80,
                colConnectors = 'black',
                boxedLabels = T,
                labSize = 2,
                shape = 16,
                colCustom = diffexpressed
                #col = c("grey","grey","#1D5B2D")
                ) + 
  theme_classic() +
  theme(legend.title=element_blank())

ggsave("Bulk/volcano_MS_2.pdf", plot = vol, device = "pdf", useDingbats = FALSE, width = 6.5, height = 4.5)
```

### Enrichment Analysis

Hover the mouse over the dot to see the term

```{r LPE GO organoid}
#| echo: false


go_up = gost(res$gene_name[res$diffexpressed =="UP"], organism = 'mmusculus', evcodes = T)

#write.csv(go_up$result[,c(11,3:11,16)],paste0("",Con1,"_vs_",Con2,"_up_fo.csv"))

gostplot(go_up,interactive = T)%>% layout(title = paste0("Enriched terms in LPE"))


```

```{r LPE terms organoid}
#| echo: false


browsable(
  tagList(
    tags$h4("Enrichment Analysis on LPE Upregulated Genes"),
    datatable(go_up$result[,c(11, 3:11,16)])
))
```

#### Reactome

```{r reactome analysis}
#| include: false
#| echo: false
up <- clusterProfiler::bitr(geneID = res[res$diffexpressed =="UP",]$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db" ) 
x <- enrichPathway(gene= up$ENTREZID, pvalueCutoff = 0.05, organism = "mouse", readable=TRUE)

```

```{r}
df <- x@result


terms <- c("GABA receptor activation","Neurotransmitter receptors and postsynaptic signal transmission","MET activates PTK2 signaling","Transport of connexons to the plasma membrane")


df <- as.data.frame(x@result)
df <- df[df$Description %in% terms , ] %>% mutate(qscore = -log(p.adjust, base = 10))

df$GeneRatio <-df$Count/255
pdf("reactome_bar.pdf",width = 7.5,height = 3)
ggplot(df, aes(x = GeneRatio, y = forcats::fct_reorder(Description, GeneRatio, .desc = FALSE), fill= qscore) ) +
  scale_fill_continuous(
    low = "#68a8cd", # Starting color for less significant (higher p.adjust)
    high = "#f16569",   # Ending color for highly significant (lower p.adjust)
    guide = guide_colorbar(reverse = FALSE) # Reverse the color key for p.adjust
  ) +
  geom_bar(stat = "identity") +
  labs(y = NULL, title = "Reactome" ) + # Set hjust = 0.5 to center the plot title
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
  axis.text.y = element_text(size = 12) )
dev.off()
```

```{r}




pdf("reactome_bar.pdf", width = 5, height = 3)
barplot(x, x = "GeneRatio", showCategory = c("GABA receptor activation","Neurotransmitter receptors and postsynaptic signal transmission","MET activates PTK2 signaling","Transport of connexons to the plasma membrane"))
dev.off()
viewPathway( c("GABA receptor activation"), organism = "mouse", readable = T)
```

```{r Reactome with foldchange}
genelist <- res[,c("log2FoldChange","gene_name")]
genelist <- distinct(genelist, gene_name, .keep_all = T)
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
entrez <- mapIds(
  org.Mm.eg.db,
  keys      =rownames(res),
  keytype   = "ENSEMBL",
  column    = "ENTREZID",
  multiVals = "first"     # <- force 1:1
)
gene_list <- res$log2FoldChange

names(gene_list) <- entrez
gene_list <-gene_list[!duplicated(names(gene_list))]

pdf("reactome_pathway.pdf", width = 8, height = 6)
viewPathway( c("GABA receptor activation"), organism = "mouse", readable = T, foldChange = gene_list) +ggtitle("GABA receptor activation")+scale_color_gradient(low = "blue", high = "red", name = "log2 fold change")+theme(plot.title = element_text(hjust = 0.5))



viewPathway( c("Neurotransmitter receptors and postsynaptic signal transmission"), organism = "mouse", readable = T, foldChange = gene_list) +ggtitle("Neurotransmitter receptors and postsynaptic signal transmission")+scale_color_gradient(low = "blue", high = "red", name = "log2 fold change")+theme(plot.title = element_text(hjust = 0.5))
viewPathway( c("MET activates PTK2 signaling"), organism = "mouse", readable = T, foldChange = gene_list) +ggtitle("MET activates PTK2 signaling")+scale_color_gradient(low = "blue", high = "red", name = "log2 fold change")+theme(plot.title = element_text(hjust = 0.5))
dev.off()
```

```{r LPE reactome organoid}
#| echo: false


browsable(
  tagList(
    tags$h4("reactome Analysis on LPE Upregulated Genes"),
    datatable(x@result )
))
```

```{r Control GO organoid}
#| echo: false
go = gost(res$gene_name[res$diffexpressed =="DOWN"], organism = 'mmusculus')



gostplot(go,interactive = T)%>% layout(title = paste0("Enriched terms in Down-regulated Genes"))
```

```{r Control term organoid}
#| echo: false

browsable(
  tagList(
    tags$h4("Enrichment Analysis on LPE Down regulated Genes"),
    datatable(go$result[,c(11, 3:11)])
))
```

```{r interesting terms plot}
#| fig-width: 6
#| fig-height: 8
clusterprofiler_result <- go_up$result %>%
  mutate(
    ID = term_id,
    Description = term_name,
    pvalue = p_value,
    p.adjust = p_value,  # Adjust p-values
    geneID = sapply(intersection, function(genes) paste(genes, collapse = ",")), # Convert to comma-separated
    Count = intersection_size# Count of genes
    #gene = res.dds[res.dds$diffexpressed == "UP", ]
  ) %>% dplyr::select(ID, Description, pvalue, p.adjust, geneID, Count)#,gene)

enrich_result <- new(
  "enrichResult",
  result = clusterprofiler_result,
  gene = res[res$diffexpressed == "UP", ]$gene_name,
  ontology = "GO")

enrich_result@result$GeneRatio <- paste(enrich_result@result$Count, length(res[res$diffexpressed =="UP",]$gene_name), sep = "/")


enrich_result@result$`-log10(p.adjust)` <- -log10(enrich_result@result$p.adjust)

terms <- c("neurogenesis", "generation of neurons", "neuron differentiation", "neuron development", "axon ensheathment", "neuron projection development", "cell morphogenesis", "gliogenesis",  "oligodendrocyte differentiation", "brain development",  "positive regulation of nervous system development", "synapse assembly", "cell morphogenesis involved in neuron differentiation", "regulation of glial cell differentiation", "neurotransmitter transport", "axonogenesis", "gamma-aminobutyric acid signaling pathway", "neural precursor cell proliferation", "neuron maturation")
enrich_result@result$GeneRatio <-enrich_result@result$Count/length(length(res[res$diffexpressed =="UP",]$gene_name))

df <- as.data.frame(enrich_result)
df <- df[df$Description %in% terms , ] %>% mutate(qscore = -log(p.adjust, base = 10))

pdf("UP_GOBP.pdf",width = 6.5,height = 6.5)
ggplot(df, aes(x = GeneRatio, y = forcats::fct_reorder(Description, GeneRatio, .desc = FALSE), fill= qscore) ) +
  scale_fill_continuous(
    low = "#68a8cd", # Starting color for less significant (higher p.adjust)
    high = "#f16569",   # Ending color for highly significant (lower p.adjust)
    guide = guide_colorbar(reverse = FALSE) # Reverse the color key for p.adjust
  ) +
  geom_bar(stat = "identity") +
  labs(y = NULL, title = "Enriched terms from LPE" ) + # Set hjust = 0.5 to center the plot title
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
  axis.text.y = element_text(size = 12) )
dev.off()

```

### Gene plot

```{r wnt pathway}


data_plot  <- t(vst(data, blind = F) %>% assay()) %>% as.data.frame()
data_plot$sample <- rownames(data_plot)
data_plot$group <- c(rep(c("Control", "LPE"), each = 4))


gene_list <-  c(
  # Wnt ligands
  "Wnt1","Wnt3a","Wnt4","Wnt5a","Wnt7a","Wnt7b",
  # Wnt receptors & co-receptors
  "Fzd3","Fzd5","Fzd7","Lrp5","Lrp6",
  # Wnt intracellular mediators
  "Dvl1","Dvl2","Ctnnb1","Axin2","Gsk3b","Apc",
  # Wnt transcriptional regulators
  "Tcf7l2","Lef1",
  # Wnt antagonists/modulators
  "Dkk1","Sfrp1","Sfrp2",
  # BMP ligands
  "Bmp2","Bmp4","Bmp7","Bmp7",
  # BMP receptors
  "Bmpr1a","Bmpr1b","Bmpr2",
  # BMP Smads
  "Smad1","Smad5","Smad9","Smad4",
  # BMP antagonists
  "Nog","Grem1","Chrd",
  # Crosstalk regulators
  "Rnf146","Pitx2","Ccnd1","Smurf1","Smurf2"
)%>% as.data.frame()

gene_list <- res[res$threshold_OE == "TRUE", ]$gene_name %>% as.data.frame() %>% rev()

gene_list$gene_id <-DATA.org[match(sapply(strsplit(gene_list$.,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

data_plot2 <- data_plot[,colnames(data_plot) %in% gene_list$gene_id] %>% t() %>% as.data.frame()
data_plot2$gene_name<-DATA.org[match(sapply(strsplit(rownames(data_plot2),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]
data_plot2$gene_name[is.na(data_plot2$gene_name)] <- rownames(data_plot2)[is.na(data_plot2$gene_name)]
rownames(data_plot2) <- data_plot2$gene_name
data_plot2$gene_name <- NULL
data_plot2 <- t(data_plot2 %>% as.data.frame()) %>% as.data.frame()


```

```{r heatmap wnt}
pdf("wnt.pdf", width= 5, height = 8)
Heatmap(t(data_plot2) %>% pheatmap:::scale_rows(),  name = "Relative expression",
        cluster_columns = F,
        cluster_rows = T,
        #top_annotation = Annotation,
        column_title = "Wnt Signalling from ChatGPT",
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right"
)
dev.off()
```

# Gene Correlation with proteomics

```{r}
data_pro <- readxl::read_excel("Bulk/Ms_Data.xlsx")
up <- UniProt.ws(taxId = 9606)
#keytypes(up)   # includes "UNIPROTKB"
#columns(up)    # look for "PROTEIN-NAMES", "GENES", etc.

res <- select(up, keys = data_pro$Accession, columns = c("UniProtKB", "gene_names", "protein_name"), keytype = "UniProtKB")
data_pro$gene_name <- res[match(data_pro$Accession, res$From),"Gene.Names"]
data_pro$gene_name_org <- data_pro$gene_name
data_pro$gene_name <- sapply(strsplit(data_pro$gene_name, " "), `[`, 1)
data_pro$gene_name2 <- NULL
write.csv(data_pro, "Bulk/MS_Data.csv")

data_plot  <- t(vst(data, blind = F) %>% assay()) %>% as.data.frame()
#data_plot <- data_plot[,order(colSums(data_plot), decreasing = T)[1:1000]]

data_plot$sample <- rownames(data_plot)
data_plot$group <- c(rep(c("Control", "LPE"), each = 4))




#gene_list <-  c("Sv2a","Nrcam","Nomo1","Grm5","Itgb8","Erbb2","Met","Atp2b1")%>% as.data.frame()
gene_list <-  data_pro$gene_name %>% as.data.frame()



#gene_list <- res[res$threshold_OE == "TRUE", ]$gene_name %>% as.data.frame() %>% rev()

gene_list$gene_id <-DATA.org[match(sapply(strsplit(gene_list$.,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

data_plot2 <- data_plot[,colnames(data_plot) %in% gene_list$gene_id] %>% t() %>% as.data.frame()
data_plot2$gene_name<-DATA.org[match(sapply(strsplit(rownames(data_plot2),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]
data_plot2$gene_name[is.na(data_plot2$gene_name)] <- rownames(data_plot2)[is.na(data_plot2$gene_name)]
rownames(data_plot2) <- data_plot2$gene_name
data_plot2$gene_name <- NULL
data_plot2 <- t(data_plot2 %>% as.data.frame()) %>% as.data.frame()

data_plot2$Condition <- rep(c("Control", "LPE"), each = 4)
data_plot2$sample <- rownames(data_plot2)
data_plot_long <- data_plot2 %>% pivot_longer(cols = -c("sample", "Condition"), names_to = "gene",values_to = "expression")

mycolors <- c("#317EC2", "#C03830")
  names(mycolors) <- c("Control", "LPE")


  links <- data.frame(
  source = c(
    "Akt1", "Akt1", "Akt1", "Akt1",
    "Erbb2", "Erbb2", "Erbb2", "Erbb2",
    "Met", "Met", "Src", "Src", "Src",
    "Grm5", "Grm5", "Stat3", "Mapk1", "Mapk1", 
    "Mapk1", "Mapk1", "Atp2b1", "Mcu", "Itgb8", "Strap"
  ),
  target = c(
    "Bdnf", "Hey1", "Hey2", "Flt1",
    "Bdnf", "Stat5a", "Prkg1", "Flt1",
    "Stat5a", "Mras", "Prkg1", "Mras", "Ptpre",
    "Gria1", "Gria2", "Bcl3", "Zbtb20", "Klf9", 
    "Pou4f1", "Ptprr", "Atp2b2", "Cacna1d", "Tgfa", "Fn1"
  ),
  value = 1 # Strength of connection (equal weight)
)

# Transform list to long format
data_plot  <- t(vst(data, blind = F) %>% assay()) %>% as.data.frame()
#data_plot <- data_plot[,order(colSums(data_plot), decreasing = T)[1:1000]]

data_plot$sample <- rownames(data_plot)
data_plot$group <- c(rep(c("Control", "LPE"), each = 4))




#gene_list <-  c("Sv2a","Nrcam","Nomo1","Grm5","Itgb8","Erbb2","Met","Atp2b1")%>% as.data.frame()
gene_list <-  links$target %>% as.data.frame()



#gene_list <- res[res$threshold_OE == "TRUE", ]$gene_name %>% as.data.frame() %>% rev()

gene_list$gene_id <-DATA.org[match(sapply(strsplit(gene_list$.,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

data_plot2 <- data_plot[,colnames(data_plot) %in% gene_list$gene_id] %>% t() %>% as.data.frame()

rownames(data_plot2) <-DATA.org[match(sapply(strsplit(rownames(data_plot2),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"] 
library(ggalluvial)
links$source <- factor(links$source, levels = unique(links$source))
links$target <- factor(links$target, levels = rownames(data_plot2))


# Use the 'links' data frame created above
sankey <- ggplot(links, aes(y = value, axis1 = source, axis2 = target)) +
  geom_alluvium(aes(fill = source), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Parent Gene", "Downstream Target"), expand = c(.05, .05)) +
  theme_minimal() +
  theme(legend.position = "none") + theme(plot.margin = margin(r = 0))+ 
  scale_y_continuous(expand = c(0, 0))


data_plot2 <- data_plot2[c(1:5, 5,6:8,8,9 ,10,10,11, 12,13,13,14,14,15:19 ),]
data_hp<- data_plot2 %>% pheatmap:::scale_rows() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, 
               names_to = "Sample", 
               values_to = "Expression")

data_hp$Gene <-  factor(data_hp$Gene, levels = rev(rownames(data_plot2)))
heatmap <- ggplot(data_hp, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_classic() +
  theme(axis.title.y = element_blank(),   # Removes the "Gene" title
    axis.text.y  = element_blank(),   # Removes the actual gene names
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    # angle = 45 tilts the text
    # vjust and hjust align the text so the end of the word sits under the tick mark
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )+ theme(plot.margin = margin(l = 0))+ 
  scale_y_discrete(expand = c(0, 0))

combined_plot <- (sankey + heatmap) + 
  plot_layout(widths = c(2, 3)) & 
  theme(
    # This ensures the actual 'boxes' of the plots align perfectly
    plot.margin = margin(t = 5, r = 0, b = 5, l = 0) 
  ) + plot_layout(guides = "collect") & 
  theme(aspect.ratio = NULL)

ggsave("Final/Srp_Bulk.pdf", device= "pdf", plot = combined_plot, width = 8, height = 6 , dpi =600)


Heatmap(data_plot2 %>% pheatmap:::scale_rows(),  name = "Relative expression",
        cluster_columns = F,
        cluster_rows = F,
        #top_annotation = Annotation,
        column_title = "Wnt Signalling from ChatGPT",
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right"
)
```

```{r heatmap wnt}
ggplot(data_plot_long, aes(x =expression, y = gene, fill = Condition)) +
  scale_fill_manual(values = mycolors)+
    stat_summary(fun = "mean", geom = "bar", alpha = 0.9,color = "black", linewidth = 0.5,)+
    #stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),size = 0.9, color = "black")+
    theme_classic() +ylab(NULL)+xlab(NULL)+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,  hjust=1, face = "bold"))+ 
    theme(legend.position="none")
```

```{r}
#data_plot  <- t(vst(data, blind = F) %>% assay()) %>% as.data.frame() %>% t()
data_plot  <- t(counts(data,normalized = T)) %>% as.data.frame() %>% t()
data_top <- data.frame(
  gene = rownames(data_plot[,c(5:8)]),
  #gene = rownames(data_plot),
  mean_count = rowMeans(as.matrix(data_plot[,c(5:8)]), na.rm = TRUE),
  #mean_count = rowMeans(as.matrix(data_plot), na.rm = TRUE),
  stringsAsFactors = FALSE) %>%
  arrange(desc(mean_count)) %>%
  #slice_head(n = 4000) %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),  # order for plotting
    #highlight = if_else(gene %in% highlight_genes, "highlight", "other")
    rank = row_number()
  )


data_top$gene_name <- DATA.org[match(sapply(strsplit(rownames(data_top),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]

#data_top %>% filter(gene_name %in% c("Sv2a","Nrcam","Nomo1","Grm5","Itgb8","Erbb2","Met","Atp2b1") )
rank <- data_top %>% filter(gene_name %in% data_pro$gene_name )

#data_top <- data_top %>% mutate(highlight = if_else(gene_name %in% c("Sv2a","Nrcam","Nomo1","Grm5","Itgb8","Erbb2","Met","Atp2b1"), "highlight", "other"))

data_top <- data_top %>% mutate(highlight = if_else(gene_name %in% data_pro$gene_name, "highlight", "other"))

rank <- inner_join(rank, data_pro, by ="gene_name")
write.csv(rank, "Bulk/MS_data_rank.csv")
#write.csv(rank, "Bulk/MS_data_rank2.csv") #rank is in LPE group, while rank2 is in all samples

ggplot(data_top, aes(x = reorder(gene_name, mean_count), y = mean_count, fill = highlight)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c(other = "grey80", highlight = "red3")) +
  labs(x = NULL, y = "Mean counts", title = "Top 1000 genes by mean count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank() 
  )+ coord_flip()

```

```{r}
genes <- downstream_targets <- list(  Ncan = c("Sox9", "Tnr", "Bcan", "Acan"),  Pgk1 = c("Ldha", "Slc2a1", "Pdk1", "Aldoa"),  Maged1 = c("Jun", "Fas", "Bax"),  Vcp = c("Atf4", "Ddit3", "Hspa5", "Xbp1"),  Itgb8 = c("Serpine1", "Smad7", "Ctgf", "Col1a1"),  Grm5 = c("Arc", "Egr1", "Fos", "Bdnf"),  Anp32a = c("Ccnd1", "Myc", "E2f1"),  Strap = c("Smad7", "Serpine1", "Id1"),  Mapk1 = c("Fos", "Egr1", "Arc", "Myc", "Elk1"),  Stat3 = c("Socs3", "Bcl2", "Ccnd1", "Gfap", "Myc"),  Akt1 = c("Ccnd1", "Bcl2", "Mdm2", "Rps6kb1"),  Mat2a = c("Dnmt1", "Dnmt3a", "Ezh2"),  Ptpa = c("Myc", "Ccnd1", "Ctnnb1"),   Nrcam = c("Egr1", "Fos", "Mapk1"),  Atp2b1 = c("Egr1", "Fos", "Nfatc1"),  Xpo5 = c("Pten", "Cdkn1a"),   # via miRNA regulation  Nomo1 = c("Lefty1", "Pitx2", "Foxh1"), ,  Serpinb1a = c("Il1b", "Tnf", "Cxcl1")
Ppme1 = c("Myc", "Ccnd1", "Ctnnb1"),  Src = c("Fos", "Myc", "Stat3", "Mmp9", "Vegfa"),  Hook3 = c("Fos", "Egr1"),  Pcmt1 = c("Hsf1", "Hmox1"),  Gstp1 = c("Nqo1", "Hmox1", "Gclc"),  Fam171a2 = c("Map2", "Syn1"),  Mcu = c("Nrf1", "Ppargc1a", "Atf4"),  Erbb2 = c("Ccnd1", "Myc", "Fos", "Vegfa"),  Sv2a = c("Arc", "Fos", "Bdnf"),  Prkag1 = c("Ppargc1a", "Cpt1a", "Slc2a4"),  Kazn = c("Itga6", "Ctnnb1"),  Uimc1 = c("Cdkn1a", "Gadd45a", "Brca1"),  Arpp21 = c("Syn1", "Camk2a"),  Srcin1 = c("Fos", "Myc", "Stat3"),  Sfn = c("Cdkn1a", "Ccnb1"),  Pomc = c("Cartpt", "Fos", "Creb1"),  Met = c("Ccnd1", "Myc", "Mmp9", "Vegfa", "Fos"))
ht_list <- NULL
a = 1
for (i in genes){

  cat(paste0("Current pathway is ", names(genes[a]), "\n"))
gene_list = i %>% as.data.frame()
gene_list$gene_id <-DATA.org[match(sapply(strsplit(i ,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

data_plot2 <- data_plot[rownames(data_plot) %in% gene_list$gene_id, ]  %>% as.data.frame()
data_plot2$gene_name<-DATA.org[match(sapply(strsplit(rownames(data_plot2),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]
data_plot2$gene_name[is.na(data_plot2$gene_name)] <- rownames(data_plot2)[is.na(data_plot2$gene_name)]
rownames(data_plot2) <- data_plot2$gene_name
data_plot2$gene_name <- NULL
data_plot2 <- t(data_plot2 %>% as.data.frame()) %>% as.data.frame()
ht= Heatmap(t(data_plot2) %>% pheatmap:::scale_rows(),  name = names(genes[a]),
        cluster_columns = F,
        #cluster_rows = T,
        #top_annotation = Annotation,
        row_title = names(genes[a]),
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 10),
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right")
ht_list <- ht_list %v% ht
a = a+1
}

pdf("Bulk/heatmaps.pdf", width = 5, height = 30)
draw(ht_list,show_heatmap_legend = FALSE)
dev.off()

```

```{r UCell}
genes <- downstream_targets <- list(  Ncan = c("Sox9", "Tnr", "Bcan", "Acan"),  Pgk1 = c("Ldha", "Slc2a1", "Pdk1", "Aldoa"),  Maged1 = c("Jun", "Fas", "Bax"),  Vcp = c("Atf4", "Ddit3", "Hspa5", "Xbp1"),  Itgb8 = c("Serpine1", "Smad7", "Ctgf", "Col1a1"),  Grm5 = c("Arc", "Egr1", "Fos", "Bdnf"),  Anp32a = c("Ccnd1", "Myc", "E2f1"),  Strap = c("Smad7", "Serpine1", "Id1"),  Mapk1 = c("Fos", "Egr1", "Arc", "Myc", "Elk1"),  Stat3 = c("Socs3", "Bcl2", "Ccnd1", "Gfap", "Myc"),  Akt1 = c("Ccnd1", "Bcl2", "Mdm2", "Rps6kb1"),  Mat2a = c("Dnmt1", "Dnmt3a", "Ezh2"),  Ptpa = c("Myc", "Ccnd1", "Ctnnb1"),   Nrcam = c("Egr1", "Fos", "Mapk1"),  Atp2b1 = c("Egr1", "Fos", "Nfatc1"),  Xpo5 = c("Pten", "Cdkn1a"),   # via miRNA regulation  Nomo1 = c("Lefty1", "Pitx2", "Foxh1"), ,  Serpinb1a = c("Il1b", "Tnf", "Cxcl1")
Ppme1 = c("Myc", "Ccnd1", "Ctnnb1"),  Src = c("Fos", "Myc", "Stat3", "Mmp9", "Vegfa"),  Hook3 = c("Fos", "Egr1"),  Pcmt1 = c("Hsf1", "Hmox1"),  Gstp1 = c("Nqo1", "Hmox1", "Gclc"),  Fam171a2 = c("Map2", "Syn1"),  Mcu = c("Nrf1", "Ppargc1a", "Atf4"),  Erbb2 = c("Ccnd1", "Myc", "Fos", "Vegfa"),  Sv2a = c("Arc", "Fos", "Bdnf"),  Prkag1 = c("Ppargc1a", "Cpt1a", "Slc2a4"),  Kazn = c("Itga6", "Ctnnb1"),  Uimc1 = c("Cdkn1a", "Gadd45a", "Brca1"),  Arpp21 = c("Syn1", "Camk2a"),  Srcin1 = c("Fos", "Myc", "Stat3"),  Sfn = c("Cdkn1a", "Ccnb1"),  Pomc = c("Cartpt", "Fos", "Creb1"),  Met = c("Ccnd1", "Myc", "Mmp9", "Vegfa", "Fos"))

ht_list <- NULL
a = 1
for (i in genes){

  cat(paste0("Current pathway is ", names(genes[a]), "\n"))
gene_list = i %>% as.data.frame()
gene_list$gene_id <-DATA.org[match(sapply(strsplit(i ,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]
ucell_scores <- ScoreSignatures_UCell(as.matrix(data_plot), features = gene_list$gene_id)

data_plot2 <- ucell_scores[,1] %>% t()
rownames(data_plot2) <- names(genes[a])
data_plot2 <- t(data_plot2 %>% as.data.frame()) %>% as.data.frame()
ht= Heatmap(t(data_plot2) ,  name = names(genes[a]),
        cluster_columns = F,
        #cluster_rows = T,
        #top_annotation = Annotation,
        row_title = names(genes[a]),
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 10),
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right")
ht_list <- ht_list %v% ht
a = a+1
}

pdf("Bulk/Score_heatmaps.pdf", width = 5, height = 30)
draw(ht_list,show_heatmap_legend = FALSE)
dev.off()

```

```{r neurgenesis}
gene_list <-  c("Fn1", "Atp1b2", "Ednrb", "Sgk1", "Mgll", "Alcam", "Kcnj10", "Sema4d", "Pmp22", "Plec", "Sirt2", "Cyfip2", "Arsb", "Gpr17", "Tnr", "S1pr1", "Cntn1", "Clu", "Sorl1", "Gfap", "Dab2", "Enpp2", "Rnd2", "Cdh11", "Kank1", "Daam2", "Prkg1", "Erbb4", "Sema5a", "Nfasc", "Atxn1", "Gpr37l1", "Epha4", "Plp1", "Vxn", "Atcay", "Ugt8a", "Prex2", "Hgf", "Plppr5", "Lrrc4c", "Fry", "Slc1a1", "Smim45", "Sema6b", "P2ry1", "Myrf", "Disp3", "Bdnf", "Brinp1", "Pls1", "Bmp6", "Cend1", "Dnm3", "Id4", "Unc5c", "Reln", "Zfp488", "Septin4", "Slc39a12", "Bmp4", "Sh3gl3", "Ptpn5", "Grid2", "Sema3d", "Slc4a10", "Rom1", "C1ql1", "Erbb3", "Csmd3", "Dok5", "Dusp15", "Tmem108", "Scarf1", "Adora2a", "Hey2", "Mag", "S100b", "Atp2b2", "Cntn6", "Cthrc1", "Neu4", "Ripor2", "Sema3e", "Rgs6", "Cntnap2", "Gabrb1", "Wnt4", "Flrt1", "Nyap2", "Spock1", "Rarb", "Gabrb2", "Plppr4", "Rassf10", "Syt3", "Tecta", "Olfm3", "Slitrk6", "Cdk5r2", "Nkx6-2", "Kcna1", "Pou4f1", "Acan", "Emx1", "Cntn2", "Dct", "Gldn", "Zfp804a", "Chrna3", "Hapln1", "Lgi4")%>% as.data.frame()

gene_list$gene_id <-DATA.org[match(sapply(strsplit(gene_list$.,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

data_plot2 <- data_plot[,colnames(data_plot) %in% gene_list$gene_id] %>% t() %>% as.data.frame()
data_plot2$gene_name<-DATA.org[match(sapply(strsplit(rownames(data_plot2),"\\."), function(x) x[1]), DATA.org$gene_id),"gene_name"]
data_plot2$gene_name[is.na(data_plot2$gene_name)] <- rownames(data_plot2)[is.na(data_plot2$gene_name)]
rownames(data_plot2) <- data_plot2$gene_name
data_plot2$gene_name <- NULL
data_plot2 <- t(data_plot2 %>% as.data.frame()) %>% as.data.frame()

```

```{r heatmap neurogenesis}
#| fig-width: 6
#| fig-height: 16
Heatmap(t(data_plot2) %>% pheatmap:::scale_rows(),  name = "Relative expression",
        cluster_columns = F,
        cluster_rows = T,
        #top_annotation = Annotation,
        column_title = "Neurogenesis GO term",
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right"
)
```

### Module score for Wnt

```{r}
gene_list <-  c(
  # Wnt ligands
  "Wnt1","Wnt3a","Wnt4","Wnt5a","Wnt7a","Wnt7b",
  # Wnt receptors & co-receptors
  "Fzd3","Fzd5","Fzd7","Lrp5","Lrp6",
  # Wnt intracellular mediators
  "Dvl1","Dvl2","Ctnnb1","Axin2","Gsk3b","Apc",
  # Wnt transcriptional regulators
  "Tcf7l2","Lef1",
  # Wnt antagonists/modulators
  "Dkk1","Sfrp1","Sfrp2",
  # BMP ligands
  "Bmp2","Bmp4","Bmp7","Bmp7",
  # BMP receptors
  "Bmpr1a","Bmpr1b","Bmpr2",
  # BMP Smads
  "Smad1","Smad5","Smad9","Smad4",
  # BMP antagonists
  "Nog","Grem1","Chrd",
  # Crosstalk regulators
  "Rnf146","Pitx2","Ccnd1","Smurf1","Smurf2"
)%>% as.data.frame()
#gene_list <-  c("Tcf7l2", "Nog", "Bmp4", "Smad9", "Wnt7b", "Wnt4", "Bmp2", "Ctnnb1", "Sfrp1", "Gsk3b", "Bmpr2", "Fzd3", "Bmpr1a", "Wnt5a", "Wnt7a")%>% as.data.frame()

gene_list <-  c("Wnt3a", "Wnt5a","Wnt7a", "Fzd1","Fzd3","Fzd6","Fzd7", "Ngn2","Neurod1","Prox1", "Atp6ap2")%>% as.data.frame()

data_plot  <- t(vst(data, blind = F) %>% assay()) %>% as.data.frame()


gene_list$gene_id <-DATA.org[match(sapply(strsplit(gene_list$.,"\\."), function(x) x[1]), DATA.org$gene_name),"gene_id"]

gene_list <- gene_list[gene_list$gene_id %in% colnames(data_plot), ] 
```

```{r}
data2 <-data_plot[,colnames(data_plot) %in% gene_list$gene_id]

data2 <-  as.data.frame(t(data2))

rownames(data2) <-tx2gene[match(sapply(strsplit(rownames(data2),"\\."), function(x) x[1]), tx2gene$GENEID),"GENENAME"]

gene_list$gene_name[ as.character(gene_list$gene_name) %in% as.character(rownames(data2))]
```

```{r}
gene_set <- list(Module = gene_list$gene_id)

# Perform GSVA to compute module scores

para <- ssgseaParam(as.matrix(t(data_plot)), gene_set)
module_scores <- gsva(para)

wnt <- as.data.frame(t(module_scores))
wnt$Samples <- rownames(wnt)
wnt$condition <- rep(c("Control", "LPE"), each = 4)

```

```{r}
pdf("../Result/DAM_score.pdf" ,width = 8, height = 6 )
ggplot(wnt)+ aes(x = condition, y = Module) + 
    #geom_line(data =  aggregate(Arg1 ~ Group, data = data, mean), color = "azure3", group = 1)+
    geom_boxplot(aes(color = condition),staplewidth = 0.2, width = 0.3)+ 
    geom_point(aes(color = condition),position=position_jitter(w = 0,h = 0)) +
    #geom_line(aes(group=Samples,color = Samples), alpha=0.5)+
  #geom_text_repel(aes(label = rownames(d)), max.overlaps = Inf) + 
    theme_classic() +
    ylab(NULL)+ xlab(NULL)+
    ggtitle("Wnt module score") +
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,  hjust=1))

ggplot(DAM2)+ aes(x = condition, y = Module, color = condition) + 
    #geom_line(data =  aggregate(Arg1 ~ Group, data = data, mean), color = "azure3", group = 1)+
    geom_boxplot(staplewidth = 0.2, width = 0.3)+ 
    geom_point(position=position_jitter(w = 0,h = 0)) +
    #geom_line(aes(group=biosample), color = "black", alpha=0.3)+
  #geom_text_repel(aes(label = rownames(d)), max.overlaps = Inf) + 
    theme_classic() +
    ylab(NULL)+ xlab(NULL)+
    ggtitle("DAM module score") +
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,  hjust=1))+ theme(legend.position="none")
dev.off()

```



```{r}
links <- data.frame(
  source = c(
    "Akt1", "Akt1", "Akt1", "Akt1",
    "Erbb2", "Erbb2", "Erbb2", "Erbb2",
    "Met", "Met", "Src", "Src", "Src",
    "Grm5", "Grm5", "Stat3", "Mapk1", "Mapk1", 
    "Mapk1", "Mapk1", "Atp2b1", "Mcu", "Itgb8", "Strap"
  ),
  target = c(
    "Bdnf", "Hey1", "Hey2", "Flt1",
    "Bdnf", "Stat5a", "Prkg1", "Flt1",
    "Stat5a", "Mras", "Prkg1", "Mras", "Ptpre",
    "Gria1", "Gria2", "Bcl3", "Zbtb20", "Klf9", 
    "Pou4f1", "Ptprr", "Atp2b2", "Cacna1d", "Tgfa", "Fn1"
  ),
  value = 1 # Strength of connection (equal weight)
)

# Transform list to long format
links <- stack(downstream_map) %>%
  rename(source = ind, target = values) %>%
  mutate(value = 1) # Assign equal weight to all connections

library(ggalluvial)

# Use the 'links' data frame created above
ggplot(links, aes(y = value, axis1 = source, axis2 = target)) +
  geom_alluvium(aes(fill = source), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Parent Gene", "Downstream Target"), expand = c(.05, .05)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Gene to Downstream Target Mapping", y = "Connection Count")


Heatmap(t(data_plot2) %>% pheatmap:::scale_rows(),  name = "Relative expression",
        cluster_columns = F,
        cluster_rows = T,
        #top_annotation = Annotation,
        column_title = "Wnt Signalling from ChatGPT",
        #left_annotation = ha,
        #column_order = meta2$SampleID,
        row_names_side = "left", row_dend_side = "right"
)
```
