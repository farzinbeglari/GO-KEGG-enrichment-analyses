# GO-KEGG-enrichment-analyses
Gene Ontology and Kyoto Encyclopedia of Genes and Genomes (filtered.counts in WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)
library(ggplot2)

# Biological Process (BP)
genes <- rownames(filtered.counts)
genes <- as.data.frame(genes)
BP <- enrichGO(gene = genes$genes,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = 'BP',
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               pAdjustMethod = 'BH')
a1 <- dotplot(BP, showCategory = 5) +
  ggtitle('Biological Process')

# Molecular Function (MF)
MF <- enrichGO(gene = genes$genes,
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = 'MF',
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjustMethod = 'BH')
a2 <- dotplot(MF, showCategory = 5) +
  ggtitle('Molecular Function')

# Cellular Component (CC)
CC <- enrichGO(gene = genes$genes,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = 'CC',
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               pAdjustMethod = 'BH')
a3 <- dotplot(CC, showCategory = 5) +
  ggtitle('Cellular Component')

# KEGG enrichment analysis
gene.id <- bitr(geneID = genes$genes,
                fromType = 'SYMBOL',
                toType = 'ENTREZID',
                OrgDb = org.Hs.eg.db)
gene.id <- unique(gene.id$ENTREZID)
KEGG <- enrichKEGG(gene.id,
                   organism = 'hsa',
                   keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.05)
a4 <- dotplot(KEGG, showCategory = 5) +
  ggtitle('KEGG enrichment analysis')
# Visualization
grid.arrange(a1, a2, a3, a4, nrow = 2, ncol = 2)
