# testing go enrichment and simplification with kappa score similarity
library(clusterProfiler)
library(org.Hs.eg.db)
data(geneList, package="DOSE")

ego3 <- gseGO(geneList     = geneList,
OrgDb        = org.Hs.eg.db,
ont          = "CC",
nPerm        = 1000,
minGSSize    = 100,
maxGSSize    = 500,
pvalueCutoff = 0.05,
verbose      = FALSE)

ego <- enrichGO(gene          = names(geneList)[abs(geneList) > 2],
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                minGSSize = 50,
                maxGSSize = 500,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

source("/Users/vitalii/Desktop/VItalii_EBI/imex_vs_uniprot/multiplot.R")
source("/Users/vitalii/Desktop/VItalii_EBI/imex_vs_uniprot/GO_analyse_simplify_plot.R")
source("/Users/vitalii/Desktop/VItalii_EBI/imex_vs_uniprot/kappa_score.R")
source("/Users/vitalii/Desktop/VItalii_EBI/imex_vs_uniprot/dotplot_modified.R")
source("/Users/vitalii/Desktop/VItalii_EBI/imex_vs_uniprot/clusterProfiler_simplify_methods_modified.R")

simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL, use_data_table = T, use_bioc_annotationdbi = T)

simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min, measure="kappa", semData = NULL, use_data_table = T, use_bioc_annotationdbi = T)

simplify(ego3, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL, use_data_table = T, use_bioc_annotationdbi = T)

simplify(ego3, cutoff=0.7, by="p.adjust", select_fun=min, measure="kappa", semData = NULL, use_data_table = T, use_bioc_annotationdbi = T)

########################################
### multi gseaplot
for(i in ego3@result$ID) assign(paste0("gseaplot",gsub(":","_",i)), gseaplot(ego3, i),.GlobalEnv)
eval(parse(text = paste0("gsea_list = c(",paste0("gseaplot",gsub(":","_",ego3@result$ID), collapse = ","),")")))
#names(gsea_list) = paste0("gseaplot",gsub(":","_",ego3@result$ID))
#multiplot(gseaplotGO_0031012,gseaplotGO_0005578,gseaplotGO_0042383,gseaplotGO_0030496,gseaplotGO_0005788, cols = 3)
gsea_list[-names(gsea_list) == "runningScore"]
multiplot(plotlist = gsea_list, cols = 5)

dotplot(ego3, colorBy = "enrichmentScore", orderBy = "p.adjust")
dotplot(ego3, colorBy = "enrichmentScore")
dotplot(ego)
enrichMap(ego3, layout = layout_with_kk, vertex.label.cex = 0.8,vertex.size = 8, rescale=T)

x1 = clusterProfiler:::get_GO_data(OrgDb = 'org.Hs.eg.db', ont = ego3@setType, keytype = ego3@keytype)
x2 = GOSemSim::godata(OrgDb = 'org.Hs.eg.db', ont = ego3@setType, keytype = ego3@keytype, computeIC = F)

GO_from_GODB = toTable(GO.db::GOTERM)

GOSemSim::load_OrgDb

unique(Reduce(c, x1$EXTID2PATHID))
x2$
      clusterProfiler:::get
goenv = clusterProfiler:::get_GO_Env()

mean((names(x1$GO2ONT)[x1$GO2ONT == ego3@setType]) %in% (x2@geneAnno$GO))

x1$PATHID2EXTID[unique(names(x1$GO2ONT)[x1$GO2ONT == ego3@setType][!(names(x1$GO2ONT)[x1$GO2ONT == ego3@setType]) %in% (x2@geneAnno$GO)])]

#source("https://raw.githubusercontent.com/GuangchuangYu/clusterProfiler/master/R/enrichGO.R")
#source("https://raw.githubusercontent.com/GuangchuangYu/clusterProfiler/master/R/go-utilities.R")

downloader::download("http://geneontology.org/gene-associations/goa_human.gaf.gz", "goa_human.gaf.gz")
R.utils::gunzip("goa_human.gaf.gz", remove = F, overwrite = T)
goa = fread("goa_human.gaf", stringsAsFactors = F,skip = 35, header = F)
goa_terms = unique(goa$V5)

mean((names(x1$GO2ONT)[x1$GO2ONT == ego3@setType]) %in% goa_terms)
mean((x2@geneAnno$GO) %in% goa_terms)

names(x1$GO2ONT)[x1$GO2ONT == ego3@setType][!(names(x1$GO2ONT)[x1$GO2ONT == ego3@setType]) %in% goa_terms]

clusterProfiler:::GSEA_internal
DOSE:::GSEA_fgsea
DOSE:::getGeneSet
DOSE:::build_Anno

clusterProfiler:::gseGO
clusterProfiler:::get_GO_data
DOSE:::build_Anno
