##' @description High level GO enrichment and GSEA functions based on clusterprofiler (analyse, simplify, plot)
##' @author Vitalii Kleshchevnikov 
#####################################
##' @author Vitalii Kleshchevnikov
GO_enrich_simplify_plot_bioc = function(protein_set, reference_protein_set, identifier_type, ontology, pAdjustMethod_ = "BH", minSetSize, maxSetSize, simplify_by = "p.adjust", simplify_fun = "min", similarity_calc_method = "kappa", similarity_cutoff = 0.7, visualize_result = "enrichMap", above_corrected_pval = 1, use_bioc_annotationdbi = T, plot_title = "", xlabel = ""){
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(Homo.sapiens)
  })
  
  ego <- enrichGO(gene = protein_set,
                  OrgDb = org.Hs.eg.db,
                  universe = reference_protein_set,
                  keytype = identifier_type,
                  ont = ontology,
                  pAdjustMethod = pAdjustMethod_,
                  minGSSize = minSetSize,
                  maxGSSize = maxSetSize,
                  pvalueCutoff = above_corrected_pval,
                  qvalueCutoff = 1)

  if(similarity_calc_method != "none"){
  # simplify output from enrichGO by removing redundancy of enriched GO terms
  ego2 = simplify(x = ego, cutoff = similarity_cutoff,
                  by = simplify_by, 
                  select_fun = eval(parse(text = simplify_fun)), 
                  measure = similarity_calc_method,
                  semData = NULL, 
                  use_bioc_annotationdbi = use_bioc_annotationdbi)
  # "x[which.max(eval(parse(text = paste0("c(",paste0(x, collapse = ","),")"))))]"
  }
  if(similarity_calc_method == "none") ego2 = ego
  
  # visualize results
  if(visualize_result == "enrichMap") plot_res = enrichMap(ego2, layout = igraph:::layout_with_kk, vertex.label.cex = 0.8,vertex.size = 5, rescale=T)
  if(visualize_result == "dotplot") plot_res = dotplot(ego2, title = plot_title, xlabel = xlabel)
  return(list(enrichment_result = ego, simplified_enrichment_result = ego2, plot = plot_res))
}
#####################################
##' @author Vitalii Kleshchevnikov
GSEA_simplify_plot_bioc = function(ranked_protein_list, identifier_type, ontology, nPerm = 1000, pAdjustMethod_ = "BH", minSetSize, maxSetSize, simplify_by = "p.adjust", simplify_fun = "min", similarity_calc_method = "kappa", similarity_cutoff = 0.7, visualize_result = "enrichMap", above_corrected_pval = 1, use_bioc_annotationdbi = T, plot_title = "", xlabel = ""){
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(Homo.sapiens)
  })
  
  ego = gseGO(geneList = ranked_protein_list, ont = ontology, org.Hs.eg.db, keytype = identifier_type, exponent = 1,
              nPerm = nPerm, minGSSize = minSetSize, maxGSSize = maxSetSize, pvalueCutoff = above_corrected_pval,
              pAdjustMethod = pAdjustMethod_, verbose = F, seed = FALSE, by = "fgsea")
  
  if(similarity_calc_method != "none"){
  # simplify output from enrichGO by removing redundancy of enriched GO terms
  ego2 = simplify(x = ego, cutoff = similarity_cutoff,
                  by = simplify_by, 
                  select_fun = eval(parse(text = simplify_fun)), 
                  measure = similarity_calc_method,
                  semData = NULL, 
                  use_data_table = T, 
                  use_bioc_annotationdbi = use_bioc_annotationdbi)
  # "x[which.max(eval(parse(text = paste0("c(",paste0(x, collapse = ","),")"))))]"
  }
  if(similarity_calc_method == "none") ego2 = ego
  
  # visualize results
  if(visualize_result == "enrichMap") plot_res = enrichMap(ego2, layout = layout_with_kk, vertex.label.cex = 0.8,vertex.size = 5, rescale=T)
  if(visualize_result == "dotplot") plot_res = dotplot(ego2, title = plot_title, colorBy = "enrichmentScore", orderBy = "p.adjust", xlabel = xlabel)
  return(list(enrichment_result = ego, simplified_enrichment_result = ego2, plot = plot_res))
}
#####################################
##' @author Vitalii Kleshchevnikov
cluster_GO_enrich_simplify_plot_bioc = function(formula, protein_groups.dt, reference_protein_set, identifier_type, ontology, pAdjustMethod_ = "BH", minSetSize, maxSetSize, simplify_by = "p.adjust", simplify_fun = "min", similarity_calc_method = "kappa", similarity_cutoff = 0.7, visualize_result = "dotplot", above_corrected_pval = 1, plot_title = ""){
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(Homo.sapiens)
  })
  
  ego = compareCluster(formula, fun = "enrichGO", data = protein_groups.dt,
                       OrgDb = org.Hs.eg.db,
                       universe = reference_protein_set,
                       keytype = identifier_type,
                       ont = ontology,
                       pAdjustMethod = pAdjustMethod_,
                       minGSSize = minSetSize,
                       maxGSSize = maxSetSize,
                       pvalueCutoff = above_corrected_pval,
                       qvalueCutoff = 1)
  if(similarity_calc_method != "none"){
  # simplify output from enrichGO by removing redundancy of enriched GO terms
  ego2 = simplify(x = ego, cutoff = similarity_cutoff,
                  by = simplify_by, 
                  select_fun = eval(parse(text = simplify_fun)), 
                  measure = similarity_calc_method,
                  semData = NULL)
  # "x[which.max(eval(parse(text = paste0("c(",paste0(x, collapse = ","),")"))))]"
  }
  if(similarity_calc_method == "none") ego2 = ego
  
  # visualize results
  if(visualize_result == "enrichMap") plot_res = enrichMap(ego2, layout = layout_with_kk, vertex.label.cex = 0.8,vertex.size = 5, rescale=T)
  if(visualize_result == "dotplot") plot_res = dotplot(ego2, title = plot_title)
  return(list(enrichment_result = ego, simplified_enrichment_result = ego2, plot = plot_res))
}
#####################################