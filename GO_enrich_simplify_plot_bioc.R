GO_enrich_simplify_plot_bioc = function(protein_set, reference_protein_set, identifier_type, ontology, pAdjustMethod_ = "BH", minSetSize, maxSetSize, simplify_by = "p.adjust", simplify_fun = "min", similarity_calc_method = "Wang", similarity_cutoff = 0.7, visualize_result = "enrichMap", above_corrected_pval = 1){
  suppressPackageStartupMessages({
    library(igraph)
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
                   pvalueCutoff = 1,
                   qvalueCutoff = above_corrected_pval)
  
  # simplify output from enrichGO by removing redundancy of enriched GO terms
  ego2 = simplify(ego, cutoff = similarity_cutoff,
                  by = simplify_by, 
                  select_fun = eval(parse(text = simplify_fun)), 
                  measure = similarity_calc_method)
  
  # visualize results
  if(visualize_result == "enrichMap") enrichMap(ego2, layout = layout_with_kk, vertex.label.cex = 0.8,vertex.size = 5, rescale=T)
  if(visualize_result == "dotplot") dotplot(ck2)
  return(list(ego, ego2))
}

cluster_GO_enrich_simplify_plot_bioc = function(protein_table, reference_protein_set, ID_by_group_formula, ontology, pAdjustMethod_ = "BH", minSetSize, maxSetSize, simplify_by = "GeneRatio", simplify_fun = "function(x) eval(parse(text = x))", similarity_calc_method = "Wang", similarity_cutoff = 0.7){
  suppressPackageStartupMessages({
    library(igraph)
    library(clusterProfiler)
    library(Homo.sapiens)
  })
  
  ck <- compareCluster(geneCluster = ID_by_group_formula, 
                       data = protein_table,
                       OrgDb         = org.Hs.eg.db,
                       universe = reference_protein_set,
                       keytype       = identifier_type,
                       ont           = ontology,
                       pAdjustMethod = pAdjustMethod_,
                       minGSSize = minSetSize,
                       maxGSSize = maxSetSize,
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       fun = "enrichGO")
  
  # simplify output from enrichGO by removing redundancy of enriched GO terms
  ck2 = simplify(ck, cutoff = similarity_cutoff, by = simplify_by, select_fun = eval(parse(text = simplify_fun)), measure = similarity_calc_method)
  
  # visualize results
  dotplot(ck2)
  #dotplot(ck2, x=~group) + ggplot2::facet_grid(~othergroup)
  #enrichMap(ego2)
  return(list(ck, ck2))
}

