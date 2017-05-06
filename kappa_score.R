# GO term similarity - kappa score

#####################################
##' @author Vitalii Kleshchevnikov
# function to calculate kappa score between two categories describing a set of elements(GO terms, KEGG pathways, genes)
# defined as described in this article:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2375021/figure/F2/
# relies on data.table for subsetting and vcd package for calculating kappa score given 2x2 incidence matrix
# the function is faster than irr::kappa2(ratings = value, weight = "unweighted")
# the function is ~25 times slower than dist()
# function takes a data.table with 2 columns corresponding to two categories to which each element in rows can belong
kappa_score = function(value){
  colnames(value) = c("x","y")
  table_ = value[,.N, by = .(x,y)]
  table_2 = matrix(,2,2)
  table_2[1,1] = ifelse(length(table_[x == 1 & y == 1, N]), table_[x == 1 & y == 1, N],0)
  table_2[1,2] = ifelse(length(table_[x == 1 & y == 0, N]), table_[x == 1 & y == 0, N],0)
  table_2[2,1] = ifelse(length(table_[x == 0 & y == 1, N]), table_[x == 0 & y == 1, N],0)
  table_2[2,2] = ifelse(length(table_[x == 0 & y == 0, N]), table_[x == 0 & y == 0, N],0)
  kappa_score = vcd::Kappa(table_2,weights = "Fleiss-Cohen")$Unweighted
  return(kappa_score)
}

#####################################

#filename = "/Users/vitalii/Downloads/goa_human.gaf"
#mapping_table = fread(filename, skip = 34)[,.(UNIPROT = V2, GO = V5, evidence = V7, ontology = V9)]
#mapping_table = mapping_table[ontology == "P", .(UNIPROT, GO)]
#mapping_table = unique(mapping_table)

#####################################
##' @author Vitalii Kleshchevnikov
# the categ_dist function to calculate categorical distance (Cohen's Kappa score) between multiple terms 
# the function is intended to measure distances between GO terms based on proteins they annotate
# more generally, the function can be used to measure categorical distances between any terms(categories) annotating objects
# objects should be provided as a first column of a data.table, terms should be provided as a second column
categ_dist = function(mapping_table, terms_to_compare = unlist(unique(mapping_table[,2,with = F])), ignore_limit = F){
  if(ncol(mapping_table) > 2) stop("table has more than 2 columns, object id column and term column")
  if(ignore_limit == F) if(length(terms_to_compare) > 1000) stop("more than 1000 terms to compare, set ignore_limit = T if you are sure to proceed")
  if(!is.data.table(mapping_table)) stop("provided mapping / annotation table may not be in the right format (wrong class: not data.table)")
  
  mapping_table = copy(unique(mapping_table))
  print(mapping_table)
  colnames(mapping_table) = c("UNIPROT", "GO")
  z2 = dcast(mapping_table[,.(UNIPROT, GO, value = 1)], UNIPROT ~ GO, fill = 0, drop = F)[,UNIPROT := NULL][,terms_to_compare, with=F]
  combinations = t(caTools::combs(colnames(z2),2))
  dist = t(sapply(as.data.table(combinations), function(x) kappa_score(z2[,c(x[1],x[2]),with = F])))
  dist = cbind(as.data.table(dist), as.data.table(t(combinations)))
  colnames(dist) = c("kappa_score", "kappa_error", "GO1", "GO2")
  dist_temp = unique(rbind(dist,dist[,.(kappa_score,kappa_error, GO1 = GO2, GO2 = GO1)]))
  
  dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2))
  rownames_dist2 = dist2[,"GO1"]
  dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2)[,GO1 := NULL])
  rownames(dist2) = rownames_dist2
  dist2 = dist2[sort(rownames(dist2)), sort(colnames(dist2))]
  diag(dist2) = 1
  return(list(similarity_matrix = dist2, kappa_score_table = dist, kappa_score_table_redundant = dist_temp))
}