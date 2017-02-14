# biogrid

filename_biogrid = paste("./biogrid/","BIOGRID-ALL-3.4.142.mitab.txt", sep = "")
biogrid = read.delim(filename_biogrid, stringsAsFactors =F)
head(biogrid,3)

library(dplyr)
biogrid_temp1 = as.character(c(biogrid$X.ID.Interactor.A, biogrid$ID.Interactor.B))
biogrid_temp2 = as.character(c(biogrid$Taxid.Interactor.A, biogrid$Taxid.Interactor.B))
biogrid_temp3 = data.frame(interactors = biogrid_temp1, taxid = biogrid_temp2, stringsAsFactors = F)
head(biogrid_temp3,3)
biogrid_human = filter(biogrid_temp3, taxid == "taxid:9606")
biogrid_human_unique = unique(biogrid_human)

x = gsub("entrez gene/locuslink:","", biogrid_human_unique$interactors)

filename = paste("./biogrid/","entrez_geneID_",SPECIES_NAME,"_BIOGRID_", Sys.Date(),".txt", sep = "")
write.table(x,filename,col.names=T,row.names=F,sep="\t",quote=F)


#15512 out of 16170 Entrez Gene (GeneID) identifiers were successfully mapped to 27450 UniProtKB IDs in the table below.
#Click here to download the 658 unmapped identifiers.
#Reviewed (15,432) Swiss-Prot
#Unreviewed (12,018) TrEMBL
#http://www.uniprot.org/uniprot/?query=yourlist:M201611118A530B6CA0138AFAA6D2B97CE8C2A924B42AA06&sort=yourlist:M201611118A530B6CA0138AFAA6D2B97CE8C2A924B42AA06&columns=yourlist(M201611118A530B6CA0138AFAA6D2B97CE8C2A924B42AA06),isomap(M201611118A530B6CA0138AFAA6D2B97CE8C2A924B42AA06),id%2Creviewed%2Corganism%2Cproteome%2Corganism-id%2Cprotein%20names%2Ccomment%28ALTERNATIVE%20PRODUCTS%29