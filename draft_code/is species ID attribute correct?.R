## Checking identifiers - is species ID attribute correct?
filename = "./biogrid/query-taxidA_9606 AND taxi-28112016_1158.txt"
biogrid_taxid_9606 = read.delim(filename, stringsAsFactors = F)
colnames(biogrid_taxid_9606)

temp1 = c(biogrid_taxid_9606$X.ID.s..interactor.A,biogrid_taxid_9606$ID.s..interactor.B)
temp1_1 = c(biogrid_taxid_9606$Taxid.interactor.A, biogrid_taxid_9606$Taxid.interactor.B)

biogrid_taxid_9606_interactors = data.frame(interactor_IDs = temp1, interactor_SPECIES_ID = temp1_1, stringsAsFactors = F)
biogrid_taxid_9606_interactors2 = unique(biogrid_taxid_9606_interactors)
table(biogrid_taxid_9606_interactors2$interactor_SPECIES_ID)
interactorIDs = gsub("entrez gene/locuslink:", "", biogrid_taxid_9606_interactors2$interactor_IDs)

filename2 = "./biogrid/query-taxidA_9606 AND taxi-28112016_1158_interactorIDs.txt"
write(interactorIDs,filename2)

## Check mentha
SPECIES_NAME = "Homo sapiens"
SPECIES_NAME = "Mus musculus"
SPECIES_NAME = "Caenorhabditis elegans"
  
source("SPECIES_NAME_TO_ID.R")
SPECIES_IDs = SPECIES_NAME_TO_ID(SPECIES_NAME)

SPECIES_ID = SPECIES_IDs$SPECIES_ID;
Proteome_ID = SPECIES_IDs$Proteome_ID;


databases <- c("mentha")
## Query PSICQUIC for interactions, get MI-TAB-2.5, save, return
source("query_PSICQUIC_for_interactions.R")
mentha = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, SPECIES_NAME = SPECIES_NAME, databases=databases)

library(dplyr)
biogrid = filter(mentha, V13 == "psi-mi:MI:0463(biogrid)")
dip = filter(mentha, V13 == "psi-mi:MI:0465(dip)")
IntAct = filter(mentha, V13 == "psi-mi:MI:0469(IntAct)")
MINT = filter(mentha, V13 == "psi-mi:MI:0471(MINT)")
matrixdb = filter(mentha, V13 == "psi-mi:MI:0917(matrixdb)")

SPECIES_ID

biogrid_interactors = unique(interactions_to_interactors(biogrid))
biogrid_interactors_9606 = filter(biogrid_interactors, interactor_SPECIES_ID == SPECIES_ID)
filenamex = paste("./biogrid/biogrid_interactors_",SPECIES_ID, ".txt")
write(biogrid_interactors_9606$interactor_IDs, filenamex)

dip_interactors = interactions_to_interactors(dip)
dip_interactors_9606 = filter(dip_interactors, interactor_SPECIES_ID == SPECIES_ID)
filenamex = paste("./biogrid/dip_interactors_",SPECIES_ID, ".txt")
write(dip_interactors_9606$interactor_IDs, filenamex)

IntAct_interactors = interactions_to_interactors(IntAct)
IntAct_interactors_9606 = filter(IntAct_interactors, interactor_SPECIES_ID == SPECIES_ID)
filenamex = paste("./biogrid/IntAct_interactors_",SPECIES_ID, ".txt")
write(IntAct_interactors_9606$interactor_IDs, filenamex)

MINT_interactors = interactions_to_interactors(MINT)
MINT_interactors_9606 = filter(MINT_interactors, interactor_SPECIES_ID == SPECIES_ID)
filenamex = paste("./biogrid/MINT_interactors_",SPECIES_ID, ".txt")
write(MINT_interactors_9606$interactor_IDs, filenamex)

matrixdb_interactors = interactions_to_interactors(matrixdb)
matrixdb_interactors_9606 = filter(matrixdb_interactors, interactor_SPECIES_ID == SPECIES_ID)
filenamex = paste("./biogrid/matrixdb_interactors_",SPECIES_ID, ".txt")
write(matrixdb_interactors_9606$interactor_IDs, filenamex)

y = c("Q6PB44 - is mouse, 10090",
"Q6P1X5 - is human, 9606",
"Q6NZY7 - is human, 9606",
"Q6N069 - is human, 9606",
"Q6NZB1 - is mouse, 10090",
"Q6NYC8 - is human, 9606",
"P22670 - is human, 9606",
"P22626 - is human, 9606",
"Q15819 - is human, 9606")

y2 = c("Q6PB44",
      "Q6P1X5",
      "Q6NZY7",
      "Q6N069",
      "Q6NZB1",
      "Q6NYC8",
      "P22670",
      "P22626",
      "Q15819")

weird = data.frame()
for(i in 1:length(y)) {
  weird = rbind(weird, biogrid[c(grep(y2[i],biogrid$V1), grep(y2[i],biogrid$V2)),])
}

weird = unique(weird)

weird_interactors = unique(interactions_to_interactors(weird))

weird_interactors_x = data.frame()
for(i in 1:length(y)) {
  weird_interactors_x = rbind(weird_interactors_x, weird_interactors[grep(y2[i],weird_interactors$interactor_IDs),])
}

weird_interactors_y = weird_interactors_x
for(i in 1:length(y)) {
  weird_interactors_y[weird_interactors_y$interactor_IDs == y2[i],5] = y[i]
}
colnames(weird_interactors_y)[5] = "actual taxonomy"

filename2 = "./biogrid/biogrid_weird_interactions.txt"
write.table(weird,filename2, col.names=T,row.names=F,sep="\t",quote=F)

filename2 = "./biogrid/biogrid_weird_interactors.txt"
write.table(weird_interactors_y,filename2, col.names=T,row.names=F,sep="\t",quote=F)


#Dirty stuff

biogrid_intersp <-unique(subset(biogrid_from_mentha,V10=="taxid:9606(Homo sapiens)"& V11!="taxid:9606(Homo sapiens)"))

biogrid_not_human_prots <- unique(select(biogrid_intersp,upac=V2,taxid=V11))

check<-unique(biogrid_not_human_prots$upac)

biogrid_human_prots <- unique(select(biogrid_intersp,upac=V1,taxid=V10))

check2<-unique(biogrid_human_prots$upac)
