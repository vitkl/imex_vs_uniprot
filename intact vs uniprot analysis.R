# intact_vs_uniprot analysis

# SPECIES_NAME = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "strain ATCC 204508", "strain K12", "Drosophila melanogaster", "Caenorhabditis elegans", "Schizosaccharomyces pombe", "Arabidopsis thaliana", "Tetrahymena thermophila", "Neurospora crassa", "Chlamydomonas reinhardtii", "Gallus gallus", "Danio rerio")

SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "strain K12", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")

## "Heterocephalus glaber" - naked mole rat - 0 interactions
## "strain K12" is E.coli
## "strain ATCC 204508 / S288c" is Baker's yeast

date = Sys.Date()
## Please specify the date for which you want to perform analysis (if not today)
date = as.Date("2016-12-01")


## ## Use all Uniprot if reviewed == 1, only Swissprot data if reviewed == 2, 
## ## TrEMBL data if reviewed == 3
reviewed = c(1, 2)

## ## Distinguish between isoforms or use only generic Uniprot IDs: TRUE / FALSE?
isoforms = c(TRUE, FALSE)

## missing_proteins = TRUE => proteins missing protein evidence are shown
## missing_proteins = FALSE => proteins for which protein evidence exists are shown
missing_proteins = TRUE

## The code below queries databases, saves and processes results, 
## gives summary on how many interactors have Uniprot IDs or belong to the SPECIES_NAME
## generates table with 0 and 1 for the combination of SPECIES_NAME, reviewed, isoforms
source("intact_vs_uniprot.R")
for (i in isoforms) {
  for (r in reviewed) {
    for (n in SPECIES_NAME) {
      intact_vs_uniprot(SPECIES_NAME = n, reviewed = r, isoforms = i, missing_proteins = TRUE, date = date)
      message(paste("- ",n," -  reviewed", r, " -  isoforms included", i))
    }
  }
}


## The code below does overlap comparisons and saves summaries
source("intact_vs_uniprot_overlap.R")
source("SPECIES_NAME_TO_ID.R")
for (r in reviewed) {
  for (i in isoforms) {
    for (n in SPECIES_NAME) {
      SPECIES_IDs = SPECIES_NAME_TO_ID(n)
      SPECIES_ID = SPECIES_IDs$SPECIES_ID
      intact_vs_uniprot_overlap(SPECIES_NAME = n, SPECIES_ID = SPECIES_ID, reviewed = r, isoforms = i, date)
      message(paste("- ",n," -  reviewed", r, " -  isoforms included", i))
    }
  }
}

## to try if works
#intact_vs_uniprot(SPECIES_NAME = "strain K12", reviewed = 1, isoforms = TRUE)

### combine the summaries for multiple species
proteome_vs_interactome_summary.all = data.frame()
whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all = data.frame()
whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all = data.frame()
for (r in reviewed) {
  for (i in isoforms) {
    for(n in SPECIES_NAME) {
      #=====================================
      filename.summary = paste("./summaries/","proteome_vs_interactome_summary_",n,"_reviewed_",r,"_isoforms_",i,"_",date,".txt", sep = "")
      if(!file.exists(filename.summary)) {
        proteome_vs_interactome_summary.all
      }
      if(file.exists(filename.summary)) {
        proteome_vs_interactome_summary_o = as.data.frame(read.delim(filename.summary, stringsAsFactors = F))
        proteome_vs_interactome_summary.all = rbind(proteome_vs_interactome_summary.all, proteome_vs_interactome_summary_o)
      }
      #=====================================
      filename.summary2 = paste("./summaries/","whole_proteome_Uniprot_vs_IMEx_vs_protein.exist_summary_",n,"_reviewed_",r,"_isoforms_",i,"_",date,".txt", sep = "")
      if(!file.exists(filename.summary2)) {
        whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all
      }
      if(file.exists(filename.summary2)) {
        whole_proteome_Uniprot_vs_IMEx_vs_protein.exist = as.data.frame(read.delim(filename.summary2, stringsAsFactors = F))
        whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all = rbind(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all, whole_proteome_Uniprot_vs_IMEx_vs_protein.exist)
      }
      #=====================================
      filename.summary3 = paste("./summaries/","whole_proteome_Uniprot_vs_IntAct_vs_protein.exist_summary_",n,"_reviewed_",r,"_isoforms_",i,"_",date,".txt", sep = "")
      if(!file.exists(filename.summary3)) {
        whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all
      }
      if(file.exists(filename.summary3)) {
        whole_proteome_Uniprot_vs_IntAct_vs_protein.exist = as.data.frame(read.delim(filename.summary3, stringsAsFactors = F))
        whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all = rbind(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all, whole_proteome_Uniprot_vs_IntAct_vs_protein.exist)
      }
  }
}
}
## save all species summary
filename.summary.all = paste("./results/","proteome_vs_interactome_summary_",date,".txt", sep = "")
write.table(proteome_vs_interactome_summary.all, filename.summary.all, col.names=T,row.names=F,sep="\t",quote=F)
filename.summary.all2 = paste("./results/","whole_proteome_Uniprot_vs_IMEx_vs_protein.exist_summary_",date,".txt", sep = "")
write.table(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all, filename.summary.all2, col.names=T,row.names=F,sep="\t",quote=F)
filename.summary.all3 = paste("./results/","whole_proteome_Uniprot_vs_IntAct_vs_protein.exist_summary_",date,".txt", sep = "")
write.table(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all, filename.summary.all3, col.names=T,row.names=F,sep="\t",quote=F)

### ======================================================================== ###

### combine the summaries for multiple species of how many interactors have non-uniprot and non-given-species identifiers
interactome_identifiers_summary.all = data.frame()
  for (i in isoforms) {
    for(n in SPECIES_NAME) {
      filename.summary = paste("./summaries/","uniprotKB_IDs_and_",n,"_interactors_summary_isoforms_",i,"_",date,".txt", sep = "")
      if(!file.exists(filename.summary)) {
        interactome_identifiers_summary.all
      }
      if(file.exists(filename.summary)) {
        interactome_identifiers_summary.o = as.data.frame(read.delim(filename.summary, stringsAsFactors = F))
        interactome_identifiers_summary.all = rbind(interactome_identifiers_summary.all, interactome_identifiers_summary.o)
      }
    }
  }
interactome_identifiers_summary.all = unique(interactome_identifiers_summary.all)
## save all species summary
filename.summary.all = paste("./results/","interactome_identifiers_summary_",date,".txt", sep = "")
write.table(interactome_identifiers_summary.all, filename.summary.all, col.names=T,row.names=F,sep="\t",quote=F)

### ======================================================================== ###
## Plotting the number of interactors which have UniprotKB identifiers and belong to the species queried
interactome_identifiers_summary.all.s = interactome_identifiers_summary.all[,c(3,4,5,6,7)]
## Transform table for plotting with ggplot2
library(reshape2)
interactome_identifiers_summary.all.s.m= melt(data = interactome_identifiers_summary.all.s,
                                              id.vars = c("SPECIES_NAME", "SPECIES_ID"),
                                              variable.name = "decription",
                                              value.name = "number")
## rename yeast and E.coli
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("strain K12", "E. coli, strain K12", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("Caenorhabditis elegans", "C. elegans", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("Drosophila melanogaster", "D. melanogaster", interactome_identifiers_summary.all.s.m$SPECIES_NAME)

# plot
library(ggplot2)
library(dplyr)

proteome_vs_interactome_plot <- ggplot(interactome_identifiers_summary.all.s.m, aes(x=SPECIES_NAME, y=number, fill=decription)) + geom_bar(width = 0.9, stat = "identity", position = "stack") + geom_label(aes(label=number), position = "stack", size = 4, label.padding = unit(0.08, "lines")) +
  ggtitle("The number of interactors which have UniprotKB identifiers and belong to the species queried")# + theme(text=element_text(size=13,  family="serif"))
proteome_vs_interactome_plot
# save plot
filename=paste("./results/", "interactome_identifiers_summary_plot_SMALL",date,".png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 12)
### ======================================================================== ###



### ======================================================================== ###
## Plotting with venn.diagram - Uniprot vs IMEx (or IntAct)
# code to be modified for different graphs
SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "strain K12", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")
reviewed_venn = 2                 # 1 or 2
IMExdatabase = "IMEx"    # "X0469.IntAct." or "IMEx"
### ======================================================================== ###
# non-modifyable code
proteome_vs_interactome_summary.all$species.name = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", proteome_vs_interactome_summary.all$species.name)
proteome_vs_interactome_summary.all$species.name = gsub("strain K12", "Escherichia coli, strain K12", proteome_vs_interactome_summary.all$species.name)
SPECIES_NAME = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", SPECIES_NAME)
SPECIES_NAME = gsub("strain K12", "Escherichia coli, strain K12", SPECIES_NAME)
library(dplyr)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all, reviewed == reviewed_venn)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, isoforms == TRUE)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, database == IMExdatabase)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all, reviewed == reviewed_venn)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, database == IMExdatabase)

library(VennDiagram)
grid.newpage()

if(reviewed_venn == 2){plotname = paste0("Proteome (SwissProt) coverage by interaction databases(",IMExdatabase,")")
databasename = "SwissProt"}
if(reviewed_venn == 1){plotname = paste0("Proteome (all UniprotKB) coverage by interaction databases(",IMExdatabase,")")
databasename = "UniprotKB"}

pushViewport(viewport(layout=grid.layout(nrow = length(SPECIES_NAME)+2, ncol=4, widths = unit(c(2/7,2/7,2/7,1/7), "npc"), 
                                         heights = unit(c(1/((length(SPECIES_NAME)+1)*4),1/((length(SPECIES_NAME)+1)*4),rep(1/(length(SPECIES_NAME)+1/2),length(SPECIES_NAME))), "npc"))))

pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(0.7, "npc"),y= unit(0.5, "npc"))
popViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("isoforms included", x = unit(0.5, "npc"),y= unit(0.5, "npc"))
popViewport()
pushViewport(viewport(layout.pos.col=3, layout.pos.row = 2))
x =grid.text("isoforms excluded", x = unit(0.5, "npc"),y= unit(0.5, "npc"))
popViewport()

for (i in 1:length(SPECIES_NAME)) {
  
  pushViewport(viewport(layout.pos.col=1, layout.pos.row = i+2))
  x =grid.text(SPECIES_NAME[i], x = unit(0.5, "npc"),y= unit(0.5, "npc"))
  popViewport()
  
  pushViewport(viewport(layout.pos.col=2, layout.pos.row = i+2))
  source("double.venn.prot.R")
  venn = double.venn.prot(proteome_vs_interactome_summary.all_isof[i,], scaled = FALSE, IMExdatabase = IMExdatabase)
  popViewport()
  
  pushViewport(viewport(layout.pos.col=3, layout.pos.row = i+2))
  source("double.venn.prot.R")
  venn = double.venn.prot(proteome_vs_interactome_summary.all_noisof[i,], scaled = FALSE, IMExdatabase = IMExdatabase)
  popViewport()
}

## code below generates filename
paste("whole_proteome_", databasename,"_vs_interactome_",IMExdatabase,"_vennD_", date, sep = "")

### ======================================================================== ###



### ======================================================================== ###
## Plotting overlaps -Uniprot -IMEx -protein evidence- without numbers
# code to be modified for different graphs
SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "strain K12", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")
IMExdatabase = "IMEx"    # "X0469.IntAct." or "IMEx"
## missing_proteins = TRUE => proteins missing protein evidence are shown
## missing_proteins = FALSE => proteins for which protein evidence exists are shown
missing_proteins = TRUE
### ======================================================================== ###
# non-modifyable code
SPECIES_NAME = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", SPECIES_NAME)
SPECIES_NAME = gsub("strain K12", "Escherichia coli, strain K12", SPECIES_NAME)

if(IMExdatabase == "IMEx"){
whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$species_name = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$species_name)
whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$species_name = gsub("strain K12", "Escherichia coli, strain K12", whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$species_name)
library(dplyr)
if(missing_proteins == TRUE){whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$C = rep("missing_protein_evidence",length(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$C))}
proteome_vs_interactome_summary.all_isof = filter(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all, reviewed == 1)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, isoforms == FALSE)
proteome_vs_interactome_summary.all_noisof = filter(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all, reviewed == 2)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
}
if(IMExdatabase == "X0469.IntAct."){
whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$species_name = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$species_name)
whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$species_name = gsub("strain K12", "Escherichia coli, strain K12", whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$species_name)
library(dplyr)
if(missing_proteins == TRUE){whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$C = rep("missing_protein_evidence",length(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all$C))}
proteome_vs_interactome_summary.all_isof = filter(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all, reviewed == 1)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, isoforms == FALSE)
proteome_vs_interactome_summary.all_noisof = filter(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist.all, reviewed == 2)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
}
library(VennDiagram)
grid.newpage()

if(missing_proteins == TRUE){plotname = paste0("Proteome coverage by interaction databases (", IMExdatabase,"), overlap with missing evidence at the protein level (from Uniprot)")
databasename = "UniprotKB"}
if(missing_proteins == FALSE){plotname = paste0("Proteome coverage by interaction databases (", IMExdatabase,"), overlap with evidence at the protein level (from Uniprot)")
databasename = "UniprotKB"}

pushViewport(viewport(layout=grid.layout(nrow = length(SPECIES_NAME)+2, ncol=4, widths = unit(c(2/7,2/7,2/7,1/7), "npc"), 
                                         heights = unit(c(1/((length(SPECIES_NAME)+1)*4),1/((length(SPECIES_NAME)+1)*4),rep(1/(length(SPECIES_NAME)+1/2),length(SPECIES_NAME))), "npc"))))

pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(0.7, "npc"),y= unit(0.5, "npc"))
popViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("UniprotKB, isoforms excluded", x = unit(0.5, "npc"),y= unit(0.5, "npc"))
popViewport()
pushViewport(viewport(layout.pos.col=3, layout.pos.row = 2))
x =grid.text("SwissProt, isoforms excluded", x = unit(0.5, "npc"),y= unit(0.5, "npc"))
popViewport()

for (i in 1:length(SPECIES_NAME)) {
  
  pushViewport(viewport(layout.pos.col=1, layout.pos.row = i+2))
  x =grid.text(SPECIES_NAME[i], x = unit(0.5, "npc"),y= unit(0.5, "npc"))
  popViewport()
  
  pushViewport(viewport(layout.pos.col=2, layout.pos.row = i+2))
  source("triple.venn.prot.R")
  venn = triple.venn.prot(proteome_vs_interactome_summary.all_isof[i,], scaled = TRUE)
  popViewport()
  
  pushViewport(viewport(layout.pos.col=3, layout.pos.row = i+2))
  source("triple.venn.prot.R")
  venn = triple.venn.prot(proteome_vs_interactome_summary.all_noisof[i,], scaled = TRUE)
  popViewport()
}

## code below generates filename
# missing_proteins = TRUE => shows proteins missing protein evidence, 
# missing_proteins = FALSE => shows proteins having protein evidence
paste("whole_proteome_vs_interactome_",IMExdatabase,"_vs_missing_protein_evidence_isoforms_excluded_vennD_", date, sep = "")

### ======================================================================== ###
library(dplyr)
proteome_vs_interactome_summary.all_noisof = filter(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all, reviewed == 2)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all.spl = split(proteome_vs_interactome_summary.all_noisof, whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all$species_name)
fisher.test.missing.proteins = function(y){
z = data.frame(with_protein_evidence = c(NA, NA), without_protein_evidence = c(NA, NA), row.names = c("interactome", "no_interactome"))
z$with_protein_evidence[1] = y$whole_proteome_Uniprot_and_IMEx_not_missing_protein_evidence
z$with_protein_evidence[2] = y$whole_proteome_Uniprot_exclusively
z$without_protein_evidence[1] = y$IMEx_and_missing_protein_evidence
z$without_protein_evidence[2] = y$whole_proteome_Uniprot_and_missing_protein_evidence_not_IMEx
z = t(z)
fisher.test(z)
chisq.test(z)
return(fisher.test(z)$conf.int)
}
print("Confidence interval of odds ratio (exact Fisher test)")
print("if a protein doesnt have evidence at the protein level it is more likely to be missing in IntAct")
sapply(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist.all.spl, fisher.test.missing.proteins)
