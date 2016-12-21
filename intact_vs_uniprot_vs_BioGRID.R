# intact_vs_uniprot_vs_BioGRID analysis

## Enter SPECIES_NAME 
SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")
# SPECIES_NAME = c("Caenorhabditis elegans")
## !! no "E.coli strain K12" in BioGRID

## ## Use all Uniprot if reviewed == 1, only Swissprot data if reviewed == 2, 
## ## TrEMBL data if reviewed == 3
reviewed = c(1, 2)
## ## Distinguish between isoforms or use only generic Uniprot IDs: TRUE / FALSE?
isoforms = c(FALSE) # not possible to distinguish isoforms for BioGRID

date = Sys.Date()
## Please specify the date for which you want to perform analysis (if not today)
date = as.Date("2016-12-01")


#============================================================================#
source("biogrid_from_mentha_vs_proteome_vs_imex.R")
for (r in reviewed) {
  for (i in isoforms) {
    for (n in SPECIES_NAME) {
        biogrid_from_mentha_vs_proteome_vs_imex(SPECIES_NAME = n, reviewed = r, isoforms = i, date = date)
    }
  }
}

whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha = data.frame()
reference_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha = data.frame()
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha = data.frame()
source("SPECIES_NAME_TO_ID.R")
for (r in reviewed) {
  for (i in isoforms) {
    for (n in SPECIES_NAME) {
      SPECIES_IDs = SPECIES_NAME_TO_ID(n)
      SPECIES_ID = SPECIES_IDs$SPECIES_ID
      filename_vs_3 = paste("./analysis/","proteome_vs_interactome_vs_BioGRID_f_", SPECIES_ID,"_reviewed_",reviewed = r,"_isoforms_",isoforms = i,"_", date,".txt", sep = "")
      biogrid_from_mentha_vs_proteome_vs_imex_f = as.data.frame(read.delim(filename_vs_3, stringsAsFactors = F))

#============================================================================#
## Calculating overlaps

source("A_vs_B_vs_C_overlap.R")
whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha_temp = A_vs_B_vs_C_overlap(biogrid_from_mentha_vs_proteome_vs_imex_f, 
                    A = "whole_proteome_Uniprot", 
                    B = "IMEx", 
                    C = "BioGRID_from_Mentha", 
                    SPECIES_NAME = n, reviewed = r, isoforms = i)
whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha = rbind(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha_temp)

whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_temp = A_vs_B_vs_C_overlap(biogrid_from_mentha_vs_proteome_vs_imex_f, 
                                                                                 A = "whole_proteome_Uniprot", 
                                                                                 B = "X0469.IntAct.", 
                                                                                 C = "BioGRID_from_Mentha", 
                                                                                 SPECIES_NAME = n, reviewed = r, isoforms = i)
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha = rbind(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha, whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_temp)

reference_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha_temp = A_vs_B_vs_C_overlap(biogrid_from_mentha_vs_proteome_vs_imex_f, 
                                                                                 A = "reference_proteome_Uniprot", 
                                                                                 B = "IMEx", 
                                                                                 C = "BioGRID_from_Mentha", 
                                                                                 SPECIES_NAME = n, reviewed = r, isoforms = i)
reference_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha = rbind(reference_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, reference_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha_temp)
    }
  }
}

### ======================================================================== ###
## Plotting overlaps - Uniprot IMEx BioGRID
SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")
IMExdatabase = "X0469.IntAct."    # "X0469.IntAct." or "IMEx"
### ======================================================================== ###
# non-modifyable code
if(IMExdatabase == "IMEx"){
whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$species_name = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$species_name)
library(dplyr)
proteome_vs_interactome_summary.all_isof = filter(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, reviewed == 1)
proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, isoforms == FALSE)
proteome_vs_interactome_summary.all_noisof = filter(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, reviewed == 2)
proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
}
if(IMExdatabase == "X0469.IntAct."){
  whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha$species_name = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha$species_name)
  library(dplyr)
  proteome_vs_interactome_summary.all_isof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha, reviewed == 1)
  proteome_vs_interactome_summary.all_isof = filter(proteome_vs_interactome_summary.all_isof, isoforms == FALSE)
  proteome_vs_interactome_summary.all_noisof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha, reviewed == 2)
  proteome_vs_interactome_summary.all_noisof = filter(proteome_vs_interactome_summary.all_noisof, isoforms == FALSE)
}
library(VennDiagram)
grid.newpage()

plotname = paste0("Proteome coverage by interaction databases, overlap between ", IMExdatabase," and BioGRID")

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

paste("whole_proteome_vs_interactome_",IMExdatabase,"_vs_BioGRID_f_vennD_", date, sep = "")

### ======================================================================== ###
### ======================================================================== ###

### combine the summaries for multiple species of how many interactors have non-uniprot identifiers
interactome_identifiers_summary.all = data.frame()
for (i in isoforms) {
  for(n in SPECIES_NAME) {
    filename.summary = paste("./summaries/","uniprotKB_IDs_and_",n,"_biogrid_from_mentha_interactors_summary", "_isoforms_",i,"_", date,".txt", sep = "")
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
filename.summary.all = paste("./results/","interactome_identifiers_BioGRID_summary_",date,".txt", sep = "")
write.table(interactome_identifiers_summary.all, filename.summary.all, col.names=T,row.names=F,sep="\t",quote=F)

#################
interactome_identifiers_summary.all.s = interactome_identifiers_summary.all[,c(3,4,5,6,7)]
## Transform table for plotting with ggplot2
library(reshape2)
interactome_identifiers_summary.all.s.m= melt(data = interactome_identifiers_summary.all.s,
                                              id.vars = c("SPECIES_NAME", "SPECIES_ID"),
                                              variable.name = "decription",
                                              value.name = "number")
interactome_identifiers_summary.all.s.m=unique(interactome_identifiers_summary.all.s.m)
## rename yeast and E.coli
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("strain ATCC 204508", "S. cerevisiae, strain S288c", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("strain K12", "E. coli, strain K12", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("Caenorhabditis elegans", "C. elegans", interactome_identifiers_summary.all.s.m$SPECIES_NAME)
interactome_identifiers_summary.all.s.m$SPECIES_NAME = gsub("Drosophila melanogaster", "D. melanogaster", interactome_identifiers_summary.all.s.m$SPECIES_NAME)

# plot
library(ggplot2)
library(dplyr)

proteome_vs_interactome_plot <- ggplot(interactome_identifiers_summary.all.s.m, aes(x=SPECIES_NAME, y=number, fill=decription,label=number)) + geom_bar(width = 0.9, stat = "identity", position = "stack") + geom_label(position = "stack", size = 4, label.padding = unit(0.08, "lines")) +
  ggtitle("The number of interactors which have UniprotKB identifiers and belong to the species queried") + theme(text=element_text(size=13,  family="serif"))
proteome_vs_interactome_plot
# save plot
filename=paste("./results/", "interactome_identifiers_summary_plot_SMALL",date,".png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 12)

### ======================================================================== ###
