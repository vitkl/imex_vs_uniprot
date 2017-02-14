

### ======================================================================== ###
## Plotting overlaps with IntAct- without numbers

SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "Drosophila melanogaster", "Caenorhabditis elegans", "Arabidopsis thaliana")
reviewed = 2
library(dplyr)
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_isof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha, reviewed == 2)
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_isof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_isof, isoforms == TRUE)
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_noisof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha, reviewed == 2)
whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_noisof = filter(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_noisof, isoforms == FALSE)

library(VennDiagram)
grid.newpage()

if(reviewed == 2){plotname = "Proteome (SwissProt) coverage by interaction databases, overlap between IntAct and BioGRID"
databasename = "SwissProt"}
if(reviewed == 1){plotname = "Proteome (all UniprotKB) coverage by interaction databases, overlap between IntAct and BioGRID"
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
  source("triple.venn.prot.intact.R")
  venn = triple.venn.prot.intact(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_isof[i,])
  popViewport()
  
  pushViewport(viewport(layout.pos.col=3, layout.pos.row = i+2))
  source("triple.venn.prot.intact.R")
  venn = triple.venn.prot.intact(whole_proteome_Uniprot_vs_IntAct_vs_BioGRID_from_Mentha_noisof[i,])
  popViewport()
}

paste("whole_proteome_vs_interactome_IntAct_vs_BioGRID_f_vennD_reviewed_",reviewed = databasename,"_", Sys.Date(), sep = "")
### ======================================================================== ###