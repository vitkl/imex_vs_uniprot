## Draws wenn diagrams for whole proteome, IMEx, BioGRID for each species in table

double.venn.prot = function(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, scaled = TRUE, IMExdatabase = "IMEx"){

if(scaled == TRUE){
venn.d = draw.pairwise.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..Uniprot., 
                          area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..interactome.available +
                            whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$interactome..but.not.in.Uniprot, 
                          cross.area = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..interactome.available, 
                          category = c("whole_proteome", 
                                       IMExdatabase), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(0.9,0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 0
)
}

if(scaled == FALSE){
  venn.d = draw.pairwise.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..Uniprot., 
                              area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..interactome.available +
                                whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$interactome..but.not.in.Uniprot, 
                              cross.area = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole.proteome..interactome.available, 
                              category = c("whole_proteome", 
                                           IMExdatabase), 
                              lty = rep("blank", 2), 
                              fill = c("blue", "red"), 
                              alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                              cat.dist = rep(0.035, 2), 
                              cat.cex = c(0.9,0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05
    )
}

}