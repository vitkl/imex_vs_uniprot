## Draws wenn diagrams for reference proteome, IMEx, BioGRID for each species in table

triple.venn.prot.ref = function(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, scaled = TRUE){
  
if(scaled == TRUE){
venn.d = draw.triple.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_reference_proteome_Uniprot, 
                          area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_IMEx, 
                          area3 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_BioGRID_from_Mentha, 
                          n12 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_IMEx, 
                          n23 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$IMEx_and_BioGRID_from_Mentha,
                          n13 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_BioGRID_from_Mentha, 
                          n123 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_IMEx_and_BioGRID_from_Mentha, 
                          category = c(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$A, 
                                       whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$B, 
                                       whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$C), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(0, 45, 160), 
                          cat.dist = rep(0.035, 3), 
                          cat.cex = c(0.9,0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 0
)
}

if(scaled == FALSE){
    venn.d = draw.triple.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_reference_proteome_Uniprot, 
                              area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_IMEx, 
                              area3 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_BioGRID_from_Mentha, 
                              n12 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_IMEx, 
                              n23 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$IMEx_and_BioGRID_from_Mentha,
                              n13 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_BioGRID_from_Mentha, 
                              n123 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$reference_proteome_Uniprot_and_IMEx_and_BioGRID_from_Mentha, 
                              category = c(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$A, 
                                           whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$B, 
                                           whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$C), 
                              lty = rep("blank", 3), 
                              fill = c("blue", "red", "green"), 
                              alpha = rep(0.5, 3), cat.pos = c(0, 45, 160), 
                              cat.dist = rep(0.035, 3), 
                              cat.cex = c(0.9,0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05
    )
  }
  
}