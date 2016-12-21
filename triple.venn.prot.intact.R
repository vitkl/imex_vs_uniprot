## Draws wenn diagrams for whole proteome, IMEx, BioGRID for each species in table

triple.venn.prot.intact = function(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha, scaled = TRUE){

if(scaled == TRUE){
venn.d = draw.triple.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_whole_proteome_Uniprot, 
                          area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_X0469.IntAct., 
                          area3 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_BioGRID_from_Mentha, 
                          n12 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_X0469.IntAct., 
                          n23 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$X0469.IntAct._and_BioGRID_from_Mentha,
                          n13 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_BioGRID_from_Mentha, 
                          n123 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_X0469.IntAct._and_BioGRID_from_Mentha, 
                          category = c(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$A, 
                                       whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$B, 
                                       whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$C), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(0, 45, 160), 
                          cat.dist = rep(0.035, 3), 
                          cat.cex = c(0.9,0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          #direct.area = TRUE,
                          cex = 0.8
)
}

if(scaled == FALSE){
    venn.d = draw.triple.venn(area1 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_whole_proteome_Uniprot, 
                              area2 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_X0469.IntAct., 
                              area3 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$count_BioGRID_from_Mentha, 
                              n12 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_X0469.IntAct., 
                              n23 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$X0469.IntAct._and_BioGRID_from_Mentha,
                              n13 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_BioGRID_from_Mentha, 
                              n123 = whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$whole_proteome_Uniprot_and_X0469.IntAct._and_BioGRID_from_Mentha, 
                              category = c(whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$A, 
                                           whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$B, 
                                           whole_proteome_Uniprot_vs_IMEx_vs_BioGRID_from_Mentha$C), 
                              lty = rep("blank", 3), 
                              fill = c("blue", "red", "green"), 
                              alpha = rep(0.5, 3), cat.pos = c(0, 45, 160), 
                              cat.dist = rep(0.035, 3), 
                              cat.cex = c(0.9,0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                              cex = 0.8
    )
}

}