## swissprot_vs_imex_interaction_properties
SPECIES_NAME = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "strain K12")

SPECIES_NAME = "Homo sapiens"

## ## Use all Uniprot if reviewed == 1, only Swissprot data if reviewed == 2, 
## ## TrEMBL data if reviewed == 3
## only reviewed = 2 is relevant for this analysis
reviewed = 2
## ## Distinguish between isoforms or use only generic Uniprot IDs: TRUE / FALSE?
## only isoforms = FALSE is relevant for this analysis
isoforms = FALSE

date = Sys.Date()
## Please specify the date for which you want to perform analysis (if not today)
date = as.Date("2016-12-01")
print(date)
##============================================================================##

source("SPECIES_NAME_TO_ID.R")
library(dplyr)
#for (r in reviewed) {
#  for (i in isoforms) {
#    for (n in SPECIES_NAME) {
r = reviewed
i = isoforms
n = SPECIES_NAME
##============================================================================##
SPECIES_IDs = SPECIES_NAME_TO_ID(n)
SPECIES_ID = SPECIES_IDs$SPECIES_ID
##============================================================================##
# interations from PSICQUIC
## detmethod:MI:0232 - transcriptional complementation assay, includes different two-hybrids
# http://www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0232
## detmethod:MI:0090 - protein complementation assay, includes transcriptional 
# complementation assay and bimolecular fluorescence complementation, enzyme complementations
# http://www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0090

## both detmethod and pmethod are required for specifying AP-MS
## detmethod:MI:0004 - affinity chromatography technology, co-IP and pulldowns
# http://www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0004
## pmethod:MI:0433 - partial identification of protein sequence
# http://www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0433

## type in exact database names (the list below is default for the function)
databases <- c("IntAct", "MINT", "bhf-ucl", "MPIDB", "MatrixDB", 
               "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo")
source("query_PSICQUIC_for_interactions.R")
twohybrids_all_interactions = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, 
                                                   SPECIES_NAME = SPECIES_NAME, 
                                                   databases = databases, date,
                                                   detmethod = "transcriptional complementation assay")
complementation_all_interactions = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, 
                                                              SPECIES_NAME = SPECIES_NAME, 
                                                              databases = databases, date,
                                                              detmethod = "protein complementation assay")
ap_ms_all_interactions = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, 
                                                         SPECIES_NAME = SPECIES_NAME, 
                                                         databases = databases, date,
                                                         detmethod = "affinity chromatography technology",
                                                         pmethod = "partial identification of protein sequence")
##============================================================================##
## Transforming data
## the function extracts interactor IDs from interactions
source("interactions_to_interactors.R")
twohybrids_all_interactors = interactions_to_interactors(twohybrids_all_interactions)
complementation_all_interactors = interactions_to_interactors(complementation_all_interactions)
ap_ms_all_interactors = interactions_to_interactors(ap_ms_all_interactions)

## filter interactor for uniprotkb only indentifiers 
## filter for SPECIES_ID only proteins
source("uniprotkb_and_SPECIES_ID_interactor_selector.R")
twohybrids_all_interactors_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(twohybrids_all_interactors, SPECIES_ID)
complementation_all_interactors_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(complementation_all_interactors, SPECIES_ID)
ap_ms_all_interactors_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(ap_ms_all_interactors, SPECIES_ID)
## In case isoform argument to the function is FALSE - REMOVE All isoform IDs (XXXXXX-X+) from IDs
if(isoforms == FALSE){
  source("isoform_id_all_remover.R")
  twohybrids_all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(twohybrids_all_interactors_SPECIES_ID_only$interactor_IDs)
  complementation_all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(complementation_all_interactors_SPECIES_ID_only$interactor_IDs)
  ap_ms_all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(ap_ms_all_interactors_SPECIES_ID_only$interactor_IDs)
}
##============================================================================##
## preparing for logic table: selecting unique and adding the column of ones
unique_twohybrids_interactors_SPECIES_ID_only = unique(cbind(twohybrids_all_interactors_SPECIES_ID_only[c("interactor_IDs")], 1))
colnames(unique_twohybrids_interactors_SPECIES_ID_only)[2] = "two_hybrid"
unique_complementation_interactors_SPECIES_ID_only = unique(cbind(complementation_all_interactors_SPECIES_ID_only[c("interactor_IDs")], 1))
colnames(unique_complementation_interactors_SPECIES_ID_only)[2] = "all_protein_complementation"
unique_ap_ms_interactors_SPECIES_ID_only = unique(cbind(ap_ms_all_interactors_SPECIES_ID_only[c("interactor_IDs")], 1))
colnames(unique_ap_ms_interactors_SPECIES_ID_only)[2] = "AP_MS"
##============================================================================##
## loading saved logic table from swissprot_vs_imex_protein_properties
filename_vs_2 = paste("./analysis/","proteome_vs_interactome_protein_properties_f_", n,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")
proteome_vs_imex_details_f = as.data.frame(read.delim(filename_vs_2, header = T, stringsAsFactors = F,quote=""))
proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx = factor(proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx, ordered =F)
##============================================================================##
## merging new results to the logic table
proteome_vs_imex_interaction_details_t1 = merge(proteome_vs_imex_details_f, 
                                   unique_twohybrids_interactors_SPECIES_ID_only, 
                                   by.x = "whole_proteome_IDs",
                                   by.y = "interactor_IDs",
                                   all.x = T, all.y = F)
proteome_vs_imex_interaction_details_t2 = merge(proteome_vs_imex_interaction_details_t1, 
                                               unique_complementation_interactors_SPECIES_ID_only, 
                                               by.x = "whole_proteome_IDs",
                                               by.y = "interactor_IDs",
                                               all.x = T, all.y = F)
proteome_vs_imex_interaction_details_f = merge(proteome_vs_imex_interaction_details_t2, 
                                               unique_ap_ms_interactors_SPECIES_ID_only, 
                                                by.x = "whole_proteome_IDs",
                                                by.y = "interactor_IDs",
                                                all.x = T, all.y = F)
proteome_vs_imex_interaction_details_f[is.na(proteome_vs_imex_interaction_details_f)] = 0
##============================================================================##
## adding some variables
proteome_vs_imex_interaction_details_f[,length(proteome_vs_imex_interaction_details_f)+1] = interaction(proteome_vs_imex_interaction_details_f$two_hybrid, proteome_vs_imex_interaction_details_f$AP_MS, sep = "_")
colnames(proteome_vs_imex_interaction_details_f)[length(proteome_vs_imex_interaction_details_f)] = paste0("two_hybrid", "_vs_","AP_MS")
levels(proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS) = c("not_two_hybrid_and_not_AP_MS", "two_hybrid_not_AP_MS","not_two_hybrid_but_AP_MS", "two_hybrid_and_AP_MS")
proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS = as.character(proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS)
proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS[proteome_vs_imex_interaction_details_f$IMEx!=1] = "not_in_IMEx"

##============================================================================##
## saving combined logic table + protein properties from Uniprot
filename = paste("./analysis/","proteome_vs_interactome_interaction_properties_f_", n,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")
write.table(proteome_vs_imex_interaction_details_f,filename,col.names=T,row.names=F,sep="\t",quote=F)
##============================================================================##
## Analysis
# read the saved table
filename = paste("./analysis/","proteome_vs_interactome_interaction_properties_f_", n,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")
proteome_vs_imex_interaction_details_f = as.data.frame(read.delim(filename, header = T, stringsAsFactors = F,quote=""))
proteome_vs_imex_interaction_details_f$whole_proteome_Uniprot_IMEx = factor(proteome_vs_imex_interaction_details_f$whole_proteome_Uniprot_IMEx, ordered =F)
proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS = factor(proteome_vs_imex_interaction_details_f$two_hybrid_vs_AP_MS, ordered =F)

## The density and the histogram of the protein mass (overlay for different levels)
library(ggplot2)
ggplot(proteome_vs_imex_interaction_details_f, aes(x = Mass, color =two_hybrid_vs_AP_MS, fill = two_hybrid_vs_AP_MS)) +geom_density(alpha =0.2) + scale_x_log10() #+facet_grid(two_hybrid_vs_AP_MS~Organism)
ggplot(proteome_vs_imex_interaction_details_f, aes(x = Mass,color = two_hybrid_vs_AP_MS)) + scale_x_log10() + geom_histogram(position = "identity", bins = 50,alpha =0.1) #+facet_grid(two_hybrid_vs_AP_MS~Organism)

## Vioilin plot of the protein mass
library(dplyr)
library(scales) 
## adding two_hybrid and AP_MS levels to show non-combinatory relationship (not 
## excluding proteins present in both)
two_hybrid = filter(proteome_vs_imex_interaction_details_f, two_hybrid == 1)
two_hybrid$two_hybrid_vs_AP_MS = "two_hybrid"
AP_MS = filter(proteome_vs_imex_interaction_details_f, AP_MS == 1)
AP_MS$two_hybrid_vs_AP_MS = "AP_MS"
for_plot = rbind(proteome_vs_imex_interaction_details_f, two_hybrid, AP_MS)
for_plot$two_hybrid_vs_AP_MS = factor(for_plot$two_hybrid_vs_AP_MS, ordered = F)
levels(for_plot$two_hybrid_vs_AP_MS)
## Calculating protein median protein mass for each group
yy = split(for_plot, for_plot$two_hybrid_vs_AP_MS)
Mass_median = sapply(yy, function(x){median(log10(x$Mass))})
## Violin plot
ggplot(for_plot, aes(y = Mass, x =two_hybrid_vs_AP_MS, fill = two_hybrid_vs_AP_MS)) +
  geom_violin(draw_quantiles = c(0.05,0.25,0.495,0.5,0.505,0.75,0.95),scale = "count", alpha =0.7) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_abline(slope = 0, intercept = Mass_median, alpha =0.1) +
  ylab("Mass, Da") #+facet_grid(two_hybrid_vs_AP_MS~Organism)

