intact_vs_uniprot = function(SPECIES_NAME, reviewed, isoforms, missing_proteins = TRUE, date = Sys.Date(), databases) {
  ## The function below downloads reference proteome list from Uniprot, saves it,
  ## finds and returns Proteome_ID and SPECIES_ID for given SPECIES_NAME

  source("SPECIES_NAME_TO_ID.R")
  SPECIES_IDs = SPECIES_NAME_TO_ID(SPECIES_NAME)
  
  SPECIES_ID = SPECIES_IDs$SPECIES_ID;
  Proteome_ID = SPECIES_IDs$Proteome_ID;
  
  ## Query PSICQUIC for interactions, get MI-TAB-2.5, save, return
  source("query_PSICQUIC_for_interactions.R")
  all_interactions = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, SPECIES_NAME = SPECIES_NAME, databases = databases, date)
  if(nrow(all_interactions) > 0){
  
  ## Download the reference proteome from Uniprot for SPECIES_ID, gunzip, read
  #source("download_reference_proteome.R")
  #reference_proteome_ftp = download_reference_proteome(SPECIES_ID, Proteome_ID)
  
  ## Query Uniprot for the reference proteome for SPECIES_ID
  source("download_reference_proteome_q.R")
  reference_proteome_query = download_reference_proteome_q(SPECIES_ID, Proteome_ID, date)
  
  ## Use only Swissprot data if reviewed == 2 
  if(reviewed == 2) {
    reference_proteome_query = dplyr::filter(reference_proteome_query, Status == "reviewed")
  }
  ## Use only TrEMBL data if reviewed == 3 
  if(reviewed == 3) {
    reference_proteome_query = dplyr::filter(reference_proteome_query, Status == "unreviewed")
  }
  if(nrow(reference_proteome_query) > 0){  
  
  ## Query Uniprot for the whole proteome for SPECIES_ID
  source("download_whole_proteome.R")
  whole_proteome_query = download_whole_proteome(SPECIES_ID, date)
  
  ## Use only Swissprot data if reviewed == 2
  if(reviewed == 2) {
    whole_proteome_query = dplyr::filter(whole_proteome_query, Status == "reviewed")
  }
  ## Use only TrEMBL data if reviewed == 3
  if(reviewed == 3) {
    whole_proteome_query = dplyr::filter(whole_proteome_query, Status == "unreviewed")
  }
  
  if(nrow(whole_proteome_query) > 0){
  
  #============================================================================#
  ## cleaning data
  
  ## In case isoform argument to the function is TRUE  - to extract isoform IDs 
  ## from queried proteomes and merge with generic IDs
  if(isoforms == TRUE){
  source("isoform_id_extractor_Uniprot.R")
  reference_proteome_isoforms = isoform_id_extractor_Uniprot(reference_proteome_query)
  reference_proteome = as.character(c(reference_proteome_isoforms, reference_proteome_query$Entry))
  whole_proteome_isoforms = isoform_id_extractor_Uniprot(whole_proteome_query)
  whole_proteome = as.character(c(whole_proteome_isoforms, whole_proteome_query$Entry))
  }
  
  ## In case isoform argument to the function is FALSE  - to extract generic IDs 
  ## from queried proteomes
  if(isoforms == FALSE){
    reference_proteome = reference_proteome_query$Entry
    whole_proteome = whole_proteome_query$Entry
  }
  
  ## the function extracts interactor IDs from interaction databases
  source("interactions_to_interactors.R")
  all_interactors = interactions_to_interactors(all_interactions)
  
  ## filter interactor for uniprotkb only indentifiers 
  ## filter for SPECIES_ID only proteins
  source("uniprotkb_and_SPECIES_ID_interactor_selector.R")
  all_interactors_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(all_interactors, SPECIES_ID)
  
  ## In case isoform argument to the function is TRUE - REMOVE "-1" XXXXXX-1 from IDs
  if(isoforms == TRUE){
    source("isoform_id_1_remover.R")
    all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_1_remover(all_interactors_SPECIES_ID_only$interactor_IDs)
    reference_proteome = isoform_id_1_remover(reference_proteome)
    whole_proteome = isoform_id_1_remover(whole_proteome)
  }
  
  ## In case isoform argument to the function is FALSE - REMOVE All isoform IDs (XXXXXX-X+) from IDs
  if(isoforms == FALSE){
    source("isoform_id_all_remover.R")
    all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(all_interactors_SPECIES_ID_only$interactor_IDs)
  }
  
  ## select unique IDs only
  reference_proteome_unique = data.frame(Reference_proteome_IDs = unique(reference_proteome),stringsAsFactors = F)
  whole_proteome_unique = data.frame(whole_proteome_IDs = unique(whole_proteome),stringsAsFactors = F)
  all_interactors_SPECIES_ID_only = unique(all_interactors_SPECIES_ID_only)
  
  ## unique interactors - excluding overlap between databases
  unique_interactors_SPECIES_ID_only = unique(all_interactors_SPECIES_ID_only[c("interactor_IDs","interactor_IDs_databases", "interactor_SPECIES_ID")])
  
  ## Print what have been loaded
  source("summary_intact_vs_uniprot1.R")
  summary_intact_vs_uniprot1 = summary_intact_vs_uniprot1(whole_proteome_unique, reference_proteome_unique, reference_proteome_ftp, all_interactors_SPECIES_ID_only)
  
  #============================================================================#
  ## proteins missing protein evidence
  library(dplyr)
  if(missing_proteins == FALSE){
  reference_proteome_existence = data.frame(Entry = reference_proteome_query$Entry, Protein.existence = reference_proteome_query$Protein.existence,stringsAsFactors = F)
  reference_proteome_existence = filter(reference_proteome_existence, Protein.existence == "Evidence at protein level")
  whole_proteome_existence = data.frame(Entry = whole_proteome_query$Entry, Protein.existence = whole_proteome_query$Protein.existence,stringsAsFactors = F)
  whole_proteome_existence = filter(whole_proteome_existence, Protein.existence == "Evidence at protein level")
  }
  if(missing_proteins == TRUE){
    reference_proteome_existence = data.frame(Entry = reference_proteome_query$Entry, Protein.existence = reference_proteome_query$Protein.existence,stringsAsFactors = F)
    reference_proteome_existence = filter(reference_proteome_existence, Protein.existence != "Evidence at protein level")
    whole_proteome_existence = data.frame(Entry = whole_proteome_query$Entry, Protein.existence = whole_proteome_query$Protein.existence,stringsAsFactors = F)
    whole_proteome_existence = filter(whole_proteome_existence, Protein.existence != "Evidence at protein level")
  }
  
  #============================================================================#
  ## Data transformations
  #============================================================================#
  ## take list of proteins in databases (IMEx and Uniprot) and convert into the logic table
  ### logic table has a 1 when protein(by rows) is present in a database (by columns)
  ### logic table has a 0 when protein(by rows) is absent in a database (by columns)
  source("proteome_vs_interactome_logic_table_generator.R")
  proteome_vs_interactome_f = proteome_vs_interactome_logic_table_generator(reference_proteome_unique, whole_proteome_unique, all_interactors_SPECIES_ID_only, unique_interactors_SPECIES_ID_only)
  #============================================================================#
  ## adding protein existence to the logic table
  if(nrow(reference_proteome_existence)>0){
  reference_proteome_prot.exist = data.frame(Entry = reference_proteome_existence$Entry, ref.protein_evidence = 1,stringsAsFactors = F)
  if(missing_proteins == FALSE){
    colnames(reference_proteome_prot.exist)[2] = "protein_evidence"
  }
  if(missing_proteins == TRUE){
    colnames(reference_proteome_prot.exist)[2] = "missing_protein_evidence"
  }
  ## Add reference_proteome - protein evidence exists - to the logic table
  reference_proteome_prot.exist.f = merge(x = proteome_vs_interactome_f, 
                                          y = reference_proteome_prot.exist, by.x="whole_proteome_IDs", 
                                          by.y="Entry", 
                                          all.x = T, all.y = T)
  }
  if(nrow(whole_proteome_existence)>0){
  whole_proteome_prot.exist = data.frame(Entry = whole_proteome_existence$Entry, protein_evidence = 1,stringsAsFactors = F)
  if(missing_proteins == FALSE){
    colnames(whole_proteome_prot.exist)[2] = "protein_evidence"
  }
  if(missing_proteins == TRUE){
    colnames(whole_proteome_prot.exist)[2] = "missing_protein_evidence"
  }
  ## Add whole_proteome - protein evidence exists - column of ones
  proteome_vs_interactome_f = merge(x = proteome_vs_interactome_f, 
                                  y = whole_proteome_prot.exist, by.x="whole_proteome_IDs", 
                                  by.y="Entry", 
                                  all.x = T, all.y = T)
  }
  ## substitute NAs for zeros in whole proteome-containing data.frame
  proteome_vs_interactome_f[is.na(proteome_vs_interactome_f)] <- 0
  #============================================================================#
  ## Saving logic table
  filename_vs_2 = paste("./analysis/","proteome_vs_interactome_f_", SPECIES_ID,"_reviewed_",reviewed,"_isoforms_",isoforms,"_", date,".txt", sep = "")
  write.table(proteome_vs_interactome_f,filename_vs_2,col.names=T,row.names=F,sep="\t",quote=F)
  
  
  #============================================================================#
  ### Summaries
  
  ## Summary of how many interactors have non-Uniprot identifiers 
  ## and how many interactors come from the species other than queried
  source("uniprotkb_and_SPECIES_ID_interactor_summary.R")
  uniprotkb_and_SPECIES_ID_interactor_summary = uniprotkb_and_SPECIES_ID_interactor_summary(all_interactors, SPECIES_ID, SPECIES_NAME)
  filename = paste("./summaries/","uniprotKB_IDs_and_",SPECIES_NAME,"_interactors_summary", "_isoforms_",isoforms,"_", date,".txt", sep = "")
  write.table(uniprotkb_and_SPECIES_ID_interactor_summary,filename,col.names=T,row.names=F,sep="\t",quote=F)
  

  print(uniprotkb_and_SPECIES_ID_interactor_summary)
  return(SPECIES_NAME)
  }
  }
  }
}