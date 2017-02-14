##============================================================================##
## swissprot_vs_imex_protein_properties

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
      ## querying Uniprot for the list of proteins (only mapped to Uniprot)
      ## downolading length, mass, SNPs, isoforms, annotation score, existence attribute
      SPECIES_IDs = SPECIES_NAME_TO_ID(n)
      SPECIES_ID = SPECIES_IDs$SPECIES_ID
      ## reading logic table and getting list of proteins
      if(SPECIES_NAME != "strain K12"){
      filename_vs_3 = paste("./analysis/","proteome_vs_interactome_vs_BioGRID_f_", SPECIES_ID,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")}
      if(SPECIES_NAME == "strain K12"){
      filename_vs_3 = paste("./analysis/","proteome_vs_interactome_f_", SPECIES_ID,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")}
      biogrid_from_mentha_vs_proteome_vs_imex_f = as.data.frame(read.delim(filename_vs_3, stringsAsFactors = F))
      IDs = dplyr::filter(biogrid_from_mentha_vs_proteome_vs_imex_f, whole_proteome_Uniprot==1)%>% dplyr::select(whole_proteome_IDs)
      filename_perl = paste0("./Data/IDs_for_PERL_",date,".txt")
      filename_out_of_perl = paste0("./Data/Swissprot_list_from_proteome_vs_interactome_vs_BioGRID_f_", SPECIES_ID,"_reviewed_",r,"_isoforms_",i,"_", date,"")
      file.create(filename_out_of_perl)
      write(IDs$whole_proteome_IDs, file = filename_perl)
      
      ## querying the list of proteins
      # system("./scripts/up_query.pl <YOUR INPUT FILE HERE> > <YOUR OUTPUT FILE HERE>")
      #system(paste0("perl ./scripts/up_query.pl ",filename_perl," > ", filename_out_of_perl))
      #uniprot_query = read.delim(filename_out_of_perl, header = F, stringsAsFactors = F)
      #str(uniprot_query,1)
      source("download_whole_proteome.R")
      whole_proteome_query = download_whole_proteome(SPECIES_ID, date = "2016-12-14")
      uniprot_query = dplyr::filter(whole_proteome_query, Status == "reviewed")
      #system.time({
      #uniprot_query = data.frame()
      #for (id in 1:length(IDs)) {
      #    url = paste0("http://www.uniprot.org/uniprot/?query=id:",IDs[id],"&format=tab&columns=id,reviewed,length,organism,organism-id,mass,database(dbSNP),comment(ALTERNATIVE%20PRODUCTS),annotation%20score,existence,protein%20names")
      #    url = paste0("http://www.uniprot.org/uniprot/?query=id:",IDs[id],"&format=tab&columns=id,reviewed,length,organism,organism-id,mass,annotation%20score,existence")
      #              uniprot_query_temp = UniProtLIST = as.data.frame(read.delim(url,stringsAsFactors = F, quote = ""))
      #    uniprot_query = rbind(uniprot_query, uniprot_query_temp)
      #}
      #})
      ##============================================================================##
      ## merging information from Uniprot to the logic table
      proteome_vs_imex_details_f = merge(filter(biogrid_from_mentha_vs_proteome_vs_imex_f, whole_proteome_Uniprot==1), 
                                         uniprot_query, 
                                         by.x = "whole_proteome_IDs",
                                         by.y = "Entry")
      proteome_vs_imex_details_f$Mass = gsub(",","",proteome_vs_imex_details_f$Mass)
      proteome_vs_imex_details_f$Mass = as.numeric(proteome_vs_imex_details_f$Mass)
      ## creating a factor variable for presence_in_Uniprot.presence_in_IMEx (1_0/1_1)
      proteome_vs_imex_details_f[,length(proteome_vs_imex_details_f)+1] = interaction(proteome_vs_imex_details_f$whole_proteome_Uniprot, proteome_vs_imex_details_f$IMEx, sep = "_")
      colnames(proteome_vs_imex_details_f)[length(proteome_vs_imex_details_f)] = paste0(colnames(proteome_vs_imex_details_f)[2], "_","IMEx")
      levels(proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx) = c("SwissProt_not_IMEX", "SwissProt_and_IMEX")
      ##============================================================================##
      ## saving combined logic table + protein properties from Uniprot
      filename_vs_2 = paste("./analysis/","proteome_vs_interactome_protein_properties_f_", n,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")
      write.table(proteome_vs_imex_details_f,filename_vs_2,col.names=T,row.names=F,sep="\t",quote=F)
      ##============================================================================##
      ## Analysis
      # read the saved table
      filename_vs_2 = paste("./analysis/","proteome_vs_interactome_protein_properties_f_", n,"_reviewed_",r,"_isoforms_",i,"_", date,".txt", sep = "")
      proteome_vs_imex_details_f = as.data.frame(read.delim(filename_vs_2, header = T, stringsAsFactors = F,quote=""))
      proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx = factor(proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx, ordered =F)
      ## the distribution of mass
            ## plotting the distribution of mass
      library(ggplot2)
      ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_density()+ scale_x_log10()
      ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) + scale_x_log10() + geom_histogram(position = "identity", bins = 50)
      ## removing olfactory receptors:
      proteome_vs_imex_details_f_minus_odor = proteome_vs_imex_details_f[-grep("Odor", proteome_vs_imex_details_f$Protein.names),]
      proteome_vs_imex_details_f_minus_odor_olf = proteome_vs_imex_details_f_minus_odor[-grep("Olfactory", proteome_vs_imex_details_f$Protein.names),]
      ## density without olfactory receptors:
      ggplot(proteome_vs_imex_details_f_minus_odor_olf, aes(x = Mass, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_density()+ scale_x_log10()
      # protein mass vs Protein.existence status
      # ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = Protein.existence, alpha =0.5)) +geom_density() + facet_wrap(~whole_proteome_Uniprot_IMEx) + scale_x_log10()
      # ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = Protein.existence, alpha =0.5)) +geom_histogram(position = "identity", bins = 50) + facet_wrap(~whole_proteome_Uniprot_IMEx) + scale_x_log10()
      # ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_density() + facet_wrap(~Protein.existence) + scale_x_log10()
      ggplot(proteome_vs_imex_details_f, aes(x = Mass, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_histogram(position = "identity", bins = 50) + facet_wrap(~Protein.existence) + scale_x_log10()
            ## testing the difference in protein mass distributions
            ##of proteins absent/present in IMEx using Wilcox rank test
      SwissProt_not_IMEX = proteome_vs_imex_details_f$Mass[proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx=="SwissProt_not_IMEX"]
      SwissProt_and_IMEX = proteome_vs_imex_details_f$Mass[proteome_vs_imex_details_f$whole_proteome_Uniprot_IMEx=="SwissProt_and_IMEX"]
      wilcox.test.r = wilcox.test(SwissProt_not_IMEX, SwissProt_and_IMEX, conf.int = T)
      wilcox.test.r
            ## is the distibution of mass or log(mass) normal?
      library(rafalib)
      mypar(2,2)
      qqnorm(log10(SwissProt_not_IMEX), main = "SwissProt_not_IMEX, log(protein mass)")
      qqline(log10(SwissProt_not_IMEX))
      qqnorm(log10(SwissProt_and_IMEX), main = "SwissProt_and_IMEX, log(protein mass)")
      qqline(log10(SwissProt_and_IMEX))
      qqnorm((SwissProt_not_IMEX), main = "SwissProt_not_IMEX, protein mass")
      qqline((SwissProt_not_IMEX))
      qqnorm((SwissProt_and_IMEX), main = "SwissProt_and_IMEX, protein mass")
      qqline((SwissProt_and_IMEX))
      
      mean(log10(SwissProt_not_IMEX))
      popsd(log10(SwissProt_not_IMEX))
      median(log10(SwissProt_not_IMEX))
      mad(log10(SwissProt_not_IMEX))
            ## it indeed is
      
            ## Monte-Carlo simulation
            ## is there a significant difference in protein mass?
      # code can be used for things other than mass
      set.seed(1)
      # sample size
      Ns <- seq(5, 350, 20)
      lengthNs = length(Ns)
      # number of simulations
      B = 100
      # function which takes samples and does Wilcox test
      simulation = function(n){  
        wilcox.test.rr = matrix(0, 1, 3)
        x = sample(SwissProt_not_IMEX, n)
        y = sample(SwissProt_and_IMEX, n)
        wilcox.test.rr[1,c(1,2)] = wilcox.test(x, y, conf.int = T)$conf.int
        wilcox.test.rr[1,3] = wilcox.test(x, y, conf.int = T)$p.value
        return(wilcox.test.rr)
      }
      wilcox.test.rr = matrix(0, lengthNs, 3)
      # looping Monte-Carlo over samples sizes ]
      for (n in 1:lengthNs) {
        z =replicate(B, simulation(Ns[n]), simplify = T)
        wilcox.test.rr[n,] = rowMeans(z)
      }
      xx =cbind(wilcox.test.rr, Ns)
      plot(xx[,4], xx[,3], ylab = "p-val", xlab = "N")
      ##============================================================================##
      ## Analysis
      ## the distribution of variants
      
      ## how to count variants
      variants= numeric()
      source("variant_id_extractor_Uniprot.R")
      for (i in 1:nrow(proteome_vs_imex_details_f)) {
        variants[i] = variant_id_extractor_Uniprot(proteome_vs_imex_details_f[i,])
      }
      ## adding variant number to the logic table, 1 is added to make counting 0s 
      ## in histogram possible
      proteome_vs_imex_details_f.var = cbind(proteome_vs_imex_details_f,variants)
      proteome_vs_imex_details_f.var = cbind(proteome_vs_imex_details_f,variants_1 = variants+1)
      ## plotting the histogram of variants
      ggplot(proteome_vs_imex_details_f.var, aes(x = variants_1, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_histogram(position = "identity", bins = 100) +xlim(0,50) +ylim(0,5800)
      ggplot(proteome_vs_imex_details_f.var, aes(x = variants_1, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_density() +xlim(0,50)
      table(proteome_vs_imex_details_f.var$variants_1)
      ## testing the difference in the natural variant distribution using Wilcox rank test
      SwissProt_not_IMEX.var = proteome_vs_imex_details_f.var$variants_1[proteome_vs_imex_details_f.var$whole_proteome_Uniprot_IMEx=="1.0"]
      SwissProt_and_IMEX.var = proteome_vs_imex_details_f.var$variants_1[proteome_vs_imex_details_f.var$whole_proteome_Uniprot_IMEx=="1.1"]
      wilcox.test(SwissProt_not_IMEX.var, SwissProt_and_IMEX.var)
      ##============================================================================##
      ## Analysis
      ## uniprot annotation score (1-5)
      proteome_vs_imex_details_f.var$Annotation = as.numeric(substr(proteome_vs_imex_details_f.var$Annotation, 1,1))
      ggplot(proteome_vs_imex_details_f.var, aes(x = Annotation, color = whole_proteome_Uniprot_IMEx, alpha =0.5)) +geom_histogram(position = "identity", bins = 100) +xlim(0,6)# +ylim(0,5800)
      ## testing the difference in the annotation score distribution using Wilcox rank test
      SwissProt_not_IMEX.ann = proteome_vs_imex_details_f.var$Annotation[proteome_vs_imex_details_f.var$whole_proteome_Uniprot_IMEx=="1.0"]
      SwissProt_and_IMEX.ann = proteome_vs_imex_details_f.var$Annotation[proteome_vs_imex_details_f.var$whole_proteome_Uniprot_IMEx=="1.1"]
      wilcox.test(SwissProt_not_IMEX.ann, SwissProt_and_IMEX.ann)
      ##============================================================================##
#    }
#  }
#}

##============================================================================##
## disordered proteins
      ## reading downloaded data (from blob:http://www.disprot.org/90b4c112-691e-4fbd-9895-cf60dc9cc68b)
filename_dis = "./Data/disordered_proteins_2016-12-06"
table(disordered_proteins$organism) = read.csv(filename_dis, stringsAsFactors = F)
Species_ids_disorder = c("Mouse", "Human", "Baker's yeast", "Escherichia coli (strain K12)")
SPECIES_NAMEs = c("Homo sapiens", "Mus musculus", "strain ATCC 204508", "strain K12")
Species_id_disorder = Species_ids_disorder[SPECIES_NAME ==SPECIES_NAMEs]
disordered_proteins_uniprotac = disordered_proteins$uniprot_accession[disordered_proteins$organism==Species_id_disorder]
## Saving uniprot accessions to check if they map - they indeed do
write(disordered_proteins_uniprotac, file = "disordered_proteins_2016-12-06_UniprotAC.txt")
## removing isoform ids
disordered_proteins_uniprotac.noisof =  data.frame(Entry = isoform_id_all_remover(disordered_proteins_uniprotac), disordered_proteins = 1)
## merging disordered protein list to the logic table
proteome_vs_imex_details_disordered_f = merge(proteome_vs_imex_details_f.var, 
                                   disordered_proteins_uniprotac.noisof, 
                                   by.x = "whole_proteome_IDs",
                                   by.y = "Entry", all = T)
proteome_vs_imex_details_disordered_f[is.na(proteome_vs_imex_details_disordered_f)] = 0

proteome_vs_imex_details_disordered_f$IMEx = proteome_vs_imex_details_disordered_f


## draft code
##============================================================================##
#url = "http://www.uniprot.org/uniprot/?query=id:Q5S007&format=tab&columns=id,reviewed,length,organism,organism-id,mass,database(dbSNP),comment(ALTERNATIVE%20PRODUCTS),annotation%20score,existence"
#UniProtLIST_filename <- paste("./Data/","Reference_proteome_speciesID","5364_","_","_", sep = "")
#downloader::download(url, destfile = UniProtLIST_filename) 
#print("Querying Uniprot for the reference proteome...loaded from uniprot")

## read the reference proteome, some columns may include quotes
#UniProtLIST = as.data.frame(read.delim(UniProtLIST_filename,stringsAsFactors = F, quote = ""))
#UniProtLIST = as.data.frame(read.delim(url,stringsAsFactors = F, quote = ""))