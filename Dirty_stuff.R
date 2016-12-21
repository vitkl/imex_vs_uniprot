## Dirty stuff









dim(IMEx)

IMEx2 = IMEx[1:8]


IMEx2[3] = IMEx[1]-IMEx[2]
IMEx2[7] = IMEx[5]-IMEx[6]

IMEx[5] = IMEx[1]-IMEx[2]

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp

barplot(IMEx[c(1,2,4)])

library(reshape2)
proteome_vs_interactome_summary.m = melt(proteome_vs_interactome_summary)
colnames(proteome_vs_interactome_summary.m) = c("Description", "interactome_database", "number")
proteome_vs_interactome_summary.m.d = as.data.frame(proteome_vs_interactome_summary.m)

proteome_vs_interactome_summary.m.d.f = dplyr::filter(proteome_vs_interactome_summary.m.d, interactome_database == "IMEx")


?geom_()


barplot(IMEx[c(1,2,4)])
?barplot()

library(UpSetR)

upset(data = proteome_vs_interactome_f,
      mainbar.y.label = "Number of proteins", 
      sets.x.label = "Number of proteins",
      order.by = c("freq"), decreasing = TRUE,
      sets = c("IMEx","reference_proteome_Uniprot"),
      mb.ratio = c(0.6, 0.4),
      point.size = 6, 
      line.size = 2, 
      queries = list(
        list(query = intersects, params = list("IMEx"), active =T, color = "#B00000"),
        list(query = intersects, params = list("IMEx","reference_proteome_Uniprot"), active =T, color = "#00A000")))




library(UpSetR)
str(proteome_vs_interactome_o)
?upset
dim(proteome_vs_interactome_o)[2]

upset(proteome_vs_interactome_f, 
      #nsets = 4, 
      #nintersects = 7,
      point.size = 6, 
      # name.size = 12, 
      line.size = 2, 
      mainbar.y.label = "Common publications", 
      sets.x.label = "Nr of publications in dataset", 
      order.by="freq",
      decreasing=FALSE,
      queries = list(
        list(query = intersects, params = list("X0469.IntAct.","whole_proteome_Uniprot"), color = "blue", active = T),
        list(query = intersects, params = list("X0469.IntAct.","reference_proteome_Uniprot"), color= "darkorange2",active = T),
        list(query = intersects, params = list("IMEx","reference_proteome_Uniprot"), color= "lightgoldenrod2",active = T),
        list(query = intersects, params = list("IMEx","whole_proteome_Uniprot"), color= "lightgoldenrod2",active = T),
        list(query = intersects, params = list("IMEx"), color= "gray70",active = T),
        list(query = intersects, params = list("reference_proteome_Uniprot"), color= "gray70",active = T),
        list(query = intersects, params = list("whole_proteome_Uniprot"), color= "gray70",active = T)))


upset(proteome_vs_interactome_o2[dim(proteome_vs_interactome_o)[2]], nsets = 7, nintersects = 30, mb.ratio = c(0.8, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

###### dirty stuff

## Interactions per database vs reference proteome
ref_proteome_vs_interactome = merge(x = reference_proteome_unique, 
                                    y = all_interactors_unique, by.x="Reference_proteome_IDs", 
                                    by.y="interactor_IDs", 
                                    all.x = T, all.y = F)
ref_proteome_vs_interactome_summary = as.data.frame(table(ref_proteome_vs_interactome$interactor_databases,useNA="ifany"))
colnames(ref_proteome_vs_interactome_summary) = c("database", "number of proteins")

## Interactions per database vs whole proteome
whole_proteome_vs_interactome = merge(x = whole_proteome_unique, 
                                      y = all_interactors_unique, by.x="whole_proteome_IDs", 
                                      by.y="interactor_IDs", 
                                      all.x = T, all.y = F)
whole_proteome_vs_interactome_summary = as.data.frame(table(whole_proteome_vs_interactome$interactor_databases,useNA="ifany"))
colnames(whole_proteome_vs_interactome_summary) = c("database", "number of proteins")

#============================================================================#

## unique interactors vs whole proteome

source("unique_vs_proteome.R")
source("proteome_vs_unique.R")
interactome_in_proteome = unique_vs_proteome(all_interactors_unique, whole_proteome_unique)
proteome_in_interactome = proteome_vs_unique(all_interactors_unique, whole_proteome_unique)

## unique interactors vs reference proteome
unique_vs_proteome(all_interactors_unique, reference_proteome_unique)
proteome_vs_unique(all_interactors_unique, reference_proteome_unique)

library(dplyr)

source("proteome_vs_interactome_merge.R")
proteome_vs_interactome_merge = proteome_vs_interactome_merge(all_interactors_unique, whole_proteome_unique)
proteome_vs_interactome_merge[proteome_vs_interactome_merge$Uniprot == "<NA>",]
head(proteome_vs_interactome_merge,1000)

filter(proteome_vs_interactome_merge, Uniprot != "Uniprot")

library(VennDiagram)
## https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
## http://www.r-graph-gallery.com/14-venn-diagramm/

y = UniProtLIST[,2]

yp = data.frame(UniProtLIST[,2],stringsAsFactors = F)
colnames(yp) <- "up_ac"

## checks w
diff_check <- unique(merge(yp,x1,by.x="up_ac",by.y="id",all.x=T,all.y=F))
table(diff_check$db,useNA="ifany")

save(difference,file = "difference")


Summary_ = matrix()

Summary_[1,1] = mean(difference)
sum(difference)
length(x)
length(unique(x))
Summary_[1,5] = length(y)
Summary_[1,6] = length(unique(y))

colnames(Summary_) = c("Fraction - Uniprot covered by IntAct",
                       "Number - Uniprot covered by IntAct"
                       "Number of proteins in IntAct")


# dirtystuffstarts here

head(reference_proteome_ftp)
ftp_isoforms<-reference_proteome_ftp[grepl("-",reference_proteome_ftp$V2),]
query_isoforms <- reference_proteome_query[grepl("-",reference_proteome_query$V2),]

x =table(all_interactors$interactor_IDs,useNA="ifany")
dim(x)


movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=TRUE, sep=";" )

require(ggplot2); require(plyr); require(gridExtra); require(grid);

between <- function(row, min, max){
  newData <- (row["ReleaseDate"] < max) & (row["ReleaseDate"] > min)
}

plot1 <- function(mydata, x){
  myplot <- (ggplot(mydata, aes_string(x= x, fill = "color"))
             + geom_histogram() + scale_fill_identity()
             + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

plot2 <- function(mydata, x, y){
  myplot <- (ggplot(data = mydata, aes_string(x=x, y=y, colour = "color"), alpha = 0.5)
             + geom_point() + scale_color_identity()
             + theme_bw() + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

attributeplots <- list(gridrows = 55,
                       plots = list(list(plot = plot1, x= "ReleaseDate",  queries = FALSE),
                                    list(plot = plot1, x= "ReleaseDate", queries = TRUE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = FALSE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = TRUE)),
                       ncols = 3)

upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.8, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      upset(movies, sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
            queries = list(list(query = intersects, params = list("Drama", "Action")),
                           list(query = between, params = list(1970, 1980), color = "red", active = TRUE)))
      
      upset(movies, attribute.plots = attributeplots,
            queries = list(list(query = between, params = list(1920, 1940)),
                           list(query = intersects, params = list("Drama"), color= "red"),
                           list(query = elements, params = list("ReleaseDate", 1990, 1991, 1992))),
            main.bar.color = "yellow")
      
      
      