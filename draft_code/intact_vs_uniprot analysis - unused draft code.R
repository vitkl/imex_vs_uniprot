# intact_vs_uniprot analysis - unused draft code
## Transform table for plotting with ggplot2
library(reshape2)
interactome_identifiers_summary.all.m= melt(data = interactome_identifiers_summary.all,
                                            id.vars = c("SPECIES_NAME", "SPECIES_ID"),
                                            variable.name = "decription",
                                            value.name = "number")
## rename yeast and E.coli
interactome_identifiers_summary.all.m$SPECIES_NAME = gsub("strain ATCC 204508", "Saccharomyces cerevisiae, strain S288c", interactome_identifiers_summary.all.m$SPECIES_NAME)
interactome_identifiers_summary.all.m$SPECIES_NAME = gsub("strain K12", "Escherichia coli, strain K12", interactome_identifiers_summary.all.m$SPECIES_NAME)

# plot
library(ggplot2)
library(dplyr)

proteome_vs_interactome_plot <- ggplot(interactome_identifiers_summary.all.m, aes(x=decription, y=number, fill=decription)) + geom_bar(width = 0.9, stat = "identity", position = "stack") + geom_label(aes(label=number), position = "stack", size = 4, label.padding = unit(0.08, "lines")) + facet_wrap(~SPECIES_NAME, scales ="free_y") +
  ggtitle("The number of interactors which have UniprotKB identifiers and belong to the species queried") + theme(text=element_text(size=13,  family="serif"))
proteome_vs_interactome_plot
# save plot
filename=paste("./results/", "interactome_identifiers_summary_plot_",Sys.Date(),".png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 9)

#################
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
  ggtitle("The number of interactors which have UniprotKB identifiers and belong to the species queried") + theme(text=element_text(size=13,  family="serif"))
proteome_vs_interactome_plot
# save plot
filename=paste("./results/", "interactome_identifiers_summary_plot_SMALL",Sys.Date(),".png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 12)

### ======================================================================== ###

# selecting whole proteome summary
proteome_vs_interactome_summary.whole = proteome_vs_interactome_summary.all[c(1,2,3,4,8,7,6)]
# selecting reference proteome summary
proteome_vs_interactome_summary.ref = proteome_vs_interactome_summary.all[c(1,2,3,4,12,11,10)]

### ======================================================================== ###
# plotting whole proteome comparison with ggplot2
## whole proteome ggplot2 - transforms the data (melt) into the format 
## acceptable by ggplot2 and plots bar plots 

library(reshape2)
proteome_vs_interactome_summary.whole.m= melt(data = proteome_vs_interactome_summary.whole,
                                              id.vars = c("species.name", "database", "reviewed", "isoforms"),
                                              variable.name = "decription",
                                              value.name = "number")

proteome_vs_interactome_summary.whole.m$species.name = gsub("strain ATCC 204508", "Saccharomyces cerevisiae, strain S288c", proteome_vs_interactome_summary.whole.m$species.name)
proteome_vs_interactome_summary.whole.m$species.name = gsub("strain K12", "Escherichia coli, strain K12", proteome_vs_interactome_summary.whole.m$species.name)

library(ggplot2)
library(ggrepel)
library(dplyr)
#proteome_vs_interactome_summary.whole.m = filter(proteome_vs_interactome_summary.whole.m, species.name == "Homo sapiens")
proteome_vs_interactome_summary.whole.m$revieweddatabase = interaction(proteome_vs_interactome_summary.whole.m$isoforms, proteome_vs_interactome_summary.whole.m$database)

proteome_vs_interactome_plot <- ggplot(proteome_vs_interactome_summary.whole.m, aes(x=revieweddatabase, y=number, fill=decription))+  geom_bar(width = 0.9, stat = "identity", position = "stack") + geom_label(aes(label=number), position = "stack", size = 4, label.padding = unit(0.08, "lines")) + facet_grid(species.name~reviewed, scales ="free_y", space ="fixed")
proteome_vs_interactome_plot

filename=paste("./results/", "all_Uniprot","whole_proteome_vs_interactome_plot",Sys.Date(),"2.png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 36)

### ======================================================================== ###
# plotting reference proteome comparison with ggplot2
## reference proteome ggplot2 - transforms the data (melt) into the format 
## acceptable by ggplot2 and plots bar plots 

library(reshape2)
proteome_vs_interactome_summary.ref.m= melt(data = proteome_vs_interactome_summary.ref,
                                            id.vars = c("species.name", "database", "reviewed", "isoforms"),
                                            variable.name = "decription",
                                            value.name = "number")

proteome_vs_interactome_summary.ref.m$species.name = gsub("strain ATCC 204508", "Saccharomyces cerevisiae, strain S288c", proteome_vs_interactome_summary.ref.m$species.name)
proteome_vs_interactome_summary.ref.m$species.name = gsub("strain K12", "Escherichia coli, strain K12", proteome_vs_interactome_summary.ref.m$species.name)

library(ggplot2)
library(ggrepel)
library(dplyr)
#proteome_vs_interactome_summary.whole.m = filter(proteome_vs_interactome_summary.whole.m, species.name == "Homo sapiens")
proteome_vs_interactome_summary.ref.m$revieweddatabase = interaction(proteome_vs_interactome_summary.ref.m$isoforms, proteome_vs_interactome_summary.ref.m$database)

proteome_vs_interactome_plot <- ggplot(proteome_vs_interactome_summary.ref.m, aes(x=revieweddatabase, y=number, fill=decription))+  geom_bar(width = 0.9, stat = "identity", position = "stack") + geom_label(aes(label=number), position = "stack", size = 4,label.padding = unit(0.08, "lines")) + facet_grid(species.name~reviewed, scales ="free_y", space ="fixed")
proteome_vs_interactome_plot

filename=paste("./results/", "all_Uniprot","reference_proteome_vs_interactome_plot",Sys.Date(),"2.png", sep = "")
ggsave(filename, proteome_vs_interactome_plot, width = 12, height = 36)
### ======================================================================== ###

SPECIES_NAME = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "strain ATCC 204508", "strain K12", "Drosophila melanogaster", "Caenorhabditis elegans", "Schizosaccharomyces pombe", "Arabidopsis thaliana", "Tetrahymena thermophila", "Neurospora crassa", "Chlamydomonas reinhardtii", "Gallus gallus", "Danio rerio")
SPECIES_NAME = SPECIES_NAME[1:6]
### ======================================================================== ###

## to make into function

database = "IMEx"
reviewed = 1
isoforms = TRUE
interactome_vs_proteome_venn(SPECIES_NAME, proteome_vs_interactome_summary.all, database = database, reviewed = reviewed, isoforms = isoforms)

interactome_vs_proteome_venn = function(SPECIES_NAME, proteome_vs_interactome_summary.all, database, reviewed, isoforms) {
  
  library(dplyr)
  library(VennDiagram)
  grid.newpage()
  
  pushViewport(viewport(layout=grid.layout(nrow = length(SPECIES_NAME), ncol=3, widths = unit(c(2/6,2/6,2/6), "npc"), heights = unit(rep(1/length(SPECIES_NAME),length(SPECIES_NAME)), "npc"))))
  
  z = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[1]) %>% filter(database == database)
  z1 = z$whole.proteome..Uniprot.
  
  for (i in 1:length(SPECIES_NAME)) {
    
    pushViewport(viewport(layout.pos.col=2, layout.pos.row = i))
    xxx = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[i]) %>% filter(database == database) %>% filter(reviewed == reviewed) %>% filter(isoforms == isoforms)
    draw.pairwise.venn(xxx$whole.proteome..Uniprot., xxx$whole.proteome..interactome.available+xxx$interactome..but.not.in.Uniprot, xxx$whole.proteome..interactome.available, category = c("Whole proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 180), cat.dist = rep(0.035, 2), cat.cex = c(0.9,0.9), scaled = TRUE, margin = 0.05)
    popViewport()
    
    pushViewport(viewport(layout.pos.col=3, layout.pos.row = i))
    xxx2 = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[i]) %>% filter(database == database) %>% filter(reviewed ==reviewed) %>% filter(isoforms == isoforms)
    draw.pairwise.venn(xxx2$reference.proteome..Uniprot., xxx2$reference.proteome..interactome.available+xxx2$interactome..but.not.in.Uniprot.1, xxx2$reference.proteome..interactome.available, category = c("Reference proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 180), cat.dist = rep(0.035, 2), cat.cex = c(0.9,0.9), scaled = TRUE, margin = 0.05)
    popViewport()
  }
}

### ======================================================================== ###
library(dplyr)
library(VennDiagram)
grid.newpage()

pushViewport(viewport(layout=grid.layout(nrow = length(SPECIES_NAME), ncol=3, widths = unit(c(2/4,1/4,1/4), "npc"), heights = unit(rep(1/length(SPECIES_NAME),length(SPECIES_NAME)), "npc"))))

z = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[1]) %>% filter(database == "IMEx")
z1 = z$whole.proteome..Uniprot.

for (i in 1:length(SPECIES_NAME)) {
  
  pushViewport(viewport(layout.pos.col=2, layout.pos.row = i))
  xxx = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[i]) %>% filter(database == "IMEx") %>% filter(reviewed ==2) %>% filter(isoforms == FALSE)
  draw.pairwise.venn(xxx$whole.proteome..Uniprot., xxx$whole.proteome..interactome.available+xxx$interactome..but.not.in.Uniprot, xxx$whole.proteome..interactome.available, category = c("Whole proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 180), cat.dist = rep(0.035, 2), cat.cex = c(0.9,0.9), scaled = TRUE, margin = 0.05)
  popViewport()
  
  pushViewport(viewport(layout.pos.col=3, layout.pos.row = i))
  xxx2 = filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[i]) %>% filter(database == "IMEx") %>% filter(reviewed ==2) %>% filter(isoforms == FALSE)
  draw.pairwise.venn(xxx2$reference.proteome..Uniprot., xxx2$reference.proteome..interactome.available+xxx2$interactome..but.not.in.Uniprot.1, xxx2$reference.proteome..interactome.available, category = c("Reference proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 180), cat.dist = rep(0.035, 2), cat.cex = c(0.9,0.9), scaled = TRUE, margin = 0.05)
  popViewport()
}
### ======================================================================== ###



png(filename = "Pairwise_Venn_diagram.png", width = 1200, height = 1200)
grid.draw(xxxx);
dev.off();

pushViewport(viewport(layout.pos.col=1, layout.pos.row = ))
xxx =filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[1]) %>% filter(database == "IMEx")
draw.pairwise.venn(xxx$whole.proteome..Uniprot., xxx$whole.proteome..interactome.available+xxx$interactome..but.not.in.Uniprot, xxx$whole.proteome..interactome.available, category = c("Whole proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
popViewport()

pushViewport(viewport(layout.pos.col=2, layout.pos.row = ))
xxx2 =filter(proteome_vs_interactome_summary.all, species.name == SPECIES_NAME[1]) %>% filter(database == "IMEx")
xxxx = draw.pairwise.venn(xxx2$reference.proteome..Uniprot., xxx2$reference.proteome..interactome.available+xxx2$interactome..but.not.in.Uniprot.1, xxx2$reference.proteome..interactome.available, category = c("Reference proteome", "interactome from IMEx"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
ggsave(filename = "trndjmndj.png")

grid.draw(xxxx)

library(ggplot2)
library(grid)
library(gridExtra)

grid.arrange(gTree(children = xxx.plot), gTree(children = xxx2.plot), ncol = 2, main = "Venn")