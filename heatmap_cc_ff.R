#!/usr/bin/env Rscript

################
### 
### R script created by Juan Salamanca Viloria. 
### 
###
#################

#    PyinKnife script to plot heatmap for the connected component.
#    Copyright (C) 2016 Juan Salamanca Viloria <juan.salamanca.viloria@gmail.com>, Elena Papaleo <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.




library(gplots)
library(lattice)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)


new.df <- read.table(args[1],quote=" ")



abr_aa <- function(x){
  i <- 0
  new_colnames <- c()
  for (aa in (substr(colnames(x),0,3) )){
    i <- i+1
  if (aa=="ALA"){
    new_colnames[i] <- "A"
  }
  else if (aa== "ARG"){
    new_colnames[i] <-  "R"
  }
  else if (aa== "ASN"){
    new_colnames[i] <-  "N"
  }
  else if (aa== "ASP"){
    new_colnames[i] <-  "D"
  }
  else if (aa== "CYS"){
    new_colnames[i] <-  "C"
  }
  else if (aa== "GLU"){
    new_colnames[i] <-  "E"
  }
  else if (aa== "GLN"){
    new_colnames[i] <-  "Q"
  }
  else if (aa=="GLY"){
    new_colnames[i] <-  "G"
  }
  else if (aa== "HIS"){
    new_colnames[i] <-   "H" 
  }
  else if (aa== "ILE"){
    new_colnames[i] <-   "I" 
  }
  else if (aa== "LEU"){
    new_colnames[i] <-   "L" 
  }
  else if (aa== "LYS"){
    new_colnames[i] <-   "K" 
  }
  else if (aa== "MET"){
    new_colnames[i] <-   "M" 
  }
  else if (aa== "PHE"){
    new_colnames[i] <-   "F" 
  }
  else if (aa== "PRO"){
    new_colnames[i] <-   "P" 
  }
  else if (aa== "SER"){
    new_colnames[i] <-   "S" 
  }
  else if (aa== "THR"){
    new_colnames[i] <-   "T" 
  }
  else if (aa== "TRP"){
    new_colnames[i] <-   "W"
  }
  else if (aa== "TYR"){
    new_colnames[i] <-   "Y" 
  }
  else if (aa== "VAL"){
    new_colnames[i] <-   "V"
  }
  }
  return(new_colnames)
}



new.df <- t(unique(new.df[rowSums(is.na(new.df)) != ncol(new.df),1:5]))

colnames(new.df) <- new.df[1,]
new.df<- new.df[-1,]


new_col<- abr_aa(new.df)
new_col <- paste(new_col,as.numeric(new.df[1,]), sep ="")
colnames(new.df)<- new_col
new.df<- new.df[-1,]
row.names(new.df) <- c("A99*ILDN", "C22*", "C36")

###### TRANSFORM

data_num <- factor(as.numeric(new.df[1,]))

tb <- table(data_num)
sorted_1 <- factor(data_num,levels = names(tb[order(tb, decreasing = TRUE)]))
levels(sorted_1)[1:5] <- as.numeric(seq(1,5))
new.df[1,] <- sorted_1

data_num2 <- factor(as.numeric(new.df[2,]))

tb2 <- table(data_num2)
sorted_2 <- factor(data_num2,levels = names(tb2[order(tb2, decreasing = TRUE)]))
levels(sorted_2)[1:5] <- as.numeric(seq(1,5))
new.df[2,] <- sorted_2

data_num3 <- factor(as.numeric(new.df[3,]))
tb3 <- table(data_num3)
sorted_3 <- factor(data_num3,levels = names(tb3[order(tb3, decreasing = TRUE)]))
levels(sorted_3)[1:5] <- as.numeric(seq(1,5))
new.df[3,] <- sorted_3
#colors 
my_palette <- c("red2","orange", "yellow","forestgreen","royalblue4", "gray84", "gray84")
pdf(args[2])

print(levelplot(t(new.df[,as.numeric(args[3]):as.numeric(args[4])]),
          col.regions= my_palette,
          xlab = "Residue",
          ylab = " ",
          main = " ",
          #key = my.ckey,
          at = c(seq(0,5), c(5.99,max(as.numeric(levels(factor(new.df)))))),
          scales=list(x=list(rot=90, cex = .7)),
          #colorkey=list(at=seq(0, 1, 0.2), labels=list(at=factor(all_hubs), labels=c("none", "a bit", "a bit more", "a lot"))),
          par.settings = list(axis.line = list(col = "transparent"),strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')),
          panel = function(...){
            panel.levelplot(...)
            #panel.grid(h = 10, v = 49, col.line = "black")
            panel.abline(v = seq(1,as.numeric(args[4])-as.numeric(args[3])+2)-0.5, h = seq(1,11)-0.5, col = "black")
          }))

dev.off()