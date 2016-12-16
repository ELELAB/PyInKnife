#!/usr/bin/env Rscript

################
### 
### R script created by Juan Salamanca Viloria. 
### 
###
#################

#    PyinKnife script to plot hubs and connected component distribution.
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




###
###Example 
### no h-bonds
###  /data/user/juan/pyinteraph_elena/final_version/all_hubs_distribution.R distribution_hubs.txt
## h-bond
### /data/user/juan/pyinteraph_elena/final_version/hubs_distribution.R distribution_hubs.txt sc-sc

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)



new.df <- read.table(args[1],quote="\"")


colnames(new.df) <- c("Position", rep("Number", each= dim(new.df)[2]-1))
#c("Position", "Number","Number","Number","Number","Number","Number","Number","Number","Number","Number", "Number")

sd2 <- function(x) sqrt(sum((x - mean(x))^2)*((length(x)-1)/ length(x)))
new.df$sd <- apply(new.df[3:12],1,sd2)

#std <- function(x) sd(x)/sqrt(length(x))
#new.df$std <- apply(new.df[,colnames(new.df)[2:4]] ,1,std)



if (length(args) > 1) {
if (args[2] == "sc-sc") {
pdf("hubs_distribution_sc.pdf")}
if (args[2] == "mc-mc"){
pdf("hubs_distribution_mc.pdf")}
if (args[2] == "mc-sc"){
pdf("hubs_distribution_mc-sc.pdf")}

} else {pdf("hubs_distribution.pdf")}
ggplot(new.df, aes(x=Position, y= Number)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Number-sd, ymax=Number+sd),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  
  xlab("Hubs") +
  ylab("Number of elements") +
  
  #scale_y_continuous(breaks=0:20*5) +
  
  theme_bw()+ theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
         axis.title=element_text(size=22,face="bold"),axis.text=element_text(size=18))
dev.off()
