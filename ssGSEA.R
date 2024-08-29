# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-05-29
#options(echo=FALSE)
source('OPAM.library.R')
source('FileCreator.R')
library(ggplot2)
library(stringr)
source("InputProcessor.R")

# code from Pablo/Diana
# Project the dataset in the space of the entire database of pathways

ssGSEA <- function(gene_list, output_filename, gene_exp_profile){
    createGMT(gene_list)
    OPAM.project.dataset.2(
    input.ds =                gene_exp_profile,
    output.ds =               paste(output_filename, ".gct", sep=""),
    gene.set.databases =      "IntermediateFiles/genes.gmt",
    gene.set.selection =      "ALL",
    sample.norm.type         = "rank",  # "rank", "log" or "log.rank"
    weight                   = 0.75,
    statistic                = "area.under.RES",
    output.score.type        = "ES",
    combine.mode             = "combine.add",  # "combine.off", "combine.replace", "combine.add"
    nperm                    =  100,
    correl.type              = "z.score")
}

ssGSEAData <- function(gene_list, exp_profile){

    b <- data.frame()
    #run ssGSEA1.0
    ssGSEA(gene_list, "IntermediateFiles/out_score", exp_profile)
    plotData <- readLines("./IntermediateFiles/out_score.gct")
    headers <- grepl("Anoymous|Name", plotData)
    plotData <- plotData[headers]
    a <- read.table(text = plotData, sep = "\t")
    a <- a[,-c(2)]
    b <- as.data.frame(t(a))
    colnames(b) <- as.character(unlist(b[1,]))
    b <- b[-c(1),]
    colnames(b)[colnames(b)=='Name'] = 'sample_id'
    return(b)
}
rename_a_column <- function(df, old, new){
    names(df)[names(df)==old] <- new
    return(df)
}






