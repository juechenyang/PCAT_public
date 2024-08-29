# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-06-07
#input a string contains multiple gene and return a gene list
multiple_gene_processor <- function(stringIn){
     # preprocess user input
    input_string = gsub("[^[:alnum:] ^[:cntrl:] ^\\. ^\\-]", "", stringIn)
    input_string <- toupper(input_string)
    input_string <- gsub("^\\s+|\\s+$", "", input_string)
    gene_list <- unlist(strsplit(input_string, "[\n ]+"))
    return(gene_list)
}

#get all specific type of files from a dir
fetch_files <- function(dir, type){
    files = list.files(pattern = paste0('\\.', type, '$'))
    if(length(files)==0){
        return(c())
    }else{
        files <- sapply(files, function(x){
            a = paste(dir, x, sep='')
            return(a)
        })
        return(files)
    }
}

#Replace all dash or invalid symbol included in gene name
single_gene_special_transfer <- function(gene){
    if(grepl("^\\d", gene, perl = T)){
        gene <- paste('X', gene, sep='')
    }else if(grepl('-', gene)){
        gene <- gsub('-','..', gene)
    }
    return(gene)
}

single_gene_special_r_transfer <- function(gene){
    if(grepl('\\.\\.', gene)){
        gene_v <- strsplit(gene, '\\.\\.')[[1]]
        # if(length(gene_v)>2){
        #     gene <- paste(gene_v[1:length(gene_v)-1], collapse='-')
        #     gene <- paste(gene, gene_v[length(gene_v)], sep='.')
        # }else{
        #     gene <- paste(gene_v[1:length(gene_v)], collapse='.')
        # }
        gene <- paste(gene_v[1:length(gene_v)], collapse='-')
    }
    if(grepl("^X\\d", gene, perl = T)){
        gene = substr(gene, 2, nchar(gene))
    }
    return(gene)
}



df_factor_to_num <- function(df){
    indx <- sapply(df, is.factor)
    df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
    return(df)
}


df_num_to_char <- function(df){
    indx <- sapply(df, is.numeric)
    df[indx] <- lapply(df[indx], function(x) as.character(x))
    return(df)
}

