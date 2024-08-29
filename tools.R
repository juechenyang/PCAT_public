# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/5/19

#rename a column name
rename_a_column <- function(df, old, new){
    names(df)[names(df)==old] <- new
    return(df)
}

#function to get a type of file under a directory
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

head_df = function(df){
    print(df[1:10,1:3])
}

keepf4 = function(x){
    a = strsplit(x, '-')[[1]]
    first_four = a[1:4]
    total <- paste(first_four, collapse='-')
    return(total)
}

firstCap = function(x){
    if(is.na(x)){
        return(x)
    }else{
        if(x=="YES"){
            x = tolower(x)
            first_letter = substring(x, 1, 1)
            return(paste0(toupper(first_letter), substring(x,2)))
        }else{
            return(x)
        }
    }
}

all_fac2char = function (df){
    i = sapply(df, is.factor)
    df[i] = lapply(df[i], as.character)
    return(df)
}

all_var2char = function (df){
    df[] = lapply(df, as.character)
    return(df)
}

comparetwo = function (a, b){
    ifelse(toupper(a)==toupper(b), "", b)
}
