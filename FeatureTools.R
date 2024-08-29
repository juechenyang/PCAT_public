# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-05-07

library(stringr)
#determine if a variable is continuous
is_continuous <- function(var){
    uniq_var <- unique(var)
    if (length(uniq_var) >= 5 ){
        return(TRUE)
    }else{
        return(FALSE)
    }
}
get_ava_features <- function(db_table, DB, mode){
    #feature filter to get the available features
    available_features <- names(db_table)
    #remove sampleID column in the clinical matrix"
    # available_features <- available_features[!str_detect(available_features, DB$clinical_id_name)]
    disqualified_features <- vector()
    potential_unknown <- c('Unknow','Other','unknown')
    if(mode=='single'){
        for(f in available_features){
            feature_vector <- db_table[,f]
            #remove na and unknown
            feature_vector <- feature_vector[!is.na(feature_vector)]
            feature_vector <- feature_vector[!(feature_vector %in% potential_unknown)]

            #find disqualified_features
            if (!is.numeric(feature_vector)){
                if (length(unique(feature_vector))>=length(rownames(db_table))*0.3 | length(unique(feature_vector)) == 1){
                    disqualified_features <- append(disqualified_features, f)

                }
            }

        }
    }else{
        for(f in available_features){
            feature_vector <- db_table[,f]
            #remove na and unknown
            feature_vector <- feature_vector[!is.na(feature_vector)]
            feature_vector <- feature_vector[!(feature_vector %in% potential_unknown)]
            #find disqualified_features
            if (!is.numeric(feature_vector)){
                if (length(unique(feature_vector))>=length(rownames(db_table))*0.2 | length(unique(feature_vector)) == 1){
                    disqualified_features <- append(disqualified_features, f)
                }
            }else{
                disqualified_features <- append(disqualified_features, f)
                available_features <- append(available_features, paste(f, '_discretized',sep=''))
                # db_table[[paste(f, '_discretized',sep='')]] <- discretize_feature(as.numeric(db_table[,f]), groups)
            }
        }
    }
    available_features <- available_features[!(available_features %in% disqualified_features)]
}

#remove invalid value in a vector
remove_invalid = function(vx){
    invalid_set = c('Unknown', 'N/A', 'Unknow', 'unknown', 'not reported', 'Not Reported', 'unknown')
    vx = vx[!vx %in% invalid_set]
    return(vx)
}

