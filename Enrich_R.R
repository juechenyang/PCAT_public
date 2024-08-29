# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-08-20
library(reticulate)

get_enrich_r_result <- function(gene_str, db){
    source_python('request_enrich_result.py')
    result <- get_enrichment_result(gene_str, db)
    #handle the error response exception
    validate(
        need(class(result)=="data.frame", 'Getting bad response from enrich r, please try to type a single letter into the signature box'),
        errorClass="ErrorOutput"
    )
    return(result)
}

