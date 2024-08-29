# 1st argument: function_type: choose whether get cor between two genes or given one genes get top numbers of correlated genes
# other arguments: from ui$input
run_single_cor <- function(method, disease, gene1, gene2=NULL, DB, log_scale=F, plot_loc, data_loc){

    #transfer the gene name into a code friendly name
    g1 = single_gene_special_transfer(gene1)
    g2 = single_gene_special_transfer(gene2)

    #get the expression data
    exp_data <-process_exp_data(DB, c(g1,g2), mode='multiple')


    #log scale converter
    # if(log_scale){
    #     exp_data[, g1] = sapply(exp_data[,g1]+0.1, log, 2)
    #     exp_data[, g2] = sapply(exp_data[,g2]+0.1, log, 2)
    # }

    #validate if the gene is available in the selected DB
    validate(
        need(g1 %in% colnames(exp_data), paste(gene1, 'is not available is this DB')),
        errorClass="ErrorOutput"
    )
    validate(
        need(g2 %in% colnames(exp_data), paste(gene2, 'is not available is this DB')),
        errorClass="ErrorOutput"
    )

    #get disease column from clinical table
    clin = QueryAll(DB$clinical_db_name)
    exp_data = merge(exp_data, clin, by=DB$biospec_id_name)


    if(disease != "pan-cancer"){
        exp_data <- exp_data[exp_data[ ,DB$disease_tag] == disease, ]
    }

    validate(
        need(length(rownames(exp_data))>=5, paste('no enough data for computing correlation')),
        errorClass="ErrorOutput"
    )

    scatter_plot = return(createScatter(exp_data, g1, g2, method=method, two_gene_mode=T,
                                        log_scale=log_scale, DB=DB, plot_loc = plot_loc, data_loc=data_loc))



    return(scatter_plot)

}

