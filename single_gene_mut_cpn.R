# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 12/10/19
 #===========================================SINGLE GENE MUTATION=================================
update_mutation_gene <- reactive({
    updateSelectizeInput(session, 'single_mutation_gene_name', server = TRUE, choices=GENE_LIST, selected='')
})
update_mutation_disease <- reactive({
    DB <- chooseDB(input$single_mutation_db)
    clin <- QueryAll(DB$clinical_db_name)
    updateSelectInput(session, 'single_mutation_disease', choices=c('pan-cancer', unique(clin[,DB$disease_tag])))
    return(list(DB=DB, clin=clin))
})

observeEvent(input$single_mut_cpn_example_gene,{
    updateSelectizeInput(session, 'single_mutation_gene_name', server = TRUE, choices=GENE_LIST, selected="CDKN2A")
})
output$mutation_box_plot <- plotly::renderPlotly({

    source('BoxPlot.R')

    shinyjs::hide("mut_cpn_plot_download")
    shinyjs::hide("mut_cpn_data_download")
    #update available genes in the user selection
    update_mutation_gene()

    #update available disease selection based on the selected DB
    obj <- update_mutation_disease()
    DB = obj$DB
    clin = obj$clin

    #validate if there is user input
    validate(
        need(input$single_mutation_gene_name!='', 'please input a gene'),
        errorClass='ValidateInput'
    )

    gene = single_gene_special_transfer(input$single_mutation_gene_name)

    #get expresison data
    exp_data <- process_exp_data(DB, gene, 'single')

    #validate if gene expression is avaliable is the selected DB
    validate(
        need(length(rownames(exp_data))!=0, "gene expression is not available in this DB"),
        errorClass="ErrorOutput"
    )

    #mutation handler
    if(input$mut_or_cpn=='mutation'){
        #Query all mutations regrding to a gene
        gene_mutations <- QueryOnRow(DB$mut_db_name, gene)

        #Validate if the selected gene has mutations in DB
        validate(
            need(length(rownames(gene_mutations))!=0, 'this gene does not have mutations in this DB'),
            errorClass="ErrorOutput"
        )

        #Sort mutation groups by specific order
        # gene_mutations$mut_group = factor(gene_mutations$mut_group, levels = c("Coding_NonSilent", "Coding_Silent", "Noncoding", 'RNA'))

        #remove duplicated sample-mutations to give only one mutation type to a single sample
        gene_mutations = gene_mutations[!duplicated(gene_mutations[c(DB$biospec_id_name)]),]

        #get additional mutation information
        mut_info = QueryAll(DB$mut_info_db)

        gene_mutations = merge(gene_mutations, mut_info, by.x=c("gene_name", "sample_id", "mut_group"), by.y=c("gene_name", "sample_id", "mut_group"), all.x=TRUE)

        #merge expression matrix with mutation matrix
        exp_mut <- merge(exp_data, gene_mutations, by = DB$biospec_id_name)


        #if sample has no mutation on the given gene, assign wild type to the mut_group
        # exp_mut$mut_group <- sapply(exp_mut$mut_group, as.character)
        # exp_mut$mut_group[is.na(exp_mut$mut_group)] <- 'Wild_Type'

        #assign to mutorcpn dataset
        mutorcpn <- exp_mut

        #assign parameter to boxplot group
        group_feature = 'mut_group'

    }
    #copy number handler
    else if(input$mut_or_cpn=='copy number'){

        #Query all copy number change regarding to a gene
        gene_cpn <- QueryOnRow(DB$cpn_db_name, input$single_mutation_gene_name)

        #change cpn value from integer type to factor
        gene_cpn$cpn<-factor(gene_cpn$cpn, levels=c("-2", "-1", "0", "1", "2"))

        #Validate if the selected gene has copy number in DB
        validate(
            need(length(rownames(gene_cpn))!=0, 'this gene is not available for copy number analysis'),
            errorClass="ErrorOutput"
        )

        #remove duplicated sample to give only one cpn value to a single sample
        gene_cpn = gene_cpn[!duplicated(gene_cpn[c(DB$biospec_id_name)]),]

        #merge expression matrix with mutation matrix
        exp_cpn <- merge(exp_data, gene_cpn, by = DB$biospec_id_name)

        #assign to mutorcpn dataset
        mutorcpn <- exp_cpn

        #assign parameter to boxplot group
        group_feature = 'cpn'
    }

    clin = QueryAll(DB$clinical_db_name)

    #merge expression&mutation dataset with clinical matrix for sub disease analysis
    # mutorcpn[[DB$disease_tag]] <- get_target_disease(mutorcpn[[DB$biospec_id_name]])

    mutorcpn = merge(mutorcpn, clin, by=DB$biospec_id_name)

    if(input$mut_or_cpn=='copy number'){
        mutorcpn$cpn = factor(mutorcpn$cpn, levels = c("-2", "-1", "0", "1", "2"))
        mutorcpn = mutorcpn[order(mutorcpn[[DB$disease_tag]],mutorcpn[['cpn']]),]
    }



    #transfer gene name for query
    gene_name = single_gene_special_transfer(input$single_mutation_gene_name)

    if(input$single_mutation_disease=='pan-cancer'){
        if(input$mut_or_cpn=='copy number'){
            final_plot = createBox(mutorcpn, gene_name, DB$disease_tag, group_feature, plot_loc = single_gene_mutcpn_plot_loc, ascending=T, data_loc = single_gene_mutcpn_data_loc,
                                   minimum_box_size=input$mut_min_box_size, drug_plot=T, log_scale=input$mut_log_trans, cpn_plot=T, DB=DB)
        }else{
            final_plot = createBox(mutorcpn, gene_name, DB$disease_tag, group_feature, plot_loc = single_gene_mutcpn_plot_loc, ascending=T, data_loc = single_gene_mutcpn_data_loc,
                                   minimum_box_size=input$mut_min_box_size, mut_plot=T, hover_string = DB$hover_string, log_scale=input$mut_log_trans, DB=DB)
        }
    }else{
        mutorcpn <- mutorcpn[mutorcpn[ ,DB$disease_tag] == input$single_mutation_disease, ]
        if(input$mut_or_cpn=='copy number'){
            final_plot = createBox(mutorcpn, gene_name, NULL, group_feature, plot_loc = single_gene_mutcpn_plot_loc, ascending=T, data_loc = single_gene_mutcpn_data_loc,
                                   minimum_box_size=input$mut_min_box_size, drug_plot=T, log_scale=input$mut_log_trans, cpn_plot=T, DB=DB)
        }else{
            final_plot = createBox(mutorcpn, gene_name, NULL, group_feature, plot_loc = single_gene_mutcpn_plot_loc, ascending=T, data_loc = single_gene_mutcpn_data_loc,
                                   minimum_box_size=input$mut_min_box_size, mut_plot=T, hover_string = DB$hover_string, log_scale=input$mut_log_trans, DB=DB)
        }
    }

    shinyjs::show("mut_cpn_plot_download")
    shinyjs::show("mut_cpn_data_download")

    return(final_plot)
})
