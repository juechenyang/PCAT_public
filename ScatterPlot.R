library(ggplot2)
source("InputProcessor.R")

lm_eqn <- function(df, y, x, method){
    R_val = format(cor(df[,x], df[,y], method = method), digits=2)
    pvalue = format(cor.test(df[,x], df[,y], method = method, exact = FALSE)$p.value, digits=2)
    n_samples = length(rownames(df))
    return(c(R_val, pvalue, method, n_samples))
}
save_scatter_Plot <- function(scatter_plot, stat, loc, xpos, ypos){
    scatter_plot <- scatter_plot+annotate(geom="text", x=xpos, y=ypos,
            label=paste("R =", stat[1], 'p-value =', stat[2], '\nmethod:', stat[3], '\n#samples',stat[4], sep=" "),
            color="red", size=3)
    ggsave(loc, scatter_plot)
}


#Regular Scatter Plot function
createScatter <- function(rawData,gene_name, feature, plot_loc=NULL,
                          data_loc=NULL, method='spearman',
                          two_gene_mode=FALSE, log_scale=F, DB){

    #removing NAs for both vectors
    rawData <- rawData[!is.na(rawData[,feature]), ]
    rawData <- rawData[!is.na(rawData[,gene_name]), ]

    ylabel = paste(single_gene_special_r_transfer(gene_name), "Expression (FPKM)")

    if(two_gene_mode){
        xlabel = paste(single_gene_special_r_transfer(feature), "Expression (FPKM)")
    }else{
        xlabel = paste(single_gene_special_r_transfer(feature))
    }

    #if log scale applied covert all expression numbers to log scale
    if(log_scale){
        if(two_gene_mode){
            rawData[,c(gene_name, feature)] = lapply(rawData[,c(gene_name, feature)]+0.1, log, 2)
            xlabel = paste(single_gene_special_r_transfer(feature), "Expression log2(FPKM+0.1)")
        }else{
            rawData[,gene_name] = sapply(rawData[,gene_name]+0.1, log, 2)
            xlabel = paste(single_gene_special_r_transfer(feature))
        }
        ylabel = paste(single_gene_special_r_transfer(gene_name), "Expression log2(FPKM+0.1)")
    }

    adjusted_feature <- rawData[,feature]
    selected_gene <- rawData[,gene_name]
    xpos <- max(adjusted_feature)*0.82222
    ypos <- max(selected_gene)*0.9

    stat <- lm_eqn(rawData, gene_name, feature, method)
    p <- ggplot(rawData, aes_string(x= feature, y= gene_name, id='sample_id'))+
    labs(
        x=xlabel,
        y=ylabel
    )+
    geom_point(shape=21) +    # Use hollow circles
    geom_smooth(method=lm)   # Add linear regression line
    if(!is.null(plot_loc)){
        save_scatter_Plot(p, stat, plot_loc, xpos, ypos)
    }
    if(!is.null(data_loc)){
        write.csv(rawData[,c(DB$biospec_id_name, DB$disease_tag, gene_name, feature)], file=data_loc, row.names=FALSE)
    }
    scatter_plot <- ggplotly(p, tooltip=c('id')) %>%
    layout(annotations = list(x = xpos, y = ypos,
    text = paste("R =", stat[1], 'p-value =', stat[2], '<br>method:', stat[3], '<br>#samples',stat[4], sep=" "), showarrow = FALSE))
    return(scatter_plot)
}

#function to create fusion scatter plot
fusionScatter = function(data, g1, g2, sample_id, cohort, DB, log_scale=F){

    #define the legend name of fusion scatter plot
    selected_sample_tag = "SELECTED SAMPLE"
    same_fusion_sample_tag = "SAMPLES HAVING QUERIED FUSION"

    #removing NAs for both vectors
    data <- data[!is.na(data[,g1]), ]
    data <- data[!is.na(data[,g2]), ]

    #assign group tags for each sample
    data$group <- ifelse(
        data[['donor_gene']] == g1 & data[['acceptor_gene']] == g2, same_fusion_sample_tag,cohort
    )
    data$group <- ifelse(
        data[[DB$biospec_id_name]] == sample_id, selected_sample_tag, data$group
    )

    #create labels for plot
    xlabel = paste(single_gene_special_r_transfer(g1), "Expression (FPKM)")
    ylabel = paste(single_gene_special_r_transfer(g2), "Expression (FPKM)")

    #apply log if log_scale is True
    if(log_scale){
        data[,c(g1, g2)] = lapply(data[,c(g1, g2)]+0.1, log, 2)
        xlabel = paste(single_gene_special_r_transfer(g1), "Expression log2(FPKM+0.1)")
        ylabel = paste(single_gene_special_r_transfer(g2), "Expression log2(FPKM+0.1)")
    }

    #keep samples unique so that it only reflects
    data$group <- factor(data$group, levels = c(selected_sample_tag, same_fusion_sample_tag, cohort))
    data = data[order(data[[DB$biospec_id_name]], data$group),]

    #remove duplicated sample id
    data = data[!duplicated(data[c(DB$biospec_id_name)]),]

    #grab expression vectors using valid gene symbol
    gene1_v = data[[g1]]
    gene2_v = data[[g2]]

    #set colors for specific group
    unique_values <- as.vector(unlist(unique(data$group)))
    unique_values <- factor(unique_values,levels=c(cohort, same_fusion_sample_tag, selected_sample_tag))
    unique_values = sort(unique_values)
    color_settings <- colorRampPalette(c("orange", "blue", "red"))(length(unique_values))
    names(color_settings) = unique_values

    #adjust the sequence to make the selected sample showing at the top layer
    data$group <- factor(data$group, levels = c(cohort, same_fusion_sample_tag, selected_sample_tag))

    #make the ggplot
    plot = ggplot(data, aes_string(x= gene1_v, y=gene2_v, sid="sample_id"))+labs(x=xlabel,y=ylabel)+
        geom_point(aes(colour = group), size = 2)+scale_color_manual(values=color_settings)

    #convert the ggplot to interactive plotly
    plotly <- plotly::ggplotly(plot, tooltip=c('sid'))%>%layout(boxmode = mode)

    #return plotly object
    return(plotly)
}




