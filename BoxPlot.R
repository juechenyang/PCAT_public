library(ggplot2)
library(EnvStats)
library(dplyr)
source("InputProcessor.R")


createBox <- function(rawData, gene_name, disease=NULL, feature=NULL, plot_loc=NULL, data_loc=NULL,
                    ascending=F, drug_plot=F, drug_mut_plot=F, minimum_box_size=1, mut_plot=F, DB,
                    hover_string=NULL, remove_jitter=F, log_scale=F, cpn_plot=F){
    if(feature==''){
        feature <- NULL
    }

    #validate if box plot is qualified to be shown
    validate(
        need(length(rownames(rawData))>0, 'no enough data for box plot'),
        errorClass="ErrorOutput"
    )


    print(paste('total rawData has', as.character(length(rownames(rawData))), 'instances', sep=' '))

    #remove NAs
    rawData[rawData == ""] = NA
    rawData <- rawData[!is.na(rawData[,feature]), ]

    #remove records like unknown or not reported
    rawData <- rawData[!grepl("Unknow|unknown|not reported|Not Reported", rawData[, feature]),]
    print(paste('after remove feature unknown and NA, rawData has', as.character(length(rownames(rawData))), 'instances', sep=' '))

    dplyr::glimpse(rawData)

    #remove records of disease contains NA if disease is NULL
    if(!is.null(disease)){
        rawData <- rawData[!grepl("Unknow|unknown|not reported|Not Reported|N/A", rawData[, disease]),]
        print(paste('after remove disease unknown and NA, rawData has', as.character(length(rownames(rawData))), 'instances', sep=' '))
    }

    #determin whether to perform log transform
    if(log_scale){
        rawData[[gene_name]] <- sapply(rawData[[gene_name]], function(x){
            log2 <- log((x+0.1), 2)
            return(log2)
        })

        ylabel = paste(single_gene_special_r_transfer(gene_name), "Expression log2(FPKM+0.1)")
    }else{
        ylabel = paste(single_gene_special_r_transfer(gene_name), "Expression (FPKM)")
    }

    if(is.null(disease)){
        x_feature<-feature
        fill_feature <- feature
    }else{
        x_feature<-disease
        fill_feature <- feature
        if(feature==disease){
            mode = ''
        }else{
            mode = 'group'
        }
    }

    #remove boxes that does not meet the minimum_box_size
    rawData <- rawData %>% dplyr::group_by_at(c(x_feature, fill_feature)) %>% filter(n() >= minimum_box_size)
    print(paste('after remove groups has very few samples, rawData has', as.character(length(rownames(rawData))), 'instances', sep=' '))
    print(paste('there are/is', as.character(length(unique(rawData[[feature]]))), 'group(s) and thoes are', paste(unique(rawData[[feature]]), collapse=','),sep=' '))

    #validate if box plot is qualified to be shown
    validate(
        need(length(rownames(rawData))>0, 'no enough data for box plot'),
        errorClass="ErrorOutput"
    )

    #determine if the plot is presented in descending or ascending order
    if (ascending){
        exp_data = rawData[[gene_name]]
    }else{
        exp_data = -rawData[[gene_name]]
    }

    #order boxes by median value
    if(!cpn_plot){
        xlabel <- x_feature
        x_feature <- reorder(rawData[[x_feature]], exp_data, FUN = median)
    }

    #color settings for mutation box plot
    #set colors for specific group
    possible_values = c("Missense_Mutation","Silent","Nonsense_Mutation","Splice_Region",
                        "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Translation_Start_Site", "Wild_Type")
    color_settings <- colorRampPalette(c("#FF00E4", "cyan", "orange", "pink", "green", "purple", "yellow",
                                         "#27AE60", "#5B2C6F", "#7E5109", "#7B241C", "#C7C7C7"))(length(possible_values))
    names(color_settings) = possible_values



    if(mut_plot){

        id = 1:length(hover_string)
        tip_v = sapply(id, function(x){
            return(paste0('d', as.character(x)))
        })

        aes_obj = aes_string(x_feature, gene_name,
            d1=hover_string[[1]],
            d2=hover_string[[2]],
            d3=hover_string[[3]],
            d4=hover_string[[4]],
            d5=hover_string[[5]]
        )

        rawData$mut_group <- factor(rawData$mut_group, levels = c("Wild_Type", "Missense_Mutation","Silent","Nonsense_Mutation","Splice_Region",
                        "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Translation_Start_Site"))

        #create box plot object for download
        box_plot <- ggplot(rawData, aes_string(x_feature, gene_name))
        box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(fill = fill_feature, y=ylabel)+
                    theme(legend.position = 'bottom',
                    axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1), axis.title.x = element_blank())+
                    geom_point(shape=19, aes_string(color=fill_feature), position=position_jitter(width=0.12))+guides(linetype = guide_legend(nrow = 2))+stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2, width = 50, height = 50)+
                    scale_color_manual(values=color_settings)
        #save the plot to plot_loc
        ggsave(plot_loc, box_plot)

        #save the data to data_loc
        if(!is.null(data_loc)){
            write.csv(rawData[,c(DB$biospec_id_name, DB$disease_tag, gene_name, feature)], file=data_loc, row.names=FALSE)
        }

        #create box plot object for display
        box_plot <- ggplot(rawData,aes_obj)
        box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(fill = fill_feature, y=ylabel)+
                    theme(legend.position = 'bottom',
                    axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1), axis.title.x = element_blank())+
                    geom_point(shape=19, aes_string(color=fill_feature), position=position_jitter(width=0.12))+guides(linetype = guide_legend(nrow = 2))+stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)+
                    scale_color_manual(values=color_settings)
        plotly_boxplot <- plotly::ggplotly(box_plot, tooltip=tip_v)
    }else if(cpn_plot){

        #set colors for specific group
        possible_values = c("-2", "-1", "0", "1", "2")
        # possible_values = c(-2, -1, 0, 1, 2)
        color_settings <- colorRampPalette(c("cyan", "orange", "pink", "purple", "red"))(length(possible_values))
        names(color_settings) = possible_values

        #order cpn boxes by ascending value of cpn
        #rawData[[feature]] <- factor(rawData[[feature]], levels = c("-2", "-1", "0", "1", "2"))
        #rawData = rawData[order(rawData[[feature]]),]sssss

        #construct box plot for download
        box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name, fill=fill_feature))
        box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(fill = fill_feature, y=ylabel)+
                    theme(legend.position = 'bottom',
                    axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1), axis.title.x = element_blank())+
                    guides(linetype = guide_legend(nrow = 2))+stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)+
                    scale_color_manual(values=color_settings)+scale_fill_manual(values=color_settings)

        #if do not remove jitter, then add jitter points
        if(!remove_jitter){
            box_plot = box_plot + geom_point(shape=21, position = position_jitterdodge(jitter.width = 0.05))
        }
        #save the plot into the file location
        ggsave(plot_loc, box_plot)

        #save the data to data_loc
        if(!is.null(data_loc)){
            write.csv(rawData[,c(DB$biospec_id_name, DB$disease_tag, gene_name, feature)], file=data_loc, row.names=FALSE)
        }

        #create box plot for display
        box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name, fill=fill_feature, sid='sample_id'))
        box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(fill = fill_feature, y=ylabel)+
                    theme(legend.position = 'bottom',
                    axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1), axis.title.x = element_blank())+
                    # geom_point(shape=21, position = position_jitterdodge(dodge.width=0.7))+guides(linetype = guide_legend(nrow = 2))+
                    # guides(linetype = guide_legend(nrow = 2))+
                    stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)+scale_fill_manual(values=color_settings)
        if(!remove_jitter){
            box_plot = box_plot + geom_point(shape=21, position = position_jitterdodge(jitter.width = 0.05))
        }
        #render the box plot into plotly for visualization
        plotly_boxplot <- plotly::ggplotly(box_plot, tooltip=c('sid'))%>%layout(boxmode = mode)
    }else{

        #drug plot settings
        if(drug_plot){

            color_settings <- colorRampPalette(c("#FF00E4", "cyan", "orange", "pink", "green", "purple", "yellow",
                                         "#27AE60", "#5B2C6F", "#7E5109", "#7B241C", "#000000"))(length(possible_values))
            names(color_settings) = possible_values


            #replace the input x feature string to a valid plot string
            replace_char = "Yang_special"
            colnames(rawData)[colnames(rawData)==xlabel] <- replace_char
            x_feature = replace_char
            fill_feature = 'mut_group'

            #construct hovering strings
            id = 1:length(hover_string)
            tip_v = sapply(id, function(x){
                return(paste0('d', as.character(x)))
            })

            aes_obj = aes_string(x_feature, gene_name,
                d1=hover_string[[1]],
                d2=hover_string[[2]],
                d3=hover_string[[3]],
                d4=hover_string[[4]],
                d5=hover_string[[5]]
            )

            #add level to xlabel
            xlabel = paste(xlabel, 'Response', 'Level')
          
            box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name))
            box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(x=xlabel, y=ylabel)+
                    theme(legend.position = 'bottom', legend.title = element_blank(), axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1))+
                    guides(linetype = guide_legend(nrow = 2))+
                    stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)

            #if mutation is selected in preclinical plot then add geom points
            if(drug_mut_plot){


                  rawData$mut_group <- factor(rawData$mut_group, levels = c("Wild_Type", "Missense_Mutation","Silent","Nonsense_Mutation","Splice_Region",
                                       "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Translation_Start_Site"))
              
                  box_plot = box_plot+geom_point(shape=19, aes_string(color=fill_feature), position=position_jitter(width=0.12))+
                    guides(linetype = guide_legend(nrow = 2))+
                    stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2, width = 50, height = 50)+
                    scale_color_manual(values=color_settings)
            }else{
                box_plot <- box_plot+geom_point(shape=19, position=position_jitter(width=0.12))
            }
            

            #save the plot into the file location
            ggsave(plot_loc, box_plot)

            if(drug_mut_plot){
                  box_plot <- ggplot(rawData, aes_obj)
                  box_plot = box_plot+geom_point(shape=19, aes_string(color=fill_feature), position=position_jitter(width=0.12))+
                    guides(linetype = guide_legend(nrow = 2))+
                    stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2, width = 50, height = 50)+
                    scale_color_manual(values=color_settings)
                tooltip = tip_v
            }else{
                box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name, sid='sample_id'))+
                geom_point(shape=19, position=position_jitter(width=0.12))
                tooltip = c('sid')
            }
            box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(x=xlabel, y=ylabel)+
                    theme(legend.position = 'bottom', legend.title = element_blank(), axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1))+
                    guides(linetype = guide_legend(nrow = 2))+
                    stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)

             #render the box plot into plotly for visualization
             plotly_boxplot <- plotly::ggplotly(box_plot, tooltip=tooltip)%>%layout(boxmode = mode)
        }else{
            #construct box plot for download
            box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name, fill=fill_feature))
            box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(x=xlabel, y=ylabel)+
                        theme(legend.position = 'bottom', axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1),
                        axis.title.x = element_blank())+ guides(linetype = guide_legend(nrow = 2))+
                        stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)

            #if do not remove jitter, then add jitter points
            if(!remove_jitter){
                box_plot = box_plot + geom_point(shape=21, position = position_jitterdodge(jitter.width = 0.05))
            }

            #save the plot into the file location
            ggsave(plot_loc, box_plot)

            #create box plot for display
            box_plot <- ggplot(rawData, aes_string(x=x_feature, y=gene_name, fill=fill_feature, sid='sample_id'))
            box_plot <- box_plot+geom_boxplot(outlier.shape = NA)+labs(x=xlabel, y=ylabel)+
                        theme(legend.position = 'bottom',
                        axis.text.x = element_text(face="bold", color="#d16a0a", angle = 60, hjust = 1), axis.title.x = element_blank())+
                        # geom_point(shape=21, position = position_jitterdodge(dodge.width=0.7))+guides(linetype = guide_legend(nrow = 2))+
                        # guides(linetype = guide_legend(nrow = 2))+
                        stat_n_text(y.pos=max(rawData[[gene_name]])*1.1, size = 2)
            if(!remove_jitter){
                box_plot = box_plot + geom_point(shape=21, position = position_jitterdodge(jitter.width = 0.05))
            }
            #render the box plot into plotly for visualization
            plotly_boxplot <- plotly::ggplotly(box_plot, tooltip=c('sid'))%>%layout(boxmode = mode)


        }




    }

    #remove outliers from plotly
    for(i in 1:length(plotly_boxplot$x$data)){
        if(!is.null(plotly_boxplot$x$data[[i]]$marker$outliercolor)){
            plotly_boxplot$x$data[[i]]$marker$opacity=0
        }
    }

    if(drug_plot){
        colnames(rawData)[colnames(rawData)=="Yang_special"] <- feature
    }else{

    }
    #save the data to data_loc
    if(!is.null(data_loc)){
        write.csv(rawData[,c(DB$biospec_id_name, DB$disease_tag, gene_name, feature)], file=data_loc, row.names=FALSE)
    }


    return(plotly_boxplot)
}


