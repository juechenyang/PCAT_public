# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 11/1/19
library(ComplexHeatmap)
library(circlize)
library(shiny)
library(corrplot)

#heapmap function for pairwise correlation
Pairwise_Heatmap = function(cor_matrix, exp_df, loc){
    col_fun = colorRampPalette(c("blue", 'white', 'red'))(20)
    p.mat <- corrplot::cor.mtest(exp_df)

    #save as pdf
    pdf(loc, height=12, width=12)
    heatmap_obj = corrplot(cor_matrix, p.mat = p.mat$p, col=col_fun, type="upper", insig = "blank", tl.cex=0.8)
    dev.off()

    #draw the heatmap obj and return the real heatmap
    heatmap_obj = corrplot(cor_matrix, p.mat = p.mat$p, col=col_fun, type="upper", insig = "blank", tl.cex=0.8,
    cl.cex = 1)
    return(heatmap_obj)


    #***********Used For Debug**************
    # heatmap_plot = draw(heatmap_obj)
    # return(invisible(NULL))
    # heatmap_obj = corrplot::corrplot(matrix, p.mat = res1$p, col = col_fun,
    # type="upper", pch.col = "white", insig = "label_sig",order="hclust", pch.cex=2, tl.pos="d",tl.srt=60)

}

#heapmap function for multi-genes analysis
Multigene_Heatmap = function(df, gene_list, feature, row_cluster,
                             column_cluster, plot_loc, data_loc, DB){

    #df: a comprehensive dataframe that include all info
    #gene_list: all valid input genes
    #feature: the selected parameter

    #make all empty values to NA
    df[df==""] = NA

    #remove NAs
    df <- df[!is.na(df[,feature]), ]

    #remove unknowns
    df <- df[!grepl('Unknow|unknown|not reported|Not Reported|Not noted', df[, feature]),]

    validate(
        need(length(rownames(df))>1, "too less data to plot"),
        errorClass="ErrorOutput"
    )

    #sort the data by feature
    df <- df[order(df[[feature]]),]

    #get the number of unique values for feature
    unique_values <- as.vector(unlist(unique(df[feature])))
    colors <- colorRampPalette(c("green", "yellow", "pink", "purple", "orange"))(length(unique_values))
    names(colors) = unique_values



    #transpose the input df so that gene name is on row side
    processed_data <-  scale(df[, gene_list])
    heatmap_matrix = data.matrix(t(processed_data))


    #save multi genes clinical data
    write.csv(df[,c(DB$biospec_id_name, DB$disease_tag, gene_list, feature)],
              file=data_loc, row.names=FALSE)

    #if feature is numeric, top annotation is a barplot
    if(!any(!grepl('^-?[0-9.]+$', df[[feature]]))){
        top_anno = HeatmapAnnotation(
            clinical_character = anno_barplot(as.numeric(df[[feature]]),
                axis_param = list(
                    side = "right"
                ),
                gp = gpar(fill = "#2194E0"), bar_width = 1,
                height = unit(3, "cm")
            ),
            show_annotation_name = FALSE
        )
    }else{
        top_anno = HeatmapAnnotation(
            clinical_character = df[[feature]],
            simple_anno_size = unit(2, "cm"),
            col = list(clinical_character = colors),
            annotation_legend_param = list(
                clinical_character = list(
                    title = feature,
                    grid_height = unit(6, "mm"),
                    grid_width = unit(6, "mm")
                    # labels_gp = gpar(fontsize = 15),
                    # title_gp = gpar(fontsize = 15, fontface = "bold"),
                    # title_gap = unit(20, "mm"),
                    # title_position = "lefttop-rot",
                )
            ),
            show_annotation_name = FALSE
        )
    }

    #define the range to show the heatmap
    col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

    #create heatmap object
    heatmap_obj = Heatmap(heatmap_matrix, name = "Expression",
        column_dend_height = unit(3, "cm"),
        row_dend_width = unit(3, "cm"),
        col = col_fun,
        cluster_rows = row_cluster,
        cluster_columns = column_cluster,
        show_column_names = FALSE,
        top_annotation = top_anno,
        heatmap_legend_param = list(
            labels_gp = gpar(fontsize = 15),
            title_gp = gpar(fontsize = 15, fontface = "bold"),
            grid_width = unit(15, "mm"),
            grid_height = unit(10, "mm"),
            title_gap = unit(30, "mm"),
            title_position = "lefttop-rot"
            # direction = "horizontal"
        )
    )

    #save as pdf for download
    pdf(plot_loc, height=9, width=12)
    ht <- draw(heatmap_obj,
        heatmap_legend_side = "left",
        annotation_legend_side = "right"
    )
    decorate_annotation("clinical_character",
        {grid.text(feature, x = unit(-2, "mm"),just = "right", gp = gpar(fontsize = 8))}
    )
    dev.off()
    #draw the heapmap object to a real heatmap
    ht <- draw(heatmap_obj,
        heatmap_legend_side = "left",
        annotation_legend_side = "right"
    )
    decorate_annotation("clinical_character",
        {grid.text(feature, x = unit(-2, "mm"),just = "right", gp = gpar(fontsize = 8))}
    )
    return(invisible(NULL))
}

#heatmap function for ssGSEA
ssGSEA_Heatmap <- function(df, gene_list, feature, loc, drug=NULL){

    #make empty entries to NAs
    df[df==""] = NA

    #remove NAs
    df <- df[!is.na(df[,feature]), ]

    #remove unknowns
    df <- df[!grepl('Unknow|unknown|not reported|Not Reported|Not noted', df[, feature]),]

    #if number of observation from a raw dataframe is less than 2, print error message
    validate(
        need(length(rownames(df))>=2, "too less data to plot"),
        errorClass="ErrorOutput"
    )

    # sort by ssGSEA scores
    df <- df[order(df[,'ssGSEA_scores']),]

    #assign colors to the unique values of the feature variable
    unique_values <- as.vector(unlist(unique(df[feature])))
    colors <- colorRampPalette(c("green", "yellow", "pink", "purple", "orange"))(length(unique_values))
    names(colors) <- unique_values

    if(!is.null(drug)){
        #define colors for drug response
        drug_response_unique_values <- c("NA", "PD1", "PD2", "SD", "PR", "CR", "MCR")
        drug_colors <- colorRampPalette(c("#cecece", "palegreen4",
                       "chartreuse3", "black", "yellow", "orange", "tomato"))(length(drug_response_unique_values))
        names(drug_colors) <- drug_response_unique_values
    }


    #scale the expression matrix
    processed_data <-  scale(df[, gene_list])

    # create plot matrix to put genes on the row side and samples on the column side
    sub = data.matrix(t(processed_data))


    #if feature is numeric, top annotation is a barplot
    if(!any(!grepl('^-?[0-9.]+$', df[[feature]]))){
        if(!is.null(drug)){
            top_anno = HeatmapAnnotation(
                ssGSEA_score = anno_barplot(as.numeric(df[,"ssGSEA_scores"]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                clinical_feature = anno_barplot(as.numeric(df[[feature]]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                drug_feature = df[,drug],
                simple_anno_size = unit(1, "cm"),
                col = list(drug_feature = drug_colors),
                annotation_legend_param = list(
                    drug_feature = list(
                        title = drug
                    )
                ),
                show_annotation_name = FALSE
            )
        }else{
            top_anno = HeatmapAnnotation(
                ssGSEA_score = anno_barplot(as.numeric(df[,"ssGSEA_scores"]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                clinical_feature = anno_barplot(as.numeric(df[[feature]]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                show_annotation_name = FALSE
            )
        }

    }else{
        if(!is.null(drug)){
            top_anno = HeatmapAnnotation(
                ssGSEA_score = anno_barplot(as.numeric(df[,"ssGSEA_scores"]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                clinical_feature = df[,feature],
                drug_feature = df[,drug],
                simple_anno_size = unit(1, "cm"),
                col = list(clinical_feature = colors, drug_feature = drug_colors),
                annotation_legend_param = list(
                    clinical_feature = list(
                        title = feature,
                        grid_height = unit(6, "mm"),
                        grid_width = unit(6, "mm")
                        # labels_gp = gpar(fontsize = 15),
                        # title_gp = gpar(fontsize = 15, fontface = "bold"),
                        # title_gap = unit(20, "mm"),
                        # title_position = "lefttop-rot",
                    ),
                    drug_feature = list(
                        title = drug
                    )
                ),
                show_annotation_name = FALSE
            )
        }else{
            top_anno = HeatmapAnnotation(
                ssGSEA_score = anno_barplot(as.numeric(df[,"ssGSEA_scores"]),
                    axis_param = list(
                        side = "right"
                    ),
                    gp = gpar(fill = "#2194E0"), bar_width = 1,
                    height = unit(3, "cm")
                ),
                clinical_feature = df[,feature],
                simple_anno_size = unit(1, "cm"),
                col = list(clinical_feature = colors),
                annotation_legend_param = list(
                    clinical_feature = list(
                        title = feature,
                        grid_height = unit(6, "mm"),
                        grid_width = unit(6, "mm")
                        # labels_gp = gpar(fontsize = 15),
                        # title_gp = gpar(fontsize = 15, fontface = "bold"),
                        # title_gap = unit(20, "mm"),
                        # title_position = "lefttop-rot",
                    )
                ),
                show_annotation_name = FALSE
            )
        }

    }

    #define the range to show the heatmap
    col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

    #create heatmap obj
    heat_map <- Heatmap(sub, name = "Expression",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        top_annotation = top_anno,
        heatmap_legend_param = list(
            labels_gp = gpar(fontsize = 15),
            title_gp = gpar(fontsize = 15, fontface = "bold"),
            grid_width = unit(15, "mm"),
            grid_height = unit(10, "mm"),
            title_gap = unit(30, "mm"),
            title_position = "lefttop-rot"
            # direction = "horizontal"
        )
    )

    #save pdf for download
    pdf(loc,height=9, width=12)
    if(is.null(drug)){
        ht <- draw(heat_map,
            column_title = "ssGSEA_score Heatmap",
            heatmap_legend_side = "left",
            annotation_legend_side = "right"
        )
        decorate_annotation("clinical_feature",
            {grid.text(feature, x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("ssGSEA_score",
            {grid.text("ssGSEA_score", x = unit(-5, "mm"),just = "right")}
        )
    }else{
        ht <- draw(heat_map,
            column_title = "ssGSEA_score Heatmap",
            heatmap_legend_side = "left",
            annotation_legend_side = "right"
        )
        decorate_annotation("clinical_feature",
            {grid.text(feature, x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("ssGSEA_score",
            {grid.text("ssGSEA_score", x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("drug_feature",
            {grid.text(drug, x = unit(-5, "mm"),just = "right")}
        )
    }

    dev.off()



    #draw the heatmap and decorate the annotation
    if(is.null(drug)){
        ht <- draw(heat_map,
            column_title = "ssGSEA_score Heatmap",
            heatmap_legend_side = "left",
            annotation_legend_side = "right"
        )
        decorate_annotation("clinical_feature",
            {grid.text(feature, x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("ssGSEA_score",
            {grid.text("ssGSEA_score", x = unit(-5, "mm"),just = "right")}
        )
    }else{
        ht <- draw(heat_map,
            column_title = "ssGSEA_score Heatmap",
            heatmap_legend_side = "left",
            annotation_legend_side = "right"
        )
        decorate_annotation("clinical_feature",
            {grid.text(feature, x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("ssGSEA_score",
            {grid.text("ssGSEA_score", x = unit(-5, "mm"),just = "right")}
        )
        decorate_annotation("drug_feature",
            {grid.text(drug, x = unit(-5, "mm"),just = "right")}
        )
    }


    return(invisible(NULL))
}
