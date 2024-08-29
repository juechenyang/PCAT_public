# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2/19/20

library(circlize)
library(grid)


drawFusionLinks = function(x, color){

    circos.link(x["donor_chr"], as.integer(x["LeftBreak"]), x["acceptor_chr"], as.integer(x["RightBreak"]), col=color)
}

createCircos = function(sample_fusion_df, gene1, gene2){

    #clear some previous figures
    circos.clear()

    #change the start point to the top center
    circos.par("start.degree" = 90)

    #initialize the genome plot
    circos_plot = circos.initializeWithIdeogram(plotType = c("ideogram", "labels"), ideogram.height = convert_height(10, "mm"))

    #draw links for all fusions
    apply(sample_fusion_df, 1, drawFusionLinks, color='black')

    #filter to dataframe to selected fusion
    sample_fusion_df = sample_fusion_df[(sample_fusion_df[['donor_gene']]==gene1) & (sample_fusion_df[['acceptor_gene']]==gene2),]

    #draw link for selected fusion
    apply(sample_fusion_df, 1, drawFusionLinks, color='red')

    #clear the circos
    circos.clear()

    #create legends of plot
    lgd_lines = ComplexHeatmap::Legend(at = c(paste0(gene1, '-', gene2), "Other fusions in this sample"), type = "lines",
    legend_gp = gpar(col = c(2,1), lwd = 2), labels_gp = gpar(fontsize = 12), title_position = "topleft", title_gp = gpar(fontsize = 15, fontface = "bold"),
                    grid_width = unit(15, "mm"), grid_height = unit(10, "mm"))

    #pack the legend
    legend_set = ComplexHeatmap::packLegend(lgd_lines)

    #draw the legend
    ComplexHeatmap::draw(legend_set, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

    #return plot object for display
    return(circos_plot)
}






