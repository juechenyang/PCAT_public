# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-05-30
createGCT <- function (data, filename, path = "./")
{
    ngene = nrow(data)
    nsample = ncol(data)
    outf = sprintf("%s/%s.gct", path, filename)
    write("#1.2", file = outf)
    write(paste(ngene, nsample, sep = "\t"), file = outf, append = T)
    tmp = c("Name", "Description", colnames(data))
    write(paste(tmp, collapse = "\t"), file = outf, append = T)
    data = cbind(rep("na", nrow(data)), data)
    write.table(data, file = outf, col.names = F, sep = "\t",
        quote = F, append = T)
}

createGMT <- function(gene_list, signature_name="Anoymous"){
    outf = sprintf("IntermediateFiles/genes.gmt")
    write(paste(signature_name, "na", paste0(gene_list, collapse="\t"), sep="\t"), file=outf)
}