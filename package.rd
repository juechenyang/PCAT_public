# this doc is used for listing dependencies of pcat web-app

#channel r

conda install -c r r-shiny r-ggplot2 r-rmysql r-survival r-dplyr r-stringr r-tidyr r-tibble r-reticulate r-dt r-feather

#channel conda-forge
conda install -c conda-forge r-ggpubr r-ggpmisc r-maxstat r-survminer r-plotly r-envstats r-shinyjs

install.packages('ComplexHeatmap')

web server ip:10.248.115.166
MySQL server ip:10.248.115.165

sudo systemctl start mysqld