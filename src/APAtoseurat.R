

library(argparse)
library(data.table)
library(dplyr)
library(Seurat)


#ArgumentParser
#==============
parser = ArgumentParser(description='Process counts file')
parser$add_argument('APA_PATH', type="character", help='path of APA counts file')
parser$add_argument('sample_name', type="character", help='sample_name')
parser$add_argument('OUTPUT', type="character", help='path of output seurat file')
args = parser$parse_args()

<<<<<<< HEAD
MIN.CELLS=0
MIN.FEATURES=0

#file opening and writing
tab = read.table(args$APA_PATH, sep="\t", header = T, row.names = 1)
myobj = Seurat::CreateSeuratObject(tab, names.field = 1, min.cells = MIN.CELLS, min.features = MIN.FEATURES, project = args$sample_name)
saveRDS(myobj, file=paste0(args$sample_name,"_seurat.RDS"))
=======

#file opening and writing
utils::read.table(args$APA_PATH, sep="\t", header = T, row.names = 1) %>% 
    Seurat::CreateSeuratObject(names.field = 1, min.cells = 0, min.features = 0, project = args$sample_name) %>%
    saveRDS(file=paste0(args$sample_name,"_seurat.RDS"))
>>>>>>> main

