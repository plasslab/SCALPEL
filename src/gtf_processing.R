#!/usr/bin/env Rscript

# =============================================================================
# Script: gtf_processing.R
# Description: Collapse isoforms based on exon position similarity and bulk expression
# Outputs: GTF with exon model (with collapsed isoforms) and unique isoform gene list
# =============================================================================

#libs
library(argparse)
library(dplyr)
library(data.table)

#ArgumentParser
parser = ArgumentParser(description='Process GTF file')
parser$add_argument('GTF_PATH', type="character", help='path of gtf file')
parser$add_argument('QF_PATH', type="character", help='path of salmon file')
parser$add_argument('DIST_TR', type="character", help='transcriptomic distance range threshold')
parser$add_argument('DIST_TR_EXON', type="character", help='transcripts inter distance threshold')
parser$add_argument('OUTPUT', type="character", help='path of output gtf file')
parser$add_argument('OUTPUT_UNIQUE', type="character", help='path of output gtf file unique isoform')
args = parser$parse_args()


GTF_PATH = args$GTF_PATH
QF_PATH = args$QF_PATH
DT_THRESHOLD = as.numeric(args$DIST_TR)
DT_EX_THRESHOLD = as.numeric(args$DIST_TR_EXON)


#Functions
#----------

groupExon_ends = function(x, distance_exon){
  #Function for grouping exon with similar ends
  #=============================================
  if (is.unsorted(x)) {
    idx <- order(x)
    x <- x[idx]
  } else idx <- integer()
  res <- cumsum(c(1L, diff(x) > distance_exon))
  res[idx] <- res
  res
}

collapseIsoforms = function(x, distance_transcriptome, distance_ex){
  #function for collapsing Isoforms
  #================================

  #get isoforms with similar ends
  strand = x$strand[1]
  x.tmp = x %>%
    dplyr::filter(exon_number==1)
  x.tmp = x.tmp %>%
    dplyr::mutate(clusters = ifelse(strand=="-",groupExon_ends(start, distance_ex),groupExon_ends(end, distance_ex)))
  simIsoforms = group_by(x.tmp, clusters) %>% dplyr::filter(n()>1)
  simIsoforms = lapply(split(simIsoforms, simIsoforms$clusters), function(x) x$transcript_name)
  x = suppressMessages(dplyr::left_join(x, distinct(x.tmp, transcript_name, clusters)))

  #filter input tab
  collapsed = data.frame(collapsed=character(0))
  #loop on similar isoforms
  for(trs in simIsoforms){
    curr = dplyr::filter(x, transcript_name %in% trs)
    max_exon = max(curr$exon_number)
    for(ind in unique(curr$exon_number)){
      current = curr %>% dplyr::filter(exon_number==ind)
      #in last exons
      if(ind==1){
        if(strand=="-"){
          if(ind==max_exon){current = dplyr::mutate(current, end=max(end))}
          current = mutate(current, end= ifelse(endR==distance_transcriptome-1,max(end),end))
          current.2 = group_by(current, clusters, end) %>%
            mutate(collapsed = paste(transcript_name, collapse = "_"), nCollapsed = n_distinct(transcript_name)) %>%
            dplyr::filter(nCollapsed!=1) %>% ungroup()
        }else{
          if(ind==max_exon){current = dplyr::mutate(current, start=min(start))}
          current = mutate(current, start= ifelse(endR==distance_transcriptome-1,min(start),start))
          current.2 = group_by(current, clusters, start) %>%
            mutate(collapsed = paste(transcript_name, collapse = "_"), nCollapsed = n_distinct(transcript_name)) %>%
            dplyr::filter(nCollapsed!=1) %>% ungroup()
        }
        #breaking rule
        if(nrow(current.2)==0){break}
        if(ind == max(curr$exon_number)){
          collapsed = rbind(collapsed, distinct(current.2, collapsed))
        }
      }else{
        if(strand=="-"){
          if(ind==max_exon){current = dplyr::mutate(current, end=max(end))}
          current.2 = group_by(current,clusters,start) %>%
            mutate(end= ifelse(n()>1 & endR==distance_transcriptome-1,max(end),end)) %>% ungroup()
        }else{
          if(ind==max_exon){current = dplyr::mutate(current, start=min(start))}
          current.2 = group_by(current,clusters,end) %>%
            mutate(start= ifelse(n()>1 & endR==distance_transcriptome-1,min(start),start)) %>% ungroup()
        }
        current.2 = group_by(current,start,end) %>%
          mutate(collapsed=NA, collapsed = paste(transcript_name, collapse="_"), nCollapsed = n_distinct(transcript_name)) %>%
          dplyr::filter(nCollapsed!=1) %>% ungroup()
        #breaking rule
        if(nrow(current.2)==0){break}
        if(ind == max(curr$exon_number)){
          collapsed = rbind(collapsed, distinct(current.2, collapsed))
        }
      }
    }
  }
  return(stats::na.omit(collapsed))
}



<<<<<<< HEAD
# File opening
#gtf..
gtf = fread(GTF_PATH, col.names=c("seqnames","start","end","width","strand","gene_id","gene_name","transcript_id","transcript_name"))

#. Read bulk quantification...
qf = fread(QF_PATH, col.names=c("gene_name","transcript_name", "bulk_TPMperc")) %>% dplyr::filter(transcript_name %in% gtf$transcript_name)

if(nrow(qf)==0){
  
  warning(paste("Any bulk quantification found for the file ", GTF_PATH))

}else{

  #. Formatting table
  gtf = gtf %>%
    arrange(seqnames,start,desc(end)) %>% 
    group_by(transcript_name) %>% 
    dplyr::mutate(exon_number = ifelse(strand=="+",n():1,1:n())) %>% 
    arrange(transcript_name) %>% 
    data.table()

  #. Calculate relative coordinates
  gtf = gtf %>% 
    arrange(transcript_name,exon_number) %>% 
    group_by(transcript_name) %>% 
=======
# 1) File opening
#----------------
gtf = fread(GTF_PATH)
qf = fread(QF_PATH, col.names=c("gene_name","transcript_name", "meanTPM", "bulk_TPMperc"))

#Formatting table (A)
gtf = gtf %>% 
    distinct(seqnames,start,end,width,strand,gene_name,transcript_name) %>%
    arrange(seqnames,start,desc(end)) %>%
    stats::na.omit() %>%
    group_by(transcript_name) %>%
    dplyr::mutate(exon_number = ifelse(strand=="+",n():1,1:n())) %>%
    arrange(transcript_name) %>%
    data.table()

#Calculate relative coordinates (B)
gtf = gtf %>%
    arrange(transcript_name,exon_number) %>%
    group_by(transcript_name) %>%
>>>>>>> main
    group_modify(~{
      starts = .x$start
      ends = .x$end
      startR = rep(NA,length(starts))
      endR = rep(NA,length(ends))
      for(i in .x$exon_number){
      if(i==1){
          startR[i] = 0
          endR[i] = (ends[i] - (starts[i])) - 1
      }else{
          startR[i] = endR[i-1] + 1
          endR[i] = startR[i] + (ends[i] - starts[i])
      }
      .x$startR = startR
      .x$endR = endR
      }
      .x
    }) %>%
<<<<<<< HEAD
  distinct(seqnames,start,end,strand,startR,endR,gene_name,transcript_name,exon_number) %>%
  #Truncate isoforms within transcriptomic distance threshold
  dplyr::filter(startR < DT_THRESHOLD) %>%
  mutate(start = ifelse(strand == "+" & endR>DT_THRESHOLD,end-(DT_THRESHOLD-startR),start),
          end = ifelse(strand== "-" & endR>DT_THRESHOLD,start+(DT_THRESHOLD-startR),end),
          endR = ifelse(endR>DT_THRESHOLD,DT_THRESHOLD-1,endR)) %>%
  arrange(seqnames,start,desc(end),transcript_name) %>%
  data.table()

  #. Identification of similar isoforms
  collapseds.tab = gtf %>%
    group_by(gene_name) %>%  
    group_modify(~{collapseIsoforms(.x, DT_THRESHOLD, DT_EX_THRESHOLD)}) %>% 
    ungroup()


  #. Format collapsed table
  gtf = left_join(gtf, collapseds.tab) %>% 
    mutate(collapsed = ifelse(is.na(collapsed),"none",collapsed))

  #writing
  distinct(gtf, gene_name, collapsed) %>% 
    dplyr::filter(collapsed!="none") %>% 
    dplyr::distinct() %>% 
    data.table::fwrite(file = paste0(gtf$seqnames[1],"_collapsed_isoforms.txt"), sep="\t", col.names = F, row.names = F)

  #. Collapsing similar transcripts
  gtf = left_join(gtf, qf) %>% 
    group_by(collapsed) %>%
    dplyr::mutate(check=ifelse(collapsed!="none",max(bulk_TPMperc),bulk_TPMperc), 
      end=ifelse(strand=="+" & collapsed!="none" & exon_number==1,max(end),end), 
      start=ifelse(strand=="-" & collapsed!="none" & exon_number==1,min(start),start)) %>% 
    dplyr::filter(bulk_TPMperc==check) %>% 
    group_by(collapsed) %>%
    dplyr::mutate(endR=ifelse(collapsed!="none" & exon_number==max(exon_number), (startR + (end-start)) - 1, endR)) %>% 
    ungroup() %>% 
    distinct(seqnames,start,end,startR,endR,strand,gene_name,transcript_name,exon_number,bulk_TPMperc,collapsed)

  #. Writing
  #all exons entries
  fwrite(gtf, file = args$OUTPUT, sep="\t", nThread=1)

  #. Writing all genes with unique isoforms
  gtf_unique = gtf %>% dplyr::select(gene_name,transcript_name) %>%
    distinct() %>% 
    group_by(gene_name) %>% 
=======
    distinct(seqnames,start,end,strand,startR,endR,gene_name,transcript_name,exon_number) %>%
    data.table()

#Truncate isoforms within transcriptomic distance threshold
gtf = gtf %>%
    dplyr::filter(startR < DT_THRESHOLD) %>%
    mutate(start = ifelse(strand == "+" & endR>DT_THRESHOLD,end-(DT_THRESHOLD-startR),start),
            end = ifelse(strand== "-" & endR>DT_THRESHOLD,start+(DT_THRESHOLD-startR),end),
            endR = ifelse(endR>DT_THRESHOLD,DT_THRESHOLD-1,endR)) %>%
    arrange(seqnames,start,desc(end),transcript_name)

#Collapse similar isoforms
collapseds.tab = gtf %>%
    group_by(gene_name) %>%
    group_modify(~{collapseIsoforms(.x, DT_THRESHOLD, DT_EX_THRESHOLD)}) %>%
    ungroup()

gtf = left_join(gtf, collapseds.tab) %>% 
    mutate(collapsed = ifelse(is.na(collapsed),"none",collapsed))

#writing
distinct(gtf, gene_name, transcript_name, collapsed) %>% 
    dplyr::filter(collapsed!="none") %>%
    distinct() %>%
fwrite(file = paste0(gtf$seqnames[1],"_collapsed_isoforms.txt"), sep="\t", col.names = F, row.names = F)

#collapsing
gtf = left_join(gtf, qf) %>% 
    group_by(collapsed) %>%
    dplyr::mutate(check=ifelse(collapsed!="none",max(bulk_TPMperc),bulk_TPMperc),
            end=ifelse(strand=="+" & collapsed!="none" & exon_number==1,max(end),end),
            start=ifelse(strand=="-" & collapsed!="none" & exon_number==1,min(start),start)) %>%
    dplyr::filter(bulk_TPMperc==check) %>%
    group_by(collapsed) %>%
    dplyr::mutate(endR=ifelse(collapsed!="none" & exon_number==max(exon_number), (startR + (end-start)) - 1, endR)) %>%
    ungroup() %>%
    distinct(seqnames,start,end,startR,endR,strand,gene_name,transcript_name,exon_number,bulk_TPMperc,collapsed)

#6) Writing
#all exons entries
fwrite(gtf, file = args$OUTPUT, sep="\t", nThread=1)

#all genes with unique isoforms
gtf_unique = gtf %>% dplyr::select(gene_name,transcript_name) %>% 
    distinct() %>%
    group_by(gene_name) %>%
>>>>>>> main
    mutate(counts = n()) %>%
    dplyr::filter(counts == 1) %>% 
    ungroup()
<<<<<<< HEAD
  
  gtf_unique = gtf_unique[,c("gene_name","counts")]
  data.table::fwrite(gtf_unique, file=args$OUTPUT_UNIQUE, sep="\t")

}
=======
gtf_unique = gtf_unique[,c("gene_name","counts")]
fwrite(gtf_unique, file=args$OUTPUT_UNIQUE, sep="\t")
>>>>>>> main
