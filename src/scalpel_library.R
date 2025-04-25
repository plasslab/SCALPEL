


# [Script for Function needed in Scalpel analysis]
# ------------------------------------------------
require(GenomicRanges)
require(Gviz)
require(ggplot2)
require(stringr)

<<<<<<< HEAD

FindIsoforms <- function (seurat_obj, group_by, split_by = NULL, assay = "RNA",
                           threshold_pval = 0.05, threshold_abund = 0.1, threshold_var=0.05,
                           extended=F) {
    
    #check presence of group.by into seurat object
    if (!(group_by %in% colnames(seurat_obj[[]]))) {
        
        stop(paste0("Error: ", group.by, "not found in Seurat object !"))
        
=======
Find_isoforms = function(seurat.obj, pval_adjusted=0.05, condition="orig.ident", assay="RNA",
                         threshold_tr_abundance = 0.10){
  # ---------------------------------------------------------------------------
  #Function to Find differentially expressed isoforms in the conditions defined
  # ---------------------------------------------------------------------------

  message("processing...")

  #1/ Find genes with at least 2 transcripts
  message("Find genes with at least 2 transcripts...")
  my.genes = data.frame(gene_tr = rownames(seurat.obj)) %>%
    tidyr::separate(col=gene_tr, into=c("gene","trs"), sep = "\\*\\*\\*", remove = F) %>%
    group_by(gene) %>%
    reframe(nb.trs = n_distinct(trs), gene_tr) %>%
    dplyr::filter(nb.trs>1)

  #2/ Get Isoform expression in defined conditions
  message("Get isoforms expression in the condition defined...")
  ALL_expression = AggregateExpression(object = subset(seurat.obj, features = my.genes$gene_tr),
                                       assays = assay, group.by = condition, verbose = T,
                                       return.seurat = F)[[assay]] %>% data.frame() %>%
    dplyr::mutate(gene_tr = rownames(.)) %>%
    tidyr::separate(col="gene_tr", into=c("gene_name","transcript_name"), remove = F, sep="\\*\\*\\*") %>%
    data.table()
  colnames(ALL_expression) = str_replace(colnames(ALL_expression),"X","")

  #filter trs based on abundance threshold & X2 testing...
  message("filter trs based on abundance threshold and X2 testing...")
  conds = unlist(unique(seurat.obj[[condition]]))
  tmp_all = pbapply::pblapply(split(ALL_expression, ALL_expression$gene_name), function(.x){
    if(nrow(.x)==1){return(NULL)}
    #get columns condition and perform normalization...
    tmp = as_tibble(apply(as_tibble(.x)[,conds], 2, function(x) round(x/sum(x),2)))
    #discard isoforms under the threshold in all the condition
    tmp = .x[which((rowSums(tmp < threshold_tr_abundance) == length(conds)) == FALSE),]
    #perform X2 test for satisfying genes if nb.isoforms > 1
    if(nrow(tmp)>1){
      #perform X2 test
      tmp.test = chisq.test(as_tibble(tmp)[,conds], correct = T)
      tmp$p_value = tmp.test$p.value
      return(tmp)
>>>>>>> main
    } else {
        #process seurat object if split.by is NULL
        if ( is.null(split_by) ) {
            message("process seurat object...")
            
            #get tab.counts (1)
            message("Aggregating counts in conditions...")
            counts.tab = Extract_aggTab(seurat_obj, group = group_by, assay)
            
            #process filtering ops by gene (2)
            message("Performing filtering and comparison test...")
            all_results = pbapply::pblapply(split(counts.tab, counts.tab$gene_name), function(geneTab) {
                #apply filtering checks & X2 test...
                if(nrow(geneTab)>1){
                    return(Perform_xTest(geneTab, threshold_abund, threshold_var))
                } else {
                    tmp = geneTab %>% dplyr::mutate(thr.abund=T, thr.var=T, p.value=NA, statistic=NA)
                    return(tmp)
                }
            }) %>% rbindlist()
            
            #calculate adjusted P.values (3)
            all_results = all_results %>%
                dplyr::mutate(p.value.adjusted = ifelse(thr.abund==F & thr.var==F & (!(is.na(p.value))), 0, NA))
            
            tmpA = dplyr::filter(all_results, p.value.adjusted==0)
            tmpA$p.value.adjusted = p.adjust(tmpA$p.value, method = "BH")
            tmpB = dplyr::filter(all_results, is.na(p.value.adjusted))
            all_results = rbind(tmpA, tmpB) %>% dplyr::arrange(p.value.adjusted, gene_name, transcript_name) %>%
                dplyr::distinct() %>%
                data.table::data.table()
            
            #apply p.value threshold (4)
            if(extended==F){
                all_results = dplyr::filter(all_results, p.value.adjusted<threshold_pval) %>%
                    select(!c("thr.abund", "thr.var"))
            } else {
                return(all_results)
            }
            
            #return
            return(all_results)
            
        } else {
            #check presence of split.by var into seurat object
            if (!(split_by %in% colnames(seurat_obj[[]]))) {
                
                stop(paste0("Error: ", group_by, "not found in Seurat object !"))
                
            } else {
                message("split seurat object and process with each object...")
                
                #splitting
                seurat.objs = Seurat::SplitObject(JoinLayers(seurat_obj), split.by = split_by)
                
                #processing
                message("Generation of aggregate counts...")
                lapply(names(seurat.objs), function(curr_obj){
                    message(curr_obj)
                    
                    #get tab.counts
                    counts.tab = Extract_aggTab(seurat.objs[[curr_obj]], group = group_by, assay, n=1) %>% 
                        dplyr::mutate(split.by = curr_obj)
                    
                })
            }
        }
    }
}

# Functions associateds..
#========================

Extract_aggTab <- function (seurat_obj, group, assay) {
    #Get aggregated counts tab
    tmp = Seurat::AggregateExpression(seurat_obj, group.by = group, assays = assay)[[1]]
    names.cols = colnames(tmp)
    tmp = tmp %>%
        tibble::as_tibble(rownames="gene.tr") %>%
        tidyr::separate(gene.tr, c("gene_name","transcript_name"), sep="\\*\\*\\*", remove = F) %>%
        dplyr::arrange(gene_name, transcript_name)
    
    #get genes with nb isoforms > n
    tmp.stats = dplyr::distinct(tmp, gene_name, transcript_name) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarise(nb.isoforms = n_distinct(transcript_name))
    
    tmp = data.table::data.table(dplyr::filter(tmp, gene_name %in% tmp.stats$gene_name))
    colnames(tmp) = c("gene.tr", "gene_name", "transcript_name", names.cols)
    return(tmp)
}

Perform_xTest <- function(tab, threshold.abund, threshold.var) {
    #process filtering ops
    tab = tab %>%
        rowwise() %>%
        mutate(
            row_sum = sum(c_across(where(is.numeric)))
        ) %>%
        ungroup() %>%
        mutate(
            thr.abund = if_all(where(is.numeric) & !matches("row_sum"), ~ (. / sum(.)) < threshold.abund),
            thr.var = all(across(where(is.numeric) & !matches("row_sum"), ~ sd(.) < 0.05))
        ) %>%
        select(-row_sum)
    
    #get numeric table
    numTab = tab %>% dplyr::filter(thr.abund==F & thr.var==F) %>%
        dplyr::select(where(is.numeric))
    
    #perform X2 test
    if(nrow(numTab)>1){
        #perform test if more than 1 isoforms following filtering ops
        chi2_res <- chisq.test(as.matrix(numTab))
        
        #check if X2 test condition are reunited (counts>=5)...
        if (all(round(chi2_res$expected,0) >= 5)) {
            #update tab
            tab = tab %>% dplyr::mutate(
                p.value = chi2_res$p.value,
                statistic = chi2_res$statistic
            )
        }else{
            #update tab
            tab = tab %>% dplyr::mutate(
                p.value = NA,
                statistic = NA
            )
        }
    } else {
        #update tab
        tab = tab %>% dplyr::mutate(
            p.value = NA,
            statistic = NA
        )
    }
    #return
    return(tab)
}




plot_relativeExp = function(seurat.obj, features_in, group.var, levels.group=NULL, assay="RNA", ...){
  #Function to plot with geom_boxplot relative isoform expression

  #extract transcripts
  A = seurat.obj[[assay]]$counts[features_in,] %>%
    data.frame() %>%
    mutate(gene_transcript=features_in) %>%
    melt() %>%
    suppressWarnings()

  seurat.obj$variable = rownames(seurat.obj@meta.data)
  A$group = dplyr::left_join(A, seurat.obj@meta.data) %>% dplyr::select(group.var) %>% as.vector() %>% unlist()
  A = arrange(A, variable)

  #Calculates relative expression
  A.tab = filter(A, value>0) %>%
    group_by(gene_transcript) %>%
    mutate(perc = value / sum(value) * 100) %>%
    arrange(gene_transcript,group) %>%
    data.table()

  #Visualization
  ggplot(A.tab, aes(group, perc, fill=gene_transcript)) +
    geom_boxplot2(width = 0.5, width.errorbar = 0.1, color="gray20") +
    theme_few(base_size = 12) +
    scale_fill_calc() +
    ggtitle("Isoform relative expression in cell types")
}



CoverPlot = function(genome_gr, gene_in, genome_sp, bamfiles, distZOOM=NULL, annot_tab=NULL,
                     transcripts_in=NULL, filter_trs=F, samtools.bin="samtools"){
  #check args
  message(genome_sp)
  if(is.null(genome_sp) || is.null(genome_gr) || is.null(bamfiles)){
    stop("Error ! Check input args : genome.sp / genome.gr / bamfiles")
  }

  #0. Ideogram
  axisTrack <- Gviz::GenomeAxisTrack(genome=genome_sp)

  #1. Build genome table & Track
  message("GenomeTrack building...")
  gtab.gr = genome_gr[genome_gr$gene_name==gene_in]
  gtab.gr$colors = "orange"
  #filter transcripts provided in case
  if(filter_trs){
    message(paste0("filtering...", transcripts_in))
    if(is.null(transcripts_in)){stop("Error! Provide list of transcripts: transcripts_in")}
    gtab.gr = gtab.gr[gtab.gr$transcript_name %in% transcripts_in]
  }
  gtab.gr$colors = as.character(gtab.gr$colors)
  gtab.gr$transcript_id = paste(gtab.gr$transcript_id, gtab.gr$transcript_name,sep=" - ")
  gtab = data.frame(gtab.gr) %>%
    dplyr::filter(type %in% c("UTR","CDS","exon")) %>%
    dplyr::distinct(seqnames,start,end,width,strand,gene_type,gene_id,
                    exon_id,transcript_id,gene_name,transcript_name,type) %>%
    dplyr::rename(Chromosome="seqnames",feature="gene_type",
                  gene="gene_id",exon="exon_id",transcript="transcript_id") %>%
    group_by(transcript) %>%
    mutate(check = ifelse("UTR" %in% type, T,F)) %>%
    dplyr::filter(!(type == "exon" & check == T)) %>%
    ungroup()
  gtab$feature = gtab$type
  gtab$feature = stringr::str_replace(gtab$feature,"UTR","utr")
  gtab$feature = stringr::str_replace(gtab$feature, "exon", "utr")
  gtab$colors = "orange"
  if(!is.null(transcripts_in)){
    message("highliting transcripts...")
    gtab = arrange(gtab, transcript, start, desc(end))
  }
  #define sets
  chrom = gtab$Chromosome[1]
  starts = min(gtab$start)
  ends = max(gtab$end)
  strands = as.character(gtab$strand[1])

  #track
  dplyr::select(gtab, Chromosome,start,end) %>%
    mutate(start = start - 25, end = end + 25) %>%
  data.table::fwrite(file="./coords.txt", sep="\t", row.names = F, col.names = F)
  Gtrack = Gviz::GeneRegionTrack(GenomicRanges::makeGRangesFromDataFrame(gtab, keep.extra.columns = T), chromosome = chrom,
                                 name = gene_in, transcriptAnnotation = "transcript",
                                 just.group="below", genome=genome_sp, fill="orange", color="black", col = "black",
                                 background.title="darkmagenta", fontsize.group=16,
                                 fill=gtab$colors)

  #1bis : Annotation Track building
  if(!is.null(annot_tab)){
    curr.tab = dplyr::filter(annot_tab, start>=starts & end<=ends)
    Atrack = Gviz::AnnotationTrack(GenomicRanges::makeGRangesFromDataFrame(curr.tab), name = "Peaks",
                             chromosome = chrom, fill="olivedrab", id=curr.tab$name,
                             background.title="black", featureAnnotation = "id", fontcolor.group="black")
  }else{
    Atrack = Gviz::AnnotationTrack(chromosome = chrom, background.title="black", name = "Peaks")
  }

  #2. Build Alignment Track
  message("Building Alignment Track...")
  Atracks.res = lapply(1:length(bamfiles), function(x){
    message(names(bamfiles)[x])
    #get coverage
    system(paste0(samtools.bin, " view -b --region-file coords.txt ", bamfiles[x], " > current.bam"))
    cov.exp = system(paste0(samtools.bin," depth -b coords.txt current.bam > current.cov"))
    cov.tab = fread("current.cov", col.names = c("seqnames","start","depth")) %>% dplyr::filter(depth>=0)
    #dataTrack
    curr.track = Gviz::DataTrack(start = cov.tab$start, width=1, data = cov.tab$depth, chromosome = chrom, genome=genome_sp,
                           type=c("hist"), background.title="coral4", name = names(bamfiles)[x], col.histogram="blue")
    return(list(max.depth=max(cov.tab$depth), track=curr.track))
  })
  Atracks = lapply(Atracks.res, function(x) x[[2]])
  YMAX = max(unlist(lapply(Atracks.res, function(x) x[[1]])))
  system("rm current.cov")
  system("rm current.bam")

  #plotting...
  message("Plotting...")
  in.tracks = c(axisTrack, Atrack,Gtrack, Atracks)

  #in case of ZOOMING option
  size.tracks = c(0.3,0.1,0.6,rep(0.9/length(bamfiles), length(bamfiles)))
  if(is.null(distZOOM)){
    plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX))
  }else{
    if(strands == "+"){
      size.tracks = c(0.1,0.4,0.6,rep(0.5/length(bamfiles), length(bamfiles)))
      plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX), from=ends-distZOOM, just.group="right")}
    if(strands == "-"){
      size.tracks = c(0.1,0.4,0.6,rep(0.5/length(bamfiles), length(bamfiles)))
      plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX), to=starts+distZOOM, just.group="left")}
  }
}
