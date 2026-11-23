# env: env/pi/SCALPEL
library(data.table)
library(argparse)
library(Matrix)
parser <- ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('DGE')
parser$add_argument('PRED')
parser$add_argument('OUTP')
args <- parser$parse_args()
preds <- if (file.size(args$PRED) == 0) {
  data.table(bc=character(), gene_name=character(), transcript_name=character(), tr_prob=numeric())
} else {
  fread(args$PRED, header=FALSE, col.names=c('bc','gene_name','transcript_name','tr_prob'), colClasses=c('character','character','character','numeric'))
}
dge <- fread(args$DGE, check.names=FALSE)
features <- if (nrow(preds)) paste0(preds$gene_name, '***', preds$transcript_name) else character()
pairs <- unique(data.table(feature=features, gene=preds$gene_name, transcript=preds$transcript_name))
if (anyDuplicated(pairs$feature)) stop('Feature name collision: ', pairs$feature[anyDuplicated(pairs$feature)])
feature_names <- unique(features)
barcodes <- unique(preds$bc)
ri <- match(preds$gene_name, dge[[1L]])
ci <- match(preds$bc, names(dge)[-1L])
counts <- numeric(nrow(preds))
valid <- which(!is.na(ri) & !is.na(ci))
for (idx in split(valid, ci[valid])) {
  col <- ci[idx[1L]]
  counts[idx] <- dge[[col + 1L]][ri[idx]]
}
counts[is.na(counts)] <- 0
entries <- data.table(i=match(features, feature_names), j=match(preds$bc, barcodes), x=preds$tr_prob * counts)
conflicts <- entries[, .(values=uniqueN(x)), by=.(i,j)][values > 1L]
if (nrow(conflicts)) stop('Conflicting duplicate feature/barcode: ', feature_names[conflicts$i[1L]], ' / ', barcodes[conflicts$j[1L]])
entries <- unique(entries, by=c('i','j'))
mat <- sparseMatrix(i=entries$i, j=entries$j, x=entries$x, dims=c(length(feature_names), length(barcodes)), dimnames=list(feature_names, barcodes))
mat <- drop0(mat)
print(c(features=nrow(mat), barcodes=ncol(mat), nonzero=length(mat@x)))
writeMM(mat, args$OUTP)
stem <- sub('\\.txt$', '', args$OUTP)
writeLines(feature_names, paste0(stem, '_features.txt'))
writeLines(barcodes, paste0(stem, '_barcodes.txt'))
