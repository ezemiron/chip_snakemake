#!/usr/bin/env Rscript
library(argparse)
suppressMessages(suppressWarnings(require(rtracklayer)))
suppressMessages(suppressWarnings(require(data.table)))



parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('-i', help='input file')
parser$add_argument('-cs', help='chrom sizes')
parser$add_argument('-bw', help='output bigwig')
parser$add_argument('-t', type='integer', help='number of cores')

argv <- parser$parse_args()

save(argv, file='test.Rdata')

setDTthreads(argv$t)

dt = fread(cmd=paste("sed 's/#//' ", argv$i, "| sed \"s/'//g\""))
dt = dt[rowSums(is.na(dt[,4:ncol(dt)]))==0,]

if ('Input' %in% colnames(dt)){
    rep_vec = grep('Input', colnames(dt)[4:ncol(dt)], invert=T, value=T)


    rep_total = colSums(dt[, ..rep_vec])
    input_total = sum(dt[,'Input'])

    scale_vec = unlist(lapply(rep_total, function(rep){min(rep, input_total)}))

    fcoc_list = lapply(rep_vec, function(rep){
        exp = dt[,..rep] / rep_total[rep] * scale_vec[rep] + 1
        inp = dt[,'Input'] / input_total * scale_vec[rep] + 1
        return(log2(exp / inp))
    })

    fcoc_dt = do.call(cbind, fcoc_list)
    dt[, mean:=rowMeans(fcoc_dt)]
    dt = dt[!is.na(mean) & chr!='chrM', ]
} else {
    dt[,mean:=rowMeans(dt[,4:ncol(dt)])]
}


gr = GRanges(dt[,chr], IRanges(dt[,start], dt[,end-1]), score=dt[,mean])

cs = read.table(argv$cs, row.names=1, col.names=c('seqnames', 'length'),
                stringsAsFactors=F)

seqlengths(gr) = cs[names(seqlengths(gr)), ]

export.bw(trim(gr), argv$bw)
