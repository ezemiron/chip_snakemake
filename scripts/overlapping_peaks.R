#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(rtracklayer)))
args = commandArgs(trailingOnly=T)
#
# args = c('hg38/peaks/K562/Schmidl2015/H3K27AC_chipmentation_r1_SRR2085875_peaks.narrowPeak',
#          'hg38/peaks/K562/Schmidl2015/H3K27AC_chipmentation_r2_SRR2085876_peaks.narrowPeak')
if (length(args) > 1){
    in_list = lapply(args, function(fn){
                         df = read.table(fn, stringsAsFactors=F)
                         if (ncol(df) > 3){
                             rownames(df) = df[,4]
                         } else {
                             rownames(df) = range(1,nrow(df))
                         }
                         return(df)
                     })

    in_gr = lapply(in_list, function(df){
                       gr = GRanges(seqnames = df[,1], IRanges(df[,2], df[,3]))
                       names(gr) = rownames(df)
                       return(gr)
                   })

    shared_gr = in_gr[[1]]
    mcols(shared_gr) = DataFrame(matrix('', nrow=length(shared_gr),
                                        ncol=length(in_gr)))
    mcols(shared_gr)[,1] = names(in_gr[[1]])
    for (i in 2:length(in_gr)){
        o = findOverlaps(shared_gr, in_gr[[i]])
        shared_gr = shared_gr[queryHits(o)]
        other_gr = in_gr[[i]][subjectHits(o)]

        mcols(shared_gr)[,i] = names(other_gr)

        start(shared_gr) = unlist(apply(cbind(start(shared_gr), start(other_gr)), 1, max))
        end(shared_gr) = unlist(apply(cbind(end(shared_gr), end(other_gr)), 1, min))
    }
    if (length(shared_gr) == 0){
        df = data.frame()
    } else{
        df = data.frame(data.frame(shared_gr)[,1:3],
                        name=paste0("peak_intersect", 1:length(shared_gr)),
                        score="0", strand=".")
        for (i in 1:length(in_gr)){
            name_vec = mcols(shared_gr)[,i]
            n = ncol(in_list[[i]])
            if (n > 6){
                df = data.frame(df, in_list[[i]][name_vec, c(4,7:n)])
            } else{
                df = data.frame(df, in_list[[i]][name_vec, 4])
            }
        }

    }
} else {
    df = read.table(args[1], stringsAsFactors=F)
}

write.table(df, col.names=F, row.names=F, quote=F, sep='\t')
