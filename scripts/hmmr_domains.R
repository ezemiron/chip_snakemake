#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(rtracklayer)))
suppressMessages(suppressWarnings(require(HMMt)))
library(argparse)

parser <- ArgumentParser(description='Process normalized bw track into HMM using HMMt')

parser$add_argument('--bw', help='bigwig track with input')
parser$add_argument('--gff', help='gff output file', nargs='?')



argv <- parser$parse_args()

save(argv, file='test.Rdata')

gr = import.bw(argv$bw)
df = data.frame(gr)[, c('seqnames', 'start', 'end', 'score')]

br = bridge(df, is.sorted=T)

sink(file="/dev/null")
fit = BaumWelchT(br$x, br$series.length)
sink()

boundstate <- which.max(fit$mu)
model <- 0 + (fit$ViterbiPath[br$nonvirtuals] == boundstate)
model <- ifelse(model == 1, "Domain", "iDomain")

gr$domains = model

result <- unlist(reduce(split(gr, ~domains)))
result$domain = names(result)

if (!is.null(argv$gff)){
	export(result, argv$gff, format='GFF3')
}

write.table(data.frame(result[order(result)])[,c('seqnames', 'start', 'end', 'domain')],
            row.names=F, col.names=F, sep='\t', quote=F)
