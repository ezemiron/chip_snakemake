#!/bin/bash
usage=("$(basename "$0") [-h] -b <...> -w <INT> -c <...> -o <...> -t <INT>

script to run sambamba window coverage and pipe output to a bigWig file

where:
    -h  show this help text
    -b  bam file
    -w  window size
    -c  chromosome size file
    -o  output file
    -t  threads
    ")

threads=1

while getopts ':h:b:w:c:o:?t:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    b) bam_file="$OPTARG"
       ;;
    w) window="$OPTARG"
       ;;
    o) out="$OPTARG"
       ;;
    c) chromsizes="$OPTARG"
       ;;
    t) threads="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done
shift $((OPTIND - 1))

tmp=$(tempfile)


sambamba depth window -w $window \
                      -t $threads \
                      $bam_file | \
     awk -vOFS='\t' '{if ($1 !~ /^#|^Process|^chrM/){print $1, $2, $3, $4}}' | \
     sort -k1,1 -k2,2n > $tmp

awk '{
    if (NR==FNR){
        sizes[$1] = $2
    } else if ($3 < sizes[$1]){
        print $0
    }}' $chromsizes $tmp > $tmp"2"


bedGraphToBigWig $tmp"2" $chromsizes $out

rm $tmp
rm $tmp"2"
