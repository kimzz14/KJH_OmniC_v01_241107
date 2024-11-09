readID=$1
threadN=$2
readDir=/archive/kimzz14/SRA_RAW/INSDC/Nibea_coibor/PRJNA827677

bwa \
    mem \
    -5 \
    -T0 \
    -a \
    -t ${threadN} \
    bwadb/ref.fa \
    ${readDir}/${readID}.fastq.gz \
    2> result/${readID}.sam.log \
    | gzip > result/${readID}.sam.gz
