readID=$1
threadN=$2
readDir=/archive/kimzz14/SRA_RAW/INSDC/Nibea_coibor/PRJNA827677

bwa \
    mem \
    -5SP \
    -T0 \
    -a \
    -t ${threadN} \
    bwadb/ref.fa \
    ${readDir}/${readID}_1.fastq.gz \
    ${readDir}/${readID}_2.fastq.gz \
    2> result/${readID}.sam.log \
    | gzip > result/${readID}.sam.gz
