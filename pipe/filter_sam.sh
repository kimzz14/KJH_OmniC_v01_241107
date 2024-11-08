prefix=$1
threadN=$2

python ./script/filter_batch.py \
    -t ${threadN} \
    -1 result/${prefix}_1.sam.gz \
    -2 result/${prefix}_2.sam.gz \
    -o result/${prefix}.uniq.bam \
    1> result/${prefix}.uniq.bam.log \
    2> result/${prefix}.uniq.bam.err