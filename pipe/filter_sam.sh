prefix=$1
threadN=$2

python ./script/filter_batch.py \
    -t ${threadN} \
    -i result/${prefix}.sam.gz \
    -o result/${prefix}.uniq.bam \
    1> result/${prefix}.uniq.bam.log \
    2> result/${prefix}.uniq.bam.err