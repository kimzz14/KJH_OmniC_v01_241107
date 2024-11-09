prefix=$1
threadN=$2

sh pipe/bwa-mem.sh ${prefix}_1 ${threadN}
sh pipe/bwa-mem.sh ${prefix}_2 ${threadN}
sh pipe/filter_sam.sh ${prefix} ${threadN}

