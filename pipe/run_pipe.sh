prefix=$1
threadN=$2

sh pipe/bwa-mem.sh    ${prefix} ${threadN}
sh pipe/filter_sam.sh ${prefix} ${threadN}

