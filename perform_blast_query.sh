headers="qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus"


   	    # headers="scomnames scomname sscinames ssciname length stitle sstart send qstart \
# qend mismatch qgi qseqid gaps staxid sblastname mismatch evalue pident score bitscore"

temp_file=$(mktemp)

blastn -db mito \
-query blast_query.fa \
-out $temp_file \
-num_threads 12 \
-outfmt "6 $headers" \
-max_target_seqs 500

if [[ -e blast_result.out ]]; then
    rm blast_result.out
fi

echo $headers | tr " " \\t > blast_result.out
cat $temp_file >> blast_result.out
