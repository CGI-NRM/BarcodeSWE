
headers="scomnames sscinames evalue length \
sstart send qstart qend mismatch qgi qseqid"

temp_file=$(mktemp)

blastn -db mito \
-query blast_query.fa \
-out $temp_file \
-num_threads 12 \
-outfmt "10 $headers" \
-max_target_seqs 500

echo $headers | tr " " , > blast_result.out
cat $temp_file >> blast_result.out
