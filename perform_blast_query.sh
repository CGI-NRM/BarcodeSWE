cd mito

headers="qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus"

if [[ -e ../blast_result.out ]]; then
    rm ../blast_result.out
fi

echo $headers | tr " " \\t > ../blast_result.out

for file in ../blast_fasta_querys/*.fasta
do
    echo $file
    species=$(printf '%s' "$file" | sed -e "s/^.*\///" -e "s/\.fasta$//")
    echo $species
    temp_file=$(mktemp)

    # -entrez_query "$species[organism]" \
    # -taxids number \ # This needs be found
    blastn -db mito \
        -query "$file" \
        -out "$temp_file" \
        -outfmt "6 $headers" \
        -num_threads 12 \
        -max_target_seqs 500 \


    cat $temp_file >> ../blast_result.out
    # break
done

cd ..

