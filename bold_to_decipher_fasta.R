library(magrittr)

generate_decipher_fasta <- function(bold) {
    fasta <- apply(bold, 1, function(row) {
        paste0(
            "> ",
            gsub("^[ ]+", "", row["sequenceID"]),
            " Root; ",
            row["phylum_name"], "; ",
            row["class_name"], "; ",
            row["order_name"], "; ",
            row["family_name"], "; ",
            # row["subfamily_name"], "; ",
            row["genus_name"], "; ",
            row["species_name"],
            "\n",
            row["nucleotides"])
    }) %>% unlist

    writeLines(fasta, con = "bold.fasta")
}

