library(magrittr)

bold = read.table(file = "bold.tsv", sep = "\t", comment.char = "", header = TRUE, na.strings = c("NA"))

total_rows <- nrow(bold)
bold <- bold[!is.na(bold[,"nucleotides"]) & bold[,"nucleotides"] != "",]

print(paste0("Removed ", total_rows - nrow(bold), " rows since they did not have any sequency on them."))

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


