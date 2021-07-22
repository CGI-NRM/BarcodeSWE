library(magrittr)

generate_blast_query <- function(species, bold_data, filter_threshold) {
    search_querys <- list()

    fs::dir_create("./blast_fasta_querys/")

    for (spe in species) {
        count <- (bold_data$species_name == spe) %>% sum
        if (count == 0) {
            this_species <- bold_data[bold_data$genus_name == strsplit(spe, " ")[[1]][1], ]
        } else if (count < filter_threshold) {
            this_species <- bold_data[bold_data$species_name == spe, ]
        } else {
            next()
        }

        max_bins <- data.frame()

        for (bin in unique(this_species$BIN)) {
            this_bin <- this_species[this_species$BIN == bin, ]
            max_size <- max(this_bin$base_number)
            max_bins <- max_bins %>% rbind(this_bin[this_bin$base_number == max_size, ])
        }

        unique_bins_max <- max_bins[match(unique(max_bins$BIN), max_bins$BIN), ]

        fa_file <- unique_bins_max %>% apply(1, function(x) {
            x <- unlist(x)
            paste0(">", spe, "\n", x["nucleotides"])
        })


        search_querys[[spe]] <- unique_bins_max



    }
}

BAGS_download <- read.table(file = "Library_checklist.tsv", sep = "\t", header = TRUE, quote = '"', comment.char = '"', na.strings = "NA")

BAGS_download_filtered <- data.frame()

for (name in unique(BAGS_download[, "species"])) {
    match_vector <- BAGS_download[,"species"] == name

    if (any(BAGS_download[match_vector, "grade"] == "A")) {
        next
    }

    BAGS_this_name <- BAGS_download[match_vector,]
    BAGS_this_name_max <- data.frame()

    for (bin in unique(BAGS_this_name$BIN)) {
        BAGS_this_bin <- BAGS_this_name[BAGS_this_name$BIN == bin,]
        BAGS_this_name_max <- BAGS_this_name_max %>% rbind(BAGS_this_bin[BAGS_this_bin$base_number == max(BAGS_this_bin$base_number),])
    }

    unique_bins_max <- BAGS_this_name_max[match(unique(BAGS_this_name_max$BIN), BAGS_this_name_max$BIN),]

    BAGS_download_filtered <- BAGS_download_filtered %>% rbind(unique_bins_max)
}

fa_file <- BAGS_download_filtered %>% apply(1, function(x) {
    x <- unlist(x)
    paste0(">", x["species"], "\n", x["sequence"])
}) %>% unlist

fileConn <- file("blast_query.fa")
writeLines(fa_file, fileConn)
close(fileConn)

