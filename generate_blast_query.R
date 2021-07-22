library(magrittr)

FASTA_QUERYS_FOLDER <- "./blast_fasta_querys/"

generate_blast_query <- function(species, bold_data, filter_threshold) {
    search_querys <- list()

    fs::dir_create(FASTA_QUERYS_FOLDER)

    lapply(list.files(FASTA_QUERYS_FOLDER), function(f) {
        file.remove(paste0(FASTA_QUERYS_FOLDER, f))
    })

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

        for (bin in unique(this_species$bin_uri)) {
            this_bin <- this_species[this_species$bin_uri == bin, ]
            max_size <- max(this_bin$base_number)
            max_bins <- max_bins %>% rbind(this_bin[this_bin$base_number == max_size, ])
        }

        unique_bins_max <- max_bins[match(unique(max_bins$bin_uri), max_bins$bin_uri), ]

        fa_file <- unique_bins_max %>% apply(1, function(x) {
            x <- unlist(x)
            paste0(">", spe, "\n", x["nucleotides"])
        }) %>% unlist

        fileConn <- file(paste0("./blast_fasta_querys/", spe, ".fasta"))
        writeLines(fa_file, fileConn)
        close(fileConn)

        search_querys[[spe]] <- unique_bins_max
    }

    search_querys
}
