library(bold)


grade <- function(data, species_column, bin_column) {
    bin_per_species_list <- list()
    species_per_bin_list <- list()

    data %>% apply(1, function(x) {
        bins <- bin_per_species_list[[x[[species_column]]]]
        if (is.null(bins)) {
            bin_per_species_list[[x[[species_column]]]] <- c()
        }

        bin_per_species_list[[x[[species_column]]]] <- c(bins, x[[bin_column]])

        species <- species_per_bin_list[[x[[bin_column]]]]
        if (is.null(species)) {
            species_per_bin_list[[x[[bin_column]]]] <- c()
        }

        species_per_bin_list[[x[[bin_column]]]] <- c(species, x[[species_column]])
    })

    species_per_bin <- lapply(species_per_bin_list, function(x) {
        x %>% unique %>% length
    })

    bin_per_species <- lapply(bin_per_species_list, function(x) {
        x %>% unique %>% length
    })
}

