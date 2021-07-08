library(magrittr)

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

