library(magrittr)

bold <- read.table(file = "bold.tsv", sep = "\t", comment.char = "", header = TRUE)

ranks <- bold %>% apply(1, function(row) {
    paste0(
        "p__",
        row["phylum_name"],
        ";c__",
        row["class_name"],
        ";o__",
        row["order_name"],
        ";f__",
        row["family_name"],
        ";g__",
        row["genus_name"],
        ";s__",
        row["species_name"])
}) %>% unlist %>% unique

taxa <- setNames(c("phylum", "class", "order", "family", "genus", "species"),
                 c("p__", "c__", "o__", "f__", "g__", "s__"))
ranks <- strsplit(ranks, ";", fix=T)
count <- 1L
groups <- "Root"
index <- -1L
level <- 0L
rank <- "rootrank"
pBar <- txtProgressBar(style=3)
for (i in seq_along(ranks)) {
    for (j in seq_along(ranks[[i]])) {
        rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
        group <- substring(ranks[[i]][j], 4)
        w <- which(groups==group & rank==rank_level)
        if (length(w) > 0) {
            parent <- match(substring(ranks[[i]][j - 1], 4), groups)
            if (j==1 || any((parent - 1L)==index[w])) {
                next # already included
            }
        }
        count <- count + 1L
        groups <- c(groups, group)
        if (j==1) {
            index <- c(index, 0)
        } else {
            parent <- match(substring(ranks[[i]][j - 1], 4), groups)
            index <- c(index, parent - 1L)
        }
        level <- c(level, j)
        rank <- c(rank, taxa[j])
    }
    setTxtProgressBar(pBar, i/length(ranks))
}

groups <- gsub("^[ ]+", "", groups)
groups <- gsub("[ ]+$", "", groups)
taxid <- paste(0:(length(index) - 1L), groups, index, level, rank, sep="*")

writeLines(taxid, con = "taxid.txt")
