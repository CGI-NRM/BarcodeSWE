library(DECIPHER)

train_decipher <- function(seqs_path, rank_path, maxGroupSize, maxIterations) {
    # TODO:
    #   Add back cat from git
    #   pipe output to shiny terminal/text
    #   Only keep message for exported 'trainingSet.Rdata'
    message <- c()
    seqs <- readDNAStringSet(seqs_path)
    taxid <- read.table(rank_path, header=FALSE, col.names=c('Index', 'Name', 'Parent', 'Level', 'Rank'), sep="*", quote="", stringsAsFactors=FALSE)

    seqs <- RemoveGaps(seqs)
    seqs <- OrientNucleotides(seqs)

    groups <- names(seqs)
    groups <- gsub("(.*)(Root;)", "\\2", groups)
    groupCounts <- table(groups)
    u_groups <- names(groupCounts) # unique groups
    length(u_groups) # number of groups

    remove <- logical(length(seqs))
    for (i in which(groupCounts > maxGroupSize)) {
        index <- which(groups == u_groups[i])
        keep = sample(length(index), maxGroupSize)
        remove[index[-keep]] <- TRUE
    }

    message <- c(message, paste0("Number of eliminated species: ", sum(remove)))

    allowGroupRemoval <- FALSE
    probSeqsPrev <- integer() # suspected problem sequences from prior iteration
    for (i in seq_len(maxIterations)) {
        # train the classifier
        trainingSet <- LearnTaxa(seqs[!remove], names(seqs)[!remove], taxid)
        # look for problem sequences
        probSeqs <- trainingSet$problemSequences$Index
        if (length(probSeqs) == 0) {
            cat("No problem sequences remaining.\n")
            break
        } else if (length(probSeqs) == length(probSeqsPrev) && all(probSeqsPrev == probSeqs)) {
            cat("Iterations converged.\n")
            break
        }
        if (i == maxIterations) {
            break
        }
        probSeqsPrev <- probSeqs
        # remove any problem sequences
        index <- which(!remove)[probSeqs]
        remove[index] <- TRUE # remove all problem sequences
        if (!allowGroupRemoval) {
            # replace any removed groups
            missing <- !(u_groups %in% groups[!remove])
            missing <- u_groups[missing]
            if (length(missing) > 0) {
                index <- index[groups[index] %in% missing]
                remove[index] <- FALSE # don't remove
            }
        }
    }

    message <- c(message, c(
        paste0("Total number of sequences eliminated: ", sum(remove)),
        paste0("Number of remaining problem sequences: ", length(probSeqs)),
        "Exported training data to 'trainingSet.Rdata'"))

    save(list = "trainingSet", file = "trainingSet.Rdata")

    list(trainingSet = trainingSet, message = message)
}
