done <- FALSE

while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "nt", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}

while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "taxdb", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}

while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "taxdb-metadata", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}
