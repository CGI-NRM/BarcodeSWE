done <- FALSE

print("Downloading nt database")
while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "nt", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}

print("Downloading taxdb database")
while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "taxdb", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}

print("Downloading taxdb metadata")
while (!done) {
    done <- TRUE
    tryCatch({
        biomartr::download.database.all(db = "taxdb-metadata", path = "./ncbi_db")
    },
    error = function(err) {
        done <<- FALSE
    })
}
