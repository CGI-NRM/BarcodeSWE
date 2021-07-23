library(magrittr)
library(DECIPHER)

source("download_bold_database.R")
source("generate_blast_query.R")
source("bold_to_decipher_fasta.R")
source("bold_to_decipher_taxid.R")
source("train_decipher.R")

options(shiny.maxRequestSize = 100 * 1024^2)

europe <- c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan", 
"Belarus", "Belgium", "Bosnia and Herzegovina", "Bosnia", "Herzegovina", 
"Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", 
"Finland", "France", "Georgia", "Germany", "Greece", "Hungary", "Iceland", 
"Ireland", "Italy", "Kazakhstan", "Kosovo", "Latvia", "Liechtenstein", 
"Lithuania", "Luxembourg", "Malta", "Moldova", "Monaco", "Montenegro", 
"Netherlands", "North Macedonia", "Macedonia", "Norway", "Poland", "Portugal", 
"Romania", "Russia", "San Marino", "Serbia", "Slovakia", "Slovenia", "Spain", 
"Sweden", "Switzerland", "Turkey", "Ukraine", "United Kingdom", "UK", "Vatican City")

ui <- shiny::fluidPage(
    shiny::titlePanel("Download BOLD and NCBI Data"),
    shiny::tags$hr(),

    shiny::tabsetPanel(id = "main_panel",
        shiny::tabPanel(title = "Load data", label = "tab_load_data",
            shiny::h3("Load BOLD Data"),
            shiny::sidebarLayout(
                sidebarPanel = shiny::sidebarPanel(width = 3,
                    shiny::selectInput(inputId = "select_download_online_or_load_local", label = "Data mode",
                        choices = c("Download Online BOLD" = "download_online_bold", "Load Local BOLD" = "load_local_bold")),
                    shiny::conditionalPanel(condition = 'input.select_download_online_or_load_local == "download_online_bold"',
                        shiny::fileInput(inputId = "file_species_list", label = "Species List", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv")),
                        shiny::checkboxInput(inputId = "checkbox_species_file_header", label = "Header", value = TRUE),
                        shiny::tags$hr(),
                        shiny::actionButton(inputId = "button_download_bold_data", label = "Download BOLD Data"),
                        shiny::tags$hr(),
                        DT::dataTableOutput(outputId = "DT_species_list")
                    ),
                    shiny::conditionalPanel(condition = 'input.select_download_online_or_load_local == "load_local_bold"',
                        shiny::fileInput(inputId = "file_bold_data", label = "BOLD Data", multiple = FALSE, accept = c(".tsv"))
                    )
                ),
                mainPanel = shiny::mainPanel(
                    DT::dataTableOutput(outputId = "DT_bold_data")
                )
            )
        ),
        shiny::tabPanel(title = "View data", label = "tab_view_data",
            shiny::uiOutput(outputId = "ui_species_viewer_holder"),
            shiny::numericInput(inputId = "numeric_keep_the_longest_n", label = "Keep only the longest N samples: ", value = 30, min = 1), 
            shiny::actionButton(inputId = "button_export_filtered_data", label = "Save Filtered Bold File"),
            DT::dataTableOutput(outputId = "DT_species_viewer")
        ),
        shiny::tabPanel(title = "Blast NCBI Data", label = "tab_blast_ncbi_data",
            shiny::actionButton(inputId = "button_generate_blast_query", label = "Generate Blast Query"),
            shiny::actionButton(inputId = "button_perform_blastn", label = "Perform BLASTN Search"),
            shiny::actionButton(inputId = "button_merge_blast_data", label = "Merge Blast Data"),
            shiny::div(id = "div_show_fa_querys"),
            DT::dataTableOutput(outputId = "DT_blast_results")
        ),
        shiny::tabPanel(title = "Merge Own Data", label = "tab_merge_own_data",
            shiny::fileInput(inputId = "file_own_fasta", label = "Fasta File", multiple = FALSE, accept = c(".fa", ".fasta", ".txt")),
        ),
        shiny::tabPanel(title = "Train DECIPHER", label = "tab_train_decipher",
            shiny::actionButton(inputId = "button_train_decipher", label = "Train DECIPHER"),
            shiny::plotOutput(outputId = "plot_trainingSet", height = 1400, width = 800),
            shiny::h4("Problematic Sequences"),
            DT::dataTableOutput(outputId = "DT_problematic_sequences")
        ),
        shiny::tabPanel(title = "Search With New Sequence", label = "tab_search_with_new_sequence",
            shiny::fileInput(inputId = "file_trainingSet", label = "DECIPHER Training Set", multiple = FALSE, accept = c("*")),
            shiny::textOutput(outputId = "text_trainingSet_message"),
            shiny::tags$hr(),
            shiny::fileInput(inputId = "file_new_unknown_fasta", label = "Fasta File", multiple = FALSE, accept = c(".fa", ".fasta", ".txt")),
            shiny::selectInput(inputId = "select_decipher_strand", label = "Strand",
                choices = c("Bottom" = "bottom", "Both" = "both", "Top" = "top"), selected = "top"),
            shiny::numericInput(inputId = "numeric_decipher_threshold", label = "Threshold", value = 60, min = 0, max = 100),
            shiny::actionButton(inputId = "button_search_for_matches", label = "Search With New Sequence"),
            shiny::plotOutput(outputId = "plot_search_match_circle", width = 500),
            shiny::plotOutput(outputId = "plot_search_match_bar", width = 500),
            shiny::div(id = "div_search_match"),
            shiny::textOutput(outputId = "text_search_match")
        )
    )
)

server <- function(input, output, session) {
    species <- NULL
    bold_data <- NULL
    filtered_bold_data <- NULL
    blast_querys <- NULL

    trainingSet <- NULL
    testSequence <- NULL

    on.exit(update_viewer_list())

    shiny::observeEvent(input$select_download_online_or_load_local, {
        changed_state <- FALSE

        if (!is.null(species)) {
            species <<- NULL
            changed_state <- TRUE
        }

        if (!is.null(bold_data) || !is.null(filtered_bold_data)) {
            bold_data <<- NULL
            filtered_bold_data <<- NULL
            changed_state <- TRUE
        }

        if (changed_state) {
            shiny::showNotification("Unloaded species list and/or bold data.", duration = 15, type = "message")
        }

        output$DT_species_list <- DT::renderDataTable({
            data.frame()
        })
        output$DT_bold_data <- DT::renderDataTable({
            data.frame()
        })
    })

    shiny::observeEvent({
        input$file_species_list
        input$checkbox_species_file_header
        1
        }, {
        shiny::req(input$file_species_list)

        tryCatch({
            species <<- read.table(file = input$file_species_list$datapath, header = input$checkbox_species_file_header, sep = ",", comment.char = "", stringsAsFactors = FALSE)[[1]]
        }, error = function(err) {
            shiny::showNotification(paste0("Could not load file.", err$message), duration = 15, type = "error")
        })

        output$DT_species_list <- DT::renderDataTable(options = list(scrollX = TRUE, paging = FALSE), rownames = FALSE, filter = "none", {
            data.frame("Species" = species)
        })
    })

    shiny::observeEvent(input$button_download_bold_data, {
        if (is.null(species)) {
            shiny::showNotification("No species list loaded.", duration = 30, type = "error")
            return()
        }

        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Downloading bold data. This could take several minutes.", value = 0)

        tryCatch({
            res <<- download_species_list(species, progress)
        }, error = function(err) {
            shiny::showNotification(paste0("Failed to download bold data: ", err$message), duration = 0, type = "error")
            res <<- NULL
        })

        if (is.null(res)) {
            bold_data <<- NULL
            filtered_bold_data <<- NULL
            update_viewer_list()
            return()
        }

        bold_data <<- res$data

        update_viewer_list()

        lapply(res$message, function(x) {
            shiny::showNotification(x, duration = 0, type = "warning")
        })

        shiny::showNotification(paste0("Downloaded: ", nrow(bold_data), " samples from bold systems."), duration = 0, type = "message")

        write.table(bold_data, file = "bold.tsv", sep = "\t", row.names = FALSE)
        shiny::showNotification("Exported bold data to 'bold.tsv'", duration = 0, type = "message")

        output$DT_bold_data <- DT::renderDataTable(options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 250, nrow(bold_data) %>% min(500) %>% max(250)), scrollX = TRUE), rownames = FALSE, filter = "top", {
            bold_data
        })

    })

    shiny::observeEvent(input$file_bold_data, {
        bold_data <<- NULL
        filtered_bold_data <<- NULL
        shiny::req(input$file_bold_data)

        tryCatch({
            bold_data <<- read.table(file = input$file_bold_data$datapath, header = TRUE, sep = "\t", comment.char = "", quote = "\"")
            species <<- bold_data$species_name %>% unique

            if ("filter_threshold" %in% colnames(bold_data)) {
                filtered_bold_data <<- data.frame(bold_data)
                updateNumericInput(session, "numeric_keep_the_longest_n", value = bold_data$filter_threshold[[1]])
            }
        }, error = function(err) {
            shiny::showNotification(paste0("Could not load local bold data: ", err$message), duration = 30, type = "error")
        })

        update_viewer_list()

        if (is.null(bold_data)) {
            return()
        }

        shiny::showNotification(paste0("Loaded: ", nrow(bold_data), " samples from local bold file."), duration = 0, type = "message")
        shiny::showNotification(paste0("Found: ", length(species), " species in local bold data."), duration = 0, type = "message")

        output$DT_bold_data <- DT::renderDataTable(options = list(pageLength = 25, lengthMenu = c(25, 50, 100, 250, nrow(bold_data) %>% min(500) %>% max(250)), scrollX = TRUE), rownames = FALSE, filter = "top", {
            bold_data
        })
    })

    shiny::observeEvent(input$numeric_keep_the_longest_n, {
        update_viewer_list()
    })

    update_viewer_list <- function() {
        if (is.null(bold_data) || is.null(species)) {
            output$ui_species_viewer_holder <- shiny::renderUI({
                shiny::h3("No bold data loaded.")
            })
            output$DT_species_viewer <- DT::renderDataTable({
                data.frame()
            })
            return()
        }

        only_europe_bold_data <- europe_filtered_bold_data(species)
        filtered_bold_data <<- data.frame()

        for (spe in species) {
            if ((only_europe_bold_data$species_name == spe) %>% sum >= input$numeric_keep_the_longest_n) {
                species_mask <- only_europe_bold_data$species_name == spe
                if (!any(species_mask)) {
                    next
                }
                sorted_order <- sort.list(only_europe_bold_data[species_mask, ]$base_number, decreasing = TRUE)
                indexes <- head(sorted_order, n = input$numeric_keep_the_longest_n)
                filtered_bold_data <<- filtered_bold_data %>% rbind(only_europe_bold_data[species_mask, ][indexes, ])
            } else {
                species_mask <- bold_data$species_name == spe
                if (!any(species_mask)) {
                    next
                }
                sorted_order <- sort.list(bold_data[species_mask, ]$base_number, decreasing = TRUE)
                indexes <- head(sorted_order, n = input$numeric_keep_the_longest_n)
                filtered_bold_data <<- filtered_bold_data %>% rbind(bold_data[species_mask, ][indexes, ])
            }
        }

        filtered_bold_data$filter_threshold <<- input$numeric_keep_the_longest_n

        output$ui_species_viewer_holder <- shiny::renderUI({})
        output$DT_species_viewer <- DT::renderDataTable(options = list(scrollX = TRUE, paging = FALSE), rownames = FALSE, filter = "top", escape = FALSE, selection = "none", {
            df <- data.frame("species" = species)
            df$total_count <- species %>% lapply(function(spe) {
                (bold_data$species_name == spe) %>% sum
            })
            df$samples_from_europe <- species %>% lapply(function(spe) {
                (only_europe_bold_data$species_name == spe) %>% sum
            })
            df$using_only_european_samples <- species %>% lapply(function(spe) {
                ((only_europe_bold_data$species_name == spe) %>% sum) >= input$numeric_keep_the_longest_n
            })
            df$kept_samples <- species %>% lapply(function(spe) {
                (filtered_bold_data$species_name == spe) %>% sum
            })

            df
        })
    }

    shiny::observeEvent(input$button_export_filtered_data, {
        if (is.null(filtered_bold_data)) {
            shiny::showNotification("There is no filtered bold data to save.", duration = 10, type = "error")
            return()
        }

        write.table(filtered_bold_data, file = "bold_filtered.tsv", sep = "\t", row.names = FALSE)
        shiny::showNotification("Exported to 'bold_filtered.tsv'", duration = 0, type = "message")
    })

    europe_filtered_bold_data <- function(keep_only_europe) {
        bold_filtered <- data.frame(bold_data)

        bold_filtered <- bold_filtered[!(bold_filtered$species_name %in% keep_only_europe) 
            | (bold_filtered$country != "" 
            & !is.na(bold_filtered$country) 
            & (bold_filtered$country %in% europe)),]

        bold_filtered
    }

    shiny::observeEvent(input$button_generate_blast_query, {
        if (is.null(filtered_bold_data)) {
            shiny::showNotification("There is no filtered bold data", duration = 10, type = "error")
            return()
        }

        shiny::removeUI(selector = "#div_show_fa_querys > *", multiple = TRUE)

        blast_querys <<- generate_blast_query(species, filtered_bold_data, input$numeric_keep_the_longest_n)
        if (length(blast_querys) == 0) {
            shiny::showNotification("No blast query needed.", duration = 15, type = "message")
            return()
        }

        shiny::showNotification("Generated blast query files", duration = 15, type = "message")

        lapply(names(blast_querys), function(spe) {
            fa_text <- blast_querys[[spe]] %>% apply(1, function(x) {
                x <- unlist(x)
                paste0(">", spe, "\n", x["nucleotides"])
            }) %>% unlist %>% paste0(collapse = "\n")

            shiny::insertUI(selector = "#div_show_fa_querys", where = "beforeEnd", ui = shiny::h4(paste0("Showing querys for: ", spe)))
            shiny::insertUI(selector = "#div_show_fa_querys", where = "beforeEnd", ui = shiny::p(fa_text))
            shiny::insertUI(selector = "#div_show_fa_querys", where = "beforeEnd", ui = shiny::tags$hr())
        })
    })

    shiny::observeEvent(input$button_train_decipher, {
        if (is.null(filtered_bold_data)) {
            shiny::showNotification("No filtered bold data.", duration = 15, type = "error")
            return()
        }

        generate_decipher_fasta(filtered_bold_data)
        shiny::showNotification("Generated fasta file from bold data. Exported to 'bold.fasta'", duration = 0, type = "message")
        generate_decipher_taxid(filtered_bold_data)
        shiny::showNotification("Generated taxid file from bold data. Exported to 'taxid.txt'", duration = 0, type = "message")

        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Training DECIPHER Algorithm", value = 0)

        res <- train_decipher(seqs_path = "bold.fasta", rank_path = "taxid.txt", maxGroupSize = Inf, maxIterations = 1, shinyProgress = progress)

        if (length(res$message) > 0) {
            lapply(res$message, function(mes) {
                shiny::showNotification(mes, duration = 0, type = "message")
            })
        }

        trainingSet <<- res$trainingSet
        save(trainingSet, file = "trainingSet.Rdata")
        output$text_trainingSet_message <- shiny::renderText("Training Set generated from data")

        output$plot_trainingSet <- shiny::renderPlot({ 
            plot(res$trainingSet) 
        }, res = 144)

        output$DT_problematic_sequences <- DT::renderDataTable(options = list(scrollX = TRUE, paging = FALSE), rownames = FALSE, filter = "none", {
            diffs <- apply(res$trainingSet$problemSequences, 1, function(row) {
                partsExpected <- strsplit(row["Expected"], ";") %>% unlist
                 partsPredicted <- strsplit(row["Predicted"], ";")  %>% unlist

                diff <- partsExpected != partsPredicted
                data.frame(expected = paste0(partsExpected[diff], collapse = ";"), paste0(predicted = partsPredicted[diff], collapse = ";"))
            })

            df <- as.data.frame(do.call(rbind, diffs))
            colnames(df) <- c("Expected", "Predicted")
            df
        })
    })

    shiny::observeEvent(input$button_perform_blastn, {
        if (is.null(blast_querys) || length(blast_querys) == 0) {
            shiny::showNotification("No Blast Querys found.", duration = 15, type = "error")
            return()
        }

        system("./perform_blast_query.sh")
        shiny::showNotification("Blast search done", duration = 0, type = "message")

        blast <- read.table(file = "blast_result.out", header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, quote = "\"")
        # TODO 
        # The 'accession' column in the blast data is named something else
        # blast <- blast[!(blast$accession %in% bold_data$genbank_accession)]

        output$DT_blast_results <- DT::renderDataTable(options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 250, nrow(blast) %>% min(1000) %>% max(250)), scrollX = TRUE), rownames = FALSE, filter = "top", {
            blast
        })

        # here the blast_result.out should probably be parsed and filtered (to remove samples already in bold)
    })

    shiny::observeEvent(input$button_merge_blast_data, {
        shiny::showNotification("This has not been implemented yet.", duration = 15, type = "warning")
        # TODO
    })

    shiny::observeEvent(input$file_own_fasta, {
        shiny::showNotification("This has not been implemented yet.", duration = 15, type = "warning")
        # Load the fasta file and append it to the merged fasta file
        # Maybe check if the sequence already exists in a 100% match
        # TODO
    })

    shiny::observeEvent(input$file_new_unknown_fasta, {
        shiny::req(input$file_new_unknown_fasta)

        testSequence <<- readDNAStringSet(input$file_new_unknown_fasta$datapath)
        testSequence <<- RemoveGaps(testSequence)
    })

    shiny::observeEvent(input$button_search_for_matches, {
        if (is.null(testSequence)) {
            shiny::showNotification("No FASTA data to serach with", duration = 15, type = "error")
            return()
        }
        if (is.null(trainingSet)) {
            shiny::showNotification("No trained DECIPHER data to search against", duration = 15, type = "error")
            return()
        }

        ids <- IdTaxa(testSequence, trainingSet, processors = NULL, threshold = input$numeric_decipher_threshold, strand = input$select_decipher_strand, type = "extended")

        output$plot_search_match_circle <- shiny::renderPlot({
            plot(ids)
        })

        output$plot_search_match_bar <- shiny::renderPlot({

        })

        matching <- lapply(ids, function (id) {
            paste0(id$taxon,
                " (",
                round(id$confidence, digits=1),
                "%)", 
                collapse="; ")
        }) %>% unlist


        shiny::removeUI(selector = "#div_search_match > *", multiple = TRUE)

        paste0(names(ids), " : ", matching) %>% lapply(function(x) {
            shiny::insertUI(selector = "#div_search_match", where = "beforeEnd", ui = shiny::h4(x)) 
        })
    })

    shiny::observeEvent(input$file_trainingSet, {
        shiny::req(input$file_trainingSet)

        load(input$file_trainingSet$datapath)
        trainingSet <<- trainingSet

        output$text_trainingSet_message <- shiny::renderText("Training Set loaded from .Rdata file")
        shiny::showNotification("Loaded training data", duration = 0, type = "message")
    })
}

shiny::shinyApp(ui = ui, server = server)
