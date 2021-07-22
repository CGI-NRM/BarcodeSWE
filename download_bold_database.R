
#'grades_checklist <- function(checklist,inputz,coordz){
#'    species <- unique(checklist[,1])
#'    x <- length(species)
#'    taxon_total = data.frame()
#'    y=x%/%300+2
#'    i <- 1
#'    while (i < y){
#'        ini <- 1+(300*(i-1))
#'        fin <- 300*i
#'        taxa <- bold_seqspec(taxon=species[ini:fin], format = "tsv")
#'        taxon_total <- rbind(taxon_total,taxa)
#'        i = i+1
#'    }
#'    taxon <- taxon_total
#'
#'    write.table(taxon_total, file = "out.tsv", sep = "\t")
#'
#'    taxon2 <- taxon[taxon$species_name!=""|is.na(taxon$species_name),]
#'    taxon2 <- taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
#'    taxon2 <- dplyr::left_join(taxon2,bps,by="species_name")
#'    taxon2$bin_uri.y=NULL
#'    names(taxon2)[names(taxon2) == "bin_uri.x"] <- "bin_uri"
#'    taxon2 <- dplyr::left_join(taxon2,spb,by="bin_uri")
#'    taxon2$species_name.y=NULL
#'    names(taxon2)[names(taxon2) == "species_name.x"] <- "species_name"
#'    taxon3 <- taxon2[taxon2$markercode=="COI-5P",]
#'    taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
#'    taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
#'    patterns_non <- c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
#'    for (p in patterns_non){
#'        taxon3$species_name <- gsub(p,"",taxon3$species_name)
#'    }
#'    patterns_fixed <- c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
#'    for (p in patterns_fixed){
#'        taxon3$species_name <- gsub(p,"",taxon3$species_name,fixed=TRUE)
#'    }
#'    taxon8 <- data.frame(taxon3$species_name,taxon3$bin_uri,taxon3$nucleotides,taxon3$country,taxon3$family_name,taxon3$order_name,taxon3$class_name,taxon3$sampleid,taxon3$processid,taxon3$species_per_bin,taxon3$bin_per_species, taxon3$lat,taxon3$lon,taxon3$genbank_accession)
#'    names(taxon8) <- c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession")
#'    taxon8$grade=NA
#'    names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession","grade")
#'    taxon9=data.frame(taxon8$species,taxon8$BIN,taxon8$sequence,taxon8$country,taxon8$grade,taxon8$family,taxon8$order,taxon8$class,taxon8$sampleid,taxon8$processid,taxon8$species_per_bin,taxon8$bin_per_species,taxon8$lat,taxon8$lon,taxon8$genbank_accession)
#'    names(taxon9)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession")
#'    taxon19 <- taxon9 %>%
#'        mutate(grade = ifelse(species_per_bin>1,"E",
#'                                                    ifelse(bin_per_species>1 & species_per_bin==1,"C",
#'                                                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
#'    dominant_grade <- "E"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "C"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "AB"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    taxon19$contains_dominant=NULL
#'    if (coordz){
#'        taxon19 <- taxon19[!(is.na(taxon19$lat)) | taxon19$country!="",]
#'    }
#'    taxon19$base_number=str_count(taxon19$sequence, pattern="[A-Z]")
#'    taxon19 <- taxon19[(taxon19$base_number>=inputz),]
#'    np=(str_count(taxon19$sequence, "N")/str_count(taxon19$sequence, "[A-Z]"))*100
#'    taxon19$n_percent=np
#'    taxon19 <- subset(taxon19,taxon19$n_percent<1)
#'    taxon19$n_percent=NULL
#'    taxon19=subset(taxon19, lengths(gregexpr("\\w+", taxon19$species)) > 1)
#'    num_species=table(taxon19$species)
#'    num_species=as.data.frame(num_species)
#'    names(num_species)=c("species","frequency_species")
#'    taxon19 <- inner_join(taxon19,num_species)
#'    taxon19 <- taxon19 %>%
#'        mutate(grade = ifelse(grade=="E","E",
#'                                                    ifelse(frequency_species<3,"D",
#'                                                                 ifelse(grade=="C","C",
#'                                                                                ifelse(grade=="AB" & frequency_species<11,"B",
#'                                                                                             ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
#'    dominant_grade <- "E"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "D"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "C"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "B"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    dominant_grade <- "A"
#'    dt <- as.data.table(taxon19)
#'    dt[, contains_dominant := any(grade == dominant_grade), by=species]
#'    dt[contains_dominant == TRUE, grade := dominant_grade]
#'    taxon19 <- setDF(dt)
#'    taxon19$contains_dominant=NULL
#'    taxon19$BIN_per_species=NULL
#'    taxon19$species_per_bin=NULL
#'    taxon19$species_frequency=NULL
#'    taxon19$species_per_bin=NULL
#'    taxon19$bin_per_species=NULL
#'    taxon19$frequency_species=NULL
#'    taxon19 <- taxon19[order(taxon19$species),]
#'    assign('taxon19',taxon19,envir=.GlobalEnv)
#'}


library(magrittr)

download_species_list <- function(species, progress = NULL) {
    message <- c()
    taxon_total = data.frame()
    for (i in 1:length(species)) {
        if (!is.null(progress)) {
            progress$inc(0, detail = paste0("Downloading: ", species[i]))
        }

        taxa <- bold::bold_seqspec(taxon=species[i], format = "tsv") %>% data.frame

        taxa[,colnames(taxon_total)[!colnames(taxon_total) %in% colnames(taxa)]] <- NA
        if (nrow(taxon_total) > 0) {
            taxon_total[,colnames(taxa)[!colnames(taxa) %in% colnames(taxon_total)]] <- NA
        }

        taxon_total <- rbind(taxon_total, taxa)

        if (!is.null(progress)) {
            progress$inc(1 / length(species))
        }

    }

    taxon <- taxon_total

    taxon <- taxon[taxon$markercode == "COI-5P", ]
    taxon$base_number <- stringr::str_count(taxon$nucleotides, pattern="[ATGC]")
    taxon <- taxon[taxon$base_number >= 300, ]
    taxon <- taxon[taxon$species_name != "" & !is.na(taxon$species_name), ]
    taxon <- taxon[taxon$nucleotides != ""  & !is.na(taxon$nucleotides), ]
    taxon <- taxon[taxon$species_name != "" & !is.na(taxon$species_name), ]
    taxon <- taxon[taxon$bin_uri != ""      & !is.na(taxon$bin_uri), ]

    found_in_bold <- species %>% lapply(function(x) {
            x %in% taxon$species_name
        }) %>% unlist

    if (any(!found_in_bold)) {
        message <- paste0("Could not find any data on '", species[!found_in_bold], "' in bold.")
        print(message)
    }

    return(list(data = taxon, message = message))
}

