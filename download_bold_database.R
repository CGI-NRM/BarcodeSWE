

grades_checklist <- function(checklist, inputz, coordz) {
  species <- unique(checklist[, 1])
  x <- length(species)
  taxon_total = data.frame()
  y = x %/% 300 + 2
  i <- 1
  while (i < y) {
    ini <- 1 + (300 * (i - 1))
    fin <- 300 * i
    taxa <- bold::bold_seqspec(taxon=species[ini:fin], format = "tsv")
    taxon_total <- rbind(taxon_total, taxa)
    i = i + 1
  }
  taxon <- taxon_total
  taxon <- taxon[taxon$species_name!=""|is.na(taxon$species_name),]
  taxon <- taxon[!(taxon$bin_uri == "" | is.na(taxon$bin_uri)), ]
  taxon <- left_join(taxon,bps,by="species_name")
  taxon$bin_uri.y=NULL
  names(taxon)[names(taxon) == "bin_uri.x"] <- "bin_uri"
  taxon <- left_join(taxon,spb,by="bin_uri")
  taxon$species_name.y=NULL
  names(taxon)[names(taxon) == "species_name.x"] <- "species_name"
  taxon <- taxon[taxon$markercode=="COI-5P",]
  taxon$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon$nucleotides)
  taxon$nucleotides=gsub("-","",taxon$nucleotides)
  patterns_non <- c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
  for (p in patterns_non){
    taxon$species_name <- gsub(p,"",taxon$species_name)
  }
  patterns_fixed <- c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
  for (p in patterns_fixed){
    taxon$species_name <- gsub(p,"",taxon$species_name,fixed=TRUE)
  }
  taxon <- data.frame(taxon$species_name,taxon$bin_uri,taxon$nucleotides,taxon$country,taxon$family_name,taxon$order_name,taxon$class_name,taxon$sampleid,taxon$processid,taxon$species_per_bin,taxon$bin_per_species, taxon$lat,taxon$lon,taxon$genbank_accession)
  names(taxon) <- c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession")
  taxon$grade=NA
  names(taxon)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession","grade")
  taxon=data.frame(taxon$species,taxon$BIN,taxon$sequence,taxon$country,taxon$grade,taxon$family,taxon$order,taxon$class,taxon$sampleid,taxon$processid,taxon$species_per_bin,taxon$bin_per_species,taxon$lat,taxon$lon,taxon$genbank_accession)
  names(taxon)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","genbank_accession")
  taxon <- taxon %>%
    mutate(grade = ifelse(species_per_bin>1,"E",
                          ifelse(bin_per_species>1 & species_per_bin==1,"C",
                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "AB"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  taxon$contains_dominant=NULL
  if (coordz) {
    taxon <- taxon[!(is.na(taxon$lat)) | taxon$country!="",]
  }
  taxon$base_number=str_count(taxon$sequence, pattern="[A-Z]")
  taxon <- taxon[(taxon$base_number>=inputz),]
  np=(str_count(taxon$sequence, "N")/str_count(taxon$sequence, "[A-Z]"))*100
  taxon$n_percent=np
  taxon <- subset(taxon,taxon$n_percent<1)
  taxon$n_percent=NULL
  taxon=subset(taxon, lengths(gregexpr("\\w+", taxon$species)) > 1)
  num_species=table(taxon$species)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon <- inner_join(taxon,num_species)
  taxon <- taxon %>%
    mutate(grade = ifelse(grade=="E","E",
                          ifelse(frequency_species<3,"D",
                                 ifelse(grade=="C","C",
                                        ifelse(grade=="AB" & frequency_species<11,"B",
                                               ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "D"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "B"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  dominant_grade <- "A"
  dt <- as.data.table(taxon)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon <- setDF(dt)
  taxon$contains_dominant=NULL
  taxon$BIN_per_species=NULL
  taxon$species_per_bin=NULL
  taxon$species_frequency=NULL
  taxon$species_per_bin=NULL
  taxon$bin_per_species=NULL
  taxon$frequency_species=NULL
  taxon <- taxon[order(taxon$species),]
  assign('taxon',taxon,envir=.GlobalEnv)
}


spec <- read.table(file = "species.txt", sep = "\t", comment.char = '"', header = TRUE)
grades_checklist(spec, 300, TRUE)
