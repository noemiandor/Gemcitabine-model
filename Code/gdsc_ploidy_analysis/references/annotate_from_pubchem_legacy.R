annotate_from_pubchem <- function(coxIn, cmap_csv) {
  stopifnot(file.exists(cmap_csv))

  isempty <- function(x) {
    length(x) == 0 || all(is.na(x))
  }

  drugs <- read.csv(cmap_csv, stringsAsFactors = FALSE)

  coxIn$drugCategory_Pubchem <- NA
  udrugs <- unique(coxIn$drugName)
  for (drug in udrugs) {
    ij <- which(coxIn$drugName == drug)
    print(paste("Scraping info from PubChem for drug", drug, "(", which(drug == udrugs), "out of", length(udrugs), ")"))

    id_PubChem <- ChemmineR::pubchemName2CID(drug)[1]
    if (is.na(id_PubChem)) {
      id_PubChem <- unique(drugs$PubChem.CID[drugs$Name == drug])
    }
    if (!isempty(id_PubChem) && !is.na(id_PubChem)) {
      url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%d/JSON", as.numeric(id_PubChem))
      page <- textreadr::read_html(url)
      page <- strsplit(page, "\n", fixed = TRUE)[[1]]

      dcI <- grep("Drug Classes", page)
      if (!isempty(dcI)) {
        dc <- page[unlist(sapply(dcI, function(x) x:(10 + x), simplify = FALSE))]
        dc <- grep("\"String\":", dc, value = TRUE)
        dc <- sub('^.*"String":\\s*"', "", dc)
        dc <- sub('".*$', "", dc)
        dc <- dc[nzchar(dc)]

        coxIn$drugCategory_Pubchem[ij] <- paste(dc, collapse = "; ")
        print(paste(drug, ": ", dc))
        next
      }

      tmpI <- grep("\"TOCHeading\": \"Mechanism of Action\"", page)
      tmpI2 <- grep("echanism of action", page)
      tmpI3 <- grep(paste(drug, "is a"), page)
      if (!isempty(tmpI)) {
        tmpI <- tmpI:(tmpI + 10)
      }
      if (!isempty(tmpI2)) {
        tmpI2 <- tmpI2:(tmpI2 + 10)
      }
      page <- unlist(page[c(tmpI, tmpI2, tmpI3)])

      iI <- sapply(c("inflammatory", "immun", "lympho", "cytokine"), function(x) !isempty(grep(x, page, ignore.case = TRUE)))
      if (any(iI)) {
        coxIn$drugCategory_Pubchem[ij] <- "(Anti-)Inflammatory"
        print(paste(drug, ": (Anti-)Inflammatory"))
      }

      cI <- sapply(
        c("platinum", "Cytotox", "Microtubul", "Mitotic", "Spindle", "Alkylating"),
        function(x) !isempty(grep(x, page, ignore.case = TRUE))
      )
      if (any(cI)) {
        coxIn$drugCategory_Pubchem[ij] <- "Cytotoxic"
        print(paste(drug, ": Cytotoxic"))
      }

      sI <- sapply(c("MAPK", "ERK kinase", "MEK1", "MEK2", "MTOR", "WNT", "EGFR", "ROS"), function(x) !isempty(grep(x, page, ignore.case = FALSE)))
      sI <- c(sI, sapply(c("kinase inhibitor"), function(x) !isempty(grep(x, page, ignore.case = TRUE))))
      if (any(sI)) {
        coxIn$drugCategory_Pubchem[ij] <- "Signaling"
        print(paste(drug, ": Signaling"))
      }

      mI <- sapply(c("Metabolic Inhibitor", "detoxification"), function(x) !isempty(grep(x, page, ignore.case = TRUE)))
      if (any(mI)) {
        coxIn$drugCategory_Pubchem[ij] <- "MetabolicInhibitor"
        print(paste(drug, ": Metabolism or Hormones"))
      }
    }
  }

  coxIn
}
