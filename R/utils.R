
#' @import bigrquery
#' @import NGCHM
#' @import magrittr
NULL

# ISB-CGC tables of interest
ISB.clinicTable <- "[isb-cgc:tcga_201510_alpha.Clinical_data]";
ISB.rnaTable <- "[isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]";

#' @export
asSqlList <- function (items) {
    sprintf ("'%s'", paste(items,collapse="','"))
}

#' @export
testColumn <- function (colname, values) {
    stopifnot (length(values) > 0);
    paste ("(",
           colname,
           if (length(values)==1) "=" else " IN (",
           asSqlList (values),
           if (length(values)==1) "" else ")",
           ")", sep="")
}

#' @export
testStudy <- function (study) testColumn ('Study', toupper(study));

#' @export
testGene <- function (genes) testColumn ('HGNC_gene_symbol', genes);

#' @export
testParticipant <- function (partids) testColumn ('ParticipantBarCode', partids);

#' @export
getStudyCohort <- function(study) {
  querySql <- sprintf ("SELECT ParticipantBarcode FROM %s WHERE %s", ISB.clinicTable, testStudy(study));
  query_exec(querySql, project=getCloudProject())$ParticipantBarcode
}



#' @export
getGeneSymbols <- function() {
  querySql <- sprintf ("
    SELECT
      HGNC_gene_symbol, original_gene_symbol, gene_id
    FROM
      %s
    WHERE
      ( original_gene_symbol IS NOT NULL
        AND HGNC_gene_symbol IS NOT NULL
        AND original_gene_symbol=HGNC_gene_symbol
        AND gene_id IS NOT NULL )
    GROUP BY
      original_gene_symbol, HGNC_gene_symbol, gene_id
    ORDER BY
      HGNC_gene_symbol", ISB.rnaTable);

  query_exec(querySql, project=getCloudProject())
}

#' @export
getReferenceGenes <- function(authority) {
  authtest <- if (length(authority)==0) "" else { paste ("AND", testColumn('authority',authority), sep=" ") };
  querySql <- sprintf ("
    SELECT
      symbol
    FROM
      [ngchm-cloud-pilot:reference.allonco]
    WHERE
      ( symbol IS NOT NULL %s )
  ", authtest);

  query_exec(querySql, project=getCloudProject())$symbol
}

#' @export
getBadGeneSymbols <- function() {
  querySql <- sprintf ("
    SELECT
      HGNC_gene_symbol, original_gene_symbol, gene_id
    FROM
      %s
    WHERE
      ( original_gene_symbol IS NOT NULL
        AND HGNC_gene_symbol IS NOT NULL
        AND original_gene_symbol!=HGNC_gene_symbol
        AND gene_id IS NOT NULL )
    GROUP BY
      original_gene_symbol, HGNC_gene_symbol, gene_id
    ORDER BY
      HGNC_gene_symbol", ISB.rnaTable);

  query_exec(querySql, project = getCloudProject())
}

#' Get gene expression data from the ISB CGC cloud
#'
#' This function downloads gene expression data from the ISB CGC cloud as
#' a gene by sample matrix.  The data will be transformed using ASINH.
#'
#' @param cohort Include samples from the specified vector of participant ids.
#' @param genes Include genes from the specified vector of HGNC symbols.
#' @return A numeric matrix of genes (rows) by samples (columns).
#' @import digest
#' @export
getExpressionData <- function (cohort, genes) {
  cohort <- sort(cohort);
  genes <- sort(genes);
  querySql <- sprintf ("
    SELECT
      HGNC_gene_symbol, gene_id, ParticipantBarCode, AliquotBarCode,
      ASINH(normalized_count) AS expression
    FROM %s
    WHERE (%s AND %s)
  ", ISB.rnaTable, testGene(genes), testParticipant(cohort));

  filename <- file.path ("data", sprintf ("q%s", digest::digest(querySql)));
  useCache <- getOption ("usecache", TRUE);
  if (useCache && file.exists(filename)) {
      ee <- new.env();
      id <- load (filename, ee);
      dm <- ee[[id]];
  } else {
      qres <- query_exec(querySql, max_pages=Inf, project = getCloudProject());

      samples <- unique(qres$AliquotBarCode)
      genes <- unique(qres$HGNC_gene_symbol)
      dm <- matrix(NA, nrow=length(genes), ncol=length(samples))
      rownames(dm) <- genes
      colnames(dm) <- samples
      dm[cbind(qres$HGNC_gene_symbol,qres$AliquotBarCode)] <- qres$expression;
      if (useCache) save (dm, file=filename);
  }
  dm
}

#' Get the specified columns from the clinical data for the specified cohort
#'
#' @param cohort Include samples from the specified vector of participant ids.
#' @param columns Include columns from the specified vector of clinical table column names.
#' @return A list of column vectors.
#' @export
getClinicalData <- function (cohort, columns) {
  cohort <- sort(unique(cohort));
  result <- list();
  useCache <- getOption ("usecache", TRUE);

  if (useCache) {
      for (column in columns) {
          filename <- file.path ("data", sprintf ("cc%s", digest::digest(list(ISB.clinicTable,cohort,column))));
	  if (file.exists(filename)) {
	      ee <- new.env();
	      id <- load (filename, ee);
	      result[[column]] <- ee[[id]];
          }
      }
  }
  columns <- setdiff (columns, names(result));
  if (length(columns) > 0) {
      querySql <- sprintf ("
        SELECT
          ParticipantBarCode,
          %s
        FROM %s
        WHERE (%s)
      ", paste(columns,collapse=","), ISB.clinicTable, testParticipant(cohort));
      qres <- query_exec(querySql, max_pages=Inf, project = getCloudProject());
      samples <- unique(qres$ParticipantBarCode);
      for (column in columns) {
          filename <- file.path ("data", sprintf ("cc%s", digest::digest(list(ISB.clinicTable,cohort,column))));
          val <- rep(NA,length(cohort));
          names(val) <- cohort;
          val[qres$ParticipantBarCode] <- qres[[column]];
          if (useCache) save (val, file=filename);
          result[[column]] <- val;
      }
  }
  result
}

#' @export
createClinicalCovariate <- function (clinicTable, chm, fullname, column, ...) {
    colLabels <- ngchmGetLabelsStr (chm@layers[[1]]@data,"column");
    vals <- clinicTable[[column]][substr(colLabels,1,12)];
    names(vals) <- colLabels;
    chmNewCovariate(fullname, vals, ...)
}

isb.env <- new.env();
covariates.file <- file.path ("data", "clinical-covariates.Rda");

initClinicalCovariates <- function() {
    if (file.exists (covariates.file)) {
	ee <- new.env();
	load (covariates.file, ee);
	isb.env$covariates <- ee$covariates;
    } else {
	isb.env$covariates <- rbind(
	    list ('Vital status', 'vital_status', NULL),
	    list ('Followup (days)', 'days_to_last_followup', NULL),
	    list ('PSA', 'psa_value', 'prad'));
	colnames(isb.env$covariates) <- c("fullName", "columnName", "studies");
    }
}

#'  Define a clinical covariate for use in NGCHMs
#'
#' @param fullName Display name of covariate
#' @param columnName Name of covariate in ISB CGC clinical table
#' @param studies Vector of TCGA studies for which covariate is valid. NULL (default) is valid for all.
#'
#' @export
#' @seealso removeClinicalCovariate
defineClinicalCovariate <- function (fullName, columnName, studies) {
    if (missing(studies)) studies <- NULL;
    newVal <- list(fullName, columnName, studies);
    if (columnName %in% isb.env$covariates[,'columnName']) {
        isb.env$covariates[which(columnName==isb.env$covariates[,'columnName']),] <- newVal;
    } else {
        isb.env$covariates <- rbind (isb.env$covariates, newVal);
    }
    with (isb.env, save (covariates, file=covariates.file))
}

#'  Remove a clinical covariate definition
#'
#' @param name Name of covariate to remove
#'
#' @export
#' @seealso defineClinicalCovariate
removeClinicalCovariate <- function (name) {
    if (name %in% isb.env$covariates[,'columnName']) {
        idx <- which(name==isb.env$covariates[,'columnName']);
        if (name %in% isb.env$covariates[,'fullName']) {
            idx2 <- which(name==isb.env$covariates[,'fullName']);
            if (idx != idx2) {
                stop (sprintf ('Covariate "%s" is ambiguous', name))
            }
        }
    } else if (name %in% isb.env$covariates[,'fullName']) {
        idx <- which(name==isb.env$covariates[,'fullName']);
    } else {
        stop (sprintf ('Unknown covariate "%s"', name))
    }
    isb.env$covariates <- isb.env$covariates[-idx,];
    with (isb.env, save (covariates, file=covariates.file))
}

#' Add clinical covariates to chm
#'
#' @param chm NG-CHM
#' @param study The TCGA study identifier(s) of the data samples
#' @param cohort The participant identifiers included in the NGCHM
#'
#' @export
addClinicalCovariates <- function (chm, study, cohort) {
    clinicTab <- getClinicalData (cohort, isb.env$covariates[,'columnName']);
    study <- tolower(study);
    cvs <- 1:nrow(isb.env$covariates) %>%
           Filter (function(row) {
               okstudies <- isb.env$covariates[row,]$studies;
               is.null(okstudies) || all(vapply(study,function(s)s %in% okstudies,TRUE))
           },.) %>%
           lapply (function(row) {
               createClinicalCovariate (clinicTab, chm, isb.env$covariates[row,]$fullName, isb.env$covariates[row,]$columnName)
           });
    chm + (chmAxis('col') + cvs)
}

#' Create a canned NG-CHM from ISB CGC cloud data
#'
#' This function creates a canned gene expression NG-CHM for the given participant cohort and gene list using data
#' from the ISB CGC cloud pilot.  If the 'viewer' option is defined, it will be used to display the NG-CHM.
#' The generated NG-CHM will contain three layers: row-centered, Z-normalized, and original.
#'
#' @param name Name the generated NG-CHM 'name'.
#' @param study TCGA study(s) of the data samples
#' @param cohort Include samples from the specified vector of participant ids.
#' @param genes Include genes from the specified vector of HGNC symbols.
#' @param caption An informative caption
#' @return The generated NG-CHM
#' @export
#'
#' @seealso getExpressionData
#' @seealso demoCHM
exprCHM <- function (name, study, cohort, genes, caption=NULL) {
    data <- getExpressionData (cohort, genes);
    cent <- data - apply (data, 1, function(x)mean(x,na.rm=TRUE));
    norm <- cent / apply (cent, 1, function(x)sd(x,na.rm=TRUE));
    chm <- chmNew (name, cent, norm, data) %>%
           chmAddAxisType ('row', 'bio.gene.hugo') %>%
           chmAddAxisType ('column', tcgaBarcodeType(colnames(data)[1]));
    if (!is.null(caption)) {
        chm <- chmAddProperty (chm, 'chm.info.caption', caption);
    }
    if (length(study)==1) {
        chm <- tcgaAddCBIOStudyId (chm, sprintf ('%s_tcga',study));
    }
    chm <- addClinicalCovariates (chm, study, cohort);
    if (length(chmListServers()) > 0) {
        chmMake (chm);
        chmInstall (chm);
        plot(chm);
    }
    chm
}

#' Create a canned NG-CHM from ISB CGC cloud data
#'
#' This function creates a canned gene expression NG-CHM for the given study and gene authority using data
#' from the ISB CGC cloud pilot.  If the 'viewer' option is defined, it will be used to display the NG-CHM.
#'
#' @param study Include samples from the specified study (or studies).
#' @param authority Include genes deemed 'cancer interesting' by the specified authority (or authorities).
#' @return The generated NG-CHM
#' @export
#'
#' @seealso exprCHM
demoCHM <- function(study='prad', authority='Vogelstein') {
    exprCHM (sprintf ('mrna-%s-%s', paste(study,collapse='+'), paste(authority,collapse='+')),
             study, getStudyCohort(study), getReferenceGenes(authority),
             caption=sprintf ('mRNA expression data for TCGA study(s) %s using genes defined by %s.', paste(study,collapse=' and '),
                              paste(authority,collapse=' and ')))
}
