
#' @import bigrquery
#' @import NGCHM
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
#' @export
getExpressionData <- function (cohort, genes) {
  querySql <- sprintf ("
    SELECT
      HGNC_gene_symbol, gene_id, ParticipantBarCode, AliquotBarCode,
      ASINH(normalized_count) AS expression
    FROM %s
    WHERE (%s AND %s)
  ", ISB.rnaTable, testGene(genes), testParticipant(cohort));
  qres <- query_exec(querySql, max_pages=Inf, project = getCloudProject());

  samples <- unique(qres$AliquotBarCode)
  genes <- unique(qres$HGNC_gene_symbol)
  dm <- matrix(NA, nrow=length(genes), ncol=length(samples))
  rownames(dm) <- genes
  colnames(dm) <- samples
  dm[cbind(qres$HGNC_gene_symbol,qres$AliquotBarCode)] <- qres$expression;
  dm
}

#' Create a canned NG-CHM from ISB CGC cloud data
#'
#' This function creates a canned gene expression NG-CHM for the given participant cohort and gene list using data
#' from the ISB CGC cloud pilot.  If the 'viewer' option is defined, it will be used to display the NG-CHM.
#' The generated NG-CHM will contain three layers: row-centered, Z-normalized, and original.
#'
#' @param name Name the generated NG-CHM 'name'.
#' @param cohort Include samples from the specified vector of participant ids.
#' @param genes Include genes from the specified vector of HGNC symbols.
#' @param caption An informative caption
#' @return The generated NG-CHM
#' @export
#'
#' @seealso getExpressionData
#' @seealso demoCHM
exprCHM <- function (name, cohort, genes, caption=NULL) {
    data <- getExpressionData (cohort, genes);
    cent <- data - apply (data, 1, function(x)mean(x,na.rm=TRUE));
    norm <- cent / apply (cent, 1, function(x)sd(x,na.rm=TRUE));
    chm <- chmNew (name, cent, norm, data);
    chm <- chmAddAxisType (chm, 'row', 'bio.gene.hugo');
    if (!is.null(caption)) {
        chm <- chmAddProperty (chm, 'chm.info.caption', caption);
    }
    chmMake (chm);
    chmInstall (chm);
    v <- getOption('viewer');
    if (!is.null(v)) v(chmGetURL(chm));
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
             getStudyCohort(study), getReferenceGenes(authority),
             caption=sprintf ('mRNA expression data for TCGA study(s) %s', paste(study,collapse=' and ')))
}

#' @export
plot.ngchmVersion2 <- function(chm) {
    v <- getOption('viewer');
    if (is.null(v)) v <- browseURL;
    v(chmGetURL(chm))
}
