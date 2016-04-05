
#' @import bigrquery
#' @import NGCHM
NULL

getCloudProject <- function() {
   getOption('cloudproject')
}

.onAttach <- function(libname, pkgname) {
   proj <- getCloudProject();
   if (length(proj)==0) {
       cat ("Option cloudproject is not set\n");
   } else {
       cat ("Option cloudproject is ", proj, "\n");
   }
}

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
testStudy <- function (study) testColumn ('Study', study);

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

#' @export
exprCHM <- function (name, cohort, genes) {
    data <- getExpressionData (cohort, genes);
    cent <- data - apply (data, 1, function(x)mean(x,na.rm=TRUE));
    norm <- data / apply (cent, 1, function(x)sd(x,na.rm=TRUE));
    chm <- chmNew (name, cent, norm, data);
    chm <- chmAddAxisType (chm, 'row', 'bio.gene.hugo');
    chmMake (chm);
    chmInstall (chm);
    getOption('viewer')(chmGetURL(chm));
    chm
}

#' @export
demoCHM <- function() {
    exprCHM ('bmbroom-isb-demo', getStudyCohort('prad'), getReferenceGenes('Vogelstein'))
}
