convertHgncToEnsembl <- function(hgncIds){
  
  library(biomaRt)
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='www.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='external_gene_name',
               values=hgncIds,
               mart=ensembl)
  return(genes)
}