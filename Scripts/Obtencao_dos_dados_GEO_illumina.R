# Pacotes necessários
library("GEOquery")
library("limma")
library('u133x3pcdf')
library('u133x3p.db')
library('affy')
library('oligo')
library('utils')
library('dplyr')

# Configurando o tempo do timeout
options(timeout = Inf)

# Baixando os dados
getGEOSuppFiles("GSE13861", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26899", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26901", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE28541", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26253", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE84426", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE84433", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE147163", baseDir = "Projeto recorrencia/Dados brutos")

# Descompactando os dados
untar("Projeto recorrencia/Dados brutos/GSE13861_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE13861")
untar("Projeto recorrencia/Dados brutos/GSE26899_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26899")
untar("Projeto recorrencia/Dados brutos/GSE26901_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26901")
untar("Projeto recorrencia/Dados brutos/GSE28541_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE28541")
untar("Projeto recorrencia/Dados brutos/GSE26253_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26253")
untar("Projeto recorrencia/Dados brutos/GSE84426_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE84426")
untar("Projeto recorrencia/Dados brutos/GSE84433_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE84433")
untar("Projeto recorrencia/Dados brutos/GSE147163_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE147163")

############################ Obtenção dos dados ###############################


obtendo_dados_il <- function(codigos) {
  for (cod in codigos){

    # Obtendo os dados
    gset <- getGEO(cod, GSEMatrix =TRUE, getGPL=FALSE)

    # Configurando os dados
    plataforma <- gset[[1]]@annotation
    if (length(gset) > 1) idx <- grep(plataforma, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # Obtendo os metadados 
    meta <- gset@phenoData@data
    meta$estudo <- cod
    meta$plataforma <- plataforma

    # Obtendo os valores de expressão
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    ex <- as.data.frame(ex)

    #Anotação
    gset_1 <-  getGEO(cod, GSEMatrix = TRUE)
    feature.data <- gset_1[[1]]@featureData@data
    feature.data <- feature.data %>% dplyr::select("ID", "ILMN_Gene")
    normalized.expr <- as.data.frame(ex)
    normalized.expr <- normalized.expr %>%
      tibble::rownames_to_column(var = 'ID') %>%
      inner_join(., feature.data, by = 'ID')

    # Lidando com genes repetidos
    normalized.expr <- as.data.frame(avereps(x=normalized.expr, ID=normalized.expr$ILMN_Gene))
    rownames(normalized.expr) <- normalized.expr$ILMN_Gene

    # Retirando as colunas que naão são mais necessárias
    normalized.expr <- normalized.expr |> relocate(ILMN_Gene) 
    normalized.expr <- normalized.expr[,-c(1, 2)]
    normalized.expr[] <- sapply(normalized.expr, as.numeric)

    # Salvando os dados
    write.csv(normalized.expr, paste0("Dados normalizados/Dados separados/expr_il_", cod, ".csv"))
    write.csv(meta, paste0("Dados clinicos/Dados separados/clin_il_", cod, ".csv"))
  }
}

obtendo_dados_il(c("GSE13861",
                   "GSE26899",
                   "GSE26901",
                   "GSE28541",
                   "GSE26253",
                   "GSE84426",
                   "GSE84433",
                   "GSE147163"))
