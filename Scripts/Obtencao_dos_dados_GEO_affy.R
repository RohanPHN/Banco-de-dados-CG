# Pacotes necessários
library("GEOquery")
library("limma")
library('u133x3pcdf')
library('u133x3p.db')
library('affy')
library('oligo')
library('utils')
library('dplyr')

# Configurand o timeout para download dos dados
options(timeout = Inf)

# Baixando os dados
getGEOSuppFiles("GSE62254", baseDir = "Dados brutos")
getGEOSuppFiles("GSE29272", baseDir = "Dados brutos")
getGEOSuppFiles("GSE15459", baseDir = "Dados brutos")
getGEOSuppFiles("GSE38749", baseDir = "Dados brutos")
getGEOSuppFiles("GSE57303", baseDir = "Dados brutos")


# Descompactando os dados
untar("Dados brutos/GSE62254_RAW.tar",  exdir = "Dados brutos/GSE62254")
untar("Dados brutos/GSE29272_RAW.tar",  exdir = "Dados brutos/GSE29272")
untar("Dados brutos/GSE15459_RAW.tar",  exdir = "Dados brutos/GSE15459")
untar("Dados brutos/GSE38749_RAW.tar",  exdir = "Dados brutos/GSE38749")
untar("Dados brutos/GSE57303_RAW.tar",  exdir = "Dados brutos/GSE57303")

############################ Obtenção dos dados ###############################

obtendo_dados_geo <- function(codigos) {
  for (cod in codigos){
    # Lendo os arquivos
    raw.data <- ReadAffy(celfile.path = paste0('Dados brutos/', cod))
    
    # Executando a normalizaççao rma
    normalized.data <- affy::rma(raw.data)
    
    # Obtendo as estimativas de expressão
    normalized.expr <- as.data.frame(exprs(normalized.data))
    
    # Ajustando o nome das colunas
    colnames(normalized.expr) <- gsub("_.*","", colnames(normalized.expr))
    
    # Mapeando as sondas para simbolos de genes
    gse <- getGEO(cod, GSEMatrix = TRUE)
    feature.data <- gse[[1]]@featureData@data
    feature.data <- feature.data[,c(1,11)]
    normalized.expr <- normalized.expr %>%
      tibble::rownames_to_column(var = 'ID') %>%
      inner_join(., feature.data, by = 'ID')
    
    # Lidando com genes repetidos
    normalized.expr <- as.data.frame(avereps(x=normalized.expr, ID=normalized.expr$`Gene Symbol`))
    rownames(normalized.expr) <- normalized.expr$`Gene Symbol`
    
    # Retirando as colunas que não são mais necessárias
    normalized.expr <- normalized.expr |> relocate(`Gene Symbol`) 
    normalized.expr <- normalized.expr[,-c(1, 2)]
    
    # Obtendo os dados clinicos
    clinico <- gse[[1]]@phenoData@data
    clinico$estudo <- cod
    
    # Salvando os dados
    write.csv(normalized.expr, paste0("Dados normalizados/Dados separados/expr_affy_", cod, ".csv"))
    write.csv(clinico, paste0("Dados clinicos/Dados separados/clin_affy_", cod, ".csv"))
  }
}

obtendo_dados_geo(c("GSE62254",
                    "GSE29272",
                    "GSE15459",
                    "GSE38749",
                    "GSE57303"))

