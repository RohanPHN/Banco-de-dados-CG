# Pacotes necessários
library("GEOquery")
library("limma")
library('u133x3pcdf')
library('u133x3p.db')
library('affy')
library('oligo')
library('utils')
library('dplyr')

# Indicando a pasta de trabalho
setwd("C:/Users/2526991/Documents/Paulo Rohan")
options(timeout = Inf)

# Baixando os dados
getGEOSuppFiles("GSE62254", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE29272", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE15459", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE38749", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE57303", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE13861", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26899", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26901", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE28541", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE26253", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE84426", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE84433", baseDir = "Projeto recorrencia/Dados brutos")
getGEOSuppFiles("GSE147163", baseDir = "Projeto recorrencia/Dados brutos")

# Descompactando os dados
untar("Projeto recorrencia/Dados brutos/GSE62254_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE62254")
untar("Projeto recorrencia/Dados brutos/GSE29272_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE29272")
untar("Projeto recorrencia/Dados brutos/GSE15459_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE15459")
untar("Projeto recorrencia/Dados brutos/GSE38749_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE38749")
untar("Projeto recorrencia/Dados brutos/GSE57303_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE57303")
untar("Projeto recorrencia/Dados brutos/GSE13861_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE13861")
untar("Projeto recorrencia/Dados brutos/GSE26899_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26899")
untar("Projeto recorrencia/Dados brutos/GSE26901_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26901")
untar("Projeto recorrencia/Dados brutos/GSE28541_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE28541")
untar("Projeto recorrencia/Dados brutos/GSE26253_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE26253")
untar("Projeto recorrencia/Dados brutos/GSE84426_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE84426")
untar("Projeto recorrencia/Dados brutos/GSE84433_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE84433")
untar("Projeto recorrencia/Dados brutos/GSE147163_RAW.tar",  exdir = "Projeto recorrencia/Dados brutos/GSE147163")

#################################################### Processamento dos dados ####################################################

#__________________________ GSE62254 __________________________# 

# Lendo os arquivos
raw.data.GSE62254 <- ReadAffy(celfile.path = 'Projeto recorrencia/Dados brutos/GSE62254')

# Executando a normalizaççao rma
normalized.data.GSE62254 <- affy::rma(raw.data.GSE62254)

# Obtendo as estimativas de expressão
normalized.expr.GSE62254 <- as.data.frame(exprs(normalized.data.GSE62254))

# Mapeando as sondas para simbolos de genes
gse.GSE62254 <- getGEO("GSE62254", GSEMatrix = TRUE)
feature.data.GSE62254 <- gse.GSE62254[[1]]@featureData@data
feature.data.GSE62254 <- feature.data.GSE62254[,c(1,11)]
normalized.expr.GSE62254 <- normalized.expr.GSE62254 %>%
  tibble::rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data.GSE62254, by = 'ID')

normalized.expr.GSE62254 <- as.data.frame(avereps(x=normalized.expr.GSE62254, ID=normalized.expr.GSE62254$`Gene Symbol`))
colnames(normalized.expr.GSE62254) <- gsub("_.*","", colnames(normalized.expr.GSE62254))
rownames(normalized.expr.GSE62254) <- normalized.expr.GSE62254$`Gene Symbol`
normalized.expr.GSE62254 <- normalized.expr.GSE62254[,-c(1, 302)]

