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

#################################################### Processamento dos dados ####################################################


