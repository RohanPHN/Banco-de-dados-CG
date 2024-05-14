# Pacotes necessários
library("TCGAbiolinks")
library("DESeq2")
library("limma")
library("edgeR")
library("biomaRt")

# Indicando a pasta de trabalho
setwd("C:/Users/2526991/Documents/Paulo Rohan")

# Lendo os dados de RNA-Seq do TCGA-STAD
dat = readRDS(file = "Dados Brutos/DOUTORADO_mRNA_TCGA.RDS")

# Matriz de expressao
rna_raw <- assay(dat)

# Dados clinicos
clinical <- colData(dat)
clinical <- clinical[,!(colnames(clinical) %in% c("treatments", "primary_site", "disease_type"))]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$sample


# Organizando os nomes das colunas da matriz
delim_fn = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}

colnames(rna_raw) <- delim_fn(x = colnames(rna_raw), n = 0, i = 4)
rna_raw <- as.data.frame(rna_raw)

# Verificando a ordem dos dados
identical(colnames(rna_raw), rownames(clinical))

# Mantendo apenas os tumores
clinical <- subset(clinical, subset = clinical$definition == "Primary solid Tumor")
rna_raw <- rna_raw[,colnames(rna_raw) %in% rownames(clinical)]
identical(colnames(rna_raw), rownames(clinical))

# Ajeitando o código ensembl na matriz bruta
rna_raw <- tibble::rownames_to_column(rna_raw, "VALUE")
rna_raw$VALUE <- gsub("\\..*","", rna_raw$VALUE)
rna_raw <- as.data.frame(avereps(x=rna_raw, ID=rna_raw$VALUE))
rownames(rna_raw) <- rna_raw$VALUE
rna_raw <- rna_raw[,-1]
rna_raw[] <- sapply(rna_raw, as.numeric)

# Eliminando genes com contagem muito baixa.
c_cpm <- cpm(rna_raw)
keep <- rowSums(c_cpm) >= 0.5
rna_raw <- rna_raw[keep,]

# Criando o objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = rna_raw,
                              colData = clinical,
                              design = ~1)

# Obtendo os símbolos do gene
ens2symbol<-function(ids){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"), )
  genes <- getBM(attributes= c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                 values=ids, mart= mart)
  return(genes)
}

genes <- ens2symbol(row.names(dds))

# Obtendo as contagens normalizadas
dds <- estimateSizeFactors(dds)
rna_norm <- counts(dds, normalized=T)

# Obtendo as contagens vst. blind = TRUE devido ao fato da análise ser exploratória
vst_counts <- vst(dds, blind=TRUE)
vst_counts <- assay(vst_counts)

##### Obtendo os nomes dos genes #####

# Contagens normalizadas dos genes que codificam proteínas
rna_norm <- as.data.frame(rna_norm)
rna_norm <- tibble::rownames_to_column(rna_norm, "ensembl_gene_id")
rna_norm <- merge(genes,rna_norm, by = "ensembl_gene_id")
rna_norm <- rna_norm[, -c(1, 3)]
rna_norm <- as.data.frame(avereps(x=rna_norm, ID=rna_norm$hgnc_symbol))
rownames(rna_norm) <- rna_norm$hgnc_symbol
rna_norm <- rna_norm[,-1]
rna_norm[] <- sapply(rna_norm, as.numeric)

# Contagens vst dos genes que codificam proteínas
rna_vst <- as.data.frame(vst_counts)
rna_vst <- tibble::rownames_to_column(rna_vst, "ensembl_gene_id")
rna_vst <- merge(genes,rna_vst, by = "ensembl_gene_id")
rna_vst <- rna_vst[, -c(1, 3)]
rna_vst <- as.data.frame(avereps(x=rna_vst, ID=rna_vst$hgnc_symbol))
rownames(rna_vst) <- rna_vst$hgnc_symbol
rna_vst <- rna_vst[,-1]
rna_vst[] <- sapply(rna_vst, as.numeric)

# Salvando os dados
write.csv(rna_norm, "Projeto recorrencia/Dados normalizados/GC_TCGA_norm.csv", row.names = T)
write.csv(rna_vst, "Projeto recorrencia/Dados normalizados/GC_TCGA_vst.csv", row.names = T)
write.csv(clinical, "Projeto recorrencia/Dados clinicos/GC_TCGA_clinical.csv", row.names = T)