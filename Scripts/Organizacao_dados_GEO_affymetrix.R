# Pacotes necessários
library("GEOquery")
library("limma")
library('u133x3pcdf')
library('u133x3p.db')
library('affy')
library('oligo')
library('utils')
library('tidyverse')
library("readxl")
library('sva')

# Carregando os dados clinicos
clin_GSE15459 <- read.csv('Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE15459.csv', row.names = 1)
clin_GSE29272 <- read.csv('Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE29272.csv', row.names = 1)
clin_GSE38749 <- read.csv('Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE38749.csv', row.names = 1)
clin_GSE57303 <- read.csv('Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE57303.csv', row.names = 1)
clin_GSE62254 <- read.csv('Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE62254.csv', row.names = 1)


# _________________________ Organizando as tabelas clinicas _________________________


# __________________________________________________________ GSE15459 __________________________________________________________ #


clin_GSE15459 <- clin_GSE15459 |> 
                filter(data_processing.2 != "THIS SAMPLE FAILED QUALITY CONTROL AND WAS EXCLUDED FROM THE PUBLISHED ANALYSIS") |>
                filter(data_processing.2 != "THIS SAMPLE FOUND NOT TO BE GASTRIC ADENOCARCINOMA AND WAS EXCLUDED FROM THE PUBLISHED ANALYSIS") |> 
                 select(titulo = title,
                        codigo_acesso = geo_accession,
                        tipo_amostra = tissue.ch1,
                        estudo,
                        plataforma = platform_id) |>
                 mutate(tipo_amostra = ifelse(tipo_amostra == "gastric tumuor",
                                              "gastric tumor tissue", 
                                              tipo_amostra))


# Lendo tabela com informações externas
GSE15459_externo <- read_excel("Dados clinicos/externo/table_GSE15459.xlsx", sheet = 1)
GSE15459_externo <- GSE15459_externo |> rename(codigo_acesso = `GSM ID`)

# Unindo as tabelas
clin_GSE154591_nova <- left_join(clin_GSE15459, 
                                 GSE15459_externo, 
                                by = "codigo_acesso")

# Selecionando as colunas
clin_GSE154591_nova <- clin_GSE154591_nova |> 
                      select(-ID,
                             -`Expression CEL file`,
                             -Subtype) |> 
                      rename(idade = Age_at_surgery,
                             genero = Gender,
                             lauren = Laurenclassification,
                             estadiamento = Stage,
                             os_status = `Outcome (1=dead)`,
                             os = `Overall.Survival (Months)**`) |> 
                      mutate(across(where(is.character), ~ na_if(., "NA"))) |> # De string NA para NA de fato.  
                      mutate(os = round(os)) |> 
                      mutate(estadiamento = as.character(estadiamento)) |> 
                      mutate(os_status = as.character(os_status))
  

# Voltando para a tabela
clin_GSE15459 <- clin_GSE154591_nova

# __________________________________________________________ GSE29272 __________________________________________________________ #                       

clin_GSE29272 <- clin_GSE29272 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         tipo_amostra = source_name_ch1,
         estudo,
         plataforma = platform_id) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "adjacent normal tissue", 
                                "normal surrounding", 
                                tipo_amostra)) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "tumor tissue", 
                                "gastric tumor tissue", 
                                tipo_amostra))


clin_GSE29272$titulo <- gsub(".*\\b(TY\\w+)\\b.*", "\\1", clin_GSE29272$titulo)

# Lendo tabela com informações externas
GSE29272_externo <- read_csv("Dados clinicos/externo/table_GSE29272.csv")
GSE29272_externo <- GSE29272_externo[-1,]
GSE29272_externo <- GSE29272_externo |> rename(titulo = ID)
GSE29272_externo$titulo <- paste0(GSE29272_externo$titulo, "T")

# Unindo as tabelas
clin_GSE29272_nova <- left_join(clin_GSE29272, 
                                GSE29272_externo, 
                                 by = "titulo")
# Selecionando as colunas
clin_GSE29272_nova <- clin_GSE29272_nova |> 
  select(-No,
         -`Cell type`,
         -`Tumor Grade`,
         -`LN metastases`,
         -`FH of UGI cancer`,
         -`Cause of death`) |> 
  rename(lauren = `Tumor subtype`,
         idade = Age,
         genero = Sex,
         estadiamento = `Tumor Stage`,
         os_status = `Surivval during follow up`,
         os = `Survival (days)`) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) |> # De string NA para NA de fato.  
  mutate(os = ifelse(os == "Unknown", NA, os)) |>
  mutate(os = as.numeric(os)) |>
  mutate(os = ifelse(is.na(os), NA, os / 30.44)) |>
  mutate(os = ifelse(!is.na(os), round(os), NA)) |>
  mutate(estadiamento = as.character(estadiamento)) |> 
  mutate(os_status = as.character(os_status))

# Voltando para a tabela
clin_GSE29272 <- clin_GSE29272_nova


# __________________________________________________________ GSE38749 __________________________________________________________ # 

clin_GSE38749 <- clin_GSE38749 |> 
  select(codigo_acesso = geo_accession,
         tipo_amostra = source_name_ch1,
         idade = age..y..ch1,
         estadiamento = tumor.stage..ajcc..ch1,
         genero = gender.ch1,
         lauren = histological.type.ch1,
         estudo,
         plataforma = platform_id,
         os = time..months..overall.survival.ch1,
         os_status = status.ch1) |>
  mutate(across(where(is.character), ~ na_if(., "NA"))) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "primary tumor", 
                                "gastric tumor tissue", 
                                tipo_amostra)) |>
  mutate(estadiamento = as.character(estadiamento)) |> 
  mutate(os_status = as.character(os_status)) |> 
  mutate(os = ifelse(!is.na(os), round(os), NA))


# __________________________________________________________ GSE57303 __________________________________________________________ # 

clin_GSE57303 <- clin_GSE57303 |> 
                 select(titulo = title,
                        codigo_acesso = geo_accession,
                        tipo_amostra = phenotype.ch1,
                        estudo,
                        plataforma = platform_id) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "gastric cancer", 
                                "gastric tumor tissue", 
                                tipo_amostra))


# __________________________________________________________ GSE62254 __________________________________________________________ # 


clin_GSE62254 <- clin_GSE62254 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         estudo,
         plataforma = platform_id) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) |>
  mutate(tipo_amostra = if_else(tipo_amostra == "Gastric tumor", 
                                "gastric tumor tissue", 
                                tipo_amostra))
         
         
# Lendo tabela com informações externas
GSE62254_externo <- read_excel("Dados clinicos/externo/table_GSE62254.xlsx", sheet = 1)
GSE62254_externo <- GSE62254_externo |> rename(codigo_acesso = GEO_ID)

# Unindo as tabelas
clin_GSE62254_nova <- left_join(clin_GSE62254, 
                                 GSE62254_externo, 
                                 by = "codigo_acesso")

# Selecionando as colunas
clin_GSE62254_nova <- clin_GSE62254_nova |> 
  select(-`Sample\r\nName`,
         -`SCRI No.`,
         -Subgroup,
         -`MLH1 IHC`,
         -ACRG.sub,
         -Pathology,
         -EBV_ish,
         -WHO_simple,
         -`tumor type`,
         -GEO,
         -pStage,
         -`# of positive node (+)`,
         -`stage(TNM)`,
         -pnode) |> 
  rename(rfs = DFS.m,
         rfs_status = Recur,
         estadiamento_t = `T`,
         lauren = Lauren,
         os_status = Death,
         os = OS.m,
         idade = age,
         genero = sex,
         estadiamento = Stage,
         estadiamento_n = N,
         estadiamento_m = M,
         topologia = `revised location`) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) |> # De string NA para NA de fato.  
  mutate(os = round(os)) |> 
  mutate(rfs = round(rfs)) |> 
  mutate(estadiamento = as.character(estadiamento)) |>
  mutate(estadiamento_t = as.character(estadiamento_t)) |> 
  mutate(estadiamento_n = as.character(estadiamento_n)) |> 
  mutate(estadiamento_m = as.character(estadiamento_m)) |> 
  mutate(os_status = as.character(os_status)) |> 
  mutate(rfs_status = as.character(rfs_status))

# Voltando para a tabela
clin_GSE62254 <- clin_GSE62254_nova

# Salvando as novas tabelas separadas
write.csv(clin_GSE15459, 'Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE15459.csv')
write.csv(clin_GSE29272, 'Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE29272.csv')
write.csv(clin_GSE38749, 'Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE38749.csv')
write.csv(clin_GSE57303, 'Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE57303.csv')
write.csv(clin_GSE62254, 'Dados clinicos/Dados separados/Affymetrix/clin_affy_GSE62254.csv')


# Juntando os dados
clin_conjunto <- bind_rows(clin_GSE15459, 
                           clin_GSE29272, 
                           clin_GSE38749, 
                           clin_GSE57303, 
                           clin_GSE62254)

# Indicando o nome das linhas
rownames(clin_conjunto) <- clin_conjunto$codigo_acesso

# Criando coerencia nas colunas
clin_conjunto <- clin_conjunto |> 
                mutate(os_status = ifelse(os_status == "Unknown", NA, os_status)) |> 
                mutate(genero = ifelse(genero %in% c("female", "Female"), "F", 
                                        ifelse(genero %in% c("male", "Male"), "M", genero))) |> 
                 mutate(topologia = ifelse(topologia %in% c("Antrum", "antrum to body"), "antrum", 
                                           ifelse(topologia == "Body", "body",
                                                  ifelse(topologia == "Cardia", "cardia",
                                                         ifelse(topologia == "Fundus", "fundus", topologia))))) |> 
                mutate(lauren = ifelse(lauren %in% c("Diffuse", "diffuse type", "PD", "PD (signet ring cell)"), "diffuse", 
                                       ifelse(lauren %in% c("Intestinal", "intestinal type", "WD"), "intestinal",
                                              ifelse(lauren %in% c("Mixed", "mixed type", "WD (+PD)"), "mixed", lauren)))) |> 
                mutate(estadiamento = ifelse(estadiamento %in% c("1", "1B", "IB"), "I",
                                             ifelse(estadiamento == "2", "II",
                                                    ifelse(estadiamento %in% c("3", "3A", "3B", "IIIa", "IIIb"), "III", 
                                                                               ifelse(estadiamento == "4", "IV", estadiamento))))) 
  

clin_conjunto <- clin_conjunto |> 
  mutate(tipo_amostra = ifelse(tipo_amostra == "gastric tumor",
                               "gastric tumor tissue", 
                               tipo_amostra)) |> 
  mutate(os = round(os))




# Salvando os dados clinicos de todas as amostras
write.csv(clin_conjunto, "Dados clinicos/Dados conjuntos/Affymetrix/clin_affy_conjunto.csv")

# _______________________________________ Dados de expressão ________________________________  #

# Carregando os dados de expressão
exp_GSE15459 <- read.csv('Dados normalizados/Dados separados/Affymetrix/expr_affy_GSE15459.csv', row.names = 1)
colnames(exp_GSE15459) <- gsub("\\..*$", "", colnames(exp_GSE15459))
exp_GSE29272 <- read.csv('Dados normalizados/Dados separados/Affymetrix/expr_affy_GSE29272.csv', row.names = 1)
exp_GSE38749 <- read.csv('Dados normalizados/Dados separados/Affymetrix/expr_affy_GSE38749.csv', row.names = 1)
exp_GSE57303 <- read.csv('Dados normalizados/Dados separados/Affymetrix/expr_affy_GSE57303.csv', row.names = 1)
exp_GSE62254 <- read.csv('Dados normalizados/Dados separados/Affymetrix/expr_affy_GSE62254.csv', row.names = 1)


# Transformando rownames em colunas
exp_GSE15459 <- exp_GSE15459 %>% tibble::rownames_to_column("rownames")
exp_GSE29272 <- exp_GSE29272 %>% tibble::rownames_to_column("rownames")
exp_GSE38749 <- exp_GSE38749 %>% tibble::rownames_to_column("rownames")
exp_GSE57303 <- exp_GSE57303 %>% tibble::rownames_to_column("rownames")
exp_GSE62254 <- exp_GSE62254 %>% tibble::rownames_to_column("rownames")


# Identificar os genes comuns
genes_comuns <- Reduce(intersect, list(exp_GSE15459$rownames, 
                                       exp_GSE29272$rownames, 
                                       exp_GSE38749$rownames,
                                       exp_GSE57303$rownames,
                                       exp_GSE62254$rownames))

# Filtrar as tabelas para manter apenas os gees comuns
exp_GSE15459 <- exp_GSE15459 %>% filter(rownames %in% genes_comuns)
exp_GSE29272 <- exp_GSE29272 %>% filter(rownames %in% genes_comuns)
exp_GSE38749 <- exp_GSE38749 %>% filter(rownames %in% genes_comuns)
exp_GSE57303 <- exp_GSE57303 %>% filter(rownames %in% genes_comuns)
exp_GSE62254 <- exp_GSE62254 %>% filter(rownames %in% genes_comuns)


# Fazendo os joins
exp_todos <- exp_GSE15459 %>%
             inner_join(exp_GSE29272, by = "rownames") %>%
             inner_join(exp_GSE38749, by = "rownames") %>%
             inner_join(exp_GSE57303, by = "rownames") %>%
             inner_join(exp_GSE62254, by = "rownames") 

# Convertendo rownames de volta para rownames da tabela
exp_todos <- exp_todos %>% tibble::column_to_rownames("rownames")
exp_todos <- exp_todos[,clin_conjunto$codigo_acesso]

# Verificando se os dados clinicos e de experimento estão na mesma ordem
identical(colnames(exp_todos), clin_conjunto$codigo_acesso)

# Salvando os dados de expressão de todas as amostras
write.csv(exp_todos, "Dados normalizados/Dados conjuntos/Affymetrix/expr_affy_conjunto.csv")
