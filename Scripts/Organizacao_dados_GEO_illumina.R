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
clin_GSE13861 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE13861.csv', row.names = 1)
clin_GSE26899 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE26899.csv', row.names = 1)
clin_GSE26901 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE26901.csv', row.names = 1)
clin_GSE28541 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE28541.csv', row.names = 1)
clin_GSE26253 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE26253.csv', row.names = 1)
clin_GSE84426 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE84426.csv', row.names = 1)
clin_GSE84433 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE84433.csv', row.names = 1)
clin_GSE147163 <- read.csv('Dados clinicos/Dados separados/Illumina/clin_il_GSE147163.csv', row.names = 1)

# _________________________ Organizando as tabelas clinicas _________________________


# _____________ GSE13861 _______________ #

clin_GSE13861 <- clin_GSE13861 |> 
                 select(titulo = title,
                        codigo_acesso = geo_accession,
                        tipo_amostra = characteristics_ch1,
                        estudo,
                        plataforma) |> 
                mutate(tipo_amostra = if_else(tipo_amostra == "normal surrounding gastric tissue", 
                                              "normal surrounding", 
                                              tipo_amostra)) |>
                mutate(tipo_amostra = if_else(tipo_amostra == "gastric adenocarcinoma", 
                                              "gastric tumor tissue", 
                                              tipo_amostra)) |> 
                filter(tipo_amostra != "GIST")

# Lendo tabela com informações externas
GSE13861_externo <- read_excel("Dados clinicos/externo/info_cg_illumina.xlsx", sheet = 1)
GSE13861_externo <- GSE13861_externo |> rename(titulo = Patients_ID)

# Obtendo apenas o código na tabela clin_GSE13861
clin_GSE13861$titulo <- gsub(".*(YG[0-9]+[A-Z]?).*", "\\1", clin_GSE13861$titulo)

# Unindo as tabelas
clin_GSE13861_nova <- left_join(clin_GSE13861, 
                                GSE13861_externo, 
                                by = "titulo")

# Selecionando as colunas
clin_GSE13861_nova <- clin_GSE13861_nova |> 
                      select(-Subgroup) |> 
                      rename(genero = Sex,
                             idade = Age,
                             topologia = Location,
                             lauren = Lauren,
                             estadiamento = AJCC6,
                             estadiamento_m = M.stage,
                             os_status = Death,
                             os = OS.m,
                             rfs_status = Recurrenc,
                             rfs = RFS.m,
                             quimio = Adjuvant.chemo) |> 
                      mutate(across(where(is.character), ~ na_if(., "NA"))) # De string NA para NA de fato.

# Voltando para a tabela
clin_GSE13861 <- clin_GSE13861_nova

# _____________ GSE26899 _______________ #                       

clin_GSE26899 <- clin_GSE26899 |> 
  select(titulo = patient.ch1,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         idade = age.ch1,
         estadiamento = ajcc.stage.ch1,
         genero = gender.ch1,
         lauren = lauren.classification.ch1,
         topologia = location.ch1,
         quimio = adjuvant.chemotherapy..1.yes..0.no..na.not.available..ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "Gastric Surrounding normal tissue", 
                                "normal surrounding", 
                                tipo_amostra)) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "Gastric tumor tissue", 
                                "gastric tumor tissue", 
                                tipo_amostra))

clin_GSE26899$estadiamento <- as.character(clin_GSE26899$estadiamento)
clin_GSE26899$quimio <- as.character(clin_GSE26899$quimio)

# _____________ GSE26901 _______________ # 

clin_GSE26901 <- clin_GSE26901 |> 
  select(titulo = patient.ch1,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         idade = age.ch1,
         estadiamento = ajcc.stage.ch1,
         genero = gender.ch1,
         lauren = lauren.classification.ch1,
         topologia = location.ch1,
         quimio = adjuvant.chemotherapy..1.yes..0.no..na.not.available..ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "Gastric cancer tissue", 
                                "gastric tumor tissue", 
                                tipo_amostra)) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "Gastric tumor tissue", 
                                "gastric tumor tissue", 
                                tipo_amostra))

# Lendo tabela com informações externas
GSE26901_externo <- read_excel("Dados clinicos/externo/info_cg_illumina.xlsx", sheet = 4)
GSE26901_externo <- GSE26901_externo |> rename(titulo = Patients_ID)

# Configurando a tabela
clin_GSE26901$titulo <- gsub("_", "", clin_GSE26901$titulo)
clin_GSE26901$titulo <- gsub("T", "", clin_GSE26901$titulo)

# Unindo as tabelas
clin_GSE26901_nova <- left_join(clin_GSE26901, 
                                GSE26901_externo, 
                                by = "titulo")

# Selecionando as colunas
clin_GSE26901_nova <- clin_GSE26901_nova |> 
  select(-Subgroup,
         -Sex,
         -Age,
         -Location,
         -Lauren,
         -AJCC.stage,
         -Adjuvant.chemo) |>
  rename(estadiamento_m = M.stage,
         os_status = Death,
         os = OS.m,
         rfs_status = Recurrence,
         rfs = RFS.m) |> 
  mutate(os = round(os)) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) # De string NA para NA de fato.

# Voltando para a tabela
clin_GSE26901 <- clin_GSE26901_nova
clin_GSE26901$quimio <- as.character(clin_GSE26901$quimio)
clin_GSE26901$estadiamento <- as.character(clin_GSE26901$estadiamento)

# _____________ GSE28541 _______________ # 

clin_GSE28541 <- clin_GSE28541 |> 
  select(titulo = sample_id.ch1,
         codigo_acesso = geo_accession,
         tipo_amostra = source_name_ch1,
         estadiamento = baseline.stage.ch1,
         genero = gender.ch1,
         estudo,
         plataforma)

# Lendo tabela com informações externas
GSE28541_externo <- read_excel("Dados clinicos/externo/info_cg_illumina.xlsx", sheet = 3)
GSE28541_externo <- GSE28541_externo |> rename(titulo = `GEO ID`)

# Unindo as tabelas
clin_GSE28541_nova <- left_join(clin_GSE28541, 
                                GSE28541_externo, 
                                by = "titulo")

# Selecionando as colunas
clin_GSE28541_nova <- clin_GSE28541_nova |> 
  select(-Subgroup,
         -Sex,
         -`Array id`,
         -Stage) |> 
  rename(idade = Age,
         os_status = Deat,
         os = OS.m,
         quimio = Chemotherapy) |> 
  mutate(os = round(os)) |> 
  mutate(across(where(is.character), ~ na_if(., "NA"))) # De string NA para NA de fato.

# Voltando para a tabela
clin_GSE28541 <- clin_GSE28541_nova

# _____________ GSE26253 _______________ # 

clin_GSE26253 <- clin_GSE26253 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         estadiamento = pathological.stage.ch1,
         rfs = recurrence.free.survival.time..month..ch1,
         status_rfs = status..0.non.recurrence..1.recurrence..ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = "gastric tumor tissue")


# _____________ GSE84426 _______________ # 

clin_GSE84426 <- clin_GSE84426 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         idade = age.ch1,
         estadiamento_t = ptstage.ch1,
         estadiamento_n = pnstage.ch1,
         os = duration.overall.survival.ch1,
         genero = Sex.ch1,
         status_os = death.ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "gastric cancer ", 
                                  "gastric tumor tissue", 
                                  tipo_amostra))

# _____________ GSE84433 _______________ # 

clin_GSE84433 <- clin_GSE84433 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         idade = age.ch1,
         estadiamento_t = ptstage.ch1,
         estadiamento_n = pnstage.ch1,
         os = duration.overall.survival.ch1,
         genero = Sex.ch1,
         status_os = death.ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "gastric cancer ", 
                                "gastric tumor tissue", 
                                tipo_amostra))

# _____________ GSE147163 _______________ # 

clin_GSE147163 <- clin_GSE147163 |> 
  select(titulo = title,
         codigo_acesso = geo_accession,
         tipo_amostra = tissue.ch1,
         estudo,
         plataforma) |> 
  mutate(tipo_amostra = if_else(tipo_amostra == "gastric cancer ", 
                                "gastric tumor tissue", 
                                tipo_amostra))

# Salvando as novas tabelas separadas
write.csv(clin_GSE13861, 'Dados clinicos/clin_il_GSE13861.csv')
write.csv(clin_GSE26899, 'Dados clinicos/clin_il_GSE26899.csv')
write.csv(clin_GSE26901, 'Dados clinicos/clin_il_GSE26901.csv')
write.csv(clin_GSE28541, 'Dados clinicos/clin_il_GSE28541.csv')
write.csv(clin_GSE26253, 'Dados clinicos/clin_il_GSE26253.csv')
write.csv(clin_GSE84426, 'Dados clinicos/clin_il_GSE84426.csv')
write.csv(clin_GSE84433, 'Dados clinicos/clin_il_GSE84433.csv')
write.csv(clin_GSE147163, 'Dados clinicos/clin_il_GSE147163.csv')

# Juntando os dados
clin_conjunto <- bind_rows(clin_GSE13861, 
                           clin_GSE26899, 
                           clin_GSE26901, 
                           clin_GSE28541, 
                           clin_GSE26253, 
                           clin_GSE84426, 
                           clin_GSE84433, 
                           clin_GSE147163)

# Indicando o nome das linhas
rownames(clin_conjunto) <- clin_conjunto$codigo_acesso

# Retirando coluna indesejada
clin_conjunto <- clin_conjunto |> 
                 select(-RadiationTherapy)

# Criando coerencia nas colunas
clin_conjunto <- clin_conjunto |> 
                 mutate(genero = ifelse(genero %in% c("female", "Female"), "F", 
                                        ifelse(genero %in% c("male", "Male"), "M", genero))) |> 
                 mutate(topologia = ifelse(topologia == "Antrum", "antrum", 
                                           ifelse(topologia == "Body", "body",
                                                  ifelse(topologia == "Cardia", "cardia",
                                                         ifelse(topologia == "Fundus", "fundus", topologia))))) |> 
                mutate(lauren = ifelse(lauren == "Diffuse", "diffuse", 
                                       ifelse(lauren == "Intestinal", "intestinal",
                                              ifelse(lauren == "Mixed", "mixed", lauren)))) |> 
                mutate(estadiamento = ifelse(estadiamento %in% c("1", "1B", "IB"), "I",
                                             ifelse(estadiamento == "2", "II",
                                                    ifelse(estadiamento %in% c("3", "3A", "3B", "IIIA", "IIIB"), "III", 
                                                                               ifelse(estadiamento == "4", "IV", estadiamento))))) |> 
                mutate(quimio = ifelse(quimio == "Yes", "1", 
                                       ifelse(quimio == "No", "0", quimio)))

# Salvando os dados clinicos de todas as amostras
write.csv(clin_conjunto, "Dados clinicos/Dados conjuntos/Illumina/clin_il_conjunto.csv")

# _______________________________________ Dados de expressão ________________________________  #

# Carregando os dados de expressão
exp_GSE13861 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE13861.csv', row.names = 1)
exp_GSE26899 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE26899.csv', row.names = 1)
exp_GSE26901 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE26901.csv', row.names = 1)
exp_GSE28541 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE28541.csv', row.names = 1)
exp_GSE26253 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE26253.csv', row.names = 1)
exp_GSE84426 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE84426.csv', row.names = 1)
exp_GSE84433 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE84433.csv', row.names = 1)
exp_GSE147163 <- read.csv('Dados normalizados/Dados separados/Illumina/expr_il_GSE147163.csv', row.names = 1)

# Transformando rownames em colunas
exp_GSE13861 <- exp_GSE13861 %>% tibble::rownames_to_column("rownames")
exp_GSE26899 <- exp_GSE26899 %>% tibble::rownames_to_column("rownames")
exp_GSE26901 <- exp_GSE26901 %>% tibble::rownames_to_column("rownames")
exp_GSE28541 <- exp_GSE28541 %>% tibble::rownames_to_column("rownames")
exp_GSE26253 <- exp_GSE26253 %>% tibble::rownames_to_column("rownames")
exp_GSE84426 <- exp_GSE84426 %>% tibble::rownames_to_column("rownames")
exp_GSE84433 <- exp_GSE84433 %>% tibble::rownames_to_column("rownames")
exp_GSE147163 <- exp_GSE147163 %>% tibble::rownames_to_column("rownames")

# Identificar os genes comuns
genes_comuns <- Reduce(intersect, list(exp_GSE13861$rownames, 
                                          exp_GSE26899$rownames, 
                                          exp_GSE26901$rownames,
                                          exp_GSE28541$rownames,
                                          exp_GSE26253$rownames,
                                          exp_GSE84426$rownames,
                                          exp_GSE84433$rownames,
                                          exp_GSE147163$rownames))

# Filtrar as tabelas para manter apenas os gees comuns
exp_GSE13861 <- exp_GSE13861 %>% filter(rownames %in% genes_comuns)
exp_GSE26899 <- exp_GSE26899 %>% filter(rownames %in% genes_comuns)
exp_GSE26901 <- exp_GSE26901 %>% filter(rownames %in% genes_comuns)
exp_GSE28541 <- exp_GSE28541 %>% filter(rownames %in% genes_comuns)
exp_GSE26253 <- exp_GSE26253 %>% filter(rownames %in% genes_comuns)
exp_GSE84426 <- exp_GSE84426 %>% filter(rownames %in% genes_comuns)
exp_GSE84433 <- exp_GSE84433 %>% filter(rownames %in% genes_comuns)
exp_GSE147163 <- exp_GSE147163 %>% filter(rownames %in% genes_comuns)


# Fazendo os joins
exp_todos <- exp_GSE13861 %>%
             inner_join(exp_GSE26899, by = "rownames") %>%
             inner_join(exp_GSE26901, by = "rownames") %>%
             inner_join(exp_GSE28541, by = "rownames") %>%
             inner_join(exp_GSE26253, by = "rownames") %>%
             inner_join(exp_GSE84426, by = "rownames") %>%
             inner_join(exp_GSE84433, by = "rownames") %>%
             inner_join(exp_GSE147163, by = "rownames")
  

# Convertendo rownames de volta para rownames da tabela
exp_todos <- exp_todos %>% tibble::column_to_rownames("rownames")
exp_todos <- exp_todos[,clin_conjunto$codigo_acesso]

# Verificando se os dados clinicos e de experimento estão na mesma ordem
identical(colnames(exp_todos), clin_conjunto$codigo_acesso)

# Salvando os dados de expressão de todas as amostras
write.csv(exp_todos, "Dados normalizados/Dados conjuntos/Illumina/expr_il_conjunto.csv")
