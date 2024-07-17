library(readxl)

#EDIT GRACIELA
setwd("/home/graciela/tblab/graciela/metadata_CES") #tblab
setwd("/home/graciela/Desktop/metadata_CES") #margarita
setwd("/Users/gracielauria/Desktop/metadata_CES") #mac
library(utils)
library(R.utils)
library(vcfR)
library(stats)
library(dplyr)
library(readxl)



read_excel_allsheets <- function(filename, tibble = FALSE, nskip=1) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) data.frame(readxl::read_excel(filename, sheet = X, skip=nskip), stringsAsFactors = F))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}



extract_df <- function(filename){
  
  df_out=data.frame()
  for (i in c(0,1,2)){
    mysheets <- read_excel_allsheets(filename, nskip=i)
    familia.list = lapply(mysheets, function(x)
      if("ADN" %in% colnames(x) && "Familia" %in% colnames(x)) x[,c("ADN","Familia")])
    familia.df = do.call("rbind", familia.list)
    if (is.null(familia.df)) next
    familia.df$Proyecto = rownames(familia.df)
    familia.df$Proyecto = gsub(".[0-9]*$", "", familia.df$Proyecto, perl = TRUE)
    df_out = rbind(df_out, familia.df)
  }
  
  return(df_out)
}

#GUR: sacar tambien la Patología
extract_df_bis <- function(filename){
  
  df_out=data.frame()
  for (i in c(0,1,2)){
    mysheets <- read_excel_allsheets(filename, nskip=i)
    familia.list = lapply(mysheets, function(x)
      if("ADN" %in% colnames(x) && "Familia" %in% colnames(x) && "Patología" %in% colnames(x)) x[,c("ADN","Familia","Patología")])
    familia.df = do.call("rbind", familia.list)
    if (is.null(familia.df)) next
    familia.df$Proyecto = rownames(familia.df)
    familia.df$Proyecto = gsub(".[0-9]*$", "", familia.df$Proyecto, perl = TRUE)
    df_out = rbind(df_out, familia.df)
  }
  
  return(df_out)
}


# Function to normalize column names by removing whitespace and special characters
normalize_colnames <- function(colnames) {
  colnames <- gsub("\\s+", "", colnames)  # Remove whitespace
  colnames <- gsub("[[:punct:]]", "", colnames)  # Remove special characters
  colnames <- gsub("á", "a", colnames)  # Replace accented characters
  colnames <- gsub("é", "e", colnames)
  colnames <- gsub("í", "i", colnames)
  colnames <- gsub("ó", "o", colnames)
  colnames <- gsub("ú", "u", colnames)
  return(tolower(colnames))  # Convert to lower case for uniformity
}

# Function to extract specific columns from the Excel sheets
extract_df_gene <- function(filename) {
  
  df_out <- data.frame()
  
  for (i in c(0, 1, 2)) {
    mysheets <- read_excel_allsheets(filename, nskip = i)
    
    familia.list <- lapply(mysheets, function(x) {
      original_colnames <- colnames(x)
      normalized_colnames <- normalize_colnames(original_colnames)  # Normalize column names
      
      required_columns <- c("adn", "familia", "patologia", "genespanel")
      matched_columns <- required_columns[required_columns %in% normalized_colnames]
      
      if (all(c("adn", "familia") %in% matched_columns)) {
        selected_data <- x[, match(matched_columns, normalized_colnames)]
        colnames(selected_data) <- original_colnames[match(matched_columns, normalized_colnames)]
        return(selected_data)
      }
    })
    
    familia.df <- do.call("rbind", familia.list)
    
    if (is.null(familia.df)) next
    
    familia.df$Proyecto <- rownames(familia.df)
    familia.df$Proyecto <- gsub(".[0-9]*$", "", familia.df$Proyecto, perl = TRUE)
    df_out <- rbind(df_out, familia.df)
  }
  
  return(df_out)
}

familia.CES5.df = extract_df_gene("Exoma Clínico_CES246_actual.xls")


# DATA LOADING - load excel with the samples (excel nvestigacion) uno a uno y en diferentes argumentos.
## Añadir otro argumento si se meten mas excels


familia.TSO.df = extract_df_bis("Exoma Clínico_TSO y anterior.xls")
familia.CES1.df = extract_df_bis("Exoma Clínico_CES55_CES122.xls")
familia.CES2.df = extract_df_bis("Exoma Clínico_CES123_CES_162.xls")
familia.CES3.df = extract_df_bis("Exoma Clínico_CES163_CES211.xls")
familia.CES4.df = extract_df_bis("Exoma Clínico_CES212_CES245.xls")
familia.CES5.df = extract_df_bis("Exoma Clínico_CES246_actual.xls")
familia.otrospacientes.df<-data.table::fread("otros_pacientes.tsv",header = TRUE)


####################### AQUI METER LOS WES/WGS Y CES DE INVESTIGACION PARA QUE ESTEN EN EL NEW DF#####
familia.WES.df = data.frame(readxl::read_excel("newfinal_df_MAF_v4.0.xlsx"))  
## quitar los casos que no me sabia yo la metadata que si que tengan la metadata en CES
#familia.WES.df <- familia.WES.df %>%
  #dplyr::filter(PROBANDUS != "yes_ces")
#remove ahora la columna esa de probandus
familia.WES.df <-familia.WES.df[,-c(4,5,7)] #quitar columna de probandus y seq type
setwd("/home/graciela/tblab/graciela/metadata_CES") #tblab
#### ahora el CES de investigacion
familia.CES.INV.df = data.frame(readxl::read_excel("Lista_casos_CESinv_y_DUDAS_CrisR.xlsx", sheet="CES_INV"))  #GUR, original la de arriba
familia.CES.INV.df<-familia.CES.INV.df[,c(1,3:5)]
#############################################################################

# EDITAR
# aqui meter el nuevo, e.g. familia.CES4.df
#GONZALO: familia = rbind(familia.TSO.df, familia.CES1.df, familia.CES2.df, familia.CES3.df, familia.CES4.df, familia.CES5.df,  otras.muestras, familia.CES.INV.df)

familia = rbind(familia.TSO.df, familia.CES1.df, familia.CES2.df, familia.CES3.df, 
                familia.CES4.df, familia.CES5.df, familia.CES.INV.df,familia.WES.df,
                familia.otrospacientes.df)

#familia = rbind(familia.TSO.df, familia.CES1.df, familia.CES2.df, familia.CES3.df, familia.CES4.df, familia.CES5.df) #sin excel investigacion
familia$Familia = gsub("\n", "|", familia$Familia)
familia = familia[!duplicated(familia),] # remove duplicates

# TAG ASSOCIATION
familia$ADN = gsub(" .*", "", familia$ADN, perl = TRUE)
familia$SAMPLE = gsub("/", "-", familia$ADN)
familia$TAG = gsub("[-_ ].*", "", familia$Familia, perl = TRUE)
# familia = familia[!is.na(familia$Familia),]
familia$tag = tolower(familia$TAG)
#tag.assoc = data.frame(read_excel(args[5]), stringsAsFactors = F, sheet="clasificadas")
#tag.assoc = data.frame(read_excel("TAG_GROUP2.xlsx"), stringsAsFactors = F, sheet="clasificadas") 

######################### GUR 1/05/2024
familia$Patología= tolower(familia$Patología) #todo minusculas en enfermedad
familia$Patología = stringi::stri_trans_general (familia$Patología, "Latin-ASCII") #quitar las tildes
familia <- familia %>%
  mutate(Patología = gsub("_", " ", Patología)) #cambiar _ por espacio

#### PREPARAR DATA FRAME FAMILIA #####
familia = familia[grep("REANÁLISIS", familia$Proyecto, invert = TRUE),] #filtra reanalisis (quita pacientes)
familia = familia[grep("Confirmaciones", familia$Proyecto, invert = TRUE),] #filtra confirmaciones
familia = familia[!is.na(familia$ADN),] #quita los que no tienen ADN
familia$Familia[is.na(familia$Familia)] <- "-"
familia$TAG[is.na(familia$TAG)] <- "-"
familia$tag[is.na(familia$Familia)] <- "-"
###pegar los datos repetidos antes (newdf de si los datos se secuenciaron 2 veces)
#### 5) quitar asteriscos de las samples que lo tienen -> COLLAPSE MÁS DUPLICATES QUE QUEDABAN DESPUES Y LOS QUE NO TIENEN DUPLICATES LIMPIAR
familia$ADN <- gsub("\\*", "", familia$ADN)
familia$SAMPLE <- gsub("\\*", "", familia$SAMPLE)

# Collapse duplicates: GONZALO -> solo se pega el proyecto 
newdf <- familia[!duplicated(familia$ADN),]
for (sample in newdf$ADN){
  for (column in colnames(newdf)){
    vector = unique(familia[familia$ADN == sample, column])
    newdf[newdf$ADN == sample, column] = paste(vector[!is.na(vector)], collapse = "|")
  }
}
newdf[newdf==""]<-NA

familia<-newdf
############## PARA LOS CASOS QUE SALIERON DE PEGAR TAG DE RP CRIS ME LOS CLASIFICO EN RP Y RP SINDROMICA: reemplazar la columna de Patología por RP-sindromica si es ese caso
dudas_RP = data.frame(readxl::read_excel("Lista_casos_CESinv_y_DUDAS_CrisR.xlsx", sheet="dudas_RP_sindromicas"))  #GUR, original la de arriba
dudas_RP <- dudas_RP %>%
  select(ADN,Etiqueta.CrisR)

familia <- familia %>%
  left_join(dudas_RP, by = "ADN") %>%
  mutate(Patología = ifelse(!is.na(Etiqueta.CrisR), Etiqueta.CrisR, Patología)) %>%
  select(-Etiqueta.CrisR) # Remove the temporary equiqueta_cris column

######## PARA LOS CASOS QUE INICIALMENTE SALIERON QUE TENIAN DOS O MAS PATOLOGIAS, REEMPLAZAR LA COLUMNA DE PATOLOGIA POR EL NOMBRE COMUN
dos_cats= data.frame(readxl::read_excel("all_CES_WES_WGS_FJD.xlsx", sheet="more_than_one_disease"))  #GUR, original la de arriba
dos_cats <- dos_cats %>%
  dplyr::filter(EDIT.PATOLOGIA!="ok") #quitar las que estaban bien
dos_cats <- dos_cats %>%
  select(ADN,EDIT.PATOLOGIA)
# Replace colons with whitespace in the 'Text' column
dos_cats <- dos_cats %>%
  mutate(EDIT.PATOLOGIA = gsub(":", " ", EDIT.PATOLOGIA))
familia <- familia %>%
  left_join(dos_cats, by = "ADN") %>%
  mutate(Patología = ifelse(!is.na(EDIT.PATOLOGIA), EDIT.PATOLOGIA, Patología)) %>%
  select(-EDIT.PATOLOGIA) # Remove the temporary equiqueta_cris column

######## PARA LOS CASOS QUE INICIALMENTE SALIERON QUE NO TENIAN MATCH, REEMPLAZAR LA COLUMNA DE PATOLOGIA POR EL NOMBRE COMUN QUE YO HE DADO A MANO
inic_no_match= data.frame(readxl::read_excel("all_CES_WES_WGS_FJD.xlsx", sheet="no_match"))  #GUR, original la de arriba
inic_no_match <- inic_no_match %>%
  select(ADN,EDIT_PATOLOGIA)
# Replace colons with whitespace in the 'Text' column
inic_no_match <- inic_no_match %>%
  mutate(EDIT_PATOLOGIA = gsub(":", " ", EDIT_PATOLOGIA))
#sustituir patologia por la info nueva que yo di
familia <- familia %>%
  left_join(inic_no_match, by = "ADN") %>%
  mutate(Patología = ifelse(!is.na(EDIT_PATOLOGIA), EDIT_PATOLOGIA, Patología)) %>%
  select(-EDIT_PATOLOGIA) # Remove the temporary equiqueta_cris column

raw_familia<-familia
familia<-raw_familia

########################### CARGAR LA INFO DE MI DICCIONARIO DE ENFERMEDADES ###############################3

#tag.assoc_ORIGINAL = data.frame(readxl::read_excel("GUR_TAG_GROUP2.xlsx", sheet="clasificadas"))  #GUR, original la de arriba
#tag.assoc = data.frame(readxl::read_excel("GUR_TAG_GROUP2.xlsx", sheet="may23"))  #GUR, original la de arriba
tag.assoc = data.frame(readxl::read_excel("all_CES_WES_WGS_FJD.xlsx", sheet="NLP_words_june12th"))  #GUR, original la de arriba
tag.assoc$tag = tolower(tag.assoc$TAG)
##quitar el . de las enfermedades que lo tengan al final: por ejemplo arritmias. = arritmias
tag.assoc$Enfermedad <- sub("\\.$", "", tag.assoc$Enfermedad)

################## 
# 
# # Replace periods with a blank space if there is no space after the period
# tag.assoc$Enfermedad <- gsub("\\.(?!\\s)", " ", tag.assoc$Enfermedad , perl = TRUE)
# 
# # Remove the period if it is followed by a space
# tag.assoc$Enfermedad <- gsub("\\.\\s", " ", tag.assoc$Enfermedad)
# 
# # Replace multiple spaces with a single space
# tag.assoc$Enfermedad  <- gsub("\\s{2,}", " ", tag.assoc$Enfermedad)

################## 

############### TODO ESTO ES POR SI NO LEO EL CLEAN EXCEL DE ENFERMEDADES ############
#si alguna palabra de la familia$patologia coinicide con alguna palabra de tag_assoc$enfermedad entonces mutate en familia el tag_asocc$super_subtipo y el tab_assoc$subtipo
#else ponle nada

#tag.assoc_unique = tag.assoc[!is.na(tag.assoc$Enfermedad),] #quitar enfermedades que son NA
# tag.assoc$Enfermedad = tolower(tag.assoc$Enfermedad) #todo minusculas en enfermedad
# tag.assoc$Enfermedad = stringi::stri_trans_general (tag.assoc$Enfermedad, "Latin-ASCII") #quitar las tildes  
tag.assoc$Enfermedad<-gsub("^\\s+|\\s+$","",tag.assoc$Enfermedad) #por si me he dejado un espacio en blanco al principio o final del string de enfermedad

tag.assoc$General = stringi::stri_trans_general (tag.assoc$General, "Latin-ASCII") #quitar las tildes   
tag.assoc$General = stringr::str_to_sentence(tag.assoc$General) #1st letter capitalized
# 
tag.assoc$Subtipo = stringi::stri_trans_general (tag.assoc$Subtipo, "Latin-ASCII") #quitar las tildes   
tag.assoc$Subtipo = stringr::str_to_sentence(tag.assoc$Subtipo) #1st letter capitalized
# 
tag.assoc$Super_Subtipo = stringi::stri_trans_general (tag.assoc$Super_Subtipo, "Latin-ASCII") #quitar las tildes   
tag.assoc$Super_Subtipo = stringr::str_to_sentence(tag.assoc$Super_Subtipo) #quitar las tildes   


########################## ORDENAR EL TAG ASSOC DEL NLP WORDS Y GUARDARLO EN UNA NUEVA PESTAÑA DEL EXCEL (NUEVAS CATEGORIAS)
# Order DataFrame by first column (Col1) alphabetically, and then by second column (Col2) alphabetically
order_tag.assoc <- tag.assoc[order(tag.assoc$General,tag.assoc$Subtipo,tag.assoc$Super_Subtipo, tag.assoc$Nombre_comun), ]
##UNIQUE DISEASES
unique_tag_general<-as.data.frame(unique(order_tag.assoc[,c("General","Subtipo")]))
unique_tag_subtipo<-as.data.frame(unique(order_tag.assoc[,c("General","Subtipo","Super_Subtipo")]))
unique_tag_super_subtipo<-as.data.frame(unique(order_tag.assoc[,c("General","Subtipo","Super_Subtipo","Nombre_comun")]))
#View(unique_tag_subtipo)
#View(unique_tag_super_subtipo)
duplicated_rows <- unique_tag_super_subtipo[duplicated(unique_tag_super_subtipo$Nombre_comun), ]

#### TRADUCCION AL INGLES ####

#ingles_cats = data.frame(readxl::read_excel("ingles_cats.xlsx"))  
ingles_cats = data.frame(readxl::read_excel("ingles_cats_3julio.xlsx"))  
tag.assoc.raw<-tag.assoc
#recuperar tag.assoc
tag.assoc<-tag.assoc.raw
# Perform left join and mutate new columns
tag.assoc <- tag.assoc %>%
  left_join(ingles_cats, by = c("General", "Subtipo", "Super_Subtipo")) %>%
  mutate(GlobalTag = ifelse(is.na(GlobalTag), 0, GlobalTag),   # Example mutation based on column X
         TypeTag = ifelse(is.na(TypeTag), 0, TypeTag),   # Example mutation based on column Y
         SubtypeTag = ifelse(is.na(SubtypeTag), 0, SubtypeTag))   # Example mutation based on column Z

tag.assoc<-tag.assoc[,c(12,14,16,4,5,6,10)]
colnames(tag.assoc)<-c("General", "Subtipo", "Super_Subtipo","Nombre_comun","Enfermedad","TAG","tag")
order_tag.assoc <- tag.assoc[order(tag.assoc$General,tag.assoc$Subtipo,tag.assoc$Super_Subtipo, tag.assoc$Nombre_comun), ]
#xlsx::write.xlsx(unique_tag_super_subtipo, "categorias_genetistas_2junio.xlsx",row.names = FALSE)
#length(unique(unique_tag_super_subtipo$Nombre_comun))
####


#compare_and_mutate_one_disease <- function(familia, tag.assoc) {
  # Initialize new columns with default value
  familia$Subtipo <- "NO MATCH"
  familia$Super_Subtipo <- "NO MATCH"
  familia$Enfermedad <- "NO MATCH"
  
  # Iterate through each row in familia
  for (i in 1:nrow(familia)) {
    # Initialize match_index
    match_index <- NULL
    # Check if any enfermedad is contained within the patologia
    for (j in 1:nrow(tag.assoc)) {
      if (grepl(paste0("\\b", tag.assoc$Enfermedad[j], "\\b"), familia$Patología[i])) {
        #if (!grepl("no", familia$Patología[i], fixed = TRUE)) {
          # If "no" is not present before the matched string -> el no este le fastidia y aumentan muchisimo los no match
          match_index <- j
          #añadir lo de la palabra "NO"
          break
        #}
      }
      # if (grepl(tag.assoc$Enfermedad[j], familia$Patología[i],fixed = TRUE)) {
      #   #tiene que machear la palabra del tab 
      #   match_index <- j
      #   break
      # }
    }
    
    # If match found, mutate corresponding rows
    if (!is.null(match_index)) {
      familia$Subtipo[i] <- tag.assoc$Subtipo[match_index]
      familia$Super_Subtipo[i] <- tag.assoc$Super_Subtipo[match_index]
      familia$Enfermedad[i] <- tag.assoc$Enfermedad[match_index]
    }
  }
  
  return(familia)
}

#compare_and_mutate_multiple_disease <- function(familia, tag.assoc) {
  # Initialize new columns with default value
  familia$Subtipo <- "NO MATCH"
  familia$Super_Subtipo <- "NO MATCH"
  familia$Enfermedad <- "NO MATCH"
  
  # Iterate through each row in familia
  for (i in 1:nrow(familia)) {
    # Initialize match_index
    match_indices <- NULL
    # Check if any enfermedad is contained within the patologia
    for (j in 1:nrow(tag.assoc)) {
      if (grepl(paste0("\\b", tag.assoc$Enfermedad[j], "\\b"), familia$Patología[i])) {
        # Check if "no" appears before the matched string
        #if (!grepl("no", familia$Patología[i], fixed = TRUE)) {
          # If "no" is not present before the matched string
          match_indices <- c(match_indices, j)
        #}
      }
    }
    
    # Initialize variables to store concatenated information
    concatenated_Super_Subtipo <- character(0)
    concatenated_Subtipo <- character(0)
    concatenated_Enfermedad <- character(0)
    seen <- character(0)
    
    # Loop through match_indices
    if (!is.null(match_indices)){
      if (length(match_indices)>1){
        for (index in match_indices) {
          # Check if the current Super_Subtipo is already in concatenated_Super_Subtipo 
          concatenated_Enfermedad <- paste(concatenated_Enfermedad, tag.assoc$Enfermedad[index], sep = " : ")
          if (!(tag.assoc$Super_Subtipo[index] %in% seen)) {
            # Concatenate information with ":" in between
            concatenated_Super_Subtipo <- paste(concatenated_Super_Subtipo, tag.assoc$Super_Subtipo[index], sep = " : ")
            concatenated_Subtipo <- paste(concatenated_Subtipo, tag.assoc$Subtipo[index], sep = " : ")
            seen <- c(seen, tag.assoc$Super_Subtipo[index])
          }
        }
        # Assign concatenated information back to familia dataframe
        familia$Super_Subtipo[i] <- concatenated_Super_Subtipo
        familia$Subtipo[i] <- concatenated_Subtipo
        familia$Enfermedad[i] <- concatenated_Enfermedad
        
      }else { #solo 1 enfermedad
        familia$Super_Subtipo[i] <- tag.assoc$Super_Subtipo[match_indices]
        familia$Subtipo[i] <- tag.assoc$Subtipo[match_indices]
        familia$Enfermedad[i] <- tag.assoc$Enfermedad[match_indices]
      }
    }

    
  }
  
  return(familia)
}

#V2 ES LA BUENA
#v2_compare_and_mutate_multiple_disease <- function(familia, tag.assoc) {
  # Initialize new columns with default value
  familia$Subtipo <- "NO MATCH"
  familia$Super_Subtipo <- "NO MATCH"
  familia$Nombre_comun <- "NO MATCH"
  familia$n_Subtipo<-0
  familia$n_Super_Subtipo<-0
  familia$n_Nombre_comun<-0
  
  
  # Iterate through each row in familia
  for (i in 1:nrow(familia)) {
    # Initialize match_index
    match_indices <- NULL
    # Check if any ENFERMEDAD is contained within the patologia -> columna ENFERMEDAD del excel
    for (j in 1:nrow(tag.assoc)) {
      if (grepl(paste0("\\b", tag.assoc$Enfermedad[j], "\\b"), familia$Patología[i])) {
        # Check if "no" appears before the matched string
        #if (!grepl("no", familia$Patología[i], fixed = TRUE)) {
        # If "no" is not present before the matched string
        match_indices <- c(match_indices, j)
        #}
      }
    }
    
    # Initialize variables to store concatenated information
    concatenated_Super_Subtipo <- character(0)
    concatenated_Subtipo <- character(0)
    concatenated_Nombre_comun<- character(0)
    seen_nombre_comun<-character(0)
    seen_super <- character(0)
    seen_sub <- character(0)
    # Loop through match_indices Y PEGAR COLUMNA NOMBRE COMUN
    if (!is.null(match_indices)){
      if (length(match_indices)>=1){
        for (index in match_indices) {
          if (!(tag.assoc$Nombre_comun[index] %in% seen_nombre_comun)) {
            # Check if the current Nombre_comun is already in concatenated_Nombre_comun
            # ojo si quiero todo con espacios en vez de sep = ":" tengo que poner sep = " : "
            concatenated_Nombre_comun <- paste(concatenated_Nombre_comun, tag.assoc$Nombre_comun[index], sep = ":")
            seen_nombre_comun<- c(seen_nombre_comun, tag.assoc$Nombre_comun[index])
            # Check if the current Super_Subtipo is already in concatenated_Super_Subtipo 
            if (!(tag.assoc$Super_Subtipo[index] %in% seen_super)) {
              # Concatenate information with ":" in between
              concatenated_Super_Subtipo <- paste(concatenated_Super_Subtipo, tag.assoc$Super_Subtipo[index], sep = ":")
              seen_super <- c(seen_super, tag.assoc$Super_Subtipo[index])
              # Check if the current Subtipo is already in concatenated_Subtipo 
              if (!(tag.assoc$Subtipo[index] %in% seen_sub)) {
                concatenated_Subtipo <- paste(concatenated_Subtipo, tag.assoc$Subtipo[index], sep = ":")
                seen_sub <- c(seen_sub, tag.assoc$Subtipo[index])
              }
            }
          }
        }
        # Assign concatenated information back to familia dataframe
        familia$Super_Subtipo[i] <- concatenated_Super_Subtipo
        familia$Subtipo[i] <- concatenated_Subtipo
        familia$Nombre_comun[i] <- concatenated_Nombre_comun
        
        #count number of "nombre comun", "subtipo", super_subtipo
        familia$n_Subtipo[i]<-stringr::str_count(concatenated_Subtipo, ":")
        familia$n_Super_Subtipo[i]<-stringr::str_count(concatenated_Super_Subtipo, ":")
        familia$n_Nombre_comun[i]<-stringr::str_count(concatenated_Nombre_comun, ":")
        
        
      }else { #solo 1 enfermedad
        # familia$Super_Subtipo[i] <- tag.assoc$Super_Subtipo[match_indices]
        # familia$Subtipo[i] <- tag.assoc$Subtipo[match_indices]
        # familia$Enfermedad[i] <- tag.assoc$Enfermedad[match_indices]
      }
    }
    
    
  }
  
  return(familia)
}


#V3 inlcluye la categoria general ES LA BUENA
v3_compare_and_mutate_multiple_disease <- function(familia, tag.assoc) {
  # Initialize new columns with default value
  familia$General <- "NO MATCH"
  familia$Subtipo <- "NO MATCH"
  familia$Super_Subtipo <- "NO MATCH"
  familia$Nombre_comun <- "NO MATCH"
  familia$n_General<-0
  familia$n_Subtipo<-0
  familia$n_Super_Subtipo<-0
  familia$n_Nombre_comun<-0
  
  
  # Iterate through each row in familia
  for (i in 1:nrow(familia)) {
    # Initialize match_index
    match_indices <- NULL
    # Check if any ENFERMEDAD is contained within the patologia -> columna ENFERMEDAD del excel
    for (j in 1:nrow(tag.assoc)) {
      if (grepl(paste0("\\b", tag.assoc$Enfermedad[j], "\\b"), familia$Patología[i])) {
        # Check if "no" appears before the matched string
        #if (!grepl("no", familia$Patología[i], fixed = TRUE)) {
        # If "no" is not present before the matched string
        match_indices <- c(match_indices, j)
        #}
      }
    }
    
    # Initialize variables to store concatenated information
    concatenated_General <- character(0)
    concatenated_Super_Subtipo <- character(0)
    concatenated_Subtipo <- character(0)
    concatenated_Nombre_comun<- character(0)
    seen_nombre_comun<-character(0)
    seen_super <- character(0)
    seen_sub <- character(0)
    seen_general <- character(0)
    # Loop through match_indices Y PEGAR COLUMNA NOMBRE COMUN
    if (!is.null(match_indices)){
      if (length(match_indices)>=1){
        for (index in match_indices) {
          if (!(tag.assoc$Nombre_comun[index] %in% seen_nombre_comun)) {
            # Check if the current Nombre_comun is already in concatenated_Nombre_comun
            # ojo si quiero todo con espacios en vez de sep = ":" tengo que poner sep = " : "
            concatenated_Nombre_comun <- paste(concatenated_Nombre_comun, tag.assoc$Nombre_comun[index], sep = ":")
            seen_nombre_comun<- c(seen_nombre_comun, tag.assoc$Nombre_comun[index])
            # Check if the current Super_Subtipo is already in concatenated_Super_Subtipo 
            if (!(tag.assoc$Super_Subtipo[index] %in% seen_super)) {
              # Concatenate information with ":" in between
              concatenated_Super_Subtipo <- paste(concatenated_Super_Subtipo, tag.assoc$Super_Subtipo[index], sep = ":")
              seen_super <- c(seen_super, tag.assoc$Super_Subtipo[index])
              # Check if the current Subtipo is already in concatenated_Subtipo 
              if (!(tag.assoc$Subtipo[index] %in% seen_sub)) {
                concatenated_Subtipo <- paste(concatenated_Subtipo, tag.assoc$Subtipo[index], sep = ":")
                seen_sub <- c(seen_sub, tag.assoc$Subtipo[index])
                if (!(tag.assoc$General[index] %in% seen_general)) {
                  concatenated_General <- paste(concatenated_General, tag.assoc$General[index], sep = ":")
                  seen_general <- c(seen_general, tag.assoc$General[index])
                }
              }
            }
          }
        }
        # Assign concatenated information back to familia dataframe
        familia$General[i] <- concatenated_General
        familia$Super_Subtipo[i] <- concatenated_Super_Subtipo
        familia$Subtipo[i] <- concatenated_Subtipo
        familia$Nombre_comun[i] <- concatenated_Nombre_comun
        
        #count number of "nombre comun", "subtipo", super_subtipo
        familia$n_General[i] <- stringr::str_count(concatenated_General, ":")
        familia$n_Subtipo[i]<-stringr::str_count(concatenated_Subtipo, ":")
        familia$n_Super_Subtipo[i]<-stringr::str_count(concatenated_Super_Subtipo, ":")
        familia$n_Nombre_comun[i]<-stringr::str_count(concatenated_Nombre_comun, ":")
        
        
      }
    }
    
    
  }
  
  return(familia)
}


# # Call the function to compare and mutate
# familia_one <- compare_and_mutate_one_disease(raw_familia, tag.assoc)
# familia<-familia_one
# 
# familia_multiple <- compare_and_mutate_multiple_disease(raw_familia, tag.assoc)
# familia<-familia_multiple
# 
# ##raw_familia o raw_newdf
# familia_multiple <- v2_compare_and_mutate_multiple_disease(raw_newdf, tag.assoc)
# familia<-familia_multiple

familia_multiple <- v3_compare_and_mutate_multiple_disease(raw_familia, tag.assoc)
familia<-familia_multiple






##################### PEGAR TAG ################
### Si no he tenido match y su patologia es exactamente alguna de estas entonces poner lo de PEGAR TAG
#amilodosis es duda porque puede ser renal o cardiaca, me puedo fijar en el tag de familia pero hay algunos que son "Varios" otro si pone NEFRO o CARDIO

info_dudosa<-c("??","¿?","ngs ces","ngs-ces","u ngs-ces","ngs-ces.","ngs - ces.","ngs","ces","ngs-ces(ver volante escaneado)",
               "validacion ns2000","ngs-ces (ver volante escaneado)","validacion hamilton",
               "ya informado","ces_ngs","nges-ces","nge-ces","nges-ces","ngs-ces","ngs-ces---","ngs-ces--","ngs-ces.",
               "repetido tso_21",
               "ngs-ces. control emqn",
               "control emqn",
               "ngs emqn",
               "ngs-cescontrol emqn",
               "ngs emqn",
               "u ngs-ces. urgente emqn",
               "exoma",
               "control positivo",
               "negativo",
               "nefro",
               "saoud")

familia$General <- ifelse(familia$General == "NO MATCH" & familia$Patología %in% info_dudosa, "PEGAR TAG", familia$General)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & familia$Patología %in% info_dudosa, "PEGAR TAG", familia$Super_Subtipo)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & familia$Patología %in% info_dudosa, "PEGAR TAG", familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & familia$Patología %in% info_dudosa, familia$Patología, familia$Nombre_comun)

#if NA patologia -> also Pegar TAG
familia$General <- ifelse(is.na(familia$Patología), "PEGAR TAG", familia$General)
familia$Super_Subtipo <- ifelse(is.na(familia$Patología), "PEGAR TAG", familia$Super_Subtipo)
familia$Subtipo <- ifelse(is.na(familia$Patología), "PEGAR TAG", familia$Subtipo)
familia$Nombre_comun <- ifelse(is.na(familia$Patología),NA, familia$Nombre_comun)


#categorias muy generales -> pegar tag -> si el string de patologia contiene: "control positivo", ->> le pego el tag
info_dudosa<-c("control positivo")
for (i in length(info_dudosa)) {
  ## Si no match pero "control_positivo" ahora aparece dentro de algun string de patologia entonces es distrofia macular -> se confunde con mody si no u otros dm de diabetes mellitus de antes, los que quedan son maculares que ya lo comprobé
  familia$General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b",info_dudosa[i], "\\b"), familia$Patología), "PEGAR TAG", familia$General)
  familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b",info_dudosa[i], "\\b"), familia$Patología), "PEGAR TAG", familia$Super_Subtipo)
  familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b",info_dudosa[i], "\\b"), familia$Patología), "PEGAR TAG", familia$Subtipo)
  familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b",info_dudosa[i], "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)
}

#### yes match nuevos

## Si no match pero "dr" ahora aparece dentro de algun string de patologia entonces es distrofia retina general -> se confunde con dr. (doctor X) por ejemplo
##poner 1 a las categorias
familia$n_General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología),1, familia$n_General)
familia$n_Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología),1, familia$n_Subtipo)
familia$n_Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología),1, familia$n_Super_Subtipo)
familia$n_Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología), 1, familia$n_Nombre_comun)

familia$General<-ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología), paste0(":",tag.assoc$General[which(tag.assoc$tag=="rp")]), familia$General)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología), paste0(":",tag.assoc$Subtipo[which(tag.assoc$tag=="rp")]), familia$Subtipo)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología), paste0(":",tag.assoc$Super_Subtipo[which(tag.assoc$tag=="rp")]), familia$Super_Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","dr", "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)

## Si no match pero "dm" ahora aparece dentro de algun string de patologia entonces es distrofia macular -> se confunde con mody si no u otros dm de diabetes mellitus de antes, los que quedan son maculares que ya lo comprobé
##poner 1 a las categorias
familia$n_General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología),1, familia$n_General)
familia$n_Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología),1, familia$n_Subtipo)
familia$n_Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología),1, familia$n_Super_Subtipo)
familia$n_Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología), 1, familia$n_Nombre_comun)

familia$General<-ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología), paste0(":",tag.assoc$General[which(tag.assoc$tag=="md")]), familia$General)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología), paste0(":",tag.assoc$Subtipo[which(tag.assoc$tag=="md")]), familia$Subtipo)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología), paste0(":",tag.assoc$Super_Subtipo[which(tag.assoc$tag=="md")]), familia$Super_Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","dm", "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)


## Si no match pero "pigmentaria" ahora aparece dentro de algun string de patologia entonces es retinosis pigmentaria -> se confunde con otras pigmentarias (de la piel)
##poner 1 a las categorias
familia$n_General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología),1, familia$n_General)
familia$n_Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología),1, familia$n_Subtipo)
familia$n_Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología),1, familia$n_Super_Subtipo)
familia$n_Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología), 1, familia$n_Nombre_comun)

familia$General<-ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología), paste0(":",tag.assoc$General[which(tag.assoc$tag=="rp")]), familia$General)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología),  paste0(":",tag.assoc$Subtipo[which(tag.assoc$tag=="rp")]), familia$Subtipo)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología),  paste0(":",tag.assoc$Super_Subtipo[which(tag.assoc$tag=="rp")]), familia$Super_Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","pigmentaria", "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)


## Si no match pero "di"ahora aparece dentro de algun string de patologia entonces es discapacidad intelectual
##poner 1 a las categorias
familia$n_General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología),1, familia$n_General)
familia$n_Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología),1, familia$n_Subtipo)
familia$n_Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología),1, familia$n_Super_Subtipo)
familia$n_Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología), 1, familia$n_Nombre_comun)

familia$General<-ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología), paste0(":",tag.assoc$General[which(tag.assoc$tag=="rm")]), familia$General)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología),  paste0(":",tag.assoc$Subtipo[which(tag.assoc$tag=="rm")]), familia$Super_Subtipo)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología),  paste0(":",tag.assoc$Super_Subtipo[which(tag.assoc$tag=="rm")]), familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","di", "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)


## Si no match pero "nf"ahora aparece dentro de algun string de patologia entonces es neurofibromatosis
##poner 1 a las categorias
familia$n_General <- ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología),1, familia$n_General)
familia$n_Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología),1, familia$n_Subtipo)
familia$n_Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología),1, familia$n_Super_Subtipo)
familia$n_Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología), 1, familia$n_Nombre_comun)

familia$General<-ifelse(familia$General == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología), paste0(":",tag.assoc$General[which(tag.assoc$tag=="neurofibromatosis")]), familia$General)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología),  paste0(":",tag.assoc$Super_Subtipo[which(tag.assoc$Enfermedad=="neurofibromatosis")]), familia$Super_Subtipo)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología),  paste0(":",tag.assoc$Subtipo[which(tag.assoc$Enfermedad=="neurofibromatosis")]), familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH" & grepl(paste0("\\b","nf", "\\b"), familia$Patología), familia$Patología, familia$Nombre_comun)


# Function to check the condition and modify the 'Tipo/subtipo' column
check_and_tag <- function(patologia,enfermedad) {
  # Check if "secuenciacion masiva simple" is present in the string
  if (grepl("secuenciacion masiva simple", patologia)) {
    # Extract the part of the string after "secuenciacion masiva simple"
    after_pattern <- sub(".*secuenciacion masiva simple", "", patologia)
    
    # Find the position of the first "-" after "secuenciacion masiva simple"
    dash_pos <- regexpr("-", after_pattern)
    
    # Check if a "-" was found
    if (dash_pos != -1) {
      # Extract everything after the first "-"
      part_after_dash <- substr(after_pattern, dash_pos + 1, nchar(after_pattern))
      # Remove leading spaces
      part_after_dash <- sub("^\\s+", "", part_after_dash)
      # Check if the first character after spaces is a letter
      if (!grepl("^[a-zA-Z]", part_after_dash)) {
        return("PEGAR TAG")
      }
    }
  }
  return(enfermedad) 
}


familia$General <- ifelse(familia$General == "NO MATCH",mapply(check_and_tag, familia$Patología, familia$General),familia$General)
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH",mapply(check_and_tag, familia$Patología, familia$Super_Subtipo),familia$Super_Subtipo)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH", mapply(check_and_tag, familia$Patología, familia$Subtipo), familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH",mapply(check_and_tag, familia$Patología, familia$Nombre_comun), familia$Nombre_comun)


ok_familia<-familia
familia<-ok_familia

########################### MAS COSAS QUE HAY QUE HACER AL DATAFRAME newdf final ############

##1) Los que sean "Control" o "Prenatal" quitarles las otras categorias (que sean unicos)
##controles son solo controles
familia <- familia %>%
  mutate(General = ifelse(grepl(":controlg", General), ":controlg", General))
familia <- familia %>%
  mutate(Subtipo = ifelse(grepl(":controlt", Subtipo), ":controlt", Subtipo))
familia <- familia %>%
  mutate(Super_Subtipo = ifelse(grepl(":scontrol", Super_Subtipo), ":scontrol", Super_Subtipo))

familia <- familia %>%
  mutate(General = ifelse(grepl(":prenatalt", Subtipo), ":sterfetalg", General)) # si el subtipo es prenatal entonces la categoria general es solo :sterfetalg
familia <- familia %>%
  mutate(Subtipo = ifelse(grepl(":prenatalt", Subtipo), ":prenatalt", Subtipo))
familia <- familia %>%
  mutate(Super_Subtipo = ifelse(grepl(":sotherprenatal", Super_Subtipo), ":sotherprenatal", Super_Subtipo))

# Replace the value in 'n_Super_Subtipo' column if 'Description' contains 'Control'
familia <- familia %>%
  mutate(n_Super_Subtipo = ifelse(grepl(":controlt|:prenatalt", Subtipo), 1, n_Super_Subtipo)) # si en subtipo detecto esto poner un 1 en ambas columnas 
familia <- familia %>%
  mutate(n_Subtipo = ifelse(grepl(":controlt|:prenatalt", Subtipo), 1, n_Subtipo))
familia <- familia %>%
  mutate(n_General = ifelse(grepl(":controlt|:prenatalt", Subtipo), 1, n_General))

##2) pegar tags para los "pegar tags" y los no "match"
###Cuando decida directamente el NO-MATCH ponerle su TAG A TODOS
familia$General <- ifelse(familia$General == "NO MATCH", "PEGAR TAG", familia$General)
familia$Subtipo <- ifelse(familia$Subtipo == "NO MATCH", "PEGAR TAG", familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$Nombre_comun == "NO MATCH", "PEGAR TAG", familia$Nombre_comun)
##tiene que ser el ultimo porque lo hago todo en base a este:
familia$Super_Subtipo <- ifelse(familia$Super_Subtipo == "NO MATCH", "PEGAR TAG", familia$Super_Subtipo)


### CAMBIAR EL PEGAR TAG por efectivamente su enfermedad y subtipo en base a la familia
indices <- which(familia$Subtipo == "PEGAR TAG")
for (i in indices) {
  # Find corresponding familia$tag that matches tag.assoc$TAG
  corresponding_tag <- familia$tag[i]
  tag_index <- which(tag.assoc$tag == corresponding_tag)

  # If a matching tag is found, substitute familia$Subtipo with tag.assoc$subtipo
  if (length(tag_index) > 0) {
    familia$General[i] <- paste0(":",tag.assoc$General[tag_index])
    familia$Subtipo[i] <- paste0(":",tag.assoc$Subtipo[tag_index])
    familia$Super_Subtipo[i] <- paste(":",tag.assoc$Super_Subtipo[tag_index])
    familia$n_General[i] <- 1
    familia$n_Subtipo[i] <- 1
    familia$n_Super_Subtipo[i] <- 1
  }
}

################ cosas de los casos rp o cohorte general de ird
##### si tenemos srp y syndroprp dejar solo srp syndromica
familia <- familia %>%
  mutate(Super_Subtipo = ifelse(grepl(":srp:ssyndromrp|:srp:ssyndromrp", Super_Subtipo), ":ssyndromrp", Super_Subtipo))


#### casos excluidos DHR de tipo "Healthy individual" ponerles que son :Control
DHR_cohorte_excluidos = data.frame(readxl::read_excel("DHR_cohorte_excluidos.xlsx"))  
DHR_healthy<-DHR_cohorte_excluidos %>%
  dplyr::filter(Exclusing_motive %in% c("Healthy_individual"))

familia$General <- ifelse(familia$ADN %in% DHR_healthy$N_DNA_BBDD, ":controlg", familia$General)
familia$Subtipo <- ifelse(familia$ADN %in% DHR_healthy$N_DNA_BBDD, ":controlt", familia$Subtipo)
familia$Nombre_comun <- ifelse(familia$ADN %in% DHR_healthy$N_DNA_BBDD, "healthy", familia$Nombre_comun)
##tiene que ser el ultimo porque lo hago todo en base a este:
familia$Super_Subtipo <- ifelse(familia$ADN %in% DHR_healthy$N_DNA_BBDD, ":scontrol", familia$Super_Subtipo)

##################### HASTA AQUI ###################

###los que no son DHR
# not_DHR<-DHR_cohorte_excluidos %>%
#   dplyr::filter(Exclusing_motive %in% c("Not_RD"))
# 
# not_DHR_fam<- familia[which(familia$ADN %in% not_DHR$N_DNA_BBDD),]


###probandus que me dicen que son DHR ponerles yo etiqueta de DHR (o distrofia macular o retinosis pigmentaria)
DHR_cohorte_probandus = data.frame(readxl::read_excel("DHR_cohorte_probandus.xlsx")) 
yes_DHR_fam<- familia[which(familia$ADN %in% DHR_cohorte_probandus$N_ADN),]


##quitar los que no tengo DHR
DHR_probandus_yo_sin_diag_DHR <- yes_DHR_fam[!grepl("irdt", yes_DHR_fam$Subtipo), ]
### a los de aqui pegarles etiqueta de ird
write.table(x = DHR_probandus_yo_sin_diag_DHR, file = "revisar_probandus_vitreoretinopatias.tsv", row.names = F, quote = F, sep = "\t")


familia$GENERAL<- ifelse(grepl("probandus_dhr", familia$Proyecto) & !grepl("eyeg", familia$CATEGORY), paste0(familia$GENERAL,":eyeg"), familia$GENERAL)






##############################---------- WRITE EL DATAFRAME AL EXCEL ---------########################


############################# PREPARAR BIEN EL DATAFRAME FINAL ######################
save<-familia
#0) cambiar unico caso iker
familia$ADN<-ifelse(familia$ADN=="3126028","0312602",familia$ADN)
familia$SAMPLE<-ifelse(familia$SAMPLE=="3126028","0312602",familia$SAMPLE)
#1) CAMBIAR NOMBRES DE CATEGORIAS IMPORTANTES
colnames(familia)[2]<-"FAMILY"
colnames(familia)[8]<-"GENERAL"
colnames(familia)[9]<-"CATEGORY"
colnames(familia)[10]<-"SUBCATEGORY"
##2) asegurarse de que no hay white spaces despues de los :
familia$GENERAL <- gsub(":\\s+", ":", familia$GENERAL)
familia$SUBCATEGORY <- gsub(":\\s+", ":", familia$SUBCATEGORY)
familia$CATEGORY <- gsub(":\\s+", ":", familia$CATEGORY)

# Remove leading and trailing whitespace from all columns ############ SUPER IMPORTANTEEEEEEE!!!!!!!!!!!!!!!!! ||||||||||||||||||||
familia[] <- lapply(familia, trimws)


### las Familias que en su columna de familia tienen varias familias: ej: ANF-0060 ONCE-1234 DC-7654 -> dejar solo la primera familia
familia <- familia %>%
  mutate(FAMILY = stringr::str_replace(FAMILY, "\\s.*", ""))

## 3) quitar los new line characters de cada field individual
# Remove newline characters from all columns, porque si no el tsv detecta 2 filas distintas en una que se supone que es la misma
familia[] <- lapply(familia, function(x) gsub("\n", " ", x))


### 4) impute rest of NA values with "-"
is_na_familia<-familia[which(is.na(familia)),]
familia[is.na(familia)] <- "-"

##5) casos que pone : en General y luego esta bien el :tumor benigno en la segunda
familia <- familia %>% 
  mutate(across(c("GENERAL"), ~ ifelse(. == ":", ":oncog", .)))


################## fin preparar bien el dataframe final ######################

# #### 4) los PEGAR TAG ponerlos en :PEGAR TAG
# familia <- familia %>% 
#   mutate(across(c("GENERAL","CATEGORY","SUBCATEGORY"), ~ ifelse(. == "PEGAR TAG", ":PEGAR TAG", .)))
# 
# 
# ####4) view and impute NA values
# ##
# df_with_na <- familia %>% filter(if_any(everything(), is.na))
# ### los - en la columna de familia cambiarlos por nombres de familia aleatorios
# 
# generate_mk_code <- function() {
#   paste0("MK-", sprintf("%04d", sample(0:9999, 1)))
# }
# familia <- familia %>%
#   mutate(across(c("FAMILY"), ~ ifelse(. == "-", generate_mk_code(), .)))
# 

###exportar como tsv para la base de datos
colnames(familia)

#exportar esto porque si no los cambios de linea de la columna de patologia no se detectan bien
tsv_familia<-familia[,c(1:2,5:10)]
colnames(tsv_familia)
df_with_na <- tsv_familia %>% filter(if_any(everything(), is.na))



##impute NA values con lo que yo quiera
#tsv_familia<-na.omit(tsv_familia) #error el hail si no
# Function to impute NA values for each column
# impute_values <- c(A = 0, B = 1, C = 6)
# impute_na <- function(df, impute_values) {
#   for (col in names(impute_values)) {
#     df[[col]][is.na(df[[col]])] <- impute_values[[col]]
#   }
#   return(df)
# }
# 
# # Apply the function to your data frame
# df_imputed <- impute_na(df, impute_values)




##all_cases
familia <- familia[order(familia$GENERAL,familia$CATEGORY,familia$SUBCATEGORY, familia$Nombre_comun), ]
##no_match
no_match<-familia%>%
  dplyr::filter(n_General==0)
no_match<-familia%>%
  dplyr::filter(CATEGORY=="PEGAR TAG")
no_match <- no_match[order(no_match$TAG,no_match$Patología), ]
View(table(no_match$TAG))
View(table(no_match$Proyecto))

#filtrar nuevos CES
# no_match_new <- no_match %>%
#   filter(stringr::str_starts(Proyecto, "CES_26"))
# writexl::write_xlsx(no_match_new,"no_match_new.xlsx")



##ok_samples
ok_samples<-familia%>%
  dplyr::filter(n_General==1)
ok_samples <- ok_samples[order(ok_samples$GENERAL,ok_samples$CATEGORY,ok_samples$SUBCATEGORY, ok_samples$Nombre_comun), ]

##more_than_one_samples
more_than_one_samples<-familia%>%
  dplyr::filter(n_General>1)
more_than_one_samples <- more_than_one_samples[order(more_than_one_samples$Proyecto,more_than_one_samples$GENERAL,more_than_one_samples$CATEGORY,more_than_one_samples$SUBCATEGORY, more_than_one_samples$n_Nombre_comun), ]
more_than_one_samples_new <- more_than_one_samples %>%
  filter(stringr::str_starts(Proyecto, "CES_26"))
writexl::write_xlsx(more_than_one_samples_new,"more_than_one_samples_new.xlsx")

named_list_df <- list(
  "all_cases" = familia,
  "no_match" = no_match,
  "ok_samples"= ok_samples,
  "more_than_one_disease"= more_than_one_samples,
  "tag_general"= unique_tag_general,
  "tag_subtipo"= unique_tag_subtipo,
  "tag_suber_subtipo"= unique_tag_super_subtipo,
  "NLP_words_june12th"=order_tag.assoc
)
writexl::write_xlsx(x=named_list_df,
                    path = "all_CES_WES_WGS_FJD.xlsx",
                    col_names=TRUE,
                    format_headers = TRUE)


write.table(x = tsv_familia, file = "~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD_9july.txt", row.names = F, quote = F, sep = "\t")
write.table(x = tsv_familia, file = "all_FJD_9july.txt", row.names = F, quote = F, sep = "\t")

#write.table(x = tsv_familia, file = "~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD.txt", row.names = F, quote = F, sep = "\t")
# View the first few lines of the file to check formatting
# Write the data.table to a tab-separated file
#data.table::fwrite(tsv_familia, "~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD.txt", sep = "\t", quote = FALSE)

head(read.table("~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD.txt", sep = "\t", header = TRUE))
df <- data.table::fread("~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD.txt")


###### descargar lista entera y comprobar los vcfs que tienen metadata
list_uam <- data.table::fread("~/tblab/graciela/metadata_CES/clean_all_FJD_UAM.txt",header = FALSE)
length(unique(list_uam$V1))
no_usar_inv = data.frame(readxl::read_excel("no_usarExomaClinico_CES_INV.xlsx"))  
ok<-list_uam[!list_uam$V1 %in% familia$SAMPLE] #not metadata available



ok$SAMPLE = gsub("-", "/", ok$V1)
extras_DHR_p1<-DHR_cohorte_probandus[which(DHR_cohorte_probandus$N_ADN %in% ok$SAMPLE),]
extras_DHR_p2<-no_usar_inv[which(no_usar_inv$...4 %in% ok$SAMPLE),]

data.table::fwrite(ok, "~/tblab/graciela/metadata_CES/not_in_metadata_all_FJD.txt", sep = "\t", quote = FALSE)



##############################---------- INSPECT EL NEW DF DIRECTAMENTE ---------########################

View(table(familia$n_General)) #grupos de enfermedades en los que esta cada ces 
View(table(familia$n_Subtipo)) #grupos de enfermedades en los que esta cada ces 
View(table(familia$n_Super_Subtipo)) #subgrupos de enfermedades en los que esta cada ces

View(table(familia$General)) 
View(table(familia$Subtipo)) 
View(table(familia$Super_Subtipo)) 

View(table(no_match$TAG))
View(table(no_match$Proyecto))


# Order DataFrame by first column (Col1) alphabetically, and then by second column (Col2) alphabetically
order_no_match <- no_match[order(no_match$TAG,no_match$Patología),]

words <- unlist(strsplit(as.character(no_match$Patología), " "))
word_freq <- View(table(words))

###################################### ni idea de que es esto #####################################################

###6) VIEW SAMPLES QUE SEAN DE LA MISMA FAMILIA cuantas personas hay de cada familia 
#-> A MANO QUITAR DE MI CARPETA QUE YA HE GUARDADO UN ARCHIVO CON LAS QUE TENGO QUE QUITAR
## QUITAR TAMBIEN LAS QUE ME HA DICHO CRIS DEL EXCEL DE INVESTIGACION
# View(table(newdf$Familia))
# df_repeated <- newdf %>%
#   group_by(Familia) %>%
#   filter(n() > 1) %>%
#   ungroup()
# df_repeated <- df_repeated[order(df_repeated$Familia,df_repeated$Proyecto),]
# 
# 
# df_repeated[is.na(df_repeated)] <- "NA"
# #xlsx::write.xlsx(df_repeated, file= "CES_varios_samples_misma_familia.xlsx") 
# #ahora solo read el excel
# 
# #xlsx::write.xlsx(df_repeated, file= "CES_varios_samples_misma_familia.xlsx") 
# 
# #### abrir archivo txt con samples repetidas:
# remove_samples_list = read.table(file="CES_SAMPLES_REMOVE_barra.txt", sep="\t", header = FALSE)
# clean_df<-newdf %>%
#   dplyr::filter(!newdf$ADN %in% remove_samples_list$V1)
# 
# clean_df_repeated <- clean_df %>%
#   group_by(Familia) %>%
#   filter(n() > 1) %>%
#   ungroup()
# clean_df_repeated <- clean_df_repeated[order(clean_df_repeated$Familia,clean_df_repeated$Proyecto),]


### EXTRA COMPARAR DATAFRAMES
no_point<-no_match
common_rows <- merge(no_match, df2)
# # Extract uncommon rows
# uncommon_rows_point <- dplyr::anti_join(point, no_point, by = c("ADN"))
# uncommon_rows_no_point <- dplyr::anti_join(no_point, point, by = c("ADN"))
# 
# uncommon_rows <- bind_rows(uncommon_rows_df1, uncommon_rows_df2)
### EXTRA COMPARAR DATAFRAMES


# newdf : ESTAN LOS SAMPLE_ids unicos JUNTOS
# familia: repetidos sample ids

familia$Categoria[is.na(familia$Categoria)] <- "Varios"
newdf$Categoria[is.na(newdf$Categoria)] <- "Varios"



write.table(x = familia, file = "mymetadatapathology_20242504.txt", row.names = F, quote = F, sep = "\t")
write.table(x = newdf, file = "mymetadatapathology_uniq_20242504.txt", row.names = F, quote = F, sep = "\t")


#write.table(x = familia, file = args[6], row.names = F, quote = F, sep = "\t")
#write.table(x = newdf, file = args[7], row.names = F, quote = F, sep = "\t")





#####################
# ionut <- read.delim("/home/gonzalo/UAMssh/fjd/MAF_FJD_v3.0/metadata/mymetadatapathology_20210315.txt", stringsAsFactors = FALSE)
# diferencia <- ionut[ionut$ADN %in% setdiff(ionut$ADN, familia$ADN), c("PROJECT", "ADN", "Familia")]
# diferencia$PROJECT = gsub("_.*$", "", diferencia$PROJECT)
# diferencia$PROJECT = gsub("-.*$", "", diferencia$PROJECT)
# colnames(diferencia) <- colnames(otros)
# write.table(diferencia, "/home/gonzalo/UAMssh/fjd/MAF_FJD_v3.0/metadata/otros_pacientes.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# 
# 
# path="/home/gonzalo/Documents/MAF_FJD_v2.0/metadata/date_2021_03_29/mymetadatapathology_2021_03_29.txt"
# aa = read.table(path, sep="\t", header = TRUE)
# 

######################################### OTROS ####################################
######### COMPROBAR CES QUE ME FALTAN #######
familia_filt<-familia[8997:nrow(familia),]
my_samples<-read.table("all_samples_CES.txt")

gur_not_have<-familia_filt[!familia_filt$SAMPLE %in% my_samples$V1,]


################## casos RP y MD

###read los casos de cris y lidia


DHR_cohorte_excluidos = data.frame(readxl::read_excel("DHR_cohorte_excluidos.xlsx"))  
DHR_cohorte_probandus = data.frame(readxl::read_excel("DHR_cohorte_probandus.xlsx")) 



DHR_probandus_yo<-familia[which(familia$ADN %in% DHR_cohorte_probandus$N_ADN),]
##quitar los que no tengo DHR
DHR_probandus_yo_sin_diag_DHR <- DHR_probandus_yo[!grepl("Distrofia de retina", DHR_probandus_yo$CATEGORY), ]


#export lista que no tengo
not_DHR_probandus_yo<-DHR_cohorte_probandus[which(!DHR_cohorte_probandus$N_ADN %in% familia$ADN),]
## de las que no tengo mirar si estan en el df_fjd de wes
df_fjd_not_DHR_probandus_yo <- merge(df_fjd, not_DHR_probandus_yo, by.x = "...3",by.y="N_ADN")
df_fjd_not_DHR_probandus_yo<-df_fjd[which(df_fjd$...3 %in% not_DHR_probandus_yo$N_ADN),]

#de los que no estan mirar si hay algun familiar en la base de datos
si_fam_not_DHR_probandus_yo <- merge(familia, not_DHR_probandus_yo, by.x = "FAMILY",by.y="N_Family")

not_DHR_probandus_yo<-DHR_cohorte_probandus[which(!DHR_cohorte_probandus$N_ADN %in% familia$ADN),]


DHR_exluidos_yo<-df[which(df$ADN %in% DHR_cohorte_excluidos$N_DNA_BBDD),]
merged_df_excluidos <- merge(df, DHR_cohorte_excluidos, by.x = "ADN",by.y="N_DNA_BBDD")

not_DHR_probandus_yo<-DHR_cohorte_probandus[!(DHR_cohorte_probandus$N_ADN %in% df$ADN),]

not_dhr<-DH

View(table(DHR_cohorte_excluidos$Exclusing_motive))

RP_MD<-order_new_df %>% 
  dplyr::filter(TAG %in% c("RP","MD"))

order_RP_MD <- RP_MD[order(RP_MD$Subtipo,RP_MD$Super_Subtipo, RP_MD$Nombre_comun), ]
View(table(order_RP_MD$Subtipo))

######################## EXTRAER NUMERO DE PACIENTES DE CADA ENFERMEDAD ############33
setwd("/home/graciela/tblab/graciela/metadata_CES") #tblab
df_all <- data.table::fread("all_FJD_25june.txt")
write.table(x = tsv_familia, file = "~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/all_FJD_26june.txt", row.names = F, quote = F, sep = "\t")
write.table(x = familia, file = "~/tblab/graciela/TODO_DBofAFs/correr_hail_tblab/metadata/con_patologia_all_FJD_9julio.txt", row.names = F, quote = F, sep = "\t")
write.table(x = familia, file = "~/tblab/graciela/metadata_CES/con_patologia_all_FJD_9julio.txt", row.names = F, quote = F, sep = "\t")

list_uam <- data.table::fread("~/tblab/graciela/metadata_CES/clean_all_FJD_UAM.txt",header = FALSE)
ingles_cats = data.frame(readxl::read_excel("ingles_cats_3julio.xlsx"))  

df_all<-tsv_familia
df_samples_db<-df_all[df_all$SAMPLE%in%list_uam$V1,] #actual samples in db


########################################### RECUENTO PACIENTES ###################################
############## GENERAL CATEGORIES #############
my_table<-as.data.frame(table(df_samples_db$GENERAL))
# Separate multiple diseases into individual rows
mytable_separated <- my_table %>%
  tidyr::separate_rows(Var1, sep = ":") %>%
  filter(Var1 != "")  
df_counts <- mytable_separated %>%
  group_by(Var1) %>%
  summarise(total_count = sum(Freq)) %>%
  arrange(desc(total_count))

global_tags<-ingles_cats[,c(1,2,3)]
global_tags<-unique(global_tags)
df_counts_general <- merge(df_counts, global_tags, by.x = "Var1",by.y="GlobalTag")

###################TYPE CATEGORIES#############################
my_table<-as.data.frame(table(df_samples_db$CATEGORY))
# Separate multiple diseases into individual rows
mytable_separated <- my_table %>%
  tidyr::separate_rows(Var1, sep = ":") %>%
  filter(Var1 != "")  
df_counts <- mytable_separated %>%
  group_by(Var1) %>%
  summarise(total_count = sum(Freq)) %>%
  arrange(desc(total_count))

general_tags<-ingles_cats[,c(1:6)]
general_tags<-unique(general_tags)

df_counts_Type <- merge(df_counts, general_tags, by.x = "Var1",by.y="TypeTag")


###################SUBTYPE CATEGORIES#############################
my_table<-as.data.frame(table(df_samples_db$SUBCATEGORY))
# Separate multiple diseases into individual rows
mytable_separated <- my_table %>%
  tidyr::separate_rows(Var1, sep = ":") %>%
  filter(Var1 != "")  
df_counts <- mytable_separated %>%
  group_by(Var1) %>%
  summarise(total_count = sum(Freq)) %>%
  arrange(desc(total_count))

subtype_tags<-ingles_cats[,c(1:9)]
subtype_tags<-unique(subtype_tags)

df_counts_Subtype <- merge(df_counts, subtype_tags, by.x = "Var1",by.y="SubtypeTag")

####### ordenar
df_counts_general <- df_counts_general[order(df_counts_general$General), ]
df_counts_Type <- df_counts_Type[order(df_counts_Type$General,df_counts_Type$Subtipo), ]
df_counts_Subtype <- df_counts_Subtype[order(df_counts_Subtype$General,df_counts_Subtype$Subtipo,df_counts_Subtype$Super_Subtipo), ]

named_list_df <- list(
  "general_pheno" = df_counts_general,
  "type_pheno" = df_counts_Type,
  "subtype_pheno"= df_counts_Subtype
)
writexl::write_xlsx(x=named_list_df,
                    path = "recuento_pacientes_9julio24.xlsx",
                    col_names=TRUE,
                    format_headers = TRUE)

ingles_cats[which(!ingles_cats$SubtypeTag %in% df_counts_Subtype$Var1),]
df_counts_Subtype[!which(df_counts_Subtype$Var1%in% ingles_cats$SubtypeTag),]
###################################################################################################

