
# Autor: Gonzalo Núñez Moreno
# Fecha: 17/03/2021


# Data loading

args = commandArgs(TRUE)

ruta_pacientes = args[1]
pacientes = read.table(ruta_pacientes, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

path_relation = args[2]
df_relation = read.table(path_relation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ruta_stad = args[3]
df_stad = read.table(ruta_stad, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df_stad$COVERED = df_stad$N_GENO - df_stad$N_MISS





# Filtering

df_35 = df_relation[df_relation$PI_HAT > 0.35,]
df_35_copy = df_35

# GUR: Remover filas que tienen todos los valores como NA (habia 3/403 filas con todo NA)
df_35 <- df_35[!apply(df_35, 1, function(x) all(is.na(x))), ]
df_35_copy = df_35


df_out = data.frame(stringsAsFactors = FALSE)
contador = 1
while (nrow(df_35) != 0) {

  # En cada iteración se calcula la suma de variantes de las relaciones
  muestras_tot = unique(c(df_35$IID1, df_35$IID2))
  cover_35 = data.frame(Var1 = muestras_tot, stringsAsFactors = FALSE)
  for (muestra in muestras_tot){
    muestras_emparentadas = c(df_35[df_35$IID1 == muestra, "IID2"], df_35[df_35$IID2 == muestra, "IID1"])
    cover_35[cover_35$Var1 == muestra, 2] = sum(df_stad[df_stad$IID %in% muestras_emparentadas, "COVERED"])
  }
  
  cover_35 = merge(cover_35, pacientes, by.x = "Var1", by.y = "SAMPLE", all.x = TRUE)
  cover_35$CATEGORY[is.na(cover_35$CATEGORY)] = "NA"
  cover_35 = cover_35[order(cover_35$V2, cover_35$CATEGORY, decreasing = TRUE),]
  
  muestra = cover_35$Var1[1]
  muestras_emparentadas = c(df_35_copy[df_35_copy$IID1 == muestra, "IID2"], df_35_copy[df_35_copy$IID2 == muestra, "IID1"])
  num_relaciones = length(muestras_emparentadas)
  if (num_relaciones > 5) {
    muestras_emparentadas = ""
    pi_hat = ""
    z0 = ""
    z1 = ""
    z2 = ""
  } else {
    muestras_emparentadas = paste(muestras_emparentadas, collapse = "|")
    pi_hat = paste(c(df_35_copy[df_35_copy$IID1 == muestra, "PI_HAT"], df_35_copy[df_35_copy$IID2 == muestra, "PI_HAT"]), collapse = "|")
    z0 = paste(c(df_35_copy[df_35_copy$IID1 == muestra, "Z0"], df_35_copy[df_35_copy$IID2 == muestra, "Z0"]), collapse = "|")
    z1 = paste(c(df_35_copy[df_35_copy$IID1 == muestra, "Z1"], df_35_copy[df_35_copy$IID2 == muestra, "Z1"]), collapse = "|")
    z2 = paste(c(df_35_copy[df_35_copy$IID1 == muestra, "Z2"], df_35_copy[df_35_copy$IID2 == muestra, "Z2"]), collapse = "|")
  }
  df_out[contador, 1] = muestra
  df_out[contador, 2] = num_relaciones
  df_out[contador, 3] = muestras_emparentadas
  df_out[contador, 4] = pi_hat
  df_out[contador, 5] = z0
  df_out[contador, 6] = z1
  df_out[contador, 7] = z2
  
  contador = contador + 1
  
  df_35 = df_35[df_35$IID1 != muestra, ]
  df_35 = df_35[df_35$IID2 != muestra, ]
}




# Output generation

colnames(df_out) = c("muestra", "num_relaciones", "muestras_emparentadas", "pi_hat", "z0", "z1", "z2")
df_out = merge(df_out, df_stad, by.x = "muestra", by.y = "IID", all.x = TRUE, all.y = FALSE)
df_out = merge(df_out, pacientes, by.x = "muestra", by.y = "SAMPLE", all.x = TRUE, all.y = FALSE)

df_out = df_out[,c("muestra", "num_relaciones", "muestras_emparentadas", "pi_hat", "z0", "z1", "z2", "N_MISS", "N_GENO", "F_MISS", "CATEGORY")]
df_out = df_out[!duplicated(df_out),]

write.table(df_out, args[4], quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(df_out$muestra, args[5], quote = FALSE, row.names = FALSE, col.names = FALSE)


