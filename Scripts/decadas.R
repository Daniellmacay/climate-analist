library(ncdf4)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)

#Establecer directorio
setwd("C:/Mi_proyecto/Data")
#Abrir archivo NetCDF
nc <- nc_open("Land_and_Ocean_Alternate_EqualArea.nc")   o
nc <- nc_open("Land_and_Ocean_EqualArea.nc")

#Extraer variables
lon <- ncvar_get(nc, "longitude")          
lat <- ncvar_get(nc, "latitude")           
temp <- ncvar_get(nc, "temperature")       
tiempo <- ncvar_get(nc, "time") 
nc_close(nc)
dim(temp)

# Convertir 'tiempo' en años (si es necesario)
tiempo_years <- floor(tiempo)  # Si 'tiempo' es en fracciones de año, toma solo el año completo

# Extraer los índices de los meses (ejemplo para enero)
enero_indices <- seq(1, length(tiempo), by = 12)
temp_eneros <- temp[, enero_indices]

# Limpiar filas con NA
na_rows <- apply(temp_eneros, 1, function(x) all(is.na(x)))
temp_matrix_clean <- temp_eneros[!na_rows, ]

# PCA para enero (todo el período)
res_pca <- PCA(temp_matrix_clean, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Scree plot con título personalizado
fviz_eig(res_pca, addlabels = TRUE) + 
  ggtitle("Enero")

# comparar varios meses juntos (estaciones)
invierno_indices <- sort(c(
  seq(12,length(tiempo) , by = 12),  # Diciembre
  seq(1, length(tiempo) , by = 12),   # Enero
  seq(2, length(tiempo) , by = 12)    # Febrero
))
# Extraer los datos de temperatura para los meses de invierno
temp_invierno <- temp[, invierno_indices]

# Limpiar filas con NA
na_rows_invierno <- apply(temp_invierno, 1, function(x) all(is.na(x)))
temp_matrix_invierno_clean <- temp_invierno[!na_rows_invierno, ]

# PCA para invierno (todo el periodo estacional)
res_pca_invierno <- PCA(temp_matrix_invierno_clean, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Scree plot con título personalizado
fviz_eig(res_pca_invierno, addlabels = TRUE) + 
  ggtitle("Invierno (Diciembre, Enero, Febrero)")

###Decadas###
# Extraer los datos para las décadas 1950-1959 y 1960-1969
enero_indices_1950s <- which(tiempo_years >= 1950 & tiempo_years < 1960)
temp_1950s <- temp[, enero_indices_1950s]

enero_indices_1960s <- which(tiempo_years >= 1960 & tiempo_years < 1970)
temp_1960s <- temp[, enero_indices_1960s]

# Limpiar datos NA para cada una de las décadas
temp_1950s_clean <- na.omit(temp_1950s)
temp_1960s_clean <- na.omit(temp_1960s)

# Graficar Scree Plot para cada década
fviz_eig(pca_1950s, addlabels = TRUE) +
  ggtitle("Scree plot Enero 1950s")

fviz_eig(pca_1950s, addlabels = TRUE) +
  ggtitle("Scree plot Enero 1960s")

# Guardar el gráfico de Scree Plot
ggsave("scree_plot_Enero_1950s.png")
ggsave("scree_plot_Enero_1960s.png")

# PCA para la década de 1950
pca_1950s <- PCA(temp_1950s_clean, scale.unit = TRUE, ncp = 5, graph = FALSE)
# PCA para la década de 1960
pca_1960s <- PCA(temp_1960s_clean, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Graficar Scree Plot para cada década
pca_1950s_plot <- fviz_pca_ind(pca_1950s, axes = c(1, 2), geom = "point", pointsize = 2) +
  labs(title = "PCA: Enero Década de 1950") +
  theme_minimal()
print(pca_1950s_plot)
pca_1960s_plot <- fviz_pca_ind(pca_1960s, axes = c(1, 2), geom = "point", pointsize = 2) +
  labs(title = "PCA: Enero Década de 1960") +
  theme_minimal()
print(pca_1960s_plot)


# Guardar el gráfico de Scree Plot para la década de 1950
ggsave("PCA_Enero_1950s.png", plot = pca_1950s_plot)
ggsave("PCA_Enero_1960s.png", plot = pca_1960s_plot)

#Heatmaps
library(ggplot2)
library(reshape2)

# Función para graficar un heatmap por década
graficar_heatmap <- function(matriz, decada) {
  df <- as.data.frame(matriz)
  df$pixel <- 1:nrow(df)
  df_long <- melt(df, id.vars = "pixel", variable.name = "anio", value.name = "temperatura")
  df_long$anio <- as.numeric(gsub("V", "", df_long$anio)) + (as.numeric(substr(decada, 1, 4)) - 1)
  
  ggplot(df_long, aes(x = anio, y = pixel, fill = temperatura)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey80") +
    labs(title = paste("Heatmap de temperaturas -", decada),
         x = "Año", y = "Pixel espacial") +
    theme_minimal()
}

# Graficar
heatmap_1950s <- graficar_heatmap(temp_1950s_clean, "1950s")
heatmap_1980s <- graficar_heatmap(temp_1980s_clean, "1980s")
heatmap_2010s <- graficar_heatmap(temp_2010s_clean, "2010s")

# Mostrar o guardar
ggsave("heatmap_1950s.png", heatmap_1950s)
ggsave("heatmap_1980s.png", heatmap_1980s)
ggsave("heatmap_2010s.png", heatmap_2010s)

#Combinacón de las tres graficas
install.packages("patchwork")
library(patchwork)

# Panel de los tres heatmaps
panel <- heatmap_1950s + heatmap_1980s + heatmap_2010s + 
  plot_layout(ncol = 3)

# Mostrar
panel

# Guardar
ggsave("heatmaps_decadas_patchwork.png", panel, width = 18, height = 6, dpi = 300)

