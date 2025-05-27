📊 Análisis de Temperaturas: PCA por Estación y Década
Este proyecto analiza los patrones de temperatura utilizando Análisis de Componentes Principales (PCA) a partir de un archivo NetCDF. Se realiza una evaluación mensual, por estación y por década (1950s, 1980s, 2000s).

---

## 📌 Índice

- [1. Extracción y apertura de datos climáticos](#1-extracción-y-apertura-de-datos-climáticos)
- [2. Análisis global por estación](#2-definicion-de-estaciones,-décadas-y-procesamiento-temporal)
- [3. Funciones para analizar PCA y gráficos](#3-análisis-PCA-por-estación)
  - [3.1 Invierno](#31-invierno)
  - [3.2 Primavera](#32-primavera)
  - [3.3 Verano](#33-verano)
  - [3.4 Otoño](#34-otoño)
- [4. Resultados](#4-resultados)

---

 1. Extracción y apertura de datos climáticos

Se trabajó con archivos `.nc` (NetCDF) que contienen datos de temperatura y coordenadas geográficas. Estos fueron extraídos y organizados para su posterior análisis.

```r
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(reshape2)
library(patchwork)

# Abrir archivo NetCDF
data <- nc_open("archivo.nc")

# Extraer variables
lon <- ncvar_get(data, "lon")
lat <- ncvar_get(data, "lat")
temp <- ncvar_get(data, "temperature")

# Cerrar archivo
nc_close(data)

# Mostrar grafico de temperaturas globales
plot(lon,lat)
```
![extracción de temperatura global](Graphics/Rplot.png)

 2. Definición de estaciones, décadas y procesamiento temporal.
Se define la asignación de meses a estaciones climáticas y las décadas que se analizarán..

```r
tiempo_years <- floor(tiempo)

meses_estaciones <- list(
  invierno = c(12, 1, 2),
  primavera = c(3, 4, 5),
  verano    = c(6, 7, 8),
  otono     = c(9, 10, 11)
)

decadas <- list(
  `1950s` = 1950:1959,
  `1980s` = 1980:1989,
  `2000s` = 2000:2009
)

get_indices <- function(months) {
  unlist(lapply(months, function(m) seq(m, length(tiempo), by = 12)))
}

```

3. Funciones para análisis PCA y gráficos
Se definen funciones para limpiar datos, hacer PCA, y generar gráficos (Scree Plot, PCA plot y Heatmap).

```r
hacer_pca <- function(temp_subset) {
  na_rows <- apply(temp_subset, 1, function(x) all(is.na(x)))
  temp_clean <- temp_subset[!na_rows, ]
  pca <- PCA(temp_clean, scale.unit = TRUE, ncp = 5, graph = FALSE)
  return(list(pca = pca, temp_clean = temp_clean))
}

graficar_scree <- function(pca_obj, titulo) {
  fviz_eig(pca_obj, addlabels = TRUE) + ggtitle(titulo)
}

graficar_pca <- function(pca_obj, titulo) {
  fviz_pca_ind(pca_obj, axes = c(1, 2), geom = "point", pointsize = 2) +
    labs(title = titulo) +
    theme_minimal()
}

graficar_heatmap <- function(matriz, decada, estacion) {
  df <- as.data.frame(matriz)
  df$pixel <- 1:nrow(df)
  df_long <- melt(df, id.vars = "pixel", variable.name = "tiempo", value.name = "temperatura")
  df_long$tiempo <- as.numeric(gsub("V", "", df_long$tiempo)) + min(decadas[[decada]]) - 1
  
  ggplot(df_long, aes(x = tiempo, y = pixel, fill = temperatura)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey80") +
    labs(title = paste("Heatmap -", toupper(estacion), "-", decada),
         x = "Año", y = "Pixel espacial") +
    theme_minimal()
}


```
4. Análisis por estación y década: PCA, Screeplots y Heatmaps
Para cada estación y década se realiza el PCA, se guardan los gráficos correspondientes y se genera el heatmap.

```r
for (estacion in names(meses_estaciones)) {
  indices_estacion <- get_indices(meses_estaciones[[estacion]])
  
  # PCA para toda la estación
  temp_estacion <- temp[, indices_estacion]
  resultado_pca <- hacer_pca(temp_estacion)
  pca_obj <- resultado_pca$pca
  
  scree_plot <- graficar_scree(pca_obj, paste("Scree Plot -", toupper(estacion)))
  pca_plot <- graficar_pca(pca_obj, paste("PCA -", toupper(estacion)))
  
  ggsave(filename = file.path("Graficos", estacion, paste0("ScreePlot_", estacion, ".png")), plot = scree_plot)
  ggsave(filename = file.path("Graficos", estacion, paste0("PCA_", estacion, ".png")), plot = pca_plot)
  
  for (decada in names(decadas)) {
    años_decada <- decadas[[decada]]
    indices_decada <- which(tiempo_years %in% años_decada)
    indices_combinados <- intersect(indices_estacion, indices_decada)
    temp_subset <- temp[, indices_combinados]
    
    if (ncol(temp_subset) > 2) {
      resultado_pca_decada <- hacer_pca(temp_subset)
      pca_obj_decada <- resultado_pca_decada$pca
      temp_clean <- resultado_pca_decada$temp_clean
      
      scree_plot_decada <- graficar_scree(pca_obj_decada, paste("Scree Plot -", toupper(estacion), decada))
      ggsave(filename = file.path("Graficos", estacion, decada, paste0("ScreePlot_", estacion, "_", decada, ".png")), plot = scree_plot_decada)
      
      pca_plot_decada <- graficar_pca(pca_obj_decada, paste("PCA -", toupper(estacion), decada))
      ggsave(filename = file.path("Graficos", estacion, decada, paste0("PCA_", estacion, "_", decada, ".png")), plot = pca_plot_decada)
      
      heatmap_plot <- graficar_heatmap(temp_clean, decada, estacion)
      ggsave(filename = file.path("Graficos", estacion, decada, paste0("Heatmap_", estacion, "_", decada, ".png")), plot = heatmap_plot)
    }
  }
}


```
Screeplot:



PCA:


Heatmap:
3.1 Invierno general

```r
invierno <- monthly_avg %>% filter(month %in% c(12, 1, 2))
mat <- spread(invierno, key = lon, value = mean_temp)
pca_inv <- prcomp(mat, scale. = TRUE)
plot(pca_inv, type = "l")
biplot(pca_inv)

```
Screeplot:
![Screeplot Invierno](Graphics/screePlot_invierno.png)


PCA:
![PCA invierno](Graphics/PCA_invierno.png)

Heatmap:
3.1.1 Análisis estacional por decada
Década 1950-1959

```r
invierno_50s <- invierno %>% filter(year %in% 1950:1959)
mat_50s <- spread(invierno_50s, key = lon, value = mean_temp)
pca_inv_50s <- prcomp(mat_50s, scale. = TRUE)
plot(pca_inv_50s, type = "l")
biplot(pca_inv_50s)
```
Década 1980-1989

```r
invierno_80s <- invierno %>% filter(year %in% 1980:1989)
mat_80s <- spread(invierno_80s, key = lon, value = mean_temp)
pca_inv_80s <- prcomp(mat_80s, scale. = TRUE)
plot(pca_inv_80s, type = "l")
biplot(pca_inv_80s)

```
Década 2000-2009

```r
invierno_00s <- invierno %>% filter(year %in% 2000:2009)
mat_00s <- spread(invierno_00s, key = lon, value = mean_temp)
pca_inv_00s <- prcomp(mat_00s, scale. = TRUE)
plot(pca_inv_00s, type = "l")
biplot(pca_inv_00s)

```

3.2 Primavera general
Código PCA:

```r
primavera <- monthly_avg %>% filter(month %in% c(3, 4, 5))
mat <- spread(primavera, key = lon, value = mean_temp)
pca_prim <- prcomp(mat, scale. = TRUE)
plot(pca_prim, type = "l")
biplot(pca_prim)

```
Screeplot:
![Screeplot Invierno](graficos/screeplot_invierno.png)


PCA:


Heatmap:

3.2.1 Análisis estacional por decada
Década 1950-1959

```r
primavera_50s <- primavera %>% filter(year %in% 1950:1959)
mat_50s <- spread(primavera_50s, key = lon, value = mean_temp)
pca_prim_50s <- prcomp(mat_50s, scale. = TRUE)
plot(pca_prim_50s, type = "l")
biplot(pca_prim_90s)

```
Década 1980-1989

```r
primavera_80s <- primavera %>% filter(year %in% 1980:1989)
mat_80s <- spread(primavera_80s, key = lon, value = mean_temp)
pca_prim_80s <- prcomp(mat_80s, scale. = TRUE)
plot(pca_prim_80s, type = "l")
biplot(pca_prim_80s)

```
Década 2000-2009

```r
primavera_00s <- primavera %>% filter(year %in% 2000:2009)
mat_00s <- spread(primavera_00s, key = lon, value = mean_temp)
pca_prim_00s <- prcomp(mat_00s, scale. = TRUE)
plot(pca_prim_00s, type = "l")
biplot(pca_prim_00s)

```

3.3 Verano general
Código PCA:

```r
verano <- monthly_avg %>% filter(month %in% c(6, 7, 8))
mat <- spread(verano, key = lon, value = mean_temp)
pca_ver <- prcomp(mat, scale. = TRUE)
plot(pca_ver, type = "l")
biplot(pca_ver)

```
Screeplot:
![Screeplot Invierno](graficos/screeplot_invierno.png)


PCA:


Heatmap:
3.2.1 Análisis estacional por decada
Década 1950-1959

```r
verano_50s <- verano %>% filter(year %in% 1950:1959)
mat_50s <- spread(verano_50s, key = lon, value = mean_temp)
pca_ver_50s <- prcomp(mat_50s, scale. = TRUE)
plot(pca_ver_50s, type = "l")
biplot(pca_ver_50s)

```
Década 1980-1989

```r
verano_80s <- verano %>% filter(year %in% 1980:1989)
mat_80s <- spread(verano_80s, key = lon, value = mean_temp)
pca_ver_80s <- prcomp(mat_80s, scale. = TRUE)
plot(pca_ver_80s, type = "l")
biplot(pca_ver_80s)

```
Década 2000-2009

```r
verano_00s <- verano %>% filter(year %in% 2000:2009)
mat_00s <- spread(verano_00s, key = lon, value = mean_temp)
pca_ver_00s <- prcomp(mat_00s, scale. = TRUE)
plot(pca_ver_00s, type = "l")
biplot(pca_ver_00s)

```
3.4 Otoño general
Código PCA:

```r
otono <- monthly_avg %>% filter(month %in% c(9, 10, 11))
mat <- spread(otono, key = lon, value = mean_temp)
pca_oto <- prcomp(mat, scale. = TRUE)
plot(pca_oto, type = "l")
biplot(pca_oto)

```
Screeplot:
![Screeplot Invierno](graficos/screeplot_invierno.png)


PCA:


Heatmap:

3.4.1 Análisis estacional por decada
Década 1950-1959

```r
otono_50s <- otono %>% filter(year %in% 1950:1959)
mat_50s <- spread(otono_50s, key = lon, value = mean_temp)
pca_oto_50s <- prcomp(mat_50s, scale. = TRUE)
plot(pca_oto_50s, type = "l")
biplot(pca_oto_50s)

```
Década 1980-1989

```r
otono_80s <- otono %>% filter(year %in% 1980:1989)
mat_80s <- spread(otono_80s, key = lon, value = mean_temp)
pca_oto_80s <- prcomp(mat_80s, scale. = TRUE)
plot(pca_oto_80s, type = "l")
biplot(pca_oto_80s)

```
Década 2000-2009

```r
otono_00s <- otono %>% filter(year %in% 2000:2009)
mat_00s <- spread(otono_00s, key = lon, value = mean_temp)
pca_oto_00s <- prcomp(mat_00s, scale. = TRUE)
plot(pca_oto_00s, type = "l")
biplot(pca_oto_00s)

```






Heatmap:

4. Resultados
El presente análisis fue desarrollado en R utilizando paquetes como ncdf4, dplyr, ggplot2, y stats para reducción de dimensionalidad (PCA). Las visualizaciones permiten explorar la variabilidad espacial y estacional de los datos de temperatura.
