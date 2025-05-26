# Análisis Climático por Estaciones

Este proyecto muestra el análisis de datos climáticos extraídos desde archivos NetCDF. A continuación se presentan los procedimientos utilizados para el procesamiento de datos, junto con los resultados gráficos divididos por estación.

---

## 📌 Índice

- [1. Extracción y apertura de datos climáticos](#1-extracción-y-apertura-de-datos-climáticos)
- [2. Análisis por estación](#2-análisis-por-estación)
  - [2.1 Invierno](#21-invierno)
  - [2.2 Primavera](#22-primavera)
  - [2.3 Verano](#23-verano)
  - [2.4 Otoño](#24-otoño)

---

## 1. Extracción y apertura de datos climáticos

Se trabajó con archivos `.nc` (NetCDF) que contienen datos de temperatura y coordenadas geográficas. Estos fueron extraídos y organizados para su posterior análisis.

```r
library(ncdf4)

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

## 2. Promedio mensual por estación
Se calcula el promedio mensual para cada estación, agrupando por coordenadas y fechas.
