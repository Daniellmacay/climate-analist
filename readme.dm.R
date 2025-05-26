# 🌡️ Análisis Climático Multidecadal por PCA

Este proyecto tiene como objetivo analizar las variaciones espaciales de temperatura superficial global utilizando técnicas de Análisis de Componentes Principales (PCA) por estación del año (invierno, primavera, verano y otoño) y por década (1950s, 1980s, 2010s).

Se incluyen los pasos de preprocesamiento de datos climáticos en formato NetCDF, extracción mensual, aplicación de PCA, y visualización mediante Scree plots, gráficos PCA y mapas de calor.

---

## 📁 1. Extracción y apertura de datos climáticos

Se trabajó con archivos `.nc` (NetCDF) que contienen datos de temperatura y coordenadas geográficas. Estos fueron extraídos y organizados para su posterior análisis.

```r
library(ncdf4)

# Abrir archivo NetCDF
data <- nc_open("data/Land_and_Ocean_Alternate_EqualArea.nc")

# Extraer variables
lon <- ncvar_get(data, "longitude")
lat <- ncvar_get(data, "latitude")
temp <- ncvar_get(data, "temperature")

# Cerrar archivo
nc_close(data)
