# Análisis PCA de Temperaturas Globales (1850 - Reciente)

Este proyecto realiza un análisis de Componentes Principales (PCA) sobre datos de temperatura media mensual global (tierra y océano) desde 1850 hasta fechas recientes, usando archivos NetCDF proporcionados por Berkeley Earth.


## 📋 Descripción

El objetivo es explorar patrones espaciales y temporales en las temperaturas globales mediante PCA, analizando diferentes décadas para identificar variaciones climáticas importantes. Se extraen datos para décadas específicas (1950s, 1960s, 2000s), y se generan gráficos como scree plots y mapas de calor para interpretar los resultados.


## 📁 Datos

Este proyecto utiliza archivos `.nc` (NetCDF) de temperatura media mensual global de superficie (tierra y océano), con datos desde 1850 hasta fechas recientes.

Debido al tamaño de los archivos (~100 MB cada uno), **no se incluyen directamente en este repositorio**. Puedes descargarlos desde el sitio oficial de Berkeley Earth:

🔗 [Berkeley Earth - Global Gridded Data](https://berkeleyearth.org/data/)

En el menú desplegable **Global Temperature Data**, selecciona:

**Global Monthly Land + Ocean**  
****Average Temperature with Air Temperatures at Sea Ice (Recommended; 1850 – Recent)**  
******Equal Area (~100 MB)**
***Average Temperature with Water Temperatures at Sea Ice (1850 – Recent)**
****Equal Area (~100 MB)**

Una vez descargados los archivos `.nc`, colócalos en la carpeta: /Data/

## Requisitos 
- R (versión >= 4.0)
- Paquetes R:
  - `ncdf4`
  - `dplyr`
  - `tidyr`
  - `FactoMineR`
  - `factoextra`
  - `ggplot2`

Puedes instalar los paquetes necesarios con:

```r
install.packages(c("ncdf4", "dplyr", "tidyr", "FactoMineR", "factoextra", "ggplot2"))
