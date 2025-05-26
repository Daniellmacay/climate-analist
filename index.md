## Extracción y apertura de datos climáticos

Se trabajó con archivos .nc (NetCDF) que contienen datos de temperatura y coordenadas geográficas.

```r
library(ncdf4)
data <- nc_open("archivo.nc")
lon <- ncvar_get(data, "lon")
lat <- ncvar_get(data, "lat")
temp <- ncvar_get(data, "temperature")
nc_close(data)
plot(lon,lat)
