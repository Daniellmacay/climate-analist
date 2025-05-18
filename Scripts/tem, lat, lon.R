# Cargar librerias necesarias
library(ncdf4)

#Establecer directorio
setwd("C:/Mi_proyecto/Data")

#brir archivo NetCDF
nc <- nc_open("Land_and_Ocean_Alternate_EqualArea.nc")   o
nc <- nc_open("Land_and_Ocean_EqualArea.nc")

#Extraer variables
lon <- ncvar_get(nc, "longitude")          
lat <- ncvar_get(nc, "latitude")           
temp <- ncvar_get(nc, "temperature")       
tiempo <- ncvar_get(nc, "time") 

#Cerrar archivo NetCDF
nc_close(nc)

#Plots
plot(lon, lat)

#Guardar plots como imagen PNG
png("Land_and_Ocean_Alternate_EqualArea.png", width = 800, height = 600)
plot(lon, lat, main = "Temperaturas", xlab = "Longitud", ylab = "Latitud")
dev.off()



