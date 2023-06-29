rm(list=ls())

# Environmental data ------------------------------------------------------
library(EnvRtype)
library(geodata)
library(raster)

# Data --------------------------------------------------------------------
pheno = readRDS('../Data/pheno.rds')

# Climatic ----------------------------------------------------------------
envdata = pheno |> reframe(lat = unique(lat), lon = unique(lon),
                           plant.dat = unique(plant.date),
                           harv.dat = unique(harv.date),
                           .by = env)

envcov = get_weather(env.id = envdata$env, lat = envdata$lat, lon = envdata$lon,
                     start.day = envdata$plant.dat,end.day = envdata$harv.dat)
envcov = processWTH(env.data = envcov, Tbase1 = 8, Tbase2 = 30, Topt1 = 37,
                    Topt2 = 45)

## Soil and elevation --------------------------------------------------------
dir.create(path = 'Rasters')

soil_world(var = 'bdod', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'clay', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'nitrogen', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'phh2o', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'sand', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'silt', depth = 15, stat = 'mean', path = 'Raster')
soil_world(var = 'soc', depth = 15, stat = 'mean', path = 'Raster')
elevation_30s(country = 'BRA', path = 'Raster')

rasts <- list.files(path = 'Raster', pattern = '.tif', full.names = T)
namerasts = c('bdod','clay','nit','phh2o','sand', 'silt','soc', 'alt')

for(i in 1:length(rasts)){
  raster_file <- raster(rasts[i])
  envcov = extract_GIS(covraster = raster_file, Latitude = 'LAT', 
                       name.out = namerasts[i], Longitude = 'LON', 
                       env.data = envcov, env.id = 'env')
}

write.csv(envcov, '../Data/envcov.csv', row.names = F)
