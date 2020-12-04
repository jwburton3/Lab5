import os
import pandas as pd
import fiona
from rasterio.plot import  show, show_hist
from shapely.geometry import Polygon
import rasterio
import glob
import numpy as np
file = glob.glob(r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/data/L5_big_elk/*.tif')
file.sort()
rasters_B3 = []
rasters_B4 = []
year = []
for raster in file:
    if 'B3' in raster:
        with rasterio.open(raster) as data:
            array = data.read(1)
            rasters_B3.append(array)
            year.append(raster[-11:-7]) 
for raster in file:
    if 'B4' in raster:
        with rasterio.open(raster) as data:
            array = data.read(1)
            rasters_B4.append(array)
fire_perimeter = rasterio.open(r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/data/fire_perimeter.tif').read(1)
final_mean = []
flat_rr = []
for b3,b4 in zip(rasters_B3, rasters_B4):
    ndvi = ((b4 - b3) / (b4 + b3)) 
    healthy = np.where(fire_perimeter ==2) 
    fire = np.where(fire_perimeter == 1)   
    ndvi_mean = ndvi[healthy].mean()    
    rr = ndvi / ndvi_mean           
    burned_mean = rr[fire].mean()   
    final_mean.append(burned_mean)
    flat = rr.flatten()            
    flat_rr.append(flat)
stacked_ratio = np.vstack(flat_rr)
trend_line = np.polyfit(range(10), stacked_ratio, 1)[0]
trend_reshape = trend_line.reshape(280, 459)
mean_coef = np.where(fire_perimeter==1, trend_reshape, np.nan)
coefficient = np.nanmean(mean_coef)
for y,r in zip(year,final_mean):
    print("In", y, "the mean recovery ratio was", r)
print("The mean coefficient of recovery across all years for the burned area is", coefficient)
lab5functions  = (r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/lab5functions.py')
from lab5functions import *
os.chdir(r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/data')
DEM_ras = rasterio.open(r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/data/bigElk_dem.tif')
DEM_read = DEM_ras.read(1)
slp, aspect = slopeAspect(DEM_read, 90)
reclass_aspect = reclassAspect(aspect)
DEM_reclass = reclassByHisto(slp, 10) 
def zonal_stats_table(zones, value_raster, csv_name):
    mean_stats = []
    max_stats = []
    min_stats = []
    count_stats = []
    std_stats = []
    zone_num = []
    for u in np.unique(zones):
        ras = np.where(zones==u, u, np.nan)
        min_stats.append(np.nanmin(ras * value_raster))
        max_stats.append(np.nanmax(ras * value_raster))
        mean_stats.append(np.nanmean(ras * value_raster))
        std_stats.append(np.nanstd(ras * value_raster))
        count_stats.append(np.where(zones == u, 1, 0).sum())
        zone_num.append(int(u))
    stats = {'Zone' : zone_num, 'MIN': min_stats, 'MAX': max_stats, 'MEAN': mean_stats, 'STD': std_stats, 'COUNT': count_stats}
    df = pd.DataFrame(stats)
    df.to_csv(csv_name)
    return df

zonal_stats_table(DEM_reclass, mean_coef, 'slp.csv')
zonal_stats_table(reclass_aspect, mean_coef, 'aspect.csv') 

with rasterio.open(r'/Users/jonathanburton/Desktop/Fall2020/Geog5092/lab5/data/bigElk_dem.tif') as dataset:
    with rasterio.open(f'//Users//jonathanburton//Desktop//Fall2020//Geog5092//lab5//datarecovery_coefficent.tif' , 'w',
                       driver='GTiff',
                       height=mean_coef.shape[0],
                       width=mean_coef.shape[1],
                       count=1,
                       dtype=mean_coef.dtype,
                       crs=dataset.crs,
                       transform=dataset.transform, 
                       nodata=dataset.nodata
                      ) as out_dataset:
        out_dataset.write(mean_coef,1)
        
conclusion = "It appears that relatively steeper slopes and a southeast to southwest aspect leads to increased vegatative recovery compared to other slope or aspects."
print(conclusion)