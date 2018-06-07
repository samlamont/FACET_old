# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:09:03 2018

@author: sam.lamont
"""
import os
import subprocess
import rasterio
from rasterio.features import shapes

# ==========================================================================
#   For points along line feature at uniform distance    
# ==========================================================================    
def taudem_gagewatershed(str_pts_path, str_d8fdr_path):
    
    # << FILE PATHS >>
    path_to_fdr, fdr_filename = os.path.split(str_d8fdr_path)     
    d8fdr = os.path.join(path_to_fdr + '\\' + fdr_filename)
    
    path_to_pts, pts_filename = os.path.split(str_pts_path)
    pts = os.path.join(path_to_pts + '\\' + pts_filename)
    
    gw = os.path.join(path_to_fdr + '\\' + fdr_filename[:-4]+'_gw.tif')
    
    gw_shp = os.path.join(path_to_fdr + '\\' + fdr_filename[:-4]+'_gw.shp')
    
    inputProc = str(4)

    # << GAGE WATERSHED >>
    cmd = 'mpiexec -n ' + inputProc + ' GageWatershed -p ' + '"' + d8fdr + '"' + ' -o ' + '"' + pts + '"'  +  ' -gw ' + '"' + gw + '"'
    
    # Submit command to operating system
    print('Running TauDEM GageWatershed...')
    os.system(cmd)
    # Capture the contents of shell command and print it to the arcgis dialog box
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    
    message = "\n"
    for line in process.stdout.readlines():
        if isinstance(line, bytes):	    # true in Python 3
            line = line.decode()
        message = message + line
    print(message)              

    # << CONVERT WATERSHED GRID TO SHAPEFILE (from https://github.com/mapbox/rasterio/blob/fb93a6425c17f25141ad308cee7263d0c491a0a9/examples/rasterio_polygonize.py)>>       
    with rasterio.open(gw) as src:
        image = src.read(1)
    
    results = (
        {'properties': {'id': v}, 'geometry': s}
        for i, (s, v) 
        in enumerate(
            shapes(image, mask=None, transform=src.transform)))

    with fiona.open(
            gw_shp, 'w', 
            driver='ESRI Shapefile',
            crs=src.crs,
            schema={'properties': [('raster_val', 'int')],
                    'geometry': 'Polygon'}) as dst:
        dst.writerecords(results)
        
    return
        
    return
    
if __name__ == '__main__':    

    print('\n<<< Start >>>\r\n')
#    start_time_0 = timeit.default_timer() 
    
    str_pts_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_net_pts.shp'
#    funcs_v2.points_along_line_features(str_diss_path, str_pts_path)
    
    str_d8fdr_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach_p.tif"
    taudem_gagewatershed(str_pts_path, str_d8fdr_path)    