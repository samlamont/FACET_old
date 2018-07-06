# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:09:03 2018

@author: sam.lamont
"""
import os
import subprocess
import timeit

import numpy as np
import fiona
import rasterio
from rasterio.features import shapes
from shapely.geometry import LineString, MultiLineString, MultiPoint, mapping, Polygon, shape
from shapely.ops import transform, split
from functools import partial
import pyproj

import geopandas as gpd
from geopandas.tools import sjoin

## ==========================================================================
##   For dissolving line features    
## ==========================================================================    
#def dissolve_line_features(str_lines_path, output_filename):
#    print('Dissolve line features...')
#
#    lst_all=[]
#    
#    with fiona.open(str_lines_path) as lines:
#        
#        crs = lines.crs
#        
#        schema = {'geometry':'MultiLineString', 'properties':{'linkno':'int:6'}}         
#        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:           
#            
#            for line in lines:
#                
#                lst_all.append(line['geometry']['coordinates'])
#                        
#            output.write({'properties':{'linkno':9999}, 'geometry':{'type':'MultiLineString','coordinates':lst_all}})                            
#        
#    return
#
## ==========================================================================
##   For points along line feature at uniform distance    
## ==========================================================================    
#def points_along_line_features(str_diss_lines_path, output_filename, p_interp_spacing_m):
#    print('Points along line features...')
#    
##    p_interp_spacing_m = 1000 # meters
#
#    with fiona.open(str_diss_lines_path) as line:
#        
#        crs = line.crs
#        line = line[0]        
#        
#        # NOTE:  Need to check whether or not to convert units here!
#        
#        # Geometry transform function based on pyproj.transform
#        project = partial(
#            pyproj.transform,
#            pyproj.Proj(init=crs['init']),  
#            pyproj.Proj(init='EPSG:26914')) # UTM14N        
#                                     
#        schema = {'geometry':'Point', 'properties':{'id':'int:6'}}         
#        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:
#                    
#            line_shply = MultiLineString(line['geometry']['coordinates'])
#            
#            line2 = transform(project, line_shply)  
#            length_m = line2.length
#           
#            length_deg = line_shply.length # units depend on crs
#            
#            p_interp_spacing = (length_deg*p_interp_spacing_m)/(length_m)
#               
#            step_lens = np.arange(0, length_deg, p_interp_spacing) # p_interp_spacing in projection units?
#    
#            for i, step in enumerate(step_lens): # lambda here instead?
#              
#                i_pt = np.array(line_shply.interpolate(step))
#    
#                output.write({'properties':{'id':i}, 'geometry':{'type':'Point','coordinates':i_pt}})                            
#        
#    return
#def split_lines_with_points(str_net_path, str_pts_path):
#    print('Splitting lines...')
#
#    # Open the delineated stream network...
#    with fiona.open(str_net_path) as lines:
#        
#        # Open points file...        
#        with fiona.open(str_pts_path) as pts:
#            
#            lst_lines=[]
#            lst_pts=[]
#            
#            # Create MultiLineString and MultiPoint layers...
#            for line in lines:
#                lst_lines.extend(line['geometry']['coordinates'])
#                
#            for pt in pts:
#                lst_pts.append(pt['geometry']['coordinates'])
#                
#            ml_lines = MultiLineString([lst_lines])
#            mp_pts = MultiPoint(lst_pts)
#            
##            test = ml_lines.intersection(mp_pts)
#            
#            print('Performing split...')
#            test = split(ml_lines, mp_pts) # THIS TAKES FOREVER
#            
#            print('pause')
#def intersect_lines_catchments(str_net_path, str_catch_path, str_catchlen_path):
 
#    # << Create a multilinelayer >>
#    # Open the delineated stream network...
#    with fiona.open(str_net_path) as lines:
#        
#        crs = lines.crs
#        
#        print('Creating MultiLineString...')
#        
#        # Open an output file for dissolved lines...THIS IS NOT REALLY NEEDED, ONLY FOR TESTING
#        schema = {'geometry':'MultiLineString', 'properties':{'linkno':'int:6'}}         
#        with fiona.open(str_diss_path, 'w', 'ESRI Shapefile', schema, crs) as output:           
#            
#            lst_all=[]
#            for line in lines:
#                
#                lst_all.append(line['geometry']['coordinates'])
#                            
#            output.write({'properties':{'linkno':9999}, 'geometry':{'type':'MultiLineString','coordinates':lst_all}})         
        
#    with fiona.open(str_catch_path, 'r') as sheds_in:
#       with fiona.open(str_net_path) as lines:
#        
#            schema = sheds_in.schema.copy()
#            crs = sheds_in.crs
#            schema['properties']['length'] = 'float'
#            
#            with fiona.open(str_catchlen_path, 'w', 'ESRI Shapefile', schema, crs) as sheds_out:
#            
#                for i, shed in enumerate(sheds_in):
#                    
#                    poly = Polygon(shed['geometry']['coordinates'][0])
#                    
#                    tot_len=0
#                    for line in lines:
#                        
#                        ls = LineString(line['geometry']['coordinates'])
#                        if poly.intersects(ls):
#                           
#                            print(i)
#                            tot_len+=ls.length
#                    
#                    shed['properties']['length'] = tot_len  
#                    sheds_out.write({'properties':shed['properties'],'geometry':mapping(shape(shed['geometry']))})
#                
#                print(i)
#                    
#                    
def intersect_validation_catchments(str_valid_path, str_catch_path):
    
    print('Creating GeoDataFrames...')
    gdf_sheds = gpd.GeoDataFrame.from_file(str_catch_path)
    gdf_valid = gpd.GeoDataFrame.from_file(str_valid_path)
    
    print('Performing Intersections...')
    intersections = gpd.overlay(gdf_sheds, gdf_valid, how='intersection')
    
    print('Write out file...')
    str_out_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_fema_intersect_TEST.shp"
    intersections.to_file(str_out_path)

    print('pause')

# ==========================================================================
#   First dissolve line features, then generate points    
# ==========================================================================    
def create_points_at_uniform_length(str_net_path, str_diss_path, str_pts_path, p_interp_spacing_m):
    
    print('Dissolve line features...')

    lst_all=[]
    
    # Open the delineated stream network...
    with fiona.open(str_net_path) as lines:
        
        crs = lines.crs
        
        # Open an output file for dissolved lines...THIS IS NOT REALLY NEEDED, ONLY FOR TESTING
        schema = {'geometry':'MultiLineString', 'properties':{'linkno':'int:6'}}         
        with fiona.open(str_diss_path, 'w', 'ESRI Shapefile', schema, crs) as output:           
            
            for line in lines:
                
                lst_all.append(line['geometry']['coordinates'])
                            
            output.write({'properties':{'linkno':9999}, 'geometry':{'type':'MultiLineString','coordinates':lst_all}})   

        print('Points along line features...')    
        # Open the dissolve network...(skip this write/read step?)
        with fiona.open(str_diss_path) as line:
            
            crs = line.crs
            line = line[0]        
            
            # NOTE:  Need to check whether or not to convert units here!
            
            # Geometry transform function based on pyproj.transform to convert length between decimal degrees and meters
            project = partial(
                pyproj.transform,
                pyproj.Proj(init=crs['init']),  
                pyproj.Proj(init='EPSG:26914')) # UTM14N For TX data;  NOTE: HARDCODED FOR NOW!      
        
            # Open output point shapefile...                                     
            schema = {'geometry':'Point', 'properties':{'id':'int:6'}}         
            with fiona.open(str_pts_path, 'w', 'ESRI Shapefile', schema, crs) as pt_output:
                        
                line_shply = MultiLineString(line['geometry']['coordinates'])
                                        
                line2 = transform(project, line_shply)  
                length_m = line2.length
               
                length_deg = line_shply.length # units depend on crs
                
                p_interp_spacing = (length_deg*p_interp_spacing_m)/(length_m) # convert from meters to dec. degrees
                   
                step_lens = np.arange(0, length_deg, p_interp_spacing) # p_interp_spacing in projection units?
                
#                # Open the output line segment shapefile...
#                schema = {'geometry':'LineString', 'properties':{'id':'int:6', 'len_m':'float'}}         
#                with fiona.open(str_segs_path, 'w', 'ESRI Shapefile', schema, crs) as ln_output:                
        
#                i_pt_prev=np.empty(0)
                for i, step in enumerate(step_lens): # lambda here instead?              
                    i_pt = np.array(line_shply.interpolate(step))
                    
#                        if i_pt_prev.any():
#                            ls = LineString([i_pt, i_pt_prev])
#                            
#                            # Write out the lines (with length)...
#                            len_m = transform(project, ls).length
#                            
#                            ln_output.write({'geometry': mapping(ls), 'properties':{'id':i, 'len_m':len_m}})
#                    
#                        i_pt_prev = i_pt
                    
                    # Write out the point...
                    pt_output.write({'properties':{'id':i}, 'geometry':{'type':'Point','coordinates':i_pt}}) 
                
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
            schema={'properties': [('id', 'int')],
                    'geometry': 'Polygon'}) as dst:
        dst.writerecords(results)
        
    return

# ==========================================================================
#                 << MAIN INTERFACE FOR FUNCTIONS >>    
# ==========================================================================    
if __name__ == '__main__':    

    print('\n<<< Start >>>\r\n')
    start_time_0 = timeit.default_timer() 
    
    # << INPUT FILE PATHS >>
    # Streamlines...
    str_net_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_10m_dem_breach_net.shp"
#    str_net_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_net.shp"   # Streamlines delineated from DEM
    
    # D8 flow direction...
    str_d8fdr_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_p.tif" # D8 flow direction grid
    
    # Catchments...
    str_catch_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_10m_dem_breach_p_uniform.shp"
    
    # Validation Polygon...
    str_valid_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_fema_a_ae_zones.shp"
    
    # << OUTPUT FILE PATHS >>
    str_segs_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_10m_dem_breach_net_segs_TEST.shp"    
    str_diss_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_10m_dem_breach_net_diss_TEST.shp"
#    str_diss_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_net_diss.shp"
    
#    str_pts_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_net_pts.shp"
    str_pts_path = r"C:\Terrain_and_Bathymetry\OWP\hand_development\sanantonio\sa_10m_dem_breach_net_pts_TEST.shp"
    
    # << SPACING >>
    p_interp_spacing_m = 1000 # meters 
    
    # << CALL FUNCTIONS FOR UNIFORM REACH POINTS AND CATCHMENTS >> 
#    intersect_lines_catchments(str_net_path, str_catch_path, str_catchlen_path)   
#    split_lines_with_points(str_net_path, str_pts_path) # NOTE:  THIS TAKES FOREVER
    
    create_points_at_uniform_length(str_net_path, str_diss_path, str_pts_path, p_interp_spacing_m)
       
    taudem_gagewatershed(str_pts_path, str_d8fdr_path)    
       
#    intersect_validation_catchments(str_valid_path, str_catch_path) 
    
        
    print('\n<<< End >>>\r\n')
    print('Run time:  {}'.format(timeit.default_timer() - start_time_0))