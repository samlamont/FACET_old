# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 10:36:01 2018

@author: Sam.Lamont
"""
import fiona
import numpy as np
import rasterio
import rasterio.mask
import timeit

# ===============================================================================
#  Delineate a FIM from the HAND grid using depth at each polygon (eg, catchment)
#
#    1. Read in catchment polygons (assume these have attributes necessary for regression?)
#    2. Calculate a HAND height (h) for each polygon based on some attribute(s)
#    3. Delineate FIM per catchment
#
# =============================================================================== 
def fim_hand_poly(str_hand_path, str_sheds_path):
    
    # Open the HAND layer...
    with rasterio.open(str(str_hand_path)) as ds_hand:  
        
        out_meta = ds_hand.meta.copy()
        arr_fim = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])  
        arr_fim[:,:] = out_meta['nodata']
        
        # Open the catchment polygon layer...
        with fiona.open(np.str(str_sheds_path), 'r') as sheds:
        
            for shed in sheds:
                
                # Get HAND height (h) here based on shed properties and regression eqn
                h = shed['properties']['MAX']
                
#                buff = mapping(shed)                            
                
                # Mask the bankpts file for each feature...
                w, out_transform = rasterio.mask.mask(ds_hand, [shed['geometry']], crop=True)  
        
#                w[(w<=h) & (w>=0.)] # Get HAND depth below h
                w[(w>h)] = out_meta['nodata'] # Assign NoData to everywhere else
#                w[(w<0.)] = out_meta['nodata'] # Assign NoData to everywhere else
                
                # Now write out the FIM for this shed...
                w = w[0]                                                        
                shp=np.shape(w)                                
                bounds = rasterio.transform.array_bounds(shp[0],shp[1],out_transform) # window bounds in x-y space (west, south, east, north)                                
                col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[3]) # upper left row and column of window?                                
                
                row_min = np.int(row_min)                            
                col_min = np.int(col_min)
                row_max = np.int(row_min + shp[0])
                col_max = np.int(col_min + shp[1])    

                arr_w = np.empty([row_max-row_min, col_max-col_min], dtype=out_meta['dtype'])
                arr_w[:,:] = arr_fim[row_min:row_max, col_min:col_max]
#                
                inds_lt = np.where(arr_fim[row_min:row_max, col_min:col_max]<w)
                arr_w[inds_lt] = w[inds_lt]

                arr_fim[row_min:row_max, col_min:col_max] = arr_w # assign the FIM window for this catchment to the total array                    
    
    # Write out the final FIM grid...
    print('Writing final FIM .tif file...')   
    str_fim_path = str_hand_path[:-4]+'_fim_uniform_femamax.tif'
    with rasterio.open(str_fim_path, "w", **out_meta) as dest:
        dest.write(arr_fim, indexes=1)    
        
    return 

# ==========================================================================
#                 << MAIN INTERFACE FOR FUNCTIONS >>    
# ==========================================================================    
if __name__ == '__main__':    

    print('\n<<< Start >>>\r\n')
    start_time_0 = timeit.default_timer() 
    
    # << FILE PATHS >>
    # Input files...
    str_hand_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_hand.tif"
#    str_hand_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\baseline\121003hand.tif"
    
    str_sheds_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\experimental\hand_breached\sa_10m_dem_breach_p_uniform_diss_femamax.shp"
#    str_sheds_path = r"B:\Terrain\HAND_Experiments\SanAntonio_121003\baseline\sanantonio_gw_catchments_max30_bd.shp"
      
    # << CALL FUNCTION >>    
    fim_hand_poly(str_hand_path, str_sheds_path)
    
        
    print('\n<<< End >>>\r\n')
    print('Run time:  {}'.format(timeit.default_timer() - start_time_0))