# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:51:45 2016

@author: sam.lamont
"""
#import pprint # neat!

import subprocess
import timeit
import numpy as np
#from numpy import array
import os
#import ntpath
from math import atan, ceil
import sys
from math import isinf, sqrt
import rasterio
import rasterio.tools.mask
from rasterio.warp import transform

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
plt.style.use('ggplot')
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

#from affine import Affine
import pandas as pd
#from scipy import ndimage
#from shapely.geometry import Point, 
from shapely.geometry import shape, mapping, LineString, MultiLineString, Point
#from jenks import jenks
from PyQt4 import QtGui, QtCore

#import scipy.io as sio

import fiona

import gospatial as gs

#plt.ion()

#from shapely import speedups
#if speedups.available:
#    speedups.enable() # enable performance enhancements written in C


# ===============================================================================
#  Utility functions
# ===============================================================================
#def open_mat_files():
#    
#    # Loads a MATLAB file (.mat) using SciPy...
#    test = sio.loadmat(r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\field_data\FieldData\smith\AtGageAB_US.mat')
#    
#    arr_data = test['AtGageAB_US'] #numpy array
#    
#    return
#def rotate(x_center, x_orig, y_center, y_orig, theta):
#    x_rot = np.cos(theta)*(x_orig-x_center) - np.sin(theta)*(y_orig-y_center) + x_center
#    y_rot = np.sin(theta)*(x_orig-x_center) + np.cos(theta)*(y_orig-y_center) + y_center    
#    return x_rot, y_rot

# ===============================================================================
#  Callback for GoSpatial tool messages
#    If a callback is not provided, it will simply print the output stream.
#    A provided callback allows for custom processing of the output stream.
# ===============================================================================
def callback(s):
    try:
        # print("{}".format(s))
        if "%" in s:
            str_array = s.split(" ")
#            label = s.replace(str_array[len(str_array)-1], "")
            progress = int(str_array[len(str_array)-1].replace("%", "").strip())
            print("Progress: {}%".format(progress))
        else:
            if "error" in s.lower():
                print("ERROR: {}".format(s))
            else:
                print("{}".format(s))
    except:
        print(s)    

# ===============================================================================
#  Create weight file for TAuDEM D8 FAC
# ===============================================================================
def create_wg_from_streamlines(str_streamlines_path, str_dem_path, str_danglepts_path):
    
    # 1. Read in the streamlines vector layer
    # 2. Feature vertices to points (dangle)
    # 3. Convert these to raster (environment variables!)

    # NOTE:  Warning!! Lots of looping ahead
    # This entire function could probably use some re-design...although it does seem pretty fast so far

    lst_coords=[]
    lst_dangles=[]

    with fiona.open(str_streamlines_path) as lines:        
        
        streamlines_crs = lines.crs
        
        for line in lines:
            lst_coords.append(line['geometry']['coordinates'])
            
        mls = MultiLineString(lst_coords)
        
        for line in lines:            
            
            pt1 = Point(line['geometry']['coordinates'][0])
            pt2 = Point(line['geometry']['coordinates'][0])              
            
            k=0
            j=0
            for ls in mls:
                
                test1 = pt1.intersects(ls)
                test2 = pt2.intersects(ls)
                
                if test1:
                    j+=1
                if test2:
                    k+=1
                    
            if j==1:
                tpl = (pt1.xy[0][0], pt1.xy[1][0])                
                lst_dangles.append(tpl)
            
            if k==1:
                tpl = (pt2.xy[0][0], pt2.xy[1][0])
                lst_dangles.append(tpl)                            
       
    # NOTE:  Or just write it directly to a grid??
    lst_x, lst_y = zip(*lst_dangles)   # Coordinate handling is a bit of a mess here!!!
                       
    # Open DEM to copy metadata...
    with rasterio.open(str_dem_path) as ds_dem:            
        out_meta = ds_dem.meta.copy() 
        
        # Construct the output array...
        arr_danglepts=np.zeros([out_meta['height'], out_meta['width']], dtype=out_meta['dtype']) 
        
        tpl_pts = transform(streamlines_crs, out_meta['crs'], lst_x, lst_y)
        lst_dangles = zip(tpl_pts[0], tpl_pts[1])
        
        for coords in lst_dangles:            
            
            col, row = ~ds_dem.affine * (coords[0], coords[1]) # BUT you have to convert coordinates from hires to dem        
            arr_danglepts[row,col] = 1.
        
    # Now open the new grid using this metadata...    
    with rasterio.open(str_danglepts_path, "w", **out_meta) as dest:
        dest.write(arr_danglepts, indexes=1)
            
            
#    # Now write out the dangle points into a new shapefile...
#    print('pause')
#    test_schema = {'geometry': 'Point', 'properties': {'linkno': 'int'}} 
#    
#    str_dangles_path= r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\dr_nhd_hires_dangles.shp'
#       
#    print('Building and Writing Dangle Point File...')
#    with fiona.open(str_dangles_path, 'w', driver='ESRI Shapefile', crs=streamlines_crs, schema=test_schema) as dangles:
#
#        # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
#        for coords in lst_dangles:
#            pt = {'type': 'Point', 'coordinates':coords}
#            prop = {'linkno': 1}
#            dangles.write({'geometry':pt, 'properties':prop})         
        
    return
    
# ===============================================================================
#  Mega-function for processing a raw DEM
#   1. Breaching and filling
#   2. TauDEM functions
# ===============================================================================        
def preprocess_dem(str_dem_path, str_danglepts_path):
    try:

        distmeth = 'v'
        statmeth = 'ave'
        inputProc = str(4)

        
        # ========== << 1. Breach Depressions with GoSpatial/Whitebox Tool >> ==========
        '''
        This tool is used to remove the sinks (i.e. topographic depressions and flat areas) from digital elevation models (DEMs) using a highly efficient and flexible breaching, or carving, method.
        Arg Name: InputDEM, type: string, Description: The input DEM name with file extension
        Arg Name: OutputFile, type: string, Description: The output filename with file extension
        Arg Name: MaxDepth, type: float64, Description: The maximum breach channel depth (-1 to ignore)
        Arg Name: MaxLength, type: int, Description: The maximum length of a breach channel (-1 to ignore)
        Arg Name: ConstrainedBreaching, type: bool, Description: Use constrained breaching?
        Arg Name: SubsequentFilling, type: bool, Description: Perform post-breach filling?    
        '''
        
        # Get the gospatial version number
        print(gs.version())   
#        print(gs.tool_help("BreachDepressions"))
        gs.set_working_dir(r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach')
    
        # Run the BreachDepressions tool, specifying the arguments.
        name = "BreachDepressions"
        args = ['dr3m_dem.tif', 'dr3m_breach.tif', '-1', '-1', 'False', 'True']        
        
        # Run the tool and check the return value
        ret = gs.run_tool(name, args, callback)
        if ret != 0:
            print("ERROR: return value={}".format(ret)) 
#            
#        # =========== << 2. Fill Pits with TauDEM (OPTION?) ==========
#        fel = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\dr3m_fel.tif'
#        cmd = 'mpiexec -n ' + inputProc + ' pitremove -z ' + '"' + str_dem_path + '"' + ' -fel ' + '"' + fel + '"'
##        if considering4way == 'true':
##            cmd = cmd + ' -4way '
##        if arcpy.Exists(maskgrid):
##            cmd = cmd + ' -depmask ' + '"' + mkgr + '"'
##        if arcpy.Exists(maskgrid) and considering4way == 'true':
##            cmd = cmd + ' -depmask ' + '"' + mkgr + '"' + ' -4way '
#        
##        arcpy.AddMessage("\nCommand Line: "+cmd)
#        print('Running TauDEM PitRemove...')
#        os.system(cmd)
#        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#        
#        message = "\n"
#        for line in process.stdout.readlines():
#            if isinstance(line, bytes):	    # true in Python 3
#                line = line.decode()
#            message = message + line
#        print(message)
            
#        # ==============  << 3. D8 FDR with TauDEM >> ================
        fel = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach.tif'
        p = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_p.tif'
        sd8 = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_sd8.tif'
                
        cmd = 'mpiexec -n ' + inputProc + ' D8FlowDir -fel ' + '"' + fel + '"' + ' -p ' + '"' + p + '"' + \
              ' -sd8 ' + '"' + sd8 + '"'
        
        # Submit command to operating system
        print('Running TauDEM D8 Flow Direction...')
        os.system(cmd)
        
        # Capture the contents of shell command and print it to the arcgis dialog box
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        
        # Get some feedback from the process to print out...
        message = "\n"
        for line in process.stdout.readlines():
            line = line.decode()
    #            if isinstance(line, bytes):	   # true in Python 3
    #                line = line.decode()
            message = message + line        
        print(message)            
            
        # ============= << 4. D8 FAC with TauDEM >> ================
#        wtgr = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_wg_rc.tif'
        wtgr = str_danglepts_path
        ad8 = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_ad8.tif'
        
        cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -wg ' + '"' + wtgr + '"' + ' -nc '
#        if arcpy.Exists(ogrfile):
#            cmd = cmd + ' -o ' + '"' + shfl + '"'
#        if arcpy.Exists(weightgrid):
#            cmd = cmd + ' -wg ' + '"' + wtgr + '"'
#        if edgecontamination == 'false':
#            cmd = cmd + ' -nc '
            
#        print("\nCommand Line:\n" + cmd)
        
        # Submit command to operating system
        print('Running TauDEM D8 FAC...')
        os.system(cmd)
        # Capture the contents of shell command and print it to the arcgis dialog box
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        
        message = "\n"
        for line in process.stdout.readlines():
            if isinstance(line, bytes):	    # true in Python 3
                line = line.decode()
            message = message + line
        print(message)            
            
        # ============= << 5. Dinf with TauDEM >> =============
        slp = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_slp.tif'
        ang = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_ang.tif'
                
        print('Running TauDEM Dinfinity...')        
        cmd = 'mpiexec -n ' + inputProc + ' DinfFlowDir -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + \
              ' -slp ' + '"' + slp + '"'
        
        # Submit command to operating system
        os.system(cmd)
        
        # Capture the contents of shell command and print it to the arcgis dialog box
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        
        # Get some feedback from the process to print out...
        message = "\n"
        for line in process.stdout.readlines():
            line = line.decode()
    #            if isinstance(line, bytes):	   # true in Python 3
    #                line = line.decode()
            message = message + line        
        print(message)            
        
        # ============= << 6. DinfDistanceDown (HAND) with TauDEM >> =============
        dd = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_breach_dd.tif'
        
        print('Running TauDEM Dinf Distance Down...')
        cmd = 'mpiexec -n ' + inputProc + ' DinfDistDown -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + \
              ' -src ' + '"' + ad8 + '"' + ' -dd ' + '"' + dd + '"' + ' -m ' + statmeth + ' ' + distmeth
    
        # Submit command to operating system
        os.system(cmd)
        
        # Get some feedback from the process to print out...
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) 
        
        message = "\n"
        for line in process.stdout.readlines():
            line = line.decode()
    #            if isinstance(line, bytes):	   # true in Python 3
    #                line = line.decode()
            message = message + line
            
        print(message)

    
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise        
    
    return


# ===============================================================================
#  NOTE:  Much of the following function was copied directly from Tarboton's 
#         TauDEM ArcToolbox .py files
# ===============================================================================
def taudem_dinfdd(fel, src, ang, dd, slp):        
    
    distmeth = 'v'
    statmeth = 'ave'
    inputProc = str(4)
    
    # ===========================================================================
    #  Run Dinf...
    print('Running TauDEM Dinfinity...')        
    cmd = 'mpiexec -n ' + inputProc + ' DinfFlowDir -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + \
          ' -slp ' + '"' + slp + '"'
    
    # Submit command to operating system
    os.system(cmd)
    
    # Capture the contents of shell command and print it to the arcgis dialog box
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    
    # Get some feedback from the process to print out...
    message = "\n"
    for line in process.stdout.readlines():
        line = line.decode()
#            if isinstance(line, bytes):	   # true in Python 3
#                line = line.decode()
        message = message + line        
    print(message)
    # ===========================================================================
    
    # ===========================================================================
    # Now run DinfDD...        NOTE: Edge contamination in on by default
    print('Running TauDEM Dinf Distance Down...')
    cmd = 'mpiexec -n ' + inputProc + ' DinfDistDown -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + \
          ' -src ' + '"' + src + '"' + ' -dd ' + '"' + dd + '"' + ' -m ' + statmeth + ' ' + distmeth

    # Submit command to operating system
    os.system(cmd)
    
    # Get some feedback from the process to print out...
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) 
    
    message = "\n"
    for line in process.stdout.readlines():
        line = line.decode()
#            if isinstance(line, bytes):	   # true in Python 3
#                line = line.decode()
        message = message + line
        
    print(message)
    return
    # ===========================================================================
       
def find_peaks_from_series(srs_in, param_window, param_thresh):

#    srs_rlocal_big = pd.rolling_mean(srs_in, window=param_big_window, center=True, min_periods=1)
#    srs_rlocal_big = pd.rolling_median(srs_in, window=param_big_window, center=True)
#    srs_rlocal_big = pd.rolling_window(srs_in, window=param_big_window, win_type='parzen', center=True, min_periods=1, mean=True)
#    srs_rlocal_big = pd.rolling_quantile(srs_in, window=param_window, quantile=0.35, center=True, min_periods=1)
#    srs_rlocal_big = srs_in.rolling(min_periods=1, window=param_window, center=True).quantile(quantile=0.35)
      
    first_diff = srs_in.diff()                    
    first_diff[first_diff<0]=0 # assign 0 to fall
    first_diff[first_diff>0]=1 # assign 1 to rise
    
    sh_first_diff = first_diff.shift(periods=-1)                    
    sub_diff = first_diff-sh_first_diff    
    
    try:
        srs_diff_vals = srs_in.loc[sub_diff.loc[sub_diff==1].index]   
#        srs_rlocal_big_vals = srs_rlocal_big.loc[sub_diff.loc[sub_diff==1].index]
#
#        # If values are at least some value times greater than the local mean value...
#        srs_out = srs_in.loc[srs_diff_vals[srs_diff_vals/srs_rlocal_big_vals > param_thresh].index] 

    except:
        pass
    
    
    return srs_diff_vals #, srs_rlocal_big #, srs_rlocal_small     
    
def interpolate(arr_in, ind_val):
    
    if ind_val == np.ceil(ind_val):
        out_val = arr_in[int(np.ceil(ind_val))]
    else:
        out_val = arr_in[np.int(ind_val)] + (ind_val - np.int(ind_val))*(arr_in[int(np.ceil(ind_val))]-arr_in[np.int(ind_val)]) # it will always be divided by 1            
    
    return out_val 

# Count the number of features in a vector file...    
def get_feature_count(str_shp_path):
    
    with fiona.open(str(str_shp_path), 'r') as features:
        i_count=len(features)
        
    return i_count

# ===============================================================================
#  Calculate floodplain width based on a HAND slice and stream line buffers of 
#   entire reach
# =============================================================================== 
def floodplain_width_reach_buffers_po(str_streamlines_path, str_fp_path, str_reachid, cell_size):
    
#    lst_output=[]
    print('Channel width from bank pixels...')
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_fpwidth.shp'
    str_outpath = head + '/' + str_outname 
    
    # Open the bankpts layer...
    with rasterio.open(str(str_fp_path)) as ds_fp:
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
            # Get the crs...
            streamlines_crs = streamlines.crs                
#                streamlines_schema = streamlines.schema.copy()
#                streamlines_schema['properties']['fp_width'] = 'float'
            schema_output = {'geometry': 'LineString', 'properties': {'width_total':'float', 'width_left':'float', 'width_right':'float'}}
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
          
                for line in streamlines:                    
                    
                    lst_tally=[]
                    
                    # Buffer each feature...
                    ls = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                                    
#                    print('linkno: {}'.format(line['properties'][str_reachid]))
                    
                    if line['properties'][str_reachid] <> 1115:
                        continue
                                   
                    # Successive buffer-mask operations to count bank pixels at certain intervals
#                        prev_val=0
                    lst_buff=range(cell_size,80,cell_size) # use maximum expected FP width here
                    for buff_dist in lst_buff:
                        
                        ls_offset_left = ls.parallel_offset(buff_dist, 'left')
                        ls_offset_rt = ls.parallel_offset(buff_dist, 'right')                 
                                    
                        try:
                            out_left, out_transform = rasterio.tools.mask.mask(ds_fp, [mapping(ls_offset_left)], crop=True)   
                            out_rt, out_transform = rasterio.tools.mask.mask(ds_fp, [mapping(ls_offset_rt)], crop=True)
                        except:
                            continue
                        
                        num_pixels_left = len(out_left[out_left==1])
                        num_pixels_rt = len(out_rt[out_rt==1])                      
                        
                        tpl_out = line['properties'][str_reachid], buff_dist, num_pixels_left, num_pixels_rt
                        lst_tally.append(tpl_out)                    
                        df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','interval_left','interval_rt'])
                        
                        
                    # Calculate weighted average                     
                    # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                    try:                        
                        weighted_avg_left=0 
                        weighted_avg_rt=0
                        n_top=2
                        for tpl in df_tally.nsmallest(n_top, 'interval_left').iloc[0:2].itertuples():
                            weighted_avg_left += tpl.buffer*(np.float(tpl.interval_left)/np.float(df_tally.nsmallest(n_top, 'interval_left').iloc[0:2].sum().interval_left))

                        for tpl in df_tally.nsmallest(n_top, 'interval_rt').iloc[0:2].itertuples():
                            weighted_avg_rt += tpl.buffer*(np.float(tpl.interval_rt)/np.float(df_tally.nsmallest(n_top, 'interval_rt').iloc[0:2].sum().interval_rt))
                            
#                                weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
                    except:
#                            weighted_avg=-9999.
                        continue
                    
                    # Write to an output file here...
                    output.write({'properties':{'width_total': weighted_avg_left+weighted_avg_rt, 'width_left': weighted_avg_left, 'width_right': weighted_avg_rt}, 'geometry':mapping(ls)})
                        
#                    # To save the output...    
#                    tpl_output = (line['properties'][str_reachid],weighted_avg)                    
#                    lst_output.append(tpl_output)
                                
                    # FOR TESTING...
    #                if line['properties'][str_reachid] == 13:                                        
    #                    print('pause')     
#                    df_linkno.plot.bar(x='buffer', y='interval', title='linkno: {} | width: {}'.format(line['properties'][str_reachid], weighted_avg))                    
                    sys.exit()
    
#    # Write to a csv...
#    pd_output = pd.DataFrame(lst_output)    
#    pd_output.to_csv('/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/channel_width_test.csv')
    
    
    return
    
# ===============================================================================
#  Calculate floodplain width based on a HAND slice and stream line buffers of 
#   entire reach
# =============================================================================== 
def floodplain_width_reach_buffers(str_streamlines_path, str_fp_path, str_reachid, cell_size):
    
#    lst_output=[]
    print('Channel width from bank pixels...')
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_fpwidth.shp'
    str_outpath = head + '/' + str_outname 
    
    # Open the bankpts layer...
    with rasterio.open(str(str_fp_path)) as ds_fp:
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
            # Get the crs...
            streamlines_crs = streamlines.crs                
            streamlines_schema = streamlines.schema.copy()
            streamlines_schema['properties']['chan_width'] = 'float'
            
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', streamlines_schema, streamlines_crs) as output:
          
                for line in streamlines:                    
                    
                    lst_tally=[]
                    
                    # Buffer each feature...
                    geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                                    
#                    print('linkno: {}'.format(line['properties'][str_reachid]))
                    
                    if line['properties'][str_reachid] <> 1115:
                        continue
                                   
                    # Successive buffer-mask operations to count bank pixels at certain intervals
                    prev_val=0
                    lst_buff=range(cell_size,80,cell_size) # use maximum expected FP width here
                    for buff_dist in lst_buff:
                        
                        geom_buff = geom.buffer(buff_dist,cap_style=2)
                        buff = mapping(geom_buff)           
                        
    #                    # Write out buffer polygon(s)...
    #                    with fiona.open('/home/sam.lamont/USGSChannelFPMetrics/testing/buffer.shp','w','ESRI Shapefile', schema) as buff_out:                     
    #                        buff_out.write({'properties': {'name': 'hey'}, 'geometry': buff})                       
    #                    sys.exit()                    
                                    
                        # Mask the bankpts file for each feature...
                        out_image, out_transform = rasterio.tools.mask.mask(ds_fp, [buff], crop=True)                    
                        num_pixels = len(out_image[out_image==1])
                        
                        # You want the number of pixels gained by each interval...                    
                        interval = num_pixels-prev_val                    
                        prev_val = num_pixels
                        
                        tpl_out = line['properties'][str_reachid], buff_dist, num_pixels, interval
                        lst_tally.append(tpl_out)                    
                        df_linkno = pd.DataFrame(lst_tally, columns=['linkno','buffer','cumulative','interval'])
                        
                        
                        
                        # Could say stop when there is no increase?
    #                    if num_pixels == 0:
    #                        break                    
                    
#                        df_linkno['interval'].plot()                        
                    
                    # Avoid division by zero                     
                    if df_linkno.interval.sum() == 0:
                        continue                    
                    
                    # Calculate weighted average                     
                    # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                    try:                        
                        weighted_avg=0 
                        n_top=2
                        for tpl in df_linkno.nlargest(n_top, 'interval').itertuples():
                            weighted_avg += tpl.buffer*(np.float(tpl.interval)/np.float(df_linkno.nlargest(n_top, 'interval').sum().interval))
                            
                        weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
                    except:
                        continue
                    
                    # Or write to the streamline shapefile attribute table here...
                    line['properties']['chan_width'] = weighted_avg
                    output.write({'properties':line['properties'], 'geometry':mapping(shape(line['geometry']))})
                        
#                    # To save the output...    
#                    tpl_output = (line['properties'][str_reachid],weighted_avg)                    
#                    lst_output.append(tpl_output)
                                
                    # FOR TESTING...
    #                if line['properties'][str_reachid] == 13:                                        
    #                    print('pause')     
#                    df_linkno.plot.bar(x='buffer', y='interval', title='linkno: {} | width: {}'.format(line['properties'][str_reachid], weighted_avg))                    
#                    sys.exit()
    
#    # Write to a csv...
#    pd_output = pd.DataFrame(lst_output)    
#    pd_output.to_csv('/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/channel_width_test.csv')
    
    
    return
    
# ===============================================================================
#  Calculate channel width based on bank pixels and stream line parallel offsets, 
#  potentially subdividing using Xn's
#  NOTE: Make this a generic metric calculator by segment?  (ie., sinuosity, etc)
# =============================================================================== 
def channel_width_bankpixels_segments_po(df_coords, str_streamlines_path, str_bankpixels_path, str_reachid, cell_size):
    
    print('Channel width from bank pixels -- segmented...')
    
#    lst_output=[]
    
#    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_chwidth_po.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')
    
#    schema_buff = {'geometry': 'Polygon', 'properties': {'buff': 'str'}}
    schema_output = {'geometry': 'LineString', 'properties': {'width_total':'float', 'width_left':'float', 'width_right':'float', 'sinuosity':'float'}}
    
    # Open the bankpts layer...
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:    
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
#            progBar.setRange(0, len(streamlines)) 
            
            # Get the crs...
            streamlines_crs = streamlines.crs                
#                streamlines_schema = streamlines.schema.copy()
#                streamlines_schema['properties']['chan_width'] = 'float'
#                streamlines_schema['properties']['sinuosity'] = 'float'
            
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
                
                for i_linkno, df_linkno in gp_coords:
                    
#                    progBar.setValue(j)
#                    j+=1                    
            
                    i_linkno = int(i_linkno)
                    
#                        if i_linkno > 0:
#                            break
                    
                    print('linkno:  {}'.format(i_linkno))
          
                    # Now loop over df_linkno
                    # Set up index array for looping...
                    arr_ind = np.arange(0, len(df_linkno.index)+1, 10) # NOTE: Change the step for resolution?
                    
                    for i, indx in enumerate(arr_ind):
                        
#                            if indx == 110:
#                                print('pause')
                        
                        if i>0: indx = indx-1
                        try:
                            arr_x = df_linkno.x.iloc[indx:arr_ind[i+1]].values
                            arr_y = df_linkno.y.iloc[indx:arr_ind[i+1]].values            
                        except:
                            break

                        # Calculate straight line distance...
                        dist_sl = np.sqrt((arr_x[0] - arr_x[len(arr_x)-1])**2 + (arr_y[0] - arr_y[len(arr_y)-1])**2)                     
              
                        # Create a line segment from endpts in df_linkno...
                        ls = LineString(zip(arr_x, arr_y))
                        
                        dist = ls.length
                        
                        sinuosity = dist/dist_sl # ratio of sinuous length to straight line length
                        
                        lst_tally=[]
                        
                        # Buffer each feature...
#                            geom = shape(ls)   # Convert to shapely geometry to operate on it
                               
                        # Successive buffer-mask operations to count bank pixels at certain intervals
#                            prev_val_left=0
#                            prev_val_rt=0
                        lst_buff=range(cell_size,30,cell_size)
                        for buff_dist in lst_buff:                                                            
                            
                            ls_offset_left = ls.parallel_offset(buff_dist, 'left')
                            ls_offset_rt = ls.parallel_offset(buff_dist, 'right')
#                                ls_offset_left
                            
#                            # Write out buffer polygon(s)...NOTE:  Could write out ind
#                            with fiona.open(r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\buffer.shp','w','ESRI Shapefile', schema_buff) as buff_out:                     
#                                buff_out.write({'properties': {'buff': 'mmmm'}, 'geometry': buff})                       
#                            sys.exit()                    
                                        
                            # Mask the bankpts file for each feature...
#                                out_image, out_transform = rasterio.tools.mask.mask(ds_bankpixels, [buff], crop=True) 
                                        
#                                print(indx, buff_dist)
                            
#                                if (indx==79) and (buff_dist==24):
#                                    print('pause')
#                                    continue
                            try:
                                out_left, out_transform = rasterio.tools.mask.mask(ds_bankpixels, [mapping(ls_offset_left)], crop=True)   
                                out_rt, out_transform = rasterio.tools.mask.mask(ds_bankpixels, [mapping(ls_offset_rt)], crop=True)
                            except:
                                continue
                            
                            num_pixels_left = len(out_left[out_left>0])
                            num_pixels_rt = len(out_rt[out_rt>0])
                            
#                                # You want the number of pixels gained by each interval...                    
#                                interval_left = num_pixels_left-prev_val_left                    
#                                prev_val_left = num_pixels_left
#                                
#                                interval_rt = num_pixels_rt-prev_val_rt                    
#                                prev_val_rt = num_pixels_rt                                
                            
                            tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
                            lst_tally.append(tpl_out)                    
                            df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','interval_left','interval_rt'])
                   
        
                        # Avoid division by zero                     
#                        if df_tally.interval.sum() == 0:
#                            continue                    
                        
                        # Calculate weighted average                     
                        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                        try:                        
                            weighted_avg_left=0 
                            weighted_avg_rt=0
                            n_top=2
                            for tpl in df_tally.nlargest(n_top, 'interval_left').iloc[0:2].itertuples():
                                weighted_avg_left += tpl.buffer*(np.float(tpl.interval_left)/np.float(df_tally.nlargest(n_top, 'interval_left').iloc[0:2].sum().interval_left))

                            for tpl in df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].itertuples():
                                weighted_avg_rt += tpl.buffer*(np.float(tpl.interval_rt)/np.float(df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].sum().interval_rt))
                                
#                                weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
                        except:
#                            weighted_avg=-9999.
                            continue
                        
                        # Write to an output file here...
                        output.write({'properties':{'width_total': weighted_avg_left+weighted_avg_rt, 'width_left': weighted_avg_left, 'width_right': weighted_avg_rt, 'sinuosity': sinuosity}, 'geometry':mapping(ls)})
                                                
#                    if i_linkno > 10:
#                        break
                    
#                    # To save the output...    
#                    tpl_output = (line['properties'][str_reachid],weighted_avg)                    
#                    lst_output.append(tpl_output)
                                
                    # FOR TESTING...
    #                if line['properties'][str_reachid] == 13:                                        
    #                    print('pause')     
#                    df_linkno.plot.bar(x='buffer', y='interval', title='linkno: {} | width: {}'.format(line['properties'][str_reachid], weighted_avg))                    
#                    sys.exit()
    
#    # Write to a csv...
#    pd_output = pd.DataFrame(lst_output)    
#    pd_output.to_csv('/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/channel_width_test.csv')
     
    return

# ===============================================================================
#  Calculate channel width based on bank pixels and stream line buffers, 
#  potentially subdividing using Xn's
#  NOTE: Make this a generic metric calculator by segment?  (ie., sinuosity, etc)
# =============================================================================== 
def channel_width_bankpixels_segments(df_coords, str_streamlines_path, str_bankpixels_path, str_reachid, cell_size):
    
    print('Channel width from bank pixels -- segmented...')
    
#    lst_output=[]
    
#    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_chwidth_buff.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')
    
#    schema_buff = {'geometry': 'Polygon', 'properties': {'buff': 'str'}}
    schema_output = {'geometry': 'LineString', 'properties': {'width':'float', 'sinuosity':'float'}}
    
    # Open the bankpts layer...
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:    
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
#            progBar.setRange(0, len(streamlines)) 
            
            # Get the crs...
            streamlines_crs = streamlines.crs                
            streamlines_schema = streamlines.schema.copy()
            streamlines_schema['properties']['chan_width'] = 'float'
            streamlines_schema['properties']['sinuosity'] = 'float'
            
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
                
                for i_linkno, df_linkno in gp_coords:
                    
#                    progBar.setValue(j)
#                    j+=1                    
            
                    i_linkno = int(i_linkno)
                    
                    print('linkno:  {}'.format(i_linkno))
          
                    # Now loop over df_linkno
                    # Set up index array for looping...
                    arr_ind = np.arange(0, len(df_linkno.index)+1, 10) # NOTE: Change the step for resolution?
                    
                    for i, indx in enumerate(arr_ind):
                        if i>0: indx = indx-1
                        try:
                            arr_x = df_linkno.x.iloc[indx:arr_ind[i+1]].values
                            arr_y = df_linkno.y.iloc[indx:arr_ind[i+1]].values            
                        except:
                            break

                        # Calculate straight line distance...
                        dist_sl = np.sqrt((arr_x[0] - arr_x[len(arr_x)-1])**2 + (arr_y[0] - arr_y[len(arr_y)-1])**2)                     
              
                        # Create a line segment from endpts in df_linkno...
                        ls = LineString(zip(arr_x, arr_y))
                        
                        dist = ls.length
                        
                        sinuosity = dist/dist_sl # ratio of sinuous length to straight line length
                        
                        lst_tally=[]
                        
                        # Buffer each feature...
                        geom = shape(ls)   # Convert to shapely geometry to operate on it
                               
                        # Successive buffer-mask operations to count bank pixels at certain intervals
                        prev_val=0
                        lst_buff=range(cell_size,30,cell_size)
                        for buff_dist in lst_buff:
                            
                            geom_buff = geom.buffer(buff_dist,cap_style=2)
                            buff = mapping(geom_buff)     
                            
#                            # Write out buffer polygon(s)...NOTE:  Could write out ind
#                            with fiona.open(r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\buffer.shp','w','ESRI Shapefile', schema_buff) as buff_out:                     
#                                buff_out.write({'properties': {'buff': 'mmmm'}, 'geometry': buff})                       
#                            sys.exit()                    
                                        
                            # Mask the bankpts file for each feature...
                            out_image, out_transform = rasterio.tools.mask.mask(ds_bankpixels, [buff], crop=True)                                       
                            
                            num_pixels = len(out_image[out_image>0])
                            
                            # You want the number of pixels gained by each interval...                    
                            interval = num_pixels-prev_val                    
                            prev_val = num_pixels
                            
                            tpl_out = i_linkno, buff_dist, num_pixels, interval
                            lst_tally.append(tpl_out)                    
                            df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','cumulative','interval'])
                   
        
                        # Avoid division by zero                     
#                        if df_tally.interval.sum() == 0:
#                            continue                    
                        
                        # Calculate weighted average                     
                        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                        try:                        
                            weighted_avg=0 
                            n_top=2
                            for tpl in df_tally.nlargest(n_top, 'interval').itertuples():
                                weighted_avg += tpl.buffer*(np.float(tpl.interval)/np.float(df_tally.nlargest(n_top, 'interval').sum().interval))
                                
                            weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
                        except:
#                            weighted_avg=-9999.
                            continue
                        
                        # Write to an output file here...
                        output.write({'properties':{'width': weighted_avg, 'sinuosity': sinuosity}, 'geometry':mapping(ls)})
                                                
#                    if i_linkno > 10:
#                        break
                    
#                    # To save the output...    
#                    tpl_output = (line['properties'][str_reachid],weighted_avg)                    
#                    lst_output.append(tpl_output)
                                
                    # FOR TESTING...
    #                if line['properties'][str_reachid] == 13:                                        
    #                    print('pause')     
#                    df_linkno.plot.bar(x='buffer', y='interval', title='linkno: {} | width: {}'.format(line['properties'][str_reachid], weighted_avg))                    
#                    sys.exit()
    
#    # Write to a csv...
#    pd_output = pd.DataFrame(lst_output)    
#    pd_output.to_csv('/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/channel_width_test.csv')
     
    return
    
# ===============================================================================
#  Calculate channel width based on bank pixels and stream line buffers, 
#  potentially subdividing using Xn's
# =============================================================================== 
def channel_width_bankpixels(str_streamlines_path, str_bankpixels_path, str_reachid, cell_size, self):
    
#    lst_output=[]
    print('Channel width from bank pixels...')
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_chwidth.shp'
    str_outpath = head + '/' + str_outname 
    
    j=0   
    progBar = self.progressBar
    progBar.setVisible(True)    
    
    # Open the bankpts layer...
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
            progBar.setRange(0, len(streamlines))        
       
            # Get the crs...
            streamlines_crs = streamlines.crs                
            streamlines_schema = streamlines.schema.copy()
            streamlines_schema['properties']['chan_width'] = 'float'
            
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', streamlines_schema, streamlines_crs) as output:
          
                for line in streamlines:    
                    
                    progBar.setValue(j)
                    j+=1
                    
                    lst_tally=[]
                    
                    # Buffer each feature...
                    geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                                    
#                    print('linkno: {}'.format(line['properties'][str_reachid]))
                    
#                    if line['properties'][str_reachid] <> 7:
#                        continue
                                   
                    # Successive buffer-mask operations to count bank pixels at certain intervals
                    prev_val=0
                    lst_buff=range(cell_size,30,cell_size)
                    for buff_dist in lst_buff:
                        
                        geom_buff = geom.buffer(buff_dist,cap_style=2)
                        buff = mapping(geom_buff)           
                        
    #                    # Write out buffer polygon(s)...
    #                    with fiona.open('/home/sam.lamont/USGSChannelFPMetrics/testing/buffer.shp','w','ESRI Shapefile', schema) as buff_out:                     
    #                        buff_out.write({'properties': {'name': 'hey'}, 'geometry': buff})                       
    #                    sys.exit()                    
                                    
                        # Mask the bankpts file for each feature...
                        out_image, out_transform = rasterio.tools.mask.mask(ds_bankpixels, [buff], crop=True)                    
                        num_pixels = len(out_image[out_image>0])
                        
                        # You want the number of pixels gained by each interval...                    
                        interval = num_pixels-prev_val                    
                        prev_val = num_pixels
                        
                        tpl_out = line['properties'][str_reachid], buff_dist, num_pixels, interval
                        lst_tally.append(tpl_out)                    
                        df_linkno = pd.DataFrame(lst_tally, columns=['linkno','buffer','cumulative','interval'])
                        
                        # Could say stop when there is no increase?
    #                    if num_pixels == 0:
    #                        break                    
    
                    # Avoid division by zero                     
                    if df_linkno.interval.sum() == 0:
                        continue                    
                    
                    # Calculate weighted average                     
                    # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                    try:                        
                        weighted_avg=0 
                        n_top=2
                        for tpl in df_linkno.nlargest(n_top, 'interval').itertuples():
                            weighted_avg += tpl.buffer*(np.float(tpl.interval)/np.float(df_linkno.nlargest(n_top, 'interval').sum().interval))
                            
                        weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
                    except:
                        continue
                    
                    # Or write to the streamline shapefile attribute table here...
                    line['properties']['chan_width'] = weighted_avg
                    output.write({'properties':line['properties'], 'geometry':mapping(shape(line['geometry']))})
                        
#                    # To save the output...    
#                    tpl_output = (line['properties'][str_reachid],weighted_avg)                    
#                    lst_output.append(tpl_output)
                                
                    # FOR TESTING...
    #                if line['properties'][str_reachid] == 13:                                        
    #                    print('pause')     
#                    df_linkno.plot.bar(x='buffer', y='interval', title='linkno: {} | width: {}'.format(line['properties'][str_reachid], weighted_avg))                    
#                    sys.exit()
    
#    # Write to a csv...
#    pd_output = pd.DataFrame(lst_output)    
#    pd_output.to_csv('/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/channel_width_test.csv')
    
    
    return
    
 # ===============================================================================
#  Analyze DEM in vertical slices successive buffers stream reaches
# ===============================================================================    
def analyze_hand_reach_buffers(str_net_path, str_hand_path, str_reachid):    
    
    cell_size = 3
    
    # Open the bankpts layer...
    with rasterio.open(str(str_hand_path)) as ds_hand:
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_net_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
        
            for line in streamlines:    
                
                test_id=62
                if line['properties'][str_reachid] <> test_id:
                    continue
                
                lst_count=[]
    
                # Buffer each feature...
                geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                
                if geom.length < 9:
                    continue
                
                fig = plt.figure()            
                ax = fig.add_subplot(111) 
                
                print('linkno: {}'.format(line['properties'][str_reachid]))
                
                # Successive buffer-mask operations to count hand pixels at certain intervals
                arr_buffs=np.arange(1,30,1)
                
                for buff_dist in arr_buffs:
                    
                    geom_buff = geom.buffer(buff_dist,cap_style=2)
                    buff = mapping(geom_buff)           
                                
                    # Mask the bankpts file for each feature...
                    w, out_transform = rasterio.tools.mask.mask(ds_hand, [buff], crop=True)
                    
                    num_pixels = w[w>=0].size
                    lst_count.append(num_pixels)    
    
                df_steps = pd.DataFrame({'count':lst_count, 'buffer':arr_buffs})
#                df_steps['count_cs'] = df_steps['count'].cumsum()
#                
#                # Slope of cumulative curve...
#                i_offset=1
#                df_steps['count_slp'] = df_steps['count'].pct_change(i_offset) 
#                
#                # Delta slope...
#                df_steps['count_delta'] = df_steps['count_slp'].diff(i_offset)
                
                # Index of max slope value...
#                test = df_steps['count_slp'].idxmax()
                
                # Height value where slope is greatest...
#                val_height = df_steps['height'].iloc[df_steps['count_delta'].idxmax()-i_offset]
#                val_cs = df_steps['count_cs'].iloc[df_steps['count_delta'].idxmax()-i_offset]
                                
#                print('xn: {} | val_height: {} | val_cs: {}\t| w_max: {:3.1f}'.format(cntr, val_height, val_cs, w_max))    
                
                # =====================================================================              
                ax.plot(df_steps['count'], df_steps['buffer'])   
#                ax.plot(val_cs, val_height, marker='*', linestyle='') 
                ax.set_title('ARCID: {}'.format(line['properties'][str_reachid]), fontsize=11) 
                # =====================================================================    
                
#                if line['properties'][str_reachid] == 5:
#                    break
    
    return
# ===============================================================================
#  Analyze DEM in vertical slices using subcatchment polygons
# =============================================================================== 
def analyze_hand_sheds(str_sheds_path, str_grid_path, i_interval):
    
#    buff_dist=10
    
    str_reachid='GRIDCODE'
           
    # Open the dem...
    with rasterio.open(str(str_grid_path)) as ds_grid:
        
        # Open the xn layer...
        with fiona.open(np.str(str_sheds_path), 'r') as sheds: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
            # Get the crs...
    #        streamlines_crs = streamlines.crs          

#            fig = plt.figure()            
#            ax = fig.add_subplot(111)  
                
            cntr = 0
            prev_linkno = 0
            for shed in sheds:    

                test_id=49
                if shed['properties'][str_reachid] <> test_id:
                    continue
                
                if shed['properties'][str_reachid] <> prev_linkno:
                    fig = plt.figure()            
                    ax = fig.add_subplot(111) 
                    cntr = 0
                
                # Buffer each feature...
#                geom = shape(shed['geometry'])   # Convert to shapely geometry to operate on it
                
#                print('linkno: {}'.format(line['properties'][str_reachid]))
       
                                     
                # Mask the bankpts file for each feature...
                w, out_transform = rasterio.tools.mask.mask(ds_grid, [shed['geometry']], crop=True)
                
                if np.size(w) < 9:
                    continue
                
#                w[w<0]=9999. # handle no data vals?
                
                # Normalize to zero...
#                min_elev = w[w>0].min()             
#                w = w - min_elev
                
                # Find max value in window...
                w_max = w.max()
                
#                # ===================================================================
#                # << 3D plot >> 
#                fig = plt.figure()            
#                ax = fig.add_subplot(111, projection='3d')
#                Y=np.arange(0,np.shape(arr_norm)[1],1)
#                X=np.arange(0,np.shape(arr_norm)[2],1)
#                X, Y = np.meshgrid(X, Y)
#                
#                arr_norm[arr_norm<0] = -0.5
#                out_image[out_image<0] = -0.5
#               
#                surf = ax.plot_surface(X, Y, out_image[0], cmap=cm.jet, linewidth=0, antialiased=True, rstride=1, cstride=1)
#                
#                # Add a color bar which maps values to colors.
#                fig.colorbar(surf, shrink=0.5, aspect=5)            
#                # ===================================================================
                                
                # Now perform vertical slices...
#                i_interval = 0.1 # m
#                i_rng = w_max # m
                i_rng=10
                arr_slices = np.arange(i_interval, i_rng, i_interval)    
                
                lst_count=[]
                # List comprehension here instead??
                for i_step in arr_slices:                    
                    lst_count.append(w[(w<i_step) & (w>=0)].size) # number of pixels greater than or equal to zero and less than the height interval

                df_steps = pd.DataFrame({'count':lst_count, 'height':arr_slices})
                
                # Slope of count...
                df_steps['count_diff'] = df_steps['count'].diff()
                
                # Slope of slope of count...
                df_steps['count_diff_2nd'] = df_steps['count_diff'].diff()
                
                # Gradient...
                df_steps['gradient'] = np.gradient(df_steps['count'])
                
                # Max slp_2nd is bank, min is FP...
                idx_bank = df_steps['count_diff_2nd'].idxmax() 
                idx_fp = df_steps['count_diff_2nd'].idxmin() 
                
#                # Rolling mean of slope...
#                df_steps['mean_slp'] = df_steps['count_slp'].rolling(window=3, center=True).mean()                
#                
#                # Cumulative of count...
#                df_steps['count_cs'] = df_steps['count'].cumsum()
#                
#                # Slope of cumulative curve...
#                i_offset=1
#                df_steps['count_pct'] = df_steps['count'].pct_change(i_offset) 
#                
#                # Delta slope...
#                df_steps['pct_delta'] = df_steps['count_pct'].diff(i_offset)
#                
#                srs_peaks = find_peaks_from_series(df_steps['mean_slp'], 6, 1.35)                
#                srs_peaks.sort_values(inplace=True, ascending=False)
                
#                srs_peaks.index[0]
                
                # Index of max slope value...
#                test = df_steps['count_slp'].idxmax()
                
#                # Enforce positive pct_change?
#                if (df_steps['pct_delta'].iloc[df_steps['pct_delta'].idxmax()-i_offset] > 0):
#                    # Height value where slope is greatest...
#                val_height = df_steps['height'].iloc[df_steps['count_pct'].idxmax()-i_offset]
#                val_count_cs = df_steps['count_cs'].iloc[df_steps['count_pct'].idxmax()-i_offset]
#    #                val_count = df_steps['count'].iloc[df_steps['count_delta'].idxmax()-i_offset]
#                else:
#                    val_height = np.nan
#                    val_count_cs = np.nan


#                print('linkno: {} | xn: {} | val_height: {} | val_cs: {}\t| w_max: {:3.1f}'.format(line['properties']['linkno'], cntr, val_height, val_count_cs, w_max))    
                
#                if shed['properties'][str_reachid] == 6 and cntr == 0:
#                    print('pause')                
                try:
                    # ===================================================================== 
                    ax.plot(df_steps['height'].iloc[idx_bank], df_steps['count'].iloc[idx_bank], marker='o', linestyle='')
                    ax.plot(df_steps['height'].iloc[idx_fp], df_steps['count'].iloc[idx_fp], marker='*', linestyle='')
    #                ax.plot(df_steps['height'], df_steps['mean_slp'], marker='.')
                    ax.plot(df_steps['height'], df_steps['count'])   
    #                ax.plot(val_count_cs, val_height, marker='*', linestyle='') 
                    ax.set_title('ARCID: {}'.format(shed['properties'][str_reachid]), fontsize=11) 
                    # =====================================================================
                except:
                    pass
                
                cntr += 1
                
                prev_linkno = shed['properties'][str_reachid]
                
#                if line['properties'][str_reachid] == test_id and cntr > 50:
#                    break
                
            ax.set_ylabel('count')
            ax.set_xlabel('height (m)') 
                             
#                if line['properties']['linkno'] > 3:
#                    sys.exit()
                    
    return
    
# ===============================================================================
#  Analyze DEM in vertical slices using 2D Xn's (buffered Xn lines)
# =============================================================================== 
def analyze_hand_2D_xns(str_xns_path, str_grid_path, i_interval):
    
    buff_dist=10
    
           
    # Open the dem...
    with rasterio.open(str(str_grid_path)) as ds_grid:
        
        out_meta = ds_grid.meta.copy()
        
        # Open the xn layer...
        with fiona.open(np.str(str_xns_path), 'r') as xns: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
            # Get the crs...
    #        streamlines_crs = streamlines.crs          

#            fig = plt.figure()            
#            ax = fig.add_subplot(111)  
                
            cntr = 0
            prev_linkno = 0
            for line in xns:    

                test_id=777
                if line['properties']['linkno'] <> test_id:
                    continue
                
                if line['properties']['linkno'] <> prev_linkno:
                    fig = plt.figure()            
                    ax = fig.add_subplot(111) 
                    cntr = 0
                
                # Buffer each feature...
                geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                
#                print('linkno: {}'.format(line['properties'][str_reachid]))
                
                # Successive buffer-mask operations to count bank pixels at certain intervals
#                for buff_dist in range(cell_size,30,cell_size):
                    
                geom_buff = geom.buffer(buff_dist,cap_style=2)
                buff = mapping(geom_buff)           
                            
                # Mask the bankpts file for each feature...
                w, out_transform = rasterio.tools.mask.mask(ds_grid, [buff], crop=True)
                
#                w[w<0]=9999. # handle no data vals?
                
                # Normalize to zero...
#                min_elev = w[w>0].min()             
#                w = w - min_elev
#                
#                # Find max value in window...
#                w_max = w.max()
                
#                # ===================================================================
#                # << 3D plot >> 
#                fig = plt.figure()            
#                ax = fig.add_subplot(111, projection='3d')
#                Y=np.arange(0,np.shape(arr_norm)[1],1)
#                X=np.arange(0,np.shape(arr_norm)[2],1)
#                X, Y = np.meshgrid(X, Y)
#                
#                arr_norm[arr_norm<0] = -0.5
#                out_image[out_image<0] = -0.5
#               
#                surf = ax.plot_surface(X, Y, out_image[0], cmap=cm.jet, linewidth=0, antialiased=True, rstride=1, cstride=1)
#                
#                # Add a color bar which maps values to colors.
#                fig.colorbar(surf, shrink=0.5, aspect=5)            
#                # ===================================================================
                                
                # Now perform vertical slices...
#                i_interval = 0.1 # m
#                i_rng = w_max # m
                i_rng = 5 # m
                arr_slices = np.arange(i_interval, i_rng, i_interval)    
                
                lst_count=[]
                # List comprehension here instead??
                for i_step in arr_slices:                    
                    lst_count.append(w[(w<i_step) & (w>=0)].size) # number of pixels greater than or equal to zero and less than the height interval

                df_steps = pd.DataFrame({'count':lst_count, 'height':arr_slices})
                
                # Slope of count...
                df_steps['count_diff'] = df_steps['count'].diff()
                
                # Slope of slope of count...
                df_steps['count_diff_2nd'] = df_steps['count_diff'].diff()
                
                # Max slp_2nd is bank, min is FP...
                idx_max = df_steps['count_diff_2nd'].idxmax() 
                idx_min = df_steps['count_diff_2nd'].idxmin() 
                
       
                try:

                    bank_idx = df_steps['height'].iloc[[idx_max, idx_min]].idxmin()
                    fp_idx = df_steps['height'].iloc[[idx_max, idx_min]].idxmax()
                    
                    # ===================================================================== 
                    ax.plot(df_steps['height'].iloc[bank_idx-1], df_steps['count'].iloc[bank_idx-1], marker='o', linestyle='')
                    ax.plot(df_steps['height'].iloc[fp_idx-1], df_steps['count'].iloc[fp_idx-1], marker='*', linestyle='')
    #                ax.plot(df_steps['height'], df_steps['mean_slp'], marker='.')
                    ax.plot(df_steps['height'], df_steps['count'])   
    #                ax.plot(val_count_cs, val_height, marker='*', linestyle='') 
                    ax.set_title('ARCID: {}'.format(line['properties']['linkno']), fontsize=11) 
                    # =====================================================================
                except:
                    pass
                
                cntr += 1
                
                prev_linkno = line['properties']['linkno']
                
                # Write back to a tif...
                out_meta.update({"driver": "GTiff",
                                 "height": w.shape[1],
                                 "width": w.shape[2],
                                 "transform": out_transform})                                     
                                 
                with rasterio.open("/home/sam.lamont/USGSChannelFPMetrics/testing/test.masked.{}.tif".format(cntr), "w", **out_meta) as dest:                                 
                    dest.write(w)
                
                if cntr>1:
                    sys.exit()
                
            ax.set_ylabel('count')
            ax.set_xlabel('height (m)') 
                    
    return

# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def analyze_hand_from_streamline_window(df_coords, str_hand_path):    

    w_height=100 # number of rows
    w_width=100  # number of columns
           
    with rasterio.open(str_hand_path) as ds_hand:
                
        out_meta = ds_hand.meta.copy()      
        
        arr_fp=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])  
   
        # Transform to pixel space
        df_coords['col'], df_coords['row'] = ~ds_hand.affine * (df_coords['x'], df_coords['y'])
        df_coords['col'], df_coords['row'] = ~ds_hand.affine * (df_coords['x'], df_coords['y'])   
        
        df_coords[['row','col']] = df_coords[['row','col']].astype(np.int32)  
        df_coords.drop_duplicates(['col','row'], inplace=True) # rounding to integer
              
        # Set slice inteval and range...
        i_interval = 0.1 # m
        i_rng = 5 # m
        
        arr_slices = np.arange(i_interval, i_rng+i_interval, i_interval)
        
        gp_coords = df_coords.groupby('linkno')
        
#        cntr=0
        for i_linkno, df in gp_coords:
            
#            if i_linkno <> 50:
#                continue            

#            fig = plt.figure()            
#            ax = fig.add_subplot(111)

#            lst_fp=[]
            lst_height=[]
            for tpl_row in df.itertuples():
                
#                if i_linkno>10:                 
#                    break
                
#                print('LINKNO: {}, w: {}, row: {}, col: {}'.format(tpl_row.linkno, tpl_row.Index, tpl_row.row, tpl_row.col))
                
                row_min = np.int(tpl_row.row - np.int(w_height/2))
                row_max = np.int(tpl_row.row + np.int(w_height/2))
                col_min = np.int(tpl_row.col - np.int(w_width/2))
                col_max = np.int(tpl_row.col + np.int(w_width/2))            
                
                # Now get the DEM specified by this window as a numpy array...
                w = ds_hand.read(1, window=((row_min, row_max),(col_min, col_max))) 
                
#                i_rng = w_max # m
                i_rng = 5 # m
                arr_slices = np.arange(i_interval, i_rng, i_interval)    
                
                # Get the number of pixels greater than or equal to zero and less than the height interval...
#                lst_count=[]
#                for i_step in arr_slices:
#                    lst_count.append(w[(w<i_step) & (w>=0)].size)
                
                lst_count = [w[(w<i_step) & (w>=0)].size for i_step in arr_slices]
                
                df_steps = pd.DataFrame({'count':lst_count, 'height':arr_slices})
                
                # Slope of count...
                df_steps['count_diff'] = df_steps['count'].diff()
                
                # Slope of slope of count...
                df_steps['count_diff_2nd'] = df_steps['count_diff'].diff()
                
                # Max slp_2nd is bank, min is FP...
                idx_max = df_steps['count_diff_2nd'].idxmax() 
                idx_min = df_steps['count_diff_2nd'].idxmin() 
                fp_idx = df_steps['height'].iloc[[idx_max, idx_min]].idxmax()

                tpl_lst=(tpl_row.Index, df_steps['height'].iloc[fp_idx-1])                                
                lst_height.append(tpl_lst)
                
                print('I. LINKNO: {}, w: {}'.format(tpl_row.linkno, tpl_row.Index))


            # Merge heights with df...
#            i_len=len(lst_height)
            df_heights=pd.DataFrame(lst_height, columns=['Index2','height'])     
            df_heights.set_index('Index2', inplace=True)
            
            df_heights['height_smoothed'] = df_heights['height'].rolling(window=2).mean()   

            # Fill in nan vals with originals...                                   
            df_heights['height_smoothed'].iloc[0]  =  df_heights['height'].iloc[0] 
                                                                
            df_h = df.join(df_heights)
                     
            # Loop over again...
            for tpl_row in df_h.itertuples():
                
    #                if i_linkno>10:                 
    #                    break
                
                print('II. LINKNO: {}, w: {}'.format(tpl_row.linkno, tpl_row.Index))
                
                if tpl_row.linkno == 102 and tpl_row.Index == 8772:
                    print('pause')
                
                row_min = np.int(tpl_row.row - np.int(w_height/2))
                row_max = np.int(tpl_row.row + np.int(w_height/2))
                col_min = np.int(tpl_row.col - np.int(w_width/2))
                col_max = np.int(tpl_row.col + np.int(w_width/2))            
                
                # Now get the DEM specified by this window as a numpy array...
                w = ds_hand.read(1, window=((row_min, row_max),(col_min, col_max))) 
                
                w[w<=tpl_row.height_smoothed] = 1
                w[w>tpl_row.height_smoothed] = 0
    
                arr_fp[row_min:row_max, col_min:col_max] = w + arr_fp[row_min:row_max, col_min:col_max] #bool_fp.astype(out_meta['dtype'])

#            if i_linkno > 100:
#                break
        
        arr_fp[arr_fp>0] = 1
        arr_fp[arr_fp<0] = 0
                                 
        print('Writing .tif file...')   
        str_fp_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMhand_TEST.tif'
        with rasterio.open(str_fp_path, "w", **out_meta) as dest:
            dest.write(arr_fp, indexes=1)
  
    return

# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def fp_from_streamline_window(df_coords, str_dem_path, str_fp_path):
    
    # Convert df_coords x-y to row-col via DEM affine
    # Loop over center row-col pairs accessing the window
    # Now loop over the linknos to get access grid by window...
    w_height=20 # number of rows
    w_width=20  # number of columns
#    j=0
    
#    lst_linknos=[]
#    lst_x1=[]
#    lst_y1=[]
#    lst_x2=[]
#    lst_y2=[]    
    
    with rasterio.open(str_dem_path) as ds_dem:
        
#        nodata_val = ds_dem.nodata
    
        # Transform to pixel space
        df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y'])
        df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y'])   
        
        df_coords[['row','col']] = df_coords[['row','col']].astype(np.int32)  
        df_coords.drop_duplicates(['col','row'], inplace=True) # rounding to integer
        
        out_meta = ds_dem.meta.copy()  
        buff=3 # cell size?
        cell_size=3
        
        arr_fp=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])         
#        arr_bankpts=out_meta['nodata']
        
#        lst_dist=[]
        
        for tpl_row in df_coords.itertuples():
            
            print(tpl_row.linkno)
            
            row_min = np.int(tpl_row.row - np.int(w_height/2))
            row_max = np.int(tpl_row.row + np.int(w_height/2))
            col_min = np.int(tpl_row.col - np.int(w_width/2))
            col_max = np.int(tpl_row.col + np.int(w_width/2))            
            
            # Now get the DEM specified by this window as a numpy array...
            w = ds_dem.read(1, window=((row_min, row_max),(col_min, col_max))) 
                            
            # Then extract the internal part of the window that contains the rotated window??
            
#            w_laplace = ndimage.filters.laplace(w)
#            w_sobel_x = ndimage.sobel(w, axis=0, mode='constant')
#            w_sobel_y = ndimage.sobel(w, axis=1, mode='constant')
#            w_sobel = np.hypot(w_sobel_x, w_sobel_y)
            
#            w[w==out_meta['nodata']] = 0
            
            if np.size(w) > 9: # make sure a window of appropriate size was returned from the DEM
                            
                # Calculate curvature...
#                w_dw_x, w_dw_y = np.gradient(w, cell_size)
#                w_ddw_x, temp1 = np.gradient(w_dw_x, cell_size)
#                temp2, w_ddw_y = np.gradient(w_dw_y, cell_size)            
#                w_curve = w_ddw_x + w_ddw_y
#                del temp1, temp2, w_dw_x, w_dw_y, w_ddw_x, w_ddw_y
                            
                Zy, Zx = np.gradient(w, cell_size)
                Zxy, Zxx = np.gradient(Zx, cell_size)
                Zyy, _ = np.gradient(Zy, cell_size)                            
                
#                # OR, gaussian curvature? (http://stackoverflow.com/questions/11317579/surface-curvature-matlab-equivalent-in-python)
#                w_curve = (Zxx * Zyy - (Zxy ** 2)) /  (1 + (Zx ** 2) + (Zy **2)) ** 2
                
                # OR, mean curvature...
                w_curve = (Zx**2 + 1)*Zyy - 2*Zx*Zy*Zxy + (Zy**2 + 1)*Zxx
                w_curve = -w_curve/(2*(Zx**2 + Zy**2 + 1)**(1.5))
                
                w_curve[w_curve<np.max(w_curve)*0.30] = out_meta['nodata']
#                w_curve = np.ma.masked_less(w_curve, np.max(w_curve)*0.30) # for plotting?
                
                arr_fp[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w_curve[buff:w_height-buff, buff:w_width-buff]
                
                if tpl_row.linkno == 13:
                    fig, ax = plt.subplots()
                    im = ax.imshow(w_curve, cmap='gist_rainbow')
                    fig.colorbar(im, orientation='vertical')
                    plt.show()     
                    break
                    
#                w_laplace = ndimage.filters.laplace(w_curve)    
#    #            w_laplace = np.ma.masked_less(w_laplace, np.max(w_laplace)*0.30)  # 0.30 for 3 m seems to work well  
#                w_laplace[w_laplace<np.max(w_laplace)*0.30] = out_meta['nodata']               
#                arr_bankpts[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w_laplace[buff:w_height-buff, buff:w_width-buff]


        arr_fp[arr_fp<=0.] = out_meta['nodata']            
        
        print('Writing bankpts .tif...')
        with rasterio.open(str_fp_path, "w", **out_meta) as dest:
            dest.write(arr_fp, indexes=1)
            
#        df_dist = pd.DataFrame(lst_dist)
    
    return
    
# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_streamline_window(df_coords, str_dem_path, str_bankpixels_path, self):
    
    print('Bank pixels from streamline windows...')
    
    # Convert df_coords x-y to row-col via DEM affine
    # Loop over center row-col pairs accessing the window
    # Now loop over the linknos to get access grid by window...
    w_height=20 # number of rows
    w_width=20  # number of columns
#    j=0
    
#    lst_linknos=[]
#    lst_x1=[]
#    lst_y1=[]
#    lst_x2=[]
#    lst_y2=[]    
    
    j=0   
    progBar = self.progressBar
    progBar.setVisible(True)   
    
    with rasterio.open(str_dem_path) as ds_dem:
        
#        nodata_val = ds_dem.nodata
    
        # Transform to pixel space
#        df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y'])
        df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y'])   
        
        df_coords[['row','col']] = df_coords[['row','col']].astype(np.int32)  
        df_coords.drop_duplicates(['col','row'], inplace=True) # rounding to integer
        
        out_meta = ds_dem.meta.copy()  
        buff=3 # cell size?
        cell_size=3
        
        arr_bankpts=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])         
#        arr_bankpts=out_meta['nodata']
        
        progBar.setRange(0, len(df_coords.index)) 
        
        for tpl_row in df_coords.itertuples():
            
            progBar.setValue(j)
            j+=1
            
#            print(tpl_row.linkno)
            
            row_min = np.int(tpl_row.row - np.int(w_height/2))
            row_max = np.int(tpl_row.row + np.int(w_height/2))
            col_min = np.int(tpl_row.col - np.int(w_width/2))
            col_max = np.int(tpl_row.col + np.int(w_width/2))            
            
            # Now get the DEM specified by this window as a numpy array...
            w = ds_dem.read(1, window=((row_min, row_max),(col_min, col_max))) 
                            
            # Then extract the internal part of the window that contains the rotated window??
            
#            w_laplace = ndimage.filters.laplace(w)
#            w_sobel_x = ndimage.sobel(w, axis=0, mode='constant')
#            w_sobel_y = ndimage.sobel(w, axis=1, mode='constant')
#            w_sobel = np.hypot(w_sobel_x, w_sobel_y)
            
#            w[w==out_meta['nodata']] = 0
            
            if np.size(w) > 9: # make sure a window of appropriate size was returned from the DEM
                            
                # Calculate curvature...
#                w_dw_x, w_dw_y = np.gradient(w, cell_size)
#                w_ddw_x, temp1 = np.gradient(w_dw_x, cell_size)
#                temp2, w_ddw_y = np.gradient(w_dw_y, cell_size)            
#                w_curve = w_ddw_x + w_ddw_y
#                del temp1, temp2, w_dw_x, w_dw_y, w_ddw_x, w_ddw_y
                            
                Zy, Zx = np.gradient(w, cell_size)
                Zxy, Zxx = np.gradient(Zx, cell_size)
                Zyy, _ = np.gradient(Zy, cell_size)                            
                
#                # OR, gaussian curvature? (http://stackoverflow.com/questions/11317579/surface-curvature-matlab-equivalent-in-python)
#                w_curve = (Zxx * Zyy - (Zxy ** 2)) /  (1 + (Zx ** 2) + (Zy **2)) ** 2
                
                # OR, mean curvature...
                w_curve = (Zx**2 + 1)*Zyy - 2*Zx*Zy*Zxy + (Zy**2 + 1)*Zxx
                w_curve = -w_curve/(2*(Zx**2 + Zy**2 + 1)**(1.5))
                
                w_curve[w_curve<np.max(w_curve)*0.30] = out_meta['nodata']
#                w_curve = np.ma.masked_less(w_curve, np.max(w_curve)*0.30) # for plotting?
                
                arr_bankpts[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w_curve[buff:w_height-buff, buff:w_width-buff]
                
#                if tpl_row.linkno == 13:
#                    fig, ax = plt.subplots()
#                    im = ax.imshow(w_curve, cmap='gist_rainbow')
#                    fig.colorbar(im, orientation='vertical')
#                    plt.show()     
#                    break
                    
#                w_laplace = ndimage.filters.laplace(w_curve)    
#    #            w_laplace = np.ma.masked_less(w_laplace, np.max(w_laplace)*0.30)  # 0.30 for 3 m seems to work well  
#                w_laplace[w_laplace<np.max(w_laplace)*0.30] = out_meta['nodata']               
#                arr_bankpts[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w_laplace[buff:w_height-buff, buff:w_width-buff]


        arr_bankpts[arr_bankpts<=0.] = out_meta['nodata']            
        
        print('Writing bankpts .tif...')
        with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
            dest.write(arr_bankpts, indexes=1)
            
#        df_dist = pd.DataFrame(lst_dist)
    
    return
    
# ===============================================================================
#  Calculate angle from vertical of left and right banks
# ===============================================================================
def find_bank_angles(tpl_bfpts, lst_total_slices, xn_len, xn_elev_n, parm_ivert, xn_ptdistance):

    try:
        # How many total slices?...
        total_slices = len(lst_total_slices)

        # Use second to last slice for bottom estimate
        if total_slices > 2:

            # Interpolate to find left position along bank...
            lf_bottom_ind = lst_total_slices[1][0]

            # LEFT BANK:  Make sure we're within bounds here
            if lf_bottom_ind == 0 or lf_bottom_ind == xn_len:
                lf_angle = 0
            else:
                x1 = lf_bottom_ind-1
                x2 = lf_bottom_ind
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                lf_bottombank_ind = interp_bank(x1,x2,y1,y2,yuk)

                if abs(lf_bottombank_ind - tpl_bfpts[1]) > 0:
                    lf_angle = atan((abs(lf_bottombank_ind - tpl_bfpts[1]))/((tpl_bfpts[3]-parm_ivert)*xn_ptdistance))*57.29578  # convert radians to degrees
                else:
                    lf_angle = 0

            # RIGHT BANK: Interpolate to find left position along bank...
            rt_bottom_ind = lst_total_slices[1][-1]

            # Make sure we're within bounds here
            if rt_bottom_ind == 0 or rt_bottom_ind == xn_len:
                rt_angle = 0
            else:
                x1 = rt_bottom_ind
                x2 = rt_bottom_ind+1
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                rt_bottombank_ind = interp_bank(x1,x2,y1,y2,yuk)

                if abs(rt_bottombank_ind - tpl_bfpts[2]) > 0:
                    rt_angle = atan((abs(rt_bottombank_ind - tpl_bfpts[2]))/((tpl_bfpts[3]-parm_ivert)*xn_ptdistance))*57.29578  # convert radians to degrees
                else:
                    rt_angle = 0

        else: # if there's only one slice, just set it to 0? or -9999?
##                            lf_angle = -9999
##                            rt_angle = -9999
            # Use bottom slice for bank angle estimate...
            lf_bottom_ind = lst_total_slices[0][0]
            rt_bottom_ind = lst_total_slices[0][-1]

            if abs(lf_bottom_ind - tpl_bfpts[1]) > 0:
                lf_angle = atan((abs(lf_bottom_ind - tpl_bfpts[1])*xn_ptdistance)/tpl_bfpts[3])*57.29578  # convert radians to degrees
            else:
                lf_angle = 0

            if abs(rt_bottom_ind - tpl_bfpts[2]) > 0:
                rt_angle = atan((abs(rt_bottom_ind - tpl_bfpts[2])*xn_ptdistance)/tpl_bfpts[3])*57.29578 # convert radians to degrees
            else:
                rt_angle = 0

        # NOTE: For now, just set any resulting negative values to -9999, until we figure out what's going on (27Mar2015, SJL)
        if lf_angle < 0:
            lf_angle = -9999.0

        if rt_angle < 0:
            rt_angle = -9999.0

        tpl_angles = (lf_angle,rt_angle)

    except Exception as e:
        print('\r\nError in find_bank_angles. Exception: {} \n'.format(e))

    return tpl_angles   
# ===============================================================================
# Search Xn outward to the right to find the first point greater than the left bank
# ===============================================================================
def search_right_gt(xnelev,prev_ind,lf_bf):

    # Search outward from end of previous slice...
    for i in range(prev_ind+1,len(xnelev),1):

        if xnelev[i] > lf_bf:
            bank_ind = i
            break
        else: # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind
# ===============================================================================
# Search Xn outward to the left to find the first point greater than the right bank
# ===============================================================================
def search_left_gt(xnelev,prev_ind,rt_bf):

    # Search outward from end of previous slice...
    for i in range(prev_ind,0,-1):

        if xnelev[i] > rt_bf:
            bank_ind = i
            break
        else: # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind
# ===============================================================================
# Interpolate to find positions of right/left bank
# ===============================================================================
def interp_bank(x1,x2,y1,y2,y_uk):

    x_uk = (((x2 - x1)*(y_uk - y1))/(y2-y1)) + x1;

    return x_uk   
# ===============================================================================
# Search for banks via slope break and vertical slices
# ===============================================================================
def find_bank_ratio_method(lst_total, ratio_threshold, xnelev_zero, slp_thresh):
    """
        Compares the length of the last gtzero slice (num of indices) vs. the previous slice

        Inputs:     lst_total   - a list of 1D array slice index values
                    ratio_threshold
                    xnelev_zero - the Xn elevation profile normalized to zero

        Output:     tpl_bfpts   - a tuple of bankfull points (left_index, right_index, height)
    """
    tpl_bfpts = () # output tuple
    num_slices = len(lst_total) - 1  # total number of slices, each at a height of param_vertstep
    xn_len = len(xnelev_zero) - 1 # length of Xn

    try:
        if num_slices > 2 and len(lst_total[num_slices-1]) > 2:
            top_area = len(lst_total[num_slices]) # should really include point distance but this will cancel out?
            below_area = len(lst_total[num_slices-1])

            # Check the ratio...
            this_ratio = float(top_area)/float(below_area)

            if (this_ratio > ratio_threshold): # USE THIS TO DRIVE THE BANK BREAK DETERMINATION INSTEAD

                # Find end indices of this and of previous slice...
                prev_lf_ind = lst_total[num_slices-1][0]
                prev_rt_ind = lst_total[num_slices-1][-1]

                this_lf_ind = lst_total[num_slices][0]
                this_rt_ind = lst_total[num_slices][-1]

                # Bottom left and right for searching...
                bottom_lf = lst_total[0][0]
                bottom_rt = lst_total[0][-1]

                # First derivative is slope...
                lf_arr = np.array(xnelev_zero[this_lf_ind:prev_lf_ind+1])
                rt_arr = np.array(xnelev_zero[prev_rt_ind-1:this_rt_ind])
                firstdiff_left = np.diff(lf_arr)
                firstdiff_right = np.diff(rt_arr)

                # Set both indices to negative 1 initially...
                rt_bank_ind = -1
                lf_bank_ind = -1

                # Look for the first occurrence of a very small slope value in both directions...
                #slp_thresh = 0.03 # ? a parameter
                for r, this_rt in enumerate(firstdiff_right):
                    if this_rt < slp_thresh:
                        rt_bank_ind = r + prev_rt_ind - 1
                        break

                # Left...reverse it first?
                firstdiff_left_rev = firstdiff_left[::-1]

                for r, this_lf in enumerate(firstdiff_left_rev):
                    if this_lf > -slp_thresh:
                        lf_bank_ind = prev_lf_ind - r
                        break

                # Make sure rt_bank_ind is not greater than total xn length?
                if prev_lf_ind > 0 and prev_rt_ind < xn_len:

                    # Find the smallest height of the two...
                    if rt_bank_ind > 0 and lf_bank_ind < 0: # only the right index exists

                        # Interpolate to find left bankfull...
                        bf_height = xnelev_zero[rt_bank_ind]
                        lf_x2 = search_left_gt(xnelev_zero,bottom_rt,bf_height)

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if(lf_y1 == lf_y2):
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(lf_x1,lf_x2,lf_y1,lf_y2,bf_height)

                            tpl_bfpts = (lfbf_ind,rt_bank_ind,bf_height)

                    elif lf_bank_ind > 0 and rt_bank_ind < 0: # only the left index exists

                        # Interpolate to find right bank index...
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(xnelev_zero,bottom_lf,bf_height)

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if (rt_y1 == rt_y2):
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(rt_x1, rt_x2, rt_y1, rt_y2, bf_height)

                            tpl_bfpts = (lf_bank_ind,rtbf_ind,bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[rt_bank_ind] < xnelev_zero[lf_bank_ind]: # right is smaller than left

                        # Interpolate to find left bankfull...
                        bf_height = xnelev_zero[rt_bank_ind]
                        lf_x2 = search_left_gt(xnelev_zero,bottom_rt,bf_height) # search all the way across?

                       # lf_x2 = search_left_lt(xnelev_zero, lf_bank_ind, bf_height) # find the index that's just smaller than bank height on left side FASTER?

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if(lf_y1 == lf_y2):
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(lf_x1,lf_x2,lf_y1,lf_y2,bf_height)

                            tpl_bfpts = (lfbf_ind,rt_bank_ind,bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[lf_bank_ind] < xnelev_zero[rt_bank_ind]: # left is smaller than right

                        # Interpolate to find right bank index...
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(xnelev_zero,bottom_lf,bf_height) # Searches all the way across channel

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if (rt_y1 == rt_y2):
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(rt_x1, rt_x2, rt_y1, rt_y2, bf_height)

                            tpl_bfpts = (lf_bank_ind,rtbf_ind,bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[lf_bank_ind] == xnelev_zero[rt_bank_ind]: # they're exactly equal
                        #print 'they are the same!'
                        bf_height = xnelev_zero[lf_bank_ind]
                        tpl_bfpts = (lf_bank_ind,rt_bank_ind,bf_height)

    except Exception as e:
        print('\r\nError in find_bank_ratio_method. Exception: {} \n'.format(e))

    return tpl_bfpts
    
# ===============================================================================    
# Check for continuity in vertical cross section slices
# ===============================================================================
def is_contiguous(gtzero_inds):
    """
        Used by analyze_elev function
    """
    if (np.max(gtzero_inds) - np.min(gtzero_inds)) == np.count_nonzero(gtzero_inds) - 1:
        # Contiguous, continue
        bool_cont = True

    else:
        # Not contiguous, trim off the extras
        bool_cont = False

    return bool_cont
    
# ===============================================================================
#    Analyze the elevation profile of each Xn and determine metrics
# ===============================================================================    
def analyze_xnelev(df_xn_elev, param_ivert, xn_ptdist, param_ratiothreshold, param_slpthreshold, nodata_val):
    """
        Input:  List of elevation values at points along all Xn's.  Each input list entry
                is a tuple of all elevation values along that Xn.
                param_ratiothreshold = 1.5
                param_slpthreshold = 0.03

        Output: Metrics - Xn number, bank locations (indices), bank height, etc... channel width, for each Xn
                Tuple: (Xn num, lf_i, rt_i, bf_height)

        Procedure:
            1. Loop over Xn list
            2. Normalize Xn to zero
            3. Loop over vertical slices using a step set by param_ivert
            4. Find contiguous sets of indices for each slice
            5. Search for slope break
            6. If slope break exists, determine "bankfull" locations along Xn
    """
    ## USE NUMPY ARRAYS VS. LISTS??
    lst_bfmetrics = [] # list of tuples to contain output

    try:
        
#        # =========== Progress Dialog ================
#        progDialog = QtGui.QProgressDialog("", "Cancel", 0, len(df_xn_elev.index), self)
#        progDialog.setGeometry(200, 80, 300, 20)
#        progDialog.setWindowTitle('Calculating metrics from Xn\'s...')
#        progDialog.setWindowModality(QtCore.Qt.WindowModal)
#        progDialog.show()
#        # ============================================  
        
#        progDialog.setRange(0, len(df_xn_elev.index))
        
#        j=0
        for tpl_row in df_xn_elev.itertuples():
            
#            progDialog.value(j)
#            j+=1
            
            this_linkno = tpl_row.linkno
            
#            print('this_linkno: {} | index: {}'.format(this_linkno, tpl_row.Index))
            
#            if tpl_row.Index == 2254:
#                print('pause')

            # A list to store the total number of indices/blocks in a Xn...
            lst_total_cnt = []

            # Normalize elevation to zero...
#            arr_elev = tpl_row.elev
#            arr_elev = arr_elev[arr_elev != np.float32(nodata_val)]
#            thisxn_norm = arr_elev - np.min(arr_elev)
            
            thisxn_norm = tpl_row.elev - np.min(tpl_row.elev)


            # Loop from zero to max(this_xn_norm) using a pre-defined vertical step (0.2 m?)...
            for this_slice in np.arange(0, np.max(thisxn_norm), param_ivert):

                # The indices of positives...
                gtzero_indices = np.nonzero((this_slice - thisxn_norm) > 0)[0] # Zero index get the first element of the returned tuple
                
                # NOTE:  DO THIS IN TERMS OF MAP COORDINATES HERE??

                # Use funtion to check if contiguous...
                if np.size(gtzero_indices) == 0: # the first loop only

                    # get the index of the zero value...
                    lst_total_cnt.append(np.where(thisxn_norm == 0)[0])
                    prev_val = lst_total_cnt[0][0]

                elif is_contiguous(gtzero_indices):

                    # Yes, it is contiguous
                    # Save it to the total count...
                    lst_total_cnt.append(gtzero_indices)
                    prev_val = gtzero_indices[0] # just need any value from the contiguous array

                else:
                    # No, it's not contiguous
                    # Find the contiguous part of the slice...
                    tpl_parts = np.array_split(gtzero_indices,np.where(np.diff(gtzero_indices)!=1)[0]+1) # splits the contiguous elements into separate tuple elements

                    # Find the one that contains an element of the previous slice?...
    ##                if prev_val in [this_arr for this_arr in tpl_parts]: # use a list comprehension here?
                    for this_arr in tpl_parts:
                        if prev_val in this_arr[:]:
                            lst_total_cnt.append(this_arr)
                            prev_val = this_arr[0]
                            break
                
                # NOTE:  DO THIS IN TERMS OF MAP COORDINATES HERE??
                tpl_bankfullpts = find_bank_ratio_method(lst_total_cnt,param_ratiothreshold,thisxn_norm,param_slpthreshold) 

                if tpl_bankfullpts: # test to see if there's anything in there?? THIS SEEMS LIKE KIND OF A HACK?

                    if tpl_bankfullpts[0] - tpl_bankfullpts[1] == 0: # BF points are so close that rounding to int gives identical values
                        # print('pause')
                        break                                   
                
                    # Add Xn number to the output...
                    #xn_elev_norm = tpl_thisxn[2] - np.min(tpl_thisxn[2]) # normalized elevation profile
                    xn_length = len(thisxn_norm)

                    #tpl_bankfullpts = (xn_cnt,) + tpl_bankfullpts + (lst_total_cnt,) + (this_linkno,) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),)

                    # Bank points tuple...
#                    tpl_bankfullpts = (tpl_thisxn[3],) + tpl_bankfullpts # add local ID
                    tpl_bankfullpts = (tpl_row.Index,) + tpl_bankfullpts

                    # Find bank angles...
                    tpl_bankangles = find_bank_angles(tpl_bankfullpts, lst_total_cnt, xn_length, thisxn_norm, param_ivert, xn_ptdist)

                    # Estimate bankfull area...
                    # (Bank height - xn_elev_norm[i])*xn_ptdist
                    # Round up on left, round down on right
                    bf_area = 0
                    lst_bf_rng = range(int(ceil(tpl_bankfullpts[1])),int(tpl_bankfullpts[2])+1,1)
                    
                    for i in lst_bf_rng:
                        bf_area += (tpl_bankfullpts[3] - thisxn_norm[i])*xn_ptdist                        

                    # Channel width...
                    ch_width = (tpl_bankfullpts[2]-tpl_bankfullpts[1])*xn_ptdist

                    # Overbank ratio...
#                    if (tpl_bankfullpts[2]-tpl_bankfullpts[1]) <= 0:
#                        overbank_ratio = -9999.0
#                    else:
#                        overbank_ratio = len(lst_total_cnt[-1])/(tpl_bankfullpts[2]-tpl_bankfullpts[1])
                        
                    try:
                        overbank_ratio = len(lst_total_cnt[-1])/(tpl_bankfullpts[2]-tpl_bankfullpts[1])
                    except:
                        overbank_ratio = -9999.0
                        
                    if bf_area == 0:
#                        print('pause')
                        bf_area=-9999.0
                        total_arearatio = -9999.0
                    else:
                        # Also try area under entire Xn length relative to BF area...
                        total_xn_area = sum(thisxn_norm*xn_ptdist)
                        try:
                            total_arearatio = (total_xn_area - bf_area)/bf_area
                        except:
                            total_arearatio = -9999.0


                    tpl_metrics = tpl_bankfullpts + (lst_total_cnt,) + (this_linkno,) + (tpl_bankfullpts[3] + np.min(tpl_row.elev),) + tpl_bankangles + (bf_area,) + (ch_width,) + (overbank_ratio,) + (total_arearatio,)

##                    # Output metrics tuple...       local ID      bank angles        bank height                        bank elevation                    xn area
##                    tpl_metrics = (this_linkno,) +  (xn_cnt,) + tpl_bankangles + (tpl_bankfullpts[3],) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),) + (bf_area,)

##                    # Output bank points tuple...
##                    tpl_bankpts = tpl_bankfullpts + (this_linkno,) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),)

                    lst_bfmetrics.append(tpl_metrics)

                    break # no need to keep slicing here, unless we want to try for FP analysis
    ##            else:
    ##                print 'no bank!'

#        progDialog.close()
    
    except Exception as e:
        print('\r\nError in analyze_xn_elev. Exception: {} \n'.format(e))
        sys.exit()

    return lst_bfmetrics    
    
# ====================================================================================
#  Calculate channel metrics based on the bankpoint slope-threshold method at each Xn,
#   writing the bank points to a shapefile 
# ====================================================================================
def chanmetrics_bankpts(df_xn_elev, str_xnsPath, str_demPath, str_bankptsPath, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh): 
       
    print('Channel metrics from bank points...')    
    
    # << BEGIN LOOP >>
    # Do the rest by looping in strides, rather than all at once, to conserve memory...(possibly using multiprocessing)
    xn_count = get_feature_count(str_xnsPath)
    
    # Striding...
    arr_strides = np.linspace(0, xn_count, xn_count/100)
    arr_strides = np.delete(arr_strides,0)  

#    progBar = self.progressBar
#    progBar.setVisible(True)
#    progBar.setRange(0, arr_strides.size)
        
    # Now loop over the linknos to get access grid by window...
    with rasterio.open(str_demPath) as ds_dem:
        
        dem_crs = ds_dem.crs
        nodata_val = ds_dem.nodata
        
        # Define the schema for the output bank points shapefile...
        schema = {'geometry': 'Point', 'properties': {'xn_num': 'int', 'linkno':'int','bank_hght':'float','bank_elev':'float',
                    'bnk_ang_1':'float','bnk_ang_2':'float','bf_area':'float','chan_width':'float','obank_rat':'float',
                    'area_ratio':'float'}} 
        
        with fiona.open(str_bankptsPath, 'w', driver='ESRI Shapefile', crs=dem_crs, schema=schema) as bankpts:
            
#            k=0
            j=0
            for indx in arr_strides:
                
#                progBar.setValue(k)
#                k+=1
                
#                print('\tIndex {} - {}/{}'.format(j,int(indx),xn_count))        
                df_xn_elev_n = df_xn_elev.iloc[j:int(indx)]
                j = int(indx)+1        
              
#                print('\tIndex 143613 - 143712/{}'.format(xn_count))
#                df_xn_elev_n = df_xn_elev.iloc[143613:143712]
            
                # << INTERPOLATE XNs >>  
                df_bank_metrics = pd.DataFrame(analyze_xnelev(df_xn_elev_n, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh, nodata_val),
                                                columns=['xn_no','left_ind','right_ind','bank_height','slices','linkno','bank_elev','lf_bank_ang','rt_bank_ang','bankful_area','chan_width','overbank_ratio','area_ratio'])
                
                df_bank_metrics.set_index('xn_no', inplace=True)
                
                df_map = pd.merge(df_xn_elev, df_bank_metrics, left_index=True, right_index=True) 

                lst_lfbank_row=[]
                lst_lfbank_col=[]
                lst_rtbank_row=[]
                lst_rtbank_col=[]
                
                for tpl_row in df_map.itertuples():
                    
                    lst_lfbank_row.append(interpolate(tpl_row.xn_row, tpl_row.left_ind))
                    lst_lfbank_col.append(interpolate(tpl_row.xn_col, tpl_row.left_ind))
                    lst_rtbank_row.append(interpolate(tpl_row.xn_row, tpl_row.right_ind))
                    lst_rtbank_col.append(interpolate(tpl_row.xn_col, tpl_row.right_ind))
                    
                df_map['lfbank_row'] = pd.Series(lst_lfbank_row).values
                df_map['lfbank_col'] = pd.Series(lst_lfbank_col).values
                df_map['rtbank_row'] = pd.Series(lst_rtbank_row).values
                df_map['rtbank_col'] = pd.Series(lst_rtbank_col).values    
                
                # Transform to pixel space...
                df_map['lfbank_x'], df_map['lfbank_y'] = ds_dem.affine * (df_map['lfbank_col'], df_map['lfbank_row'])
                df_map['rtbank_x'], df_map['rtbank_y'] = ds_dem.affine * (df_map['rtbank_col'], df_map['rtbank_row'])
                        
                # Write bankpts shapefile...
#                print('Writing to bank points shapefile...')
                for tpl_row in df_map.itertuples():
                    
#                    progDialog.setValue(k)
#                    k+=1
                    
                    tpl_left = (tpl_row.lfbank_x, tpl_row.lfbank_y)
                    tpl_right = (tpl_row.rtbank_x, tpl_row.rtbank_y)            
                    
                    # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                    lf_pt = {'type': 'Point', 'coordinates':tpl_left}
                    rt_pt = {'type': 'Point', 'coordinates':tpl_right}
                   
                    prop_lf = {'xn_num':int(tpl_row.Index),'linkno':int(tpl_row.linkno_x),'bank_hght':tpl_row.bank_height,'bank_elev':tpl_row.bank_elev,
                            'bnk_ang_1':tpl_row.lf_bank_ang,'bf_area':tpl_row.bankful_area,'bnk_ang_2':-9999.,
                            'chan_width':tpl_row.chan_width,'obank_rat':tpl_row.overbank_ratio,'area_ratio':tpl_row.area_ratio}

                    prop_rt = {'xn_num':int(tpl_row.Index),'linkno':int(tpl_row.linkno_x),'bank_hght':tpl_row.bank_height,'bank_elev':tpl_row.bank_elev,
                            'bnk_ang_2':tpl_row.rt_bank_ang,'bf_area':tpl_row.bankful_area,'bnk_ang_1':-9999.,
                            'chan_width':tpl_row.chan_width,'obank_rat':tpl_row.overbank_ratio,'area_ratio':tpl_row.area_ratio}
                            
                    bankpts.write({'geometry': lf_pt, 'properties':prop_lf})  
                    bankpts.write({'geometry': rt_pt, 'properties':prop_rt})
                    
#                sys.exit() # for testing
   
# ===================================================================================
#  Read an existing Xn file, calculate xy bounds for each linkno and read the DEM
#  according to that window
# ===================================================================================
def get_stats(group):
    return {'min': group.min(), 'max': group.max()}
    
#def transform(row):
#    return row['a']
    
def read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path):    
    
    print('Reading and interpolating elevation along Xn\'s...')
    
    lst_linknos=[]
    lst_x1=[]
    lst_y1=[]
    lst_x2=[]
    lst_y2=[]

#    start_time = timeit.default_timer()
    # First get all linknos...
    with fiona.open(np.str(str_xns_path), 'r') as xn_shp: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)

#        j=0  
#        progBar = self.progressBar
#        progBar.setVisible(True)
#        progBar.setRange(0, len(xn_shp))

#            # ==================================
#            progdialog = QtGui.QProgressDialog("Extracting elevation along cross sections...", "Cancel", 0, len(xn_shp))
#            progdialog.setWindowTitle("FACET")
#            progdialog.setWindowModality(QtCore.Qt.WindowModal)
#            progdialog.resize(350, 110)
#            progdialog.show()       
#            # ==================================    
         
        # Read each feature line...
        for line in xn_shp:
            lst_linknos.append(line['properties']['linkno'])
            lst_x1.append(line['geometry']['coordinates'][0][0])
            lst_y1.append(line['geometry']['coordinates'][0][1])      
            lst_x2.append(line['geometry']['coordinates'][1][0])
            lst_y2.append(line['geometry']['coordinates'][1][1]) 
            
    df_coords = pd.DataFrame({'linkno':lst_linknos, 'x1':lst_x1, 'y1':lst_y1, 'x2':lst_x2, 'y2':lst_y2})
         
    # Now loop over the linknos to get access grid by window...
    with rasterio.open(str(str_dem_path)) as ds_dem:
        
        nodata_val = ds_dem.nodata
    
        # Transform to pixel space
        df_coords['col1'], df_coords['row1'] = ~ds_dem.affine * (df_coords['x1'], df_coords['y1'])
        df_coords['col2'], df_coords['row2'] = ~ds_dem.affine * (df_coords['x2'], df_coords['y2'])
        
        ### OR...
        gp_coords = df_coords.groupby('linkno')
        
        lst_all_zi=[]
        j=0
        
        for linkno, df_linkno in gp_coords:
            
#            print(linkno)
            
            row_min = df_linkno[['row1','row2']].min(axis=0).min()
            row_max = df_linkno[['row1','row2']].max(axis=0).max()
            col_min = df_linkno[['col1','col2']].min(axis=0).min()
            col_max = df_linkno[['col1','col2']].max(axis=0).max()            
            
            # Now get the DEM specified by this window as a numpy array...
            w = ds_dem.read(1, window=((row_min, row_max+1),(col_min, col_max+1))) 
                        
            # NOW loop over each Xn...                                   
            for tpl_xn in df_linkno.itertuples():
                
#                progBar.setValue(j)
                j+=1    
                
#                    # ====================================== 
#                    QtCore.QCoreApplication.processEvents()
#                    if progdialog.wasCanceled():
#                        break     
#                               
#                    progdialog.setValue(j)
#                    # ======================================                     
                
                xn_len = int(np.hypot(tpl_xn.col2-tpl_xn.col1, tpl_xn.row2-tpl_xn.row1))
                
                lst_xnrow = np.linspace(tpl_xn.row1-row_min, tpl_xn.row2-row_min, xn_len)
                lst_xncol = np.linspace(tpl_xn.col1-col_min, tpl_xn.col2-col_min, xn_len)
                
                #xnptdist = xn_len/len(lst_xnrow) #this is always 1 cell or equivalent to cell_size in meters/feet?
                            
                arr_zi = w[lst_xnrow.astype(np.int), lst_xncol.astype(np.int)]   # nearest-neighbor                
#                lst_zi = ndimage.map_coordinates(w, np.vstack((lst_xnrow, lst_xncol)), order=1, mode='nearest') # use this for additonal interpolation options (ie, cubic, bilinear, etc)
                
#                plt.plot(np.arange(len(zi)), zi)   
                
                # Remove possible no data values...
                arr_zi = arr_zi[arr_zi > np.float32(nodata_val)]
                                
                # Convert these from window row/col to raster row/col for bankpt use...
                for i, xnrow in enumerate(lst_xnrow):
                    lst_xnrow[i] = lst_xnrow[i] + row_min
                    lst_xncol[i] = lst_xncol[i] + col_min

                tpl_out = (linkno, arr_zi, lst_xnrow, lst_xncol)
                lst_all_zi.append(tpl_out)
#                i += 1
           
#    print('\tTotal Xn\'s:  {}'.format(i))    
#    print('\tTime interpolating elevation along Xn\'s:  ' + str(timeit.default_timer() - start_time))

    return pd.DataFrame(lst_all_zi, columns=['linkno','elev','xn_row','xn_col'])
    
# ===================================================================================
#  Build the Xns for all reaches and write to shapefile
# ===================================================================================
def write_xns_shp(df_coords, streamlines_crs, str_xns_path, bool_isvalley, p_xngap, p_fitlength, p_xnlength):
    """
        Builds Xns from x-y pairs representing shapely interpolations along a reach

        Input: a list of tuples (row, col, linkno) for a reach
        
        Output: list of tuples of lists describing the Xn's along a reach (row, col)
        
    """
    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True)
#    progBar.setRange(0, pd.unique(df_coords['linkno']).size) # need length on unique linkno values
    
#        # ==================================
#        progdialog = QtGui.QProgressDialog("Building and Writing Cross Section File...", "Cancel", 0, pd.unique(df_coords['linkno']).size)
#        progdialog.setWindowTitle("FACET")
#        progdialog.setWindowModality(QtCore.Qt.WindowModal)
#        progdialog.resize(350, 110)
#        progdialog.show()       
#        # ==================================        
    
    slopeCutoffVertical = 20 # just a threshold determining when to call a Xn vertical (if it's above this, make it vertical. Otherwise things get whacky?)

    lst_xnrowcols = [] # the final output, a list of tuples of XY coordinate pairs for all Xn's for this reach

    XnCntr = 0
    m_init = 0
    
    gp_coords = df_coords.groupby('linkno')  
           
    # Create the Xn shapefile for writing...
#    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int', 'endpt1_x':'float', 'endpt1_y':'float', 'endpt2_x':'float', 'endpt2_y':'float'}} 
    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int'}} 
    
    print('Building and Writing Cross Section File...')
    with fiona.open(str_xns_path, 'w', driver='ESRI Shapefile', crs=streamlines_crs, schema=test_schema) as chan_xns:
        
        for i_linkno, df_linkno in gp_coords:
            
            i_linkno = int(i_linkno)
            
#                if i_linkno==67:
#                    print('pause')

#            progBar.setValue(j)
            j+=1 

#                # ====================================== 
#                QtCore.QCoreApplication.processEvents()
#                if progdialog.wasCanceled():
#                    break     
#                           
#                progdialog.setValue(j)
#                # ======================================    
                
            # NOTE:  Define Xn length (p_xnlength) -- and other parameters? -- relative to stream order
            i_order = df_linkno.order.iloc[0]
            if i_order == 1:
                p_xnlength=20
            elif i_order == 2:
                p_xnlength=23
            elif i_order == 3:
                p_xnlength=40
            elif i_order == 4:
                p_xnlength=50 
            elif i_order == 5:
                p_xnlength=60                    
    
            reach_len = len(df_linkno['x'])
        
            if reach_len <= p_xngap:
#                print('Less than!')
                continue # skip it for now
        
            # Loop along the reach at the specified intervals...(Xn loop)
            for i in range( p_xngap, reach_len-p_xngap, p_xngap ):                
        
                lstThisSegmentRows = []
                lstThisSegmentCols = []
        
                if p_fitlength > i or i + p_fitlength >= reach_len: # if i + paramFitLength > reach_len
                    fitLength = p_xngap
                else:
                    fitLength = p_fitlength
        
                lstThisSegmentRows.append(df_linkno['y'].iloc[i+fitLength])
                lstThisSegmentRows.append(df_linkno['y'].iloc[i-fitLength])
                lstThisSegmentCols.append(df_linkno['x'].iloc[i+fitLength])
                lstThisSegmentCols.append(df_linkno['x'].iloc[i-fitLength])
        
                midPtRow = df_linkno['y'].iloc[i]
                midPtCol = df_linkno['x'].iloc[i]
        
        
                if (max(lstThisSegmentCols) - min(lstThisSegmentCols) < 3):
                    m_init = 9999.0
                elif (max(lstThisSegmentRows) - min(lstThisSegmentRows) < 3):
                    m_init = 0.0001
                else:
                    LinFit = np.polyfit(lstThisSegmentCols,lstThisSegmentRows,1) # NOTE: could just use basic math here instead?!
                    m_init = LinFit[0]
        
        
                # Check for zero or infinite slope...
                if m_init == 0:
                    m_init = 0.0001
                elif isinf(m_init):
                    m_init = 9999.0
        
                # Find the orthogonal slope...
                m_ortho = -1/m_init
                          
                xn_steps = [-float(p_xnlength),float(p_xnlength)] # just the end points
        
                lst_xy=[]
                for r in xn_steps:
        
                    # Make sure it's not too close to vertical...
                    # NOTE X-Y vs. Row-Col here...
                    if (abs(m_ortho) > slopeCutoffVertical):
#                        fit_col_ortho.append(midPtCol)
#                        fit_row_ortho.append(midPtRow+r)
                        tpl_xy = (midPtCol, midPtRow+r)
                        
                    else:
                        fit_col_ortho = (midPtCol + (float(r)/(sqrt(1 + m_ortho**2))))
#                        fit_row_ortho.append((m_ortho)*(fit_col_ortho[-1]-midPtCol) + midPtRow)
                        
                        tpl_xy = float(((midPtCol + (float(r)/(sqrt(1 + m_ortho**2)))))), float(((m_ortho)*(fit_col_ortho-midPtCol) + midPtRow))
#                        tpl_xy = ((midPtCol + (float(r)/(sqrt(1 + m_ortho**2)))), ((m_ortho)*(fit_col_ortho-midPtCol) + midPtRow))
                    
                    lst_xy.append(tpl_xy)
                         
                XnCntr = XnCntr + 1
                
                # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                line = {'type': 'LineString', 'coordinates':lst_xy}
#                prop = {'linkno': i_linkno, 'endpt1_x':lst_xy[0][0], 'endpt1_y':lst_xy[0][1], 'endpt2_x':lst_xy[1][0], 'endpt2_y':lst_xy[1][1]}
                prop = {'linkno': i_linkno}
                chan_xns.write({'geometry': line, 'properties':prop}) 
                
#                if XnCntr > 10:
#                    break
        
    return lst_xnrowcols
    
# ===================================================================================
#  Build Xn's based on vector features
# ===================================================================================
def get_stream_coords_from_features(str_streams_filepath, cell_size, str_reachid, str_orderid):
        
#        try:
    lst_x=[]
    lst_y=[]
    lst_linkno=[]
    lst_order=[]
    
    print('Getting stream coords from features...')
 
    p_interp_spacing = cell_size #3 # larger numbers would simulate a more smoothed reach | NOTE: Hardcode this = grid resolution?
    j=0 # prog bar
    
    # Open the streamlines shapefile...
    with fiona.open(str(str_streams_filepath), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
   
        # Get the crs...
        streamlines_crs = streamlines.crs  
#        str_proj4 = crs.to_string(streamlines.crs)         

#            progBar.setRange(0,len(streamlines))

        # ==================================
        progdialog = QtGui.QProgressDialog("Getting stream coords from features...", "Cancel", 0, len(streamlines))
        progdialog.setWindowTitle("FACET")
        progdialog.setWindowModality(QtCore.Qt.WindowModal)
        progdialog.resize(350, 110)
        progdialog.show()       
        # ==================================
        
        for line in streamlines:
            
           # ====================================== 
           QtCore.QCoreApplication.processEvents()
           if progdialog.wasCanceled():
               break     
                       
           progdialog.setValue(j)
           # ======================================
           
           j+=1
#               self.emit(QtCore.SIGNAL("update(int)"), int(100*len(streamlines)/j)) 
           
           i_linkno = line['properties'][str_reachid]           
           i_order = line['properties'][str_orderid]
           
           print('{} | {}'.format(i_linkno, j))
           
           line_shply = LineString(line['geometry']['coordinates'])
           
           length = line_shply.length # units depend on crs
           
           if length > 9:
               
#               int_pts = np.linspace(3, length, length/p_interp_spacing)
               
               int_pts = np.arange(0, length, p_interp_spacing) # p_interp_spacing in projection units?

               for i in int_pts: # lambda here instead?
                  
                   i_pt = np.array(line_shply.interpolate(i))
                   
                   lst_x.append(i_pt[0])
                   lst_y.append(i_pt[1])
                   lst_linkno.append(i_linkno)   
                   lst_order.append(i_order)
           
               df_coords = pd.DataFrame({'x':lst_x, 'y':lst_y, 'linkno':lst_linkno, 'order':lst_order})
#               df_coords.set_index('linkno', inplace=True)  
               
#           if j == 100:
#               print('pause')
               
        df_coords.drop_duplicates(['x','y'], inplace=True) # duplicates due to interpolation (?)
                       
#        except Exception as e:
##            print('Error!: {}'.format(e.strerror))
#            QtGui.QMessageBox.information(self, 'Error!', '{}'.format(e.strerror))
        
    return df_coords, streamlines_crs # A list of lists      


    
    
    