# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:51:45 2016

@author: sam.lamont
"""
#import pprint # neat!

import math
import subprocess
import timeit
import numpy as np
from numpy import asarray

#from scipy.stats import gaussian_kde # TEST
#from scipy.optimize import curve_fit # TEST
from scipy import signal

import os
#import ntpath
from math import atan, ceil
import sys
from math import isinf, sqrt
import rasterio
import rasterio.mask
from rasterio.warp import transform
from rasterio.features import shapes
import rasterio.features

#import matplotlib
#matplotlib.use('TkAgg')

#import matplotlib.pyplot as plt
#plt.style.use('ggplot')
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

#from affine import Affine
import pandas as pd
#from scipy import ndimage
#from shapely.geometry import Point, 
from shapely.geometry import shape, mapping, LineString, MultiLineString, Point, MultiPoint
from shapely.ops import split
#from jenks import jenks
#from PyQt4 import QtGui, QtCore

#import scipy.io as sio

import fiona
from fiona import collection

#import gospatial as gs

np.seterr(over='raise')

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

# lstThisSegmentRows, lstThisSegmentCols, midpt_x, midpt_y, p_fpxnlen

# ================================================================================
#   For wavelet curvature calculation (Chandana)
# ================================================================================
def gauss_kern(sigma):    
    """ Returns a normalized 2D gauss kernel array for convolutions """ 

    sigma=int(sigma)

    x, y =np.mgrid[-5*sigma:5*sigma, -5*sigma:5*sigma]

    g2x = (1-x**2/sigma**2)*np.exp(-(x**2+y**2)/2/sigma**2)*1/np.sqrt(2*np.pi*sigma**2)/ 4 / sigma * np.exp(float(0.5));
    g2y = (1-y**2/sigma**2)*np.exp(-(x**2+y**2)/2/sigma**2)*1/np.sqrt(2*np.pi*sigma**2)/ 4 / sigma * np.exp(float(0.5));
  
    return g2x, g2y
    
# ================================================================================
#   For 2D cross sectional measurement
# ================================================================================    
def build_xns(lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow, p_xnlength):
    
    slopeCutoffVertical = 20 # another check
        
    # Find initial slope...
    if (abs(lstThisSegmentCols[0] - lstThisSegmentCols[-1]) < 3):
        m_init = 9999.0
    elif (abs(lstThisSegmentRows[0] - lstThisSegmentRows[-1]) < 3):
        m_init = 0.0001
    else:
        m_init = (lstThisSegmentRows[0] - lstThisSegmentRows[-1])/(lstThisSegmentCols[0] - lstThisSegmentCols[-1]) 
        
        
#    if (max(lstThisSegmentCols) - min(lstThisSegmentCols) < 3):
#        m_init = 9999.0
#    elif (max(lstThisSegmentRows) - min(lstThisSegmentRows) < 3):
#        m_init = 0.0001
#    else:
#        LinFit = np.polyfit(lstThisSegmentCols,lstThisSegmentRows,1) # NOTE: could just use basic math here instead?!
#        m_init = LinFit[0]

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
            tpl_xy = (midPtCol, midPtRow+r)
            
        else:
            fit_col_ortho = (midPtCol + (float(r)/(sqrt(1 + m_ortho**2))))                       
            tpl_xy = float(((midPtCol + (float(r)/(sqrt(1 + m_ortho**2)))))), float(((m_ortho)*(fit_col_ortho-midPtCol) + midPtRow))
        
        lst_xy.append(tpl_xy)    # A list of two tuple endpts          
    
    return lst_xy 

def get_cell_size(str_grid_path):
        
    with rasterio.open(str(str_grid_path)) as ds_grid:
        cs_x, cs_y = ds_grid.res    
        
    return cs_x
    
# ==========================================================================
#   For dissolving line features    
# ==========================================================================    
def dissolve_line_features(str_lines_path, output_filename):
    print('Dissolve line features...')

    lst_all=[]
    
    with fiona.open(str_lines_path) as lines:
        
        crs = lines.crs
        
        schema = {'geometry':'MultiLineString', 'properties':{'linkno':'int:6'}}         
        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:           
            
            for line in lines:
                
                lst_all.append(line['geometry']['coordinates'])
                        
            output.write({'properties':{'linkno':9999}, 'geometry':{'type':'MultiLineString','coordinates':lst_all}})                            
        
    return

# ==========================================================================
#   For points along line feature at uniform distance    
# ==========================================================================    
def points_along_line_features(str_diss_lines_path, output_filename):
    print('Points along line features...')
    
    p_interp_spacing = 1000

    with fiona.open(str_diss_lines_path) as line:
        
        crs = line.crs
        line = line[0]        
                        
        schema = {'geometry':'Point', 'properties':{'id':'int:6'}}         
        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:
                    
            line_shply = MultiLineString(line['geometry']['coordinates'])
           
            length = line_shply.length # units depend on crs
               
            step_lens = np.arange(0, length, p_interp_spacing) # p_interp_spacing in projection units?
    
            for i, step in enumerate(step_lens): # lambda here instead?
              
                i_pt = np.array(line_shply.interpolate(step))
    
                output.write({'properties':{'id':i}, 'geometry':{'type':'Point','coordinates':i_pt}})                            
        
    return
    
    
# ==========================================================================
#   For points along line feature at uniform distance    
# ==========================================================================    
def taudem_gagewatershed(str_pts_path, str_d8fdr_path):
    print('Points along line features...')
    
    inputProc = str(4)
    
    str_output_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_test_sheds.tif'
    
    cmd = 'mpiexec -n ' + inputProc + ' GageWatershed -p ' + '"' + str_d8fdr_path + '"' + ' -o ' + '"' + str_pts_path + '"'  +  ' -gw ' + '"' + str_output_path + '"'
    
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
        
    return
    
# ==========================================================================
#   For clipping features  
# ==========================================================================         
def clip_features(str_lines_path, output_filename, str_dem_path):
    print('Clipping streamlines to site DEM...')
#    # Build the output file name...
#    path_to_dem, dem_filename = os.path.split(str_dem_path)
#    output_filename = path_to_dem + '\\' + dem_filename[:-4]+'_nhdhires.shp'
    
    # Polygonize the raster DEM with rasterio...    
    with rasterio.open(str(str_dem_path)) as ds_dem:
        arr_dem = ds_dem.read(1)
  
    arr_dem[arr_dem>0] = 100
    mask = arr_dem == 100
   
    results = (
        {'properties': {'test': v}, 'geometry': s}
        for i, (s, v) 
        in enumerate(
            shapes(arr_dem, mask=mask, transform=ds_dem.transform)))
    
    poly = next(results)   
    poly_shp = shape(poly['geometry'])    

    # Now clip/intersect the streamlines file with results via Shapely...
    with fiona.open(str_lines_path) as lines:
        
        # Get the subset that falls within bounds of polygon...        
        subset = lines.filter(bbox=shape(poly['geometry']).bounds)
        lines_schema = lines.schema
        lines_crs = lines.crs
        
        with fiona.open(output_filename, 'w', 'ESRI Shapefile', lines_schema, lines_crs) as dst:
        
            for line in subset:
    
                if shape(line['geometry']).within(poly_shp):
                    dst.write(line)

    return output_filename

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

    print('Creating weight grid from streamlines...')

    lst_coords=[]
    lst_pts=[]
    lst_x=[]
    lst_y=[]

    with fiona.open(str_streamlines_path) as lines:        
        
        streamlines_crs = lines.crs # to use in the output grid
        
        # Get separate lists of start and end points...
        for line in lines:
            if line['geometry']['type'] == 'LineString': # Make sure it's a LineString
                # Add endpts...
                lst_coords.append(line['geometry']['coordinates'][-1])
                # Add startpts...
                lst_pts.append(line['geometry']['coordinates'][0])
                
        # If a start point is not also in the endpt list, it's first order...        
        for pt in lst_pts:            
            if pt not in lst_coords:
                lst_x.append(pt[0])
                lst_y.append(pt[1])
                                   
    # Open DEM to copy metadata and write a Weight Grid (WG)...
    with rasterio.open(str_dem_path) as ds_dem:            
        out_meta = ds_dem.meta.copy() 
        
        # Construct the output array...
        arr_danglepts=np.zeros([out_meta['height'], out_meta['width']], dtype=out_meta['dtype']) 
        
        tpl_pts = transform(streamlines_crs, out_meta['crs'], lst_x, lst_y)
        lst_dangles = zip(tpl_pts[0], tpl_pts[1])

        for coords in lst_dangles:                        
            col, row = ~ds_dem.transform * (coords[0], coords[1]) # BUT you have to convert coordinates from hires to dem  
            try:
                arr_danglepts[int(row),int(col)] = 1.
            except:
                continue
    #            # Could save these points here for later use building the streamlines??
    #            tpl_pts=(row,col)
    #            lst_finalpts_rowcols.append(tpl_pts)
        
    # Now write the new grid using this metadata...    
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

## ===============================================================================
##  Get row/col coords of start pts from weight grid
## ===============================================================================
#def get_rowcol_from_wg(str_danglepts_path):
#    
#    print('Reading weight grid...')
#    with rasterio.open(str_danglepts_path) as ds_wg:
#        wg_rast = ds_wg.read(1) # numpy array
##        wg_crs = ds_wg.crs
##        wg_affine = ds_wg.affine
#        
#    print('Getting indices...')
#    row_cols = np.where(wg_rast>0)    
#    
#    return row_cols

## ===============================================================================
##  For calling GoSpatial funcs
## ===============================================================================
def default_callback(str):
    print(str)
    
def run_gospatial_whiteboxtool(tool_name, args, exe_path, exe_name, wd, callback = default_callback):
    try:
        os.chdir(exe_path)
        cmd = []
        cmd.append("." + os.path.sep + exe_name)
        if len(wd) > 0:
            cmd.append("-cwd=\"{}\"".format(wd))

        cmd.append("-run={}".format(tool_name))
        args_str = ""
        for s in args:
            args_str += s.replace("\"", "") + ";"
        args_str = args_str[:-1]
        cmd.append("-args=\"{}\"".format(args_str))

        ps = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True)

        while True:
            line = ps.stdout.readline()
            if line != '':
                callback(line.strip())
            else:
                break

        return 0
    except Exception as e:
        print(e)
        return 1
        
def run_rust_whiteboxtool(tool_name, args, exe_path, exe_name, wd, callback = default_callback):
    try:
        os.chdir(exe_path)
        args2 = []
        args2.append("." + os.path.sep + exe_name)
        args2.append("--run=\"{}\"".format(tool_name))

        if wd.strip() != "":
            args2.append("--wd=\"{}\"".format(wd))

        for arg in args:
            args2.append(arg)

        # args_str = args_str[:-1]
        # a.append("--args=\"{}\"".format(args_str))

#        if self.verbose:
#            args2.append("-v")

        proc = subprocess.Popen(args2, shell=False, stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True)

        while True:
            line = proc.stdout.readline()
            sys.stdout.flush()
            if line != '':
                callback(line.strip())
#                if not self.cancel_op:
#                    callback(line.strip())
#                else:
#                    self.cancel_op = False
#                    proc.terminate()
#                    return 2

            else:
                break

        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as err:
        callback(str(err))
        return 1       
        
# ===============================================================================
#  Mega-function for processing a raw DEM
#   1. Breaching and filling
#   2. TauDEM functions
# ===============================================================================        
def preprocess_dem(str_dem_path, str_streamlines_path, str_mpi_path, str_taudem_path, str_whitebox_path, run_whitebox, run_wg, run_taudem):
    try:
        
        # Split DEM path and filename...  # NOT OS INDEPENDENT??
        path_to_dem, dem_filename = os.path.split(str_dem_path)           

        inputProc = str(4) # number of cores to use for TauDEM processes
               
        # << Define all filenames here >>
        str_danglepts_path = path_to_dem + '\\' + dem_filename[:-4]+'_wg.tif'
        
        dem_filename_tif = dem_filename[:-4]+'.tif'
        breach_filename_dep = dem_filename[:-4]+'_breach.dep'
        breach_filename_tif = dem_filename[:-4]+'_breach.tif'
        
        str_dem_path_tif = path_to_dem + '\\' + dem_filename[:-4]+'.tif'
        
#        fel = path_to_dem + '\\' + dem_filename[:-4]+'_breach.tif'
        fel_breach = os.path.join(path_to_dem + '/' + breach_filename_tif)
        p = os.path.join(path_to_dem + '/' + breach_filename_tif[:-4]+'_p.tif')
        sd8 = os.path.join(path_to_dem + '/' + breach_filename_tif[:-4]+'_sd8.tif')
        
        ad8_wg = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_ad8_wg.tif')
        wtgr = os.path.join(str_danglepts_path)
        ad8_no_wg = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_ad8_no_wg.tif')
        ord_g = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_ord_g.tif')
        tree = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_tree')
        coord =os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_coord')
        net = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_net.shp')
        w = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_w.tif')
        slp = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_slp.tif')
        ang = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_ang.tif')
        dd = os.path.join(path_to_dem + '\\' + breach_filename_tif[:-4]+'_hand.tif')
        
#        # ==================== TauDEM Paths =========================
        # Hardcode paths from user input...
        mpipath = os.path.join(str_mpi_path) #r'"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"'
        d8flowdir = '"' + str_taudem_path + '\D8FlowDir.exe"' # r'"C:\Program Files\TauDEM\TauDEM5Exe\D8FlowDir.exe"'
        areaD8 = '"' + str_taudem_path + '\AreaD8.exe"'
        streamnet = '"' + str_taudem_path + '\StreamNet.exe"'
        dinfflowdir = '"' + str_taudem_path + '\DinfFlowDir.exe"'
        dinfdistdown = '"' + str_taudem_path + '\DinfDistDown.exe"'        
               
                
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
        # =================== << Whitebox Functions >> =====================
        if run_whitebox:

            print('Whitebox Path:  ' + 'r' + '"' + str_whitebox_path + '"')

            str_whitebox_dir, str_whitebox_exe = os.path.split(str_whitebox_path) 
        
            # << Run the BreachDepressions tool, specifying the arguments >>
            name = "BreachDepressions"        
            args = [dem_filename, breach_filename_dep, '-1', '-1', 'True', 'True'] # GoSpatial verion. NOTE:  Make these last four variables accessible to user?
#            args = ['--dem='+dem_filename, '-o='+breach_filename_dep] # Rust version
            
            ret = run_gospatial_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
#            ret = run_rust_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
            if ret != 0:
                print("ERROR: return value={}".format(ret)) 
                
#             << Convert .dep to .tif here? >>  NOTE:  Only for DRB hack when using .dep files
            name = "WhiteBox2GeoTiff"
            args = [breach_filename_dep, breach_filename_tif] 
            
            ret = run_gospatial_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
            if ret != 0:
                print("ERROR: return value={}".format(ret))   
                
#            # Rust version...
#            name = "ConvertRasterFormat"
#            args = ['--input='+breach_filename_dep, '-o='+breach_filename_tif]
#            ret = run_rust_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)                
                
#            # << Also convert original DEM .dep to a .tif >>
#            name = "WhiteBox2GeoTiff"
#            args = [dem_filename, dem_filename_tif] 
#            
#            ret = run_gospatial_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
#            if ret != 0:
#                print("ERROR: return value={}".format(ret))    
                
            # Rust version...
#            name = "ConvertRasterFormat"
#            args = ['--input='+dem_filename, '-o='+dem_filename_tif]
#            ret = run_rust_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)                
        
        # TEST.  NOTE: Need to update spatial reference on .tifs here!
#        with rasterio.open(path_to_dem + '\\' + breach_filename_tif, 'r') as dem_tif: # , crs='EPSG:26918'
#            print(dem_tif.crs)
#            pass
#            
        if run_wg:
            create_wg_from_streamlines(str_streamlines_path, str_dem_path_tif, str_danglepts_path)            
            
        if run_taudem:

            # Testing...
#            mpipath = 'r' +  '"' + mpipath + '"'
#            fel = 'r' +  '"' + fel + '"'
            
#            print('mpipath: ' + mpipath)            
#            print('fel: ' + fel_breach)
#            print(' ')
            
#            # ==============  << 1. Pit Filling with TauDEM >> ================ 
#            cmd = 'mpiexec' + ' -n ' + inputProc + ' PitRemove -z ' + '"' + str_dem_path_tif + '"' + ' -fel ' + '"' + fel_pitremove + '"'
#            
#            # Submit command to operating system
#            print('Running TauDEM PitRemove...')
#            os.system(cmd)
#            
#            # Capture the contents of shell command and print it to the arcgis dialog box
#            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#            
#            # Get some feedback from the process to print out...
#            message = "\n"
#            for line in process.stdout.readlines():
#                line = line.decode()
#                if isinstance(line, bytes):	   # true in Python 3
#                    line = line.decode()
#                message = message + line        
#            print(message)                
            
 
            # ==============  << 2. D8 FDR with TauDEM >> ================        
#            cmd = '"' + mpipath + '"' + ' -n ' + inputProc + ' ' + d8flowdir + ' -fel ' + '"' + str_dem_path + '"' + ' -p ' + '"' + p + '"' + \
#                  ' -sd8 ' + '"' + sd8 + '"'
                  
            cmd = 'mpiexec' + ' -n ' + inputProc + ' D8FlowDir -fel ' + '"' + dem_filename + '"' + ' -p ' + '"' + p + '"' + \
                  ' -sd8 ' + '"' + sd8 + '"'                  
                          
            # Submit command to operating system
            print('Running TauDEM D8 Flow Direction...')
            os.system(cmd)
#            
            # Capture the contents of shell command and print it to the arcgis dialog box
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            
            # Get some feedback from the process to print out...
            message = "\n"
            for line in process.stdout.readlines():
                line = line.decode()
                if isinstance(line, bytes):	   # true in Python 3
                    line = line.decode()
                message = message + line        
            print(message)    
                
#    #        # ============= << 3.a AD8 with weight grid >> ================        
#    #        cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_wg + '"'  + ' -wg ' + '"' + wtgr + '"'  + ' -nc '
#            cmd = 'mpiexec' + ' -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_wg + '"'  + ' -wg ' + '"' + wtgr + '"'  + ' -nc '
#            
#            # Submit command to operating system
#            print('Running TauDEM D8 FAC (with weight grid)...')
#            os.system(cmd)
#            # Capture the contents of shell command and print it to the arcgis dialog box
#            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#            
#            message = "\n"
#            for line in process.stdout.readlines():
#                if isinstance(line, bytes):	    # true in Python 3
#                    line = line.decode()
#                message = message + line
#            print(message)   
#             
#            # ============= << 3.b AD8 no weight grid >> ================
#    #        cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_no_wg + '"'  +  ' -nc '
#            cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_no_wg + '"'  +  ' -nc '
#            
#            # Submit command to operating system
#            print('Running TauDEM D8 FAC (no weights)...')
#            os.system(cmd)
#            # Capture the contents of shell command and print it to the arcgis dialog box
#            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#            
#            message = "\n"
#            for line in process.stdout.readlines():
#                if isinstance(line, bytes):	    # true in Python 3
#                    line = line.decode()
#                message = message + line
#            print(message)            
#            
#            # ============= << 4 StreamReachandWatershed with TauDEM >> ================      
#            cmd = 'mpiexec -n ' + inputProc + ' StreamNet -fel ' + '"' + fel_pitremove + '"' + ' -p ' + '"' + p + '"' + \
#                  ' -ad8 ' + '"' + ad8_no_wg + '"' + ' -src ' + '"' + ad8_wg + '"' + ' -ord ' + '"' + ord_g + '"' + ' -tree ' + \
#                  '"' + tree + '"' + ' -coord ' + '"' + coord + '"' + ' -net ' + '"' + net + '"' + ' -w ' + '"' + w + \
#                  '"'        
#            # Submit command to operating system
#            print('Running TauDEM Stream Reach and Watershed...')
#            os.system(cmd)
#            
#            # Capture the contents of shell command and print it to the arcgis dialog box
#            process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#            
#            message = "\n"
#            for line in process.stdout.readlines():
#                if isinstance(line, bytes):	    # true in Python 3
#                    line = line.decode()
#                message = message + line        
#            print(message) 
#            
#            # Let's get rid of some output that we are not currently using...
#            try:
#    #            os.remove(w)
#                os.remove(coord)
#                os.remove(tree)
#                os.remove(ord_g)
#            except:
#                print('Warning: Problem removing files!')
#                pass
#                
#            # ============= << 5. Dinf with TauDEM >> =============                
#            print('Running TauDEM Dinfinity...')        
#            cmd = 'mpiexec -n ' + inputProc + ' DinfFlowDir -fel ' + '"' + fel_pitremove + '"' + ' -ang ' + '"' + ang + '"' + \
#                  ' -slp ' + '"' + slp + '"'
#            
#            # Submit command to operating system
#            os.system(cmd)
#            
#            # Capture the contents of shell command and print it to the arcgis dialog box
#            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#            
#            # Get some feedback from the process to print out...
#            message = "\n"
#            for line in process.stdout.readlines():
#                line = line.decode()
#        #            if isinstance(line, bytes):	   # true in Python 3
#        #                line = line.decode()
#                message = message + line        
#            print(message)            
#            
#            # ============= << 6. DinfDistanceDown (HAND) with TauDEM >> =============
#            distmeth = 'v'
#            statmeth = 'ave'
#            
#            # Use original DEM here...
#            print('Running TauDEM Dinf Distance Down...') # Use Breached or Raw DEM here?? Currently using Raw
#            cmd = 'mpiexec -n ' + inputProc + ' DinfDistDown -fel ' + '"' + str_dem_path_tif + '"' + ' -ang ' + '"' + ang + '"' + \
#                  ' -src ' + '"' + ad8_wg + '"' + ' -dd ' + '"' + dd + '"' + ' -m ' + statmeth + ' ' + distmeth
#        
#            # Submit command to operating system
#            os.system(cmd)
#            
#            # Get some feedback from the process to print out...
#            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) 
#            
#            message = "\n"
#            for line in process.stdout.readlines():
#                line = line.decode()
#        #            if isinstance(line, bytes):	   # true in Python 3
#        #                line = line.decode()
#                message = message + line
#                
#            print(message)

    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise        
    
    return  
    
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

# =================================================================================
#  Calculate channel width based on bank pixels and stream line parallel offsets, 
#  potentially subdividing using Xn's
#  NOTE: Make this a generic metric calculator by segment?  (ie., sinuosity, etc)
# =================================================================================
def channel_and_fp_width_bankpixels_segments_po_2Dfpxns(df_coords, str_streamlines_path, str_bankpixels_path, str_reachid, cell_size, p_fpxnlen, str_hand_path, parm_ivert):
    
    print('Channel width from bank pixels -- segmented reaches...')

#    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 

    # For Testing...
    bank_height_hand = -9999.
    chan_width_hand = -9999.
    
    reach_buff_dist = 30 # Half the total buffer width for calculating HAND-based reach geometry
    fp_height = 2.0 # How to get this?  By reach?
#    p_fitlength=30  # For FP Xn's
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_ch_fp_width_po_hand.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')
    
    schema_buff = {'geometry': 'Polygon', 'properties': {'buff': 'str'}}
    schema_output = {'geometry': 'LineString', 'properties': {'linkno':'int','ch_wid_total':'float', 'ch_wid_1':'float', 'ch_wid_2':'float', 'dist_sl':'float', 'dist':'float', 'sinuosity':'float', 'fp_width':'float', 'bh_hand':'float','cw_hand':'float','ba_hand':'float'}}
    
    # Access the floodplain grid...
    with rasterio.open(str(str_hand_path)) as ds_hand:
        
        out_meta = ds_hand.meta.copy()
        arr_out_pixels = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])                                  
        
        # Access the bank pixel layer...
        with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:    
            
            # Access the streamlines layer...
            with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
           
    #            progBar.setRange(0, len(streamlines)) 
                
                # Get the crs...
                streamlines_crs = streamlines.crs                
                
#                # For writing out the buffers...
#                with fiona.open(r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\buff_test_ls.shp','w','ESRI Shapefile', schema_buff) as buff_out:                     
                    
                # Open another file to write the width...
                with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
                    
                    for i_linkno, df_linkno in gp_coords:
                        
    #                    progBar.setValue(j)
    #                    j+=1                    
                
                        i_linkno = int(i_linkno)
#                        max_indx = len(df_linkno.index) - 1
                        
#                            if i_linkno != 116:
#                                continue
                        
                        print('linkno:  {}'.format(i_linkno))
              
                        # << Analysis by reach segments >>
                        # Set up index array to split up df_linkno into segments (these dictate the reach segment length)...
                        # NOTE:  Reach might not be long enough to break up
                        i_step=30 # is this same as fit_length defined above??
                        arr_ind = np.arange(i_step, len(df_linkno.index)+1, i_step) # NOTE: Change the step for resolution?                        
                        lst_dfsegs = np.split(df_linkno, arr_ind)                        
                        
                        for i_seg, df_seg in enumerate(lst_dfsegs): # looping over each reach segment
                            
                            arr_x = df_seg.x.values
                            arr_y = df_seg.y.values

                            try:
                                # Calculate straight line distance...
                                dist_sl = np.sqrt((arr_x[0] - arr_x[-1])**2 + (arr_y[0] - arr_y[-1])**2)                     
                            except:
                                print('Error calculated straight line distance')
                                dist_sl = 0.
                            
                            try:
                                # Create a line segment from endpts in df_seg...
                                ls = LineString(zip(arr_x, arr_y))
                                
                                dist = ls.length                                    
                                sinuosity = dist/dist_sl # ratio of sinuous length to straight line length
                                
                                lst_tally=[]                            

                                # Successive buffer-mask operations to count bank pixels at certain intervals
                                lst_buff=range(cell_size,30,cell_size)

                                for buff_dist in lst_buff:                                                            
                                    
                                    try:
                                        # Watch out for potential geometry errors here...
                                        ls_offset_left = ls.parallel_offset(buff_dist, 'left')
                                        ls_offset_rt = ls.parallel_offset(buff_dist, 'right')   
                                    except:
                                        print('Error performing offset buffer')
                                    
                                    # Buffer errors can result from complicated line geometry... 
                                    try:
                                        out_left, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_left)], crop=True)  
                                    except:
                                        print('Left offset error')
                                        out_left=np.array([0])
                                        
                                    try:
                                        out_rt, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_rt)], crop=True)
                                    except:
                                        print('Right offset error')
                                        out_rt=np.array([0])
                                        
                                    num_pixels_left = len(out_left[out_left>0])
                                    num_pixels_rt = len(out_rt[out_rt>0])
                                    
                                    # You want the number of pixels gained by each interval...                    
                                    tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
                                    lst_tally.append(tpl_out)                    
                                    df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','interval_left','interval_rt'])

                            except:
                                print('buffer exception pause')
                                # There may be errors associated with the buffer geometry.  Just skip it if so?
                                continue
                       
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
                                                                                                                
                            except Exception as e:
                                weighted_avg_left=-9999.
                                weighted_avg_rt=-9999.
                                fp_width=-9999.
                                print('Error calculating weighted average of buffer distance.  Maybe no bank pixels found? Exception: {} \n'.format(e))
#                                continue                                                                                                                        
                            
#                            if i_seg == 33:
#                                print('pause')                                
                            
                            # << CHANNEL PROPERTIES VIA HAND >>                            
#                                midpt_indx = int(len(arr_x)/2)
#                                                         
#                                # midPt x and y...
#                                midpt_x = df_seg.x.iloc[midpt_indx]
#                                midpt_y = df_seg.y.iloc[midpt_indx]                              
#                                
#                                # Send it the endpts of what you to draw a perpendicular line to...
#                                lst_xy = build_xns(list(arr_y), list(arr_x), midpt_x, midpt_y, p_fpxnlen)
#                                
#                                try:            
#                                    # Turn the list of tuples into a linestring
#                                    fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
#                                except:
#                                    print('Error converting Xn endpts to LineString')
#                                    continue
                                                                
#                                # 2D cross section buffer...         
#                                buff_len=dist_sl/1.85 # have of the line segment straight line distance
#                                geom_fpls_buff = fp_ls.buffer(fp_buff_dist, cap_style=2)
#                                buff = mapping(geom_fpls_buff)
                            
                            # Standard buffer on line segment...
                            geom_ls_buff = ls.buffer(reach_buff_dist, cap_style=3)
                            buff = mapping(geom_ls_buff)
                            
#                            if i_seg==5:
#                                # Write out buffer polygon(s)...NOTE:  Could write out ind
#                                with fiona.open(r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\buff_test.shp','w','ESRI Shapefile', schema_buff) as buff_out:                     
#                            buff_out.write({'properties': {'buff': 'mmmm'}, 'geometry': buff})                       
#                                break                              
                            
                            # Mask the bankpts file for each feature...
                            out_image, out_transform = rasterio.mask.mask(ds_hand, [buff], crop=True)
                            
                            # << Related to mapping the floodplain based on HAND height >>
#                                # Count the number of pixels in the buffered Xn...
#                                num_pixels = out_image[(out_image<=fp_height)&(out_image>=0.)].size
#                                 
#                                # Calculate area of FP pixels...
#                                area_pixels = num_pixels*(ds_hand.res[0]**2) # get grid resolution               
#                                
#                                # Calculate width by stretching it along the length of the 2D Xn...
#                                fp_width = area_pixels/(buff_len*2)       
                            fp_width=0 # For testing purposes
#        
#                                # Subtract channel width from fp width...
#                                fp_width = fp_width - (weighted_avg_left+weighted_avg_rt)
#                                
#                                if fp_width<0.: fp_width = 0 # don't have negatives              

                            # << CALL HAND VERTICAL SLICE ANALYSIS HERE >>                                                       
                            bank_height_hand, chan_width_hand, bank_ang_hand = analyze_hand_poly(out_image, dist_sl, reach_buff_dist, ds_hand.res[0], parm_ivert)
                            
                            # Write back the in-channel pixels to a tif...                                                
                            out_image = out_image[0]                                                        
                            shp=np.shape(out_image)                                
                            bounds = rasterio.transform.array_bounds(shp[0],shp[1],out_transform) # window bounds in x-y space (west, south, east, north)                                
                            col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[3]) # upper left row and column of window?                                
                            
                            row_min = np.int(row_min)                            
                            col_min = np.int(col_min)
                            row_max = np.int(row_min + shp[0])
                            col_max = np.int(col_min + shp[1])
                          
                            out_image[out_image<0] = 9999. # no data values?
#                            out_image[out_image<=bank_height_hand] = 1.
                            out_image[out_image>bank_height_hand] = 0.
            
                            # Re-assing channel pixels...
#                            arr_out_pixels[row_min:row_max, col_min:col_max] = out_image + arr_out_pixels[row_min:row_max, col_min:col_max] #bool_fp.astype(out_meta['dtype'])
                            
                            # Or take the average between new pixels and previous pixels...
                            arr_out_pixels[row_min:row_max, col_min:col_max] = arr_out_pixels[row_min:row_max, col_min:col_max] + (out_image + arr_out_pixels[row_min:row_max, col_min:col_max])/2
                            
                            # Write to the output shapefile here...
                            output.write({'properties':{'linkno':i_linkno,'ch_wid_total': weighted_avg_left+weighted_avg_rt,'ch_wid_1': weighted_avg_left,'ch_wid_2': weighted_avg_rt, 'dist_sl':dist_sl,'dist':dist,'sinuosity': sinuosity,'fp_width':fp_width,'bh_hand':bank_height_hand,'cw_hand':chan_width_hand,'ba_hand':bank_ang_hand}, 'geometry':mapping(ls)})                            
                            
    print('Writing .tif file...')   
    str_fp_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_test_pixels.tif'
    with rasterio.open(str_fp_path, "w", **out_meta) as dest:
        dest.write(arr_out_pixels, indexes=1)                              

    return   

# ===============================================================================
#  Delineate a FIM from the HAND grid using depth at each polygon (eg, catchment)
#
#    1. Read in catchment polygons (assume these have attributes necessary for regression?)
#    2. Calculate a HAND height (h) for each polygon based on some attribute(s)
#    3. Delineate FIM per catchment
#
# =============================================================================== 
def fim_hand_poly(str_hand_path, str_sheds_path, str_reachid):
    
    # Open the HAND layer...
    with rasterio.open(str(str_hand_path)) as ds_hand:  
        
        out_meta = ds_hand.meta.copy()
        arr_fim = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])  
        arr_fim[:,:] = out_meta['nodata']
        
        lst_h=[]
        lst_linkno=[]
        lst_prov=[]
        
        # Open the catchment polygon layer...
        with fiona.open(np.str(str_sheds_path), 'r') as sheds:
        
            for shed in sheds:
                
                # Get the linkno...
                linkno = shed['properties']['gridcode']
                
                # Get the Province...
                prov = shed['properties']['PROVINCE']
                
                # Get the Drainage Area in km^2...
                da_km2 = shed['properties']['DSContArea']/1000000
                
                if (prov == 'COASTAL PLAIN' and da_km2 >= 33 and da_km2 <= 2792):
                    h = 1.65
                elif (prov == 'PIEDMONT' and da_km2 >= 20 and da_km2 <= 1604):
                    h = (np.log10(da_km2)*0.471 + 0.523)**2
                elif (prov == 'VALLEY AND RIDGE' and da_km2 >= 16 and da_km2 <= 1748):
                    h = (np.log10(da_km2)*0.471 + 0.375)**2
                elif (prov == 'APPALACHIAN PLATEAUS' and da_km2 >= 52 and da_km2 <= 285):
                    h = (np.log10(da_km2)*0.471 + 0.041)**2
                else:
                    lst_h.append(-9999)
                    lst_linkno.append(linkno)
                    lst_prov.append(prov)
                    continue # skip this catchment
                
                lst_h.append(h)    
                lst_linkno.append(linkno)
                lst_prov.append(prov)
                
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
    str_fim_path = str_hand_path[:-4]+'_fim.tif'
    with rasterio.open(str_fim_path, "w", **out_meta) as dest:
        dest.write(arr_fim, indexes=1) 
        
    # Write HAND heights to csv...
    str_csv_path = str_hand_path[:-4]+'_fim_h.csv'
    df_h = pd.DataFrame({str_reachid:lst_linkno, 'prov':lst_prov, 'h':lst_h})
    df_h.to_csv(str_csv_path)
        
    return 

# ===============================================================================
#  Reach characteristics from HAND
# =============================================================================== 
def reach_characteristics_hand(str_sheds_path, str_hand_path, str_slp_path):
    '''
    Calculate the reach geometry metrics necessary for synthetic
    rating curve derivation from HAND grids.
    
    Returns: TO DO
    '''    
    ## Define the vertical slice array:
    i_interval=0.2 # vertical step height NOTE:  Make sure you're considering units here
    i_rng=10 # maximum height NOTE:  What should this be??
    arr_slices = np.arange(i_interval, i_rng, i_interval)      
    
    # Open the HAND layer:
    with rasterio.open(str(str_hand_path)) as ds_hand:  
        
#        out_meta = ds_hand.meta.copy()
        res = ds_hand.res[0]
#        arr_fim = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])  
#        arr_fim[:,:] = out_meta['nodata']

        # Open the slope layer:
        with rasterio.open(str(str_slp_path)) as ds_slp:        
        
            ## Open the catchment polygon layer:
            with fiona.open(np.str(str_sheds_path), 'r') as sheds:
            
                ## Loop over catchments:
                for shed in sheds:
                    
                    ## Get the linkno:
                    linkno = shed['properties']['gridcode']    
        
                    ## Mask the HAND grid for each catchment polygon:
                    w_hand, out_transform = rasterio.mask.mask(ds_hand, [shed['geometry']], crop=True)   
                    
                    ## ALSO mask the slp grid:
                    w_slp, out_tranform = rasterio.mask.mask(ds_slp, [shed['geometry']], crop=True)  
                    
#                    lst_count=[]
                    lst_props=[]
                    
                    ## Get reach length here from the attribute table?
                    length = shed['properties']['gridcode'] ## JOIN ATTRIBUTES FROM NET??? 
                    
    #                tot_len = reach_buff_len + 2*reach_buff_width # for rounded cap style; else tot_len = reach_buff_len (square)
                    
                    # List comprehension here instead??
                    for i, i_step in enumerate(arr_slices): 
                        
                        i_hand_w = w_hand[(w_hand<=i_step) & (w_hand>=0.)]
                
                        i_num_pixels = i_hand_w.size
                        
                        ## Create a boolean array to use for selecting the slope pixels for this step:
                        i_hand_bool = (i_hand_w > 0)
                        
                        i_slp_w = w_slp[i_hand_bool]
                           
#                        lst_count.append(i_num_pixels) # number of pixels greater than or equal to zero and less than the height interval
                            
                        ## Surface area of inundated zone:
                        area_surf = i_num_pixels*(res**2)
                        
                        ## Channel bed area of inundated zone:
                        area_bed = area_surf*(1 + i_slp_w**2)**0.5
                        
                        ## Volume of this slice:
                        volume = area_surf*i_interval # cubic meters correct?
                        
                        ## Cross-sectional area:
                        area_xn = volume/length
                        
                        ## Wetted perimeter:
                        wetted_p = area_bed/length
                        
                        ## Top width:
                        top_width = area_surf/length
                        
                        ## Hydraulic radius:
                        hyd_radius = area_xn/wetted_p
                        
                        ## Add to list:
                        lst_props.append((linkno, area_surf, area_bed, volume, area_xn, wetted_p, top_width, hyd_radius))
                                
    ## Create final output dataframe:
    df_props = pd.DataFrame(lst_props, columns=['linkno', 'area_surf', 'area_bed', 'volume',
                                                'area_xn', 'wetted_p', 'top_width', 'hyd_radius'])
    ## Write it to csv:
    df_props.to_csv('df_props.csv')
    
    return
   
# ===============================================================================
#  Analyze DEM in vertical slices using an individual polygon
# =============================================================================== 
def analyze_hand_poly(w, reach_buff_len, reach_buff_width, res, i_interval):
                        
    i_rng=10
    arr_slices = np.arange(i_interval, i_rng, i_interval)    
    
    lst_count=[]
    lst_width=[]
    
    tot_len = reach_buff_len + 2*reach_buff_width # for rounded cap style; else tot_len = reach_buff_len (square)
    
    # List comprehension here instead??
    for i_step in arr_slices: 

        num_pixels = w[(w<=i_step) & (w>=0.)].size
           
        lst_count.append(num_pixels) # number of pixels greater than or equal to zero and less than the height interval
            
        # Calculate area of FP pixels...
        area_pixels = num_pixels*(res**2)
        
        # Calculate width by stretching it along the length of the 2D Xn...
        lst_width.append(area_pixels/(tot_len))                 

    df_steps = pd.DataFrame({'count':lst_count, 'height':arr_slices, 'width':lst_width})
    
#    df_steps.plot(x='width',y='height', marker='.')            
    
    # Slope of width...
    df_steps['width_diff'] = df_steps['width'].diff()
    
    # Slope of slope of count...
    df_steps['width_diff_2nd'] = df_steps['width_diff'].diff()  # Max val for width_diff_2nd is the bank?
    
    # Find the top three maximum diff_2nd values and select the one with the lowest height?
    df_top3 = df_steps.nlargest(3, columns='width_diff_2nd')
    
    bank_height = df_steps['height'].iloc[df_top3['height'].idxmin()]
    chan_width = df_steps['width'].iloc[df_top3['height'].idxmin()]
    bank_ang = np.arctan(bank_height/chan_width)
    
#    lst_area[df_steps.index[df_steps.height==1.2][0]]
    
#    sys.exit()
    
#    print('bank height: {}'.format(bank_height))
#    print('channel width: {}'.format(chan_width))
#    
                    
    return bank_height, chan_width, bank_ang

# ===============================================================================
#  Calculates channel width and sinuosity using parallel offset buffering
# ===============================================================================    
def channel_width_from_bank_pixels(df_coords, str_streamlines_path, str_bankpixels_path, str_reachid, cell_size, i_step, max_buff):
    
    print('Channel width from bank pixels -- segmented reaches...')

    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 

    # Successive buffer-mask operations to count bank pixels at certain intervals
    lst_buff=range(cell_size,max_buff,cell_size)    
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_ch_width.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')
    
#    schema_buff = {'geometry': 'Polygon', 'properties': {'buff': 'str'}}
    schema_output = {'geometry': 'LineString', 'properties': {'linkno':'int','ch_wid_total':'float', 'ch_wid_1':'float', 'ch_wid_2':'float', 'dist_sl':'float', 'dist':'float', 'sinuosity':'float'}}                                  
        
    # Access the bank pixel layer...
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:    
        
        # Access the streamlines layer...
        with fiona.open(str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
      
#            progBar.setRange(0, len(streamlines)) 
            
            # Get the crs...
            streamlines_crs = streamlines.crs                
            
#                # For writing out the buffers...
#                with fiona.open(r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\buff_test_ls.shp','w','ESRI Shapefile', schema_buff) as buff_out:                     
                
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
                
                for i_linkno, df_linkno in gp_coords:
                    
#                    progBar.setValue(j)
                    j+=1                    
            
                    i_linkno = int(i_linkno)
#                        max_indx = len(df_linkno.index) - 1
                    
#                    if i_linkno != 1368: continue                
                    
                    print('linkno:  {}'.format(i_linkno))
          
                    # << Analysis by reach segments >>
                    # Set up index array to split up df_linkno into segments (these dictate the reach segment length)...
                    # NOTE:  Reach might not be long enough to break up
                    arr_ind = np.arange(i_step, len(df_linkno.index)+1, i_step) # NOTE: Change the step for resolution?                        
                    lst_dfsegs = np.split(df_linkno, arr_ind)                        
                    
                    for i_seg, df_seg in enumerate(lst_dfsegs): # looping over each reach segment
                        
                        arr_x = df_seg.x.values
                        arr_y = df_seg.y.values

                        try:
                            # Create a line segment from endpts in df_seg...
                            ls = LineString(zip(arr_x, arr_y))                            
                        except:
                            print('Cannot create a LineString using these points, skipping')
                            continue 
                        
                        try:
                            # Calculate straight line distance...
                            dist_sl = np.sqrt((arr_x[0] - arr_x[-1])**2 + (arr_y[0] - arr_y[-1])**2)                     
                        except:
                            print('Error calculated straight line distance')
                            dist_sl = -9999.
                            
                        dist = ls.length                                    
                        sinuosity = dist/dist_sl # ratio of sinuous length to straight line length
                        
                        lst_tally=[]                            

                        for buff_dist in lst_buff:                                                            
                            
                            try:
                                # Watch out for potential geometry errors here...
                                ls_offset_left = ls.parallel_offset(buff_dist, 'left')
                                ls_offset_rt = ls.parallel_offset(buff_dist, 'right')   
                            except:
                                print('Error performing offset buffer')
                            
                            # Buffer errors can result from complicated line geometry... 
                            try:
                                out_left, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_left)], crop=True)  
                            except:
                                print('Left offset error')
                                out_left=np.array([0])
                                
                            try:
                                out_rt, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_rt)], crop=True)
                            except:
                                print('Right offset error')
                                out_rt=np.array([0])
                                
                            num_pixels_left = len(out_left[out_left>0.])
                            num_pixels_rt = len(out_rt[out_rt>0.])
                            
                            # You want the number of pixels gained by each interval...                    
                            tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
                            lst_tally.append(tpl_out)                    
                            df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','interval_left','interval_rt'])


                        
                        # Calculate weighted average                     
                        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                        weighted_avg_left=0 
                        weighted_avg_rt=0
                        n_top=2   
                        
                        try:                        
                            for tpl in df_tally.nlargest(n_top, 'interval_left').iloc[0:2].itertuples():
                                weighted_avg_left += tpl.buffer*(np.float(tpl.interval_left)/np.float(df_tally.nlargest(n_top, 'interval_left').iloc[0:2].sum().interval_left))                        
                        except Exception as e:
                            weighted_avg_left=max_buff
                            print('Left width set to max. Exception: {} \n'.format(e))

                        try:
                            for tpl in df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].itertuples():
                                weighted_avg_rt += tpl.buffer*(np.float(tpl.interval_rt)/np.float(df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].sum().interval_rt))                                                                                                            
                        except Exception as e:
                            weighted_avg_rt=max_buff
                            print('Right width set to max. Exception: {} \n'.format(e))
                       
                        # Write to the output shapefile here...
                        output.write({'properties':{'linkno':i_linkno,'ch_wid_total': weighted_avg_left+weighted_avg_rt,'ch_wid_1': weighted_avg_left,'ch_wid_2': weighted_avg_rt, 'dist_sl':dist_sl,'dist':dist,'sinuosity': sinuosity}, 'geometry':mapping(ls)})                            
                                    
#                    if j > 50: break
    return

# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size, use_wavelet_method):
    
    print('Bank pixels from curvature windows...')
    
    # Convert df_coords x-y to row-col via DEM affine
    # Loop over center row-col pairs accessing the window
    # Now loop over the linknos to get access grid by window...
    
    # << PARAMETERS >>
    cell_size=int(cell_size)
    
    # 3 m...
    w_height=20 # number of rows
    w_width=20  # number of columns
    buff=3 # number of cells    
    curve_thresh=0.30 # good for 3m DEM?    
    
    # 10 m...
#    w_height=20 # number of rows
#    w_width=20  # number of columns
#    buff=3 # number of cells    
#    curve_thresh=0.30 # good for 3m and 10m DEM?     
    
#    j=0
    
#    lst_linknos=[]
#    lst_x1=[]
#    lst_y1=[]
#    lst_x2=[]
#    lst_y2=[]    
    
    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True)   
    
    try:
        
        sigma=1.0 # select the scale sigma=1 --> NOTE:  FIND THIS ON THE FLY? (within each window?)
        g2x1,g2y1 = gauss_kern(sigma) 
        
        with rasterio.open(str_dem_path) as ds_dem:
        
            # Transform to pixel space
            df_coords['col'], df_coords['row'] = ~ds_dem.transform * (df_coords['x'], df_coords['y'])   
            
            df_coords[['row','col']] = df_coords[['row','col']].astype(np.int32)  
            df_coords.drop_duplicates(['col','row'], inplace=True) # rounding to integer
            total_len = len(df_coords.index)
            
            out_meta = ds_dem.meta.copy()      
            out_meta['dtype'] = rasterio.uint8 # no need for float32 for bankpixels to save size of output
            out_meta['compress'] = 'lzw'
            
            arr_bankpts=np.zeros([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])         
            
    #        progBar.setRange(0, len(df_coords.index)) 
            
            for tpl_row in df_coords.itertuples():
                
#                if tpl_row.order != 6:
#                    continue
#                
                if tpl_row.order == 5:
                    w_height=30 # number of rows
                    w_width=30  # number of columns
                if tpl_row.order >= 6:
                    w_height=70
                    w_width=70

#                if tpl_row.linkno != 1368: continue
                
    #            progBar.setValue(j)
                j+=1
                
                print('{} | {} -- {}'.format(tpl_row.linkno, j, total_len))
                
                row_min = np.int(tpl_row.row - np.int(w_height/2))
                row_max = np.int(tpl_row.row + np.int(w_height/2))
                col_min = np.int(tpl_row.col - np.int(w_width/2))
                col_max = np.int(tpl_row.col + np.int(w_width/2))            
                
                # Now get the DEM specified by this window as a numpy array...
                w = ds_dem.read(1, window=((row_min, row_max),(col_min, col_max))) 
                                
                # Then extract the internal part of the window that contains the rotated window??
                w[w>9999999.0] = 0.0 # NoData values may have been corrupted by preprocessing?
                w[w<-9999999.0] = 0.0
                
                if np.size(w) > 9: # make sure a window of appropriate size was returned from the DEM
                
                    if use_wavelet_method:
                        # === Wavelet Curvature from Chandana ===
                        gradfx1 = signal.convolve2d(w, g2x1, boundary='symm', mode='same') 
                        gradfy1 = signal.convolve2d(w, g2y1, boundary='symm', mode='same')     
                        
                        w_curve=gradfx1+gradfy1                    
    
                        # Pick out bankpts...                    
                        w_curve[w_curve<np.max(w_curve)*curve_thresh] = 0.
                    
#                    plt.imshow(lum_img, cmap="hot") #, vmin=-2, vmax=2)   
#                    plt.colorbar()    
#                    fig = plt.figure()
#                    n, bins, patches=plt.hist(w_curve.flatten(), bins=256, range=(-5, 5), fc='k', ec='k')    # get the histogram
#                    plt.xlabel('Curvature')
#                    plt.ylabel('Probability')
#                    plt.title('Histogram')                    
#                    plt.show()      
#                    
#                    if j > 4:
#                        sys.exit()
                    # =======================================
                    
                    else:
                        # === Mean Curvature ====================                                
                        Zy, Zx = np.gradient(w, cell_size)
                        Zxy, Zxx = np.gradient(Zx, cell_size)
                        Zyy, _ = np.gradient(Zy, cell_size)                                                
                  
                        try:
                            w_curve = (Zx**2 + 1)*Zyy - 2*Zx*Zy*Zxy + (Zy**2 + 1)*Zxx
                            w_curve = -w_curve/(2*(Zx**2 + Zy**2 + 1)**(1.5))
                        except:
                            print('Error calculating Curvature in window...skipping')
                            continue
                                                
                        w_curve[w_curve<np.max(w_curve)*curve_thresh] = 0.
                        # =======================================
                    
                    w_curve[w_curve<-99999999.] = 0.
                    w_curve[w_curve>99999999.] = 0.
                    
                    w_curve[w_curve>0.] = 1.

                    # Note:  This assumes that the w_curve window is the specified size, which is not always the case for edge reaches...
                    arr_bankpts[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w_curve[buff:w_height-buff, buff:w_width-buff]
                             
                    out_meta['nodata'] = 0.
                    
            print('Writing bank pixels .tif...')
            with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
                dest.write(arr_bankpts.astype(rasterio.uint8), indexes=1)
#                dest.write(arr_bankpts, indexes=1)
                               
    except Exception as e:
        print('\r\nError in bankpixels_from_curvature_window. Exception: {} \n'.format(e))        
            
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
#            this_order = tpl_row.strmord
            
#            print('this_linkno: {} | index: {}'.format(this_linkno, tpl_row.Index))
            
#            if tpl_row.Index == 2254:
#                print('pause')

            # A list to store the total number of indices/blocks in a Xn...
            lst_total_cnt = []

            arr_elev = tpl_row.elev
            arr_elev = arr_elev[arr_elev != np.float32(nodata_val)]
                
            # Normalize elevation to zero...
            thisxn_norm = arr_elev - np.min(arr_elev)                
            thisxn_norm = tpl_row.elev - np.min(tpl_row.elev)    
            # Below is if you're using the breached DEM...
#            if this_order < 5:
#                thisxn_norm = arr_elev - np.min(arr_elev)                
#                thisxn_norm = tpl_row.elev - np.min(tpl_row.elev)
#            else:
#                # for order>5...(THIS ASSUMES YOU'RE USING THE BREACHED DEM, may not be necessary otherwise)
#                thisxn_norm = tpl_row.elev - np.partition(tpl_row.elev, 2)[2]
#                # then any negatives make zero...
#                thisxn_norm[thisxn_norm<0]=0

#            p_vert_zero=0.2
            # Loop from zero to max(this_xn_norm) using a pre-defined vertical step (0.2 m?)...
            for this_slice in np.arange(0., np.max(thisxn_norm), param_ivert):

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
        pass

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
    
    # NOTE:  Use MultiProcessing here over arr_strides??

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
                df_map['lfbank_x'], df_map['lfbank_y'] = ds_dem.transform * (df_map['lfbank_col'], df_map['lfbank_row']) 
                df_map['rtbank_x'], df_map['rtbank_y'] = ds_dem.transform * (df_map['rtbank_col'], df_map['rtbank_row'])
                        
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
    
    min_nodata_thresh = -99999.0
    max_nodata_thresh = 99999.0     
    
    print('Reading and interpolating elevation along Xn\'s...')
    
    lst_linknos=[]
    lst_x1=[]
    lst_y1=[]
    lst_x2=[]
    lst_y2=[]
    lst_strmord=[]

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
            lst_strmord.append(line['properties']['strmord'])
            
            
    df_coords = pd.DataFrame({'linkno':lst_linknos, 'x1':lst_x1, 'y1':lst_y1, 'x2':lst_x2, 'y2':lst_y2, 'strmord':lst_strmord})
         
    # Now loop over the linknos to get access grid by window...
    with rasterio.open(str(str_dem_path)) as ds_dem:
        
        nodata_val = ds_dem.nodata # NODATA val must be defined for this to return anything
    
        # Transform to pixel space
        df_coords['col1'], df_coords['row1'] = ~ds_dem.transform * (df_coords['x1'], df_coords['y1'])
        df_coords['col2'], df_coords['row2'] = ~ds_dem.transform * (df_coords['x2'], df_coords['y2'])
        
        ### OR...
        gp_coords = df_coords.groupby('linkno')
        
        lst_all_zi=[]
        j=0
        
        for linkno, df_linkno in gp_coords:
            
#            print(linkno)
            
#            if linkno != 120:
#                continue
            
            row_min = int(df_linkno[['row1','row2']].min(axis=0).min())
            row_max = int(df_linkno[['row1','row2']].max(axis=0).max())
            col_min = int(df_linkno[['col1','col2']].min(axis=0).min())
            col_max = int(df_linkno[['col1','col2']].max(axis=0).max())    
            strmord = int(df_linkno.strmord.iloc[0])
            
            # Now get the DEM specified by this window as a numpy array...
            w = ds_dem.read(1, window=((row_min, row_max+1),(col_min, col_max+1))) 
            
            w_min = np.min(w)
            w_max = np.max(w)                      
            
            if w_min < min_nodata_thresh:
                nodata_val = w_min
            elif w_max > max_nodata_thresh:
                nodata_val = w_max                
                        
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
                try:            
                    arr_zi = w[lst_xnrow.astype(np.int), lst_xncol.astype(np.int)]   # nearest-neighbor                
                except:
                    continue # Just skip this Xn altogether?  Could be smarter (eg, are the indices > len(w)?)
#                lst_zi = ndimage.map_coordinates(w, np.vstack((lst_xnrow, lst_xncol)), order=1, mode='nearest') # use this for additonal interpolation options (ie, cubic, bilinear, etc)
                
#                plt.plot(np.arange(len(zi)), zi)   
                
                # Remove possible no data values...NOTE:  They may not be defined in the original file
                arr_zi = arr_zi[arr_zi != np.float32(nodata_val)]
                
                if arr_zi.size < 5: continue # if it only has less than 5 elevation measurements along this Xn, skip it
                                
                # Convert these from window row/col to raster row/col for bankpt use...
                for i, xnrow in enumerate(lst_xnrow):
                    lst_xnrow[i] = lst_xnrow[i] + row_min
                    lst_xncol[i] = lst_xncol[i] + col_min

                tpl_out = (linkno, arr_zi, lst_xnrow, lst_xncol, strmord)
                lst_all_zi.append(tpl_out)
#                i += 1
           
#    print('\tTotal Xn\'s:  {}'.format(i))    
#    print('\tTime interpolating elevation along Xn\'s:  ' + str(timeit.default_timer() - start_time))

    return pd.DataFrame(lst_all_zi, columns=['linkno','elev','xn_row','xn_col','strmord'])
    
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
    
#    slopeCutoffVertical = 20 # just a threshold determining when to call a Xn vertical (if it's above this, make it vertical. Otherwise things get whacky?)

    lst_xnrowcols = [] # the final output, a list of tuples of XY coordinate pairs for all Xn's for this reach

    XnCntr = 0
#    m_init = 0
    
    gp_coords = df_coords.groupby('linkno')  
           
    # Create the Xn shapefile for writing...
#    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int', 'endpt1_x':'float', 'endpt1_y':'float', 'endpt2_x':'float', 'endpt2_y':'float'}} 
    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int', 'strmord':'int'}} 
    
    print('Building and Writing Cross Section File...')
    with fiona.open(str_xns_path, 'w', driver='ESRI Shapefile', crs=streamlines_crs, schema=test_schema) as chan_xns:
        
        for i_linkno, df_linkno in gp_coords:
            
            i_linkno = int(i_linkno)
            i_order = int(df_linkno.order.iloc[0])
            
#            if i_linkno != 177:
#                continue

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
            if not(bool_isvalley):
                
#                if i_order != 6: continue
                
                if i_order == 1:
                    p_xnlength=20
                    p_fitlength = 3
                elif i_order == 2:
                    p_xnlength=23
                    p_fitlength = 6
                elif i_order == 3:
                    p_xnlength=40
                    p_fitlength = 9
                elif i_order == 4:
                    p_xnlength=60 
                    p_fitlength = 12
                elif i_order == 5:
                    p_xnlength=80  
                    p_fitlength = 15
                elif i_order >= 6:
                    p_xnlength=250 
                    p_fitlength = 20
#                elif i_order == 7:
#                    p_xnlength=130 
#                    p_fitlength = 21
#                elif i_order == 8:
#                    p_xnlength=150  
#                    p_fitlength = 24
    
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
                
                # Send it the endpts of what you to draw a perpendicular line to...
                lst_xy = build_xns(lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow, p_xnlength) # returns a list of two endpoints
        
        
#                if (max(lstThisSegmentCols) - min(lstThisSegmentCols) < 3):
#                    m_init = 9999.0
#                elif (max(lstThisSegmentRows) - min(lstThisSegmentRows) < 3):
#                    m_init = 0.0001
#                else:
#                    LinFit = np.polyfit(lstThisSegmentCols,lstThisSegmentRows,1) # NOTE: could just use basic math here instead?!
#                    m_init = LinFit[0]
#        
#                # Check for zero or infinite slope...
#                if m_init == 0:
#                    m_init = 0.0001
#                elif isinf(m_init):
#                    m_init = 9999.0
#        
#                # Find the orthogonal slope...
#                m_ortho = -1/m_init
#                          
#                xn_steps = [-float(p_xnlength),float(p_xnlength)] # just the end points
#        
#                lst_xy=[]
#                for r in xn_steps:
#        
#                    # Make sure it's not too close to vertical...
#                    # NOTE X-Y vs. Row-Col here...
#                    if (abs(m_ortho) > slopeCutoffVertical):
#                        tpl_xy = (midPtCol, midPtRow+r)
#                        
#                    else:
#                        fit_col_ortho = (midPtCol + (float(r)/(sqrt(1 + m_ortho**2))))                       
#                        tpl_xy = float(((midPtCol + (float(r)/(sqrt(1 + m_ortho**2)))))), float(((m_ortho)*(fit_col_ortho-midPtCol) + midPtRow))
#                    
#                    lst_xy.append(tpl_xy)
                         
                XnCntr = XnCntr + 1
                
                # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                line = {'type': 'LineString', 'coordinates':lst_xy}
#                prop = {'linkno': i_linkno, 'endpt1_x':lst_xy[0][0], 'endpt1_y':lst_xy[0][1], 'endpt2_x':lst_xy[1][0], 'endpt2_y':lst_xy[1][1]}
                prop = {'linkno': i_linkno, 'strmord': i_order}
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
 
    p_interp_spacing = int(cell_size) #3 # larger numbers would simulate a more smoothed reach | NOTE: Hardcode this = grid resolution?
    j=0 # prog bar
    
    # Open the streamlines shapefile...
    with fiona.open(str(str_streams_filepath), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
   
        # Get the crs...
        streamlines_crs = streamlines.crs  
#        str_proj4 = crs.to_string(streamlines.crs)         

#            progBar.setRange(0,len(streamlines))

#        # ==================================
#        progdialog = QtGui.QProgressDialog("Getting stream coords from features...", "Cancel", 0, len(streamlines))
#        progdialog.setWindowTitle("FACET")
#        progdialog.setWindowModality(QtCore.Qt.WindowModal)
#        progdialog.resize(350, 110)
#        progdialog.show()       
#        # ==================================
        
        for line in streamlines:
            
#           # ====================================== 
#           QtCore.QCoreApplication.processEvents()
#           if progdialog.wasCanceled():
#               break     
#                       
#           progdialog.setValue(j)
#           # ======================================
           
           j+=1
#               self.emit(QtCore.SIGNAL("update(int)"), int(100*len(streamlines)/j)) 
           
#           print('{} | {}'.format(i_linkno, j))
           line_shply = LineString(line['geometry']['coordinates'])
          
           length = line_shply.length # units depend on crs
                      
           if length > 9:       
               
               i_linkno = line['properties'][str_reachid]           
               i_order = line['properties'][str_orderid]
               
#               if i_linkno != 1368: continue
              
               # Smoothing higher order reaches via Shapely...
               if i_order <= 3:
                   line_shply = line_shply.simplify(5.0, preserve_topology=False)
               elif i_order == 4:
                   line_shply = line_shply.simplify(10.0, preserve_topology=False)
               elif i_order == 5:
                   line_shply = line_shply.simplify(20.0, preserve_topology=False)
               elif i_order >= 6:
                   line_shply = line_shply.simplify(30.0, preserve_topology=False)                
       
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


    
    
    