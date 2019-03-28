# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:51:45 2016

@author: sam.lamont
"""

# import math
import subprocess
# import timeit
import numpy as np
from numpy import asarray
import logging
from scipy.ndimage import label  # TEST TEST

# from scipy.stats import gaussian_kde # TEST
# from scipy.optimize import curve_fit # TEST
from scipy import signal
from scipy.ndimage import percentile_filter

import os
# import ntpath
from math import atan, ceil, floor
import sys
from math import isinf, sqrt
import rasterio
import rasterio.mask
from rasterio.warp import transform
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import shapes
from rasterio import features

from functools import partial
import multiprocessing as mp

# import whitebox

# wbt = whitebox.WhiteboxTools()


# import matplotlib
# matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

plt.style.use('ggplot')
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm

# from affine import Affine
import pandas as pd
import geopandas as gpd
# from scipy import ndimage
# from shapely.geometry import Point,
from shapely.geometry import shape, mapping, LineString, MultiLineString, Point, MultiPoint
from shapely.ops import split
# from jenks import jenks
# from PyQt4 import QtGui, QtCore

# import scipy.io as sio

import fiona

# from fiona import collection
# from fiona.crs import from_epsg
# import jenkspy

# import gospatial as gs

np.seterr(over='raise')


# plt.ion()

# from shapely import speedups
# if speedups.available:
#    speedups.enable() # enable performance enhancements written in C

# LOGGER SETUP #
# logging.basicConfig(filename='example.log', filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
# does this need to be set? level=logging.DEBUG,
# logging.debug('This message should go to the log file')
# logging.info('So should this')
# logging.warning('And this, too')
# END LOGGER SETUP #

# ===============================================================================
#  Utility functions
# ===============================================================================
# def open_mat_files():
#    
#    # Loads a MATLAB file (.mat) using SciPy...
#    test = sio.loadmat(r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\field_data\FieldData\smith\AtGageAB_US.mat')
#    
#    arr_data = test['AtGageAB_US'] #numpy array
#    
#    return
# def rotate(x_center, x_orig, y_center, y_orig, theta):
#    x_rot = np.cos(theta)*(x_orig-x_center) - np.sin(theta)*(y_orig-y_center) + x_center
#    y_rot = np.sin(theta)*(x_orig-x_center) + np.cos(theta)*(y_orig-y_center) + y_center    
#    return x_rot, y_rot

# lstThisSegmentRows, lstThisSegmentCols, midpt_x, midpt_y, p_fpxnlen

def rasterize_gdf(gdf, str_ingrid_path, str_tempgrid, str_outgrid, huc_poly):
    '''
    Thanks to:  https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
    '''
    with rasterio.open(str_ingrid_path) as rst:
        meta = rst.meta.copy()
    meta.update(compress='lzw')
    meta.update(dtype=rasterio.int32)
    meta.update(nodata=0)

    with rasterio.open(str_tempgrid, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf['LINKNO']))

        arr_burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)

        out.write_band(1, arr_burned)

    #        w_pts, trans_pts = rasterio.mask.mask(ds_pts, huc_mask, crop=True, nodata=nodata, all_touched=True)
    #        out_mask, out_trans=rasterio.mask.mask(out, [mapping(huc_poly)], crop=True, nodata=0) #, all_touched=True)
    #        out.write_band(1, out_mask)
    #
    #        with rasterio.open(str_outgrid, 'w', **meta) as ds_out:
    #            ds_out.write_band(1, out_mask[0])

    return


def compress_grids(str_in_grid, str_out_grid, logger):
    '''
    :param str_in_grid:
    :param str_out_grid:
    :param logger:
    :return:
    '''
    # Access the interp pts .tif file:
    with rasterio.open(str_in_grid) as ds_in:
        meta_in = ds_in.meta.copy()
        meta_in.update(compress='lzw')
        meta_in.update(dtype=rasterio.float32)

        arr = ds_in.read(1)

        logger.info('Applying percentile filter...')
        arr = percentile_filter(arr, 50., size=3)

    logger.info('Compressing and saving as regular TIFF...')
    with rasterio.Env(GDAL_CACHEMAX=256, GDAL_NUM_THREADS='ALL_CPUS', BIGTIFF='NO'):
        with rasterio.open(str_out_grid, 'w', tiled=True, blockxsize=512, blockysize=512, **meta_in) as ds_out:
            ds_out.write(arr.astype(rasterio.float32), indexes=1)

        # ================================================================================


#   For wavelet curvature calculation (Chandana)
# ================================================================================
def gauss_kern(sigma):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    sigma = int(sigma)

    x, y = np.mgrid[-5 * sigma:5 * sigma, -5 * sigma:5 * sigma]

    g2x = (1 - x ** 2 / sigma ** 2) * np.exp(-(x ** 2 + y ** 2) / 2 / sigma ** 2) * 1 / np.sqrt(
        2 * np.pi * sigma ** 2) / 4 / sigma * np.exp(float(0.5));
    g2y = (1 - y ** 2 / sigma ** 2) * np.exp(-(x ** 2 + y ** 2) / 2 / sigma ** 2) * 1 / np.sqrt(
        2 * np.pi * sigma ** 2) / 4 / sigma * np.exp(float(0.5));

    return g2x, g2y


# ================================================================================
#   For 2D cross sectional measurement
# ================================================================================    
def build_xns(lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow, p_xnlength):
    slopeCutoffVertical = 20  # another check

    # Find initial slope...
    if (abs(lstThisSegmentCols[0] - lstThisSegmentCols[-1]) < 3):
        m_init = 9999.0
    elif (abs(lstThisSegmentRows[0] - lstThisSegmentRows[-1]) < 3):
        m_init = 0.0001
    else:
        m_init = (lstThisSegmentRows[0] - lstThisSegmentRows[-1]) / (lstThisSegmentCols[0] - lstThisSegmentCols[-1])

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
    m_ortho = -1 / m_init

    xn_steps = [-float(p_xnlength), float(p_xnlength)]  # just the end points

    lst_xy = []
    for r in xn_steps:

        # Make sure it's not too close to vertical...
        # NOTE X-Y vs. Row-Col here...
        if (abs(m_ortho) > slopeCutoffVertical):
            tpl_xy = (midPtCol, midPtRow + r)

        else:
            fit_col_ortho = (midPtCol + (float(r) / (sqrt(1 + m_ortho ** 2))))
            tpl_xy = float(((midPtCol + (float(r) / (sqrt(1 + m_ortho ** 2)))))), float(
                ((m_ortho) * (fit_col_ortho - midPtCol) + midPtRow))

        lst_xy.append(tpl_xy)  # A list of two tuple endpts

    return lst_xy


def get_xn_length_by_order(i_order, bool_isvalley):
    # Settings for channel cross-sections:
    if not (bool_isvalley):
        if i_order == 1:
            p_xnlength = 20
            p_fitlength = 3
        elif i_order == 2:
            p_xnlength = 23
            p_fitlength = 6
        elif i_order == 3:
            p_xnlength = 40
            p_fitlength = 9
        elif i_order == 4:
            p_xnlength = 60
            p_fitlength = 12
        elif i_order == 5:
            p_xnlength = 80
            p_fitlength = 15
        elif i_order >= 6:
            p_xnlength = 150
            p_fitlength = 20
            # Settings for floodplain cross-sections:
    elif bool_isvalley:

        if i_order == 1:
            p_xnlength = 50
            p_fitlength = 5
        elif i_order == 2:
            p_xnlength = 75
            p_fitlength = 8
        elif i_order == 3:
            p_xnlength = 100
            p_fitlength = 12
        elif i_order == 4:
            p_xnlength = 150
            p_fitlength = 20
        elif i_order == 5:
            p_xnlength = 200
            p_fitlength = 30
        elif i_order >= 6:
            p_xnlength = 500
            p_fitlength = 40

    return p_xnlength, p_fitlength


def get_cell_size(str_grid_path):
    with rasterio.open(str(str_grid_path)) as ds_grid:
        cs_x, cs_y = ds_grid.res

    return cs_x


def rugosity(arr, res, logger):
    '''
    Actual 3D area divided by 2D planar area gives a measure
    of terrain complexity or roughness
    '''
    try:
        area3d = ((res ** 2) * (1 + np.gradient(arr) ** 2) ** 0.5).sum()  # actual surface area
        area2d = len(arr) * res ** 2  # planar surface area
        rug = area3d / area2d
    except:
        logger.info(f'Error in rugosity. arr.shape: {arr.shape}')
        return -9999.

    return rug


# ==========================================================================
#   Reproject a grid layer using rasterio 
# ==========================================================================    
def define_grid_projection(str_source_grid, dst_crs, dst_file, logger):
    logger.info('Defining grid projection...')
    with rasterio.open(str_source_grid, 'r') as src:
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
        })
        arr_src = src.read(1)

    with rasterio.open(dst_file, 'w', **kwargs) as dst:
        dst.write(arr_src, indexes=1)


# ==========================================================================
#   Reproject a grid layer using rasterio 
# ==========================================================================    
def reproject_grid_layer(str_source_grid, dst_crs, dst_file, resolution, logger):
    """ Resolution is a pixel value as a tuple"""
    logger.info('Reprojecting grid layer...')
    with rasterio.open(str_source_grid) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(dst_file, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear)


# ==========================================================================
#   Reproject a vector layer using geopandas 
# ==========================================================================
def reproject_vector_layer(str_path_to_file, str_target_proj4, logger):
    logger.info('Reprojecting vector layer...')

    gdf = gpd.read_file(str_path_to_file)
    gdf = gdf.to_crs(str_target_proj4)

    str_out_path = str_path_to_file[:-4] + '_proj.shp'
    gdf.to_file(str_out_path)

    return str_out_path


# ==========================================================================
#   For dissolving line features    
# ==========================================================================    
def dissolve_line_features(str_lines_path, output_filename, logger):
    logger.info('Dissolve line features...')

    lst_all = []

    with fiona.open(str_lines_path) as lines:
        crs = lines.crs

        schema = {'geometry': 'MultiLineString', 'properties': {'linkno': 'int:6'}}
        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:
            for line in lines:
                lst_all.append(line['geometry']['coordinates'])

            output.write(
                {'properties': {'linkno': 9999}, 'geometry': {'type': 'MultiLineString', 'coordinates': lst_all}})

    return


# ==========================================================================
#   For points along line feature at uniform distance    
# ==========================================================================    
def points_along_line_features(str_diss_lines_path, output_filename, logger):
    logger.info('Points along line features...')

    p_interp_spacing = 1000

    with fiona.open(str_diss_lines_path) as line:
        crs = line.crs
        line = line[0]

        schema = {'geometry': 'Point', 'properties': {'id': 'int:6'}}
        with fiona.open(output_filename, 'w', 'ESRI Shapefile', schema, crs) as output:
            line_shply = MultiLineString(line['geometry']['coordinates'])

            length = line_shply.length  # units depend on crs

            step_lens = np.arange(0, length, p_interp_spacing)  # p_interp_spacing in projection units?

            for i, step in enumerate(step_lens):  # lambda here instead?

                i_pt = np.array(line_shply.interpolate(step))

                output.write({'properties': {'id': i}, 'geometry': {'type': 'Point', 'coordinates': i_pt}})

    return


# ==========================================================================
#   For points along line feature at uniform distance    
# ==========================================================================    
def taudem_gagewatershed(str_pts_path, str_d8fdr_path, logger):
    logger.info('Points along line features...')

    inputProc = str(4)

    str_output_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_test_sheds.tif'

    cmd = 'mpiexec -n ' + inputProc + ' GageWatershed -p ' + '"' + str_d8fdr_path + '"' + ' -o ' + '"' + str_pts_path + '"' + ' -gw ' + '"' + str_output_path + '"'

    # Submit command to operating system
    logger.info('Running TauDEM GageWatershed...')
    os.system(cmd)
    # Capture the contents of shell command and logger.info it to the arcgis dialog box
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    message = "\n"
    for line in process.stdout.readlines():
        if isinstance(line, bytes):  # true in Python 3
            line = line.decode()
        message = message + line
    logger.info(message)

    return


# ==========================================================================
#   For clipping features  
# ==========================================================================         
def clip_features_using_grid(str_lines_path, output_filename, str_dem_path, logger):
    logger.info('Clipping streamlines to site DEM...')
    #    # Build the output file name...
    #    path_to_dem, dem_filename = os.path.split(str_dem_path)
    #    output_filename = path_to_dem + '\\' + dem_filename[:-4]+'_nhdhires.shp'

    # Polygonize the raster DEM with rasterio...    
    with rasterio.open(str(str_dem_path)) as ds_dem:
        arr_dem = ds_dem.read(1)

    arr_dem[arr_dem > 0] = 100
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
#    If a callback is not provided, it will simply logger.info the output stream.
#    A provided callback allows for custom processing of the output stream.
# ===============================================================================
def callback(s, logger):
    try:
        # logger.info("{}".format(s))
        if "%" in s:
            str_array = s.split(" ")
            #            label = s.replace(str_array[len(str_array)-1], "")
            progress = int(str_array[len(str_array) - 1].replace("%", "").strip())
            logger.info("Progress: {}%".format(progress))
        else:
            if "error" in s.lower():
                logger.info("ERROR: {}".format(s))
            else:
                logger.info("{}".format(s))
    except:
        logger.info(s)

    # ===============================================================================


#  Create weight file for TAuDEM D8 FAC
# ===============================================================================
def create_wg_from_streamlines(str_streamlines_path, str_dem_path, str_danglepts_path, logger):
    logger.info('Creating weight grid from streamlines...')

    lst_coords = []
    lst_pts = []
    lst_x = []
    lst_y = []

    with fiona.open(str_streamlines_path) as lines:

        streamlines_crs = lines.crs  # to use in the output grid

        # Get separate lists of start and end points...
        for line in lines:
            if line['geometry']['type'] == 'LineString':  # Make sure it's a LineString
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
        out_meta.update(compress='lzw')
        out_meta.update(dtype=rasterio.int16)
        out_meta.update(nodata=-9999)
        out_meta.update(crs=lines.crs)  # shouldn't be necessary

        # Construct the output array...
        arr_danglepts = np.zeros([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])

        tpl_pts = transform(streamlines_crs, out_meta['crs'], lst_x, lst_y)
        lst_dangles = zip(tpl_pts[0], tpl_pts[1])

        for coords in lst_dangles:
            col, row = ~ds_dem.transform * (
            coords[0], coords[1])  # BUT you have to convert coordinates from hires to dem
            try:
                arr_danglepts[int(row), int(col)] = 1
            except:
                continue
    #            # Could save these points here for later use building the streamlines??
    #            tpl_pts=(row,col)
    #            lst_finalpts_rowcols.append(tpl_pts)

    # Now write the new grid using this metadata...
    with rasterio.open(str_danglepts_path, "w", **out_meta) as dest:
        dest.write(arr_danglepts, indexes=1)

    return


# ===============================================================================
#  Get row/col coords of start pts from weight grid
# ===============================================================================
# def get_rowcol_from_wg(str_danglepts_path):
#    
#    logger.info('Reading weight grid...')
#    with rasterio.open(str_danglepts_path) as ds_wg:
#        wg_rast = ds_wg.read(1) # numpy array
#        wg_crs = ds_wg.crs
#        wg_affine = ds_wg.affine
#        
#    logger.info('Getting indices...')
#    row_cols = np.where(wg_rast>0)    
#    
#    return row_cols

# ===============================================================================
#  For calling GoSpatial funcs
# ===============================================================================
def default_callback(str, logger):
    logger.info(str)


def run_gospatial_whiteboxtool(tool_name, args, exe_path, exe_name, wd, logger, callback=default_callback):
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

        ps = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1,
                              universal_newlines=True)

        while True:
            line = ps.stdout.readline()
            if line != '':
                callback(line.strip())
            else:
                break

        return 0
    except Exception as e:
        logger.info(e)
        return 1


def run_rust_whiteboxtool(tool_name, args, exe_path, exe_name, wd, callback=default_callback):
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
def preprocess_dem(str_dem_path, str_nhdhires_path, dst_crs, str_mpi_path, str_taudem_path, str_whitebox_path,
                   run_whitebox, run_wg, run_taudem, logger):
    try:

        # Split DEM path and filename...  # NOT OS INDEPENDENT??
        path_to_dem, dem_filename = os.path.split(str_dem_path)

        # Set this for a custom output folder:
        #        path_to_dem = '/home/sam/drb_preprocess_2018.08.28'

        inputProc = str(2)  # number of cores to use for TauDEM processes

        # << Define all filenames here >>
        str_danglepts_path = os.path.join(path_to_dem, dem_filename[:-4] + '_wg.tif')

        dem_filename_tif = dem_filename[:-4] + '.tif'
        breach_filename_dep = dem_filename[:-4] + '_breach.dep'
        breach_filename_tif = dem_filename[:-4] + '_breach.tif'
        breach_filepath_tif = os.path.join(path_to_dem, breach_filename_tif)
        breach_filepath_tif_proj = breach_filepath_tif[:-4] + '_proj.tif'

        str_dem_path_tif = os.path.join(path_to_dem, dem_filename[:-4] + '.tif')

        #        fel = path_to_dem + '\\' + dem_filename[:-4]+'_breach.tif'

        fel_pitremove = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_fel.tif')
        p = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_p.tif')
        sd8 = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_sd8.tif')

        ad8_wg = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_ad8_wg.tif')
        #        wtgr = os.path.join(str_danglepts_path)
        ad8_no_wg = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_ad8_no_wg.tif')
        ord_g = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_ord_g.tif')
        tree = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_tree')
        coord = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_coord')
        net = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_net.shp')
        w = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_w.tif')
        slp = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_slp.tif')
        ang = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_ang.tif')
        dd = os.path.join(path_to_dem, breach_filename_tif[:-4] + '_hand.tif')

        #        # ==================== TauDEM Paths =========================
        #        # Hardcode paths from user input...
        #        mpipath = os.path.join(str_mpi_path) #r'"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"'
        #        d8flowdir = '"' + str_taudem_path + '\D8FlowDir.exe"' # r'"C:\Program Files\TauDEM\TauDEM5Exe\D8FlowDir.exe"'
        #        areaD8 = '"' + str_taudem_path + '\AreaD8.exe"'
        #        streamnet = '"' + str_taudem_path + '\StreamNet.exe"'
        #        dinfflowdir = '"' + str_taudem_path + '\DinfFlowDir.exe"'
        #        dinfdistdown = '"' + str_taudem_path + '\DinfDistDown.exe"'

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
            # Whitebox python:
            wbt.set_working_dir(path_to_dem)
            #            wbt.feature_preserving_denoise(dem_filename_tif, "smoothed.tif", filter=9)
            wbt.breach_depressions(dem_filename, breach_filename_tif)

        #            logger.info('Whitebox .exe path:  ' + 'r' + '"' + str_whitebox_path + '"')
        #            str_whitebox_dir, str_whitebox_exe = os.path.split(str_whitebox_path)
        #
        #            # << Run the BreachDepressions tool, specifying the arguments >>
        #            name = "BreachDepressions"
        #            args = [dem_filename_tif, breach_filename_dep, '-1', '-1', 'True', 'True'] # GoSpatial verion. NOTE:  Make these last four variables accessible to user?
        #            args = ['--dem='+dem_filename, '-o='+breach_filename_dep] # Rust version
        #
        #            ret = run_gospatial_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
        #            ret = run_rust_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
        #            if ret != 0:
        #                logger.info("ERROR: return value={}".format(ret))
        #
        #            #  << Convert .dep to .tif here? >>  NOTE:  Only for DRB hack when using .dep files
        #            name = "WhiteBox2GeoTiff"
        #            args = [breach_filename_dep, breach_filename_tif]
        #
        #            ret = run_gospatial_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)
        #            if ret != 0:
        #                logger.info("ERROR: return value={}".format(ret))
        #
        #            # The converted TIFF file is saved without a crs, so save a projected version:
        #            define_grid_projection(breach_filepath_tif, dst_crs, breach_filepath_tif_proj)
        #
        #            # Remove native Whitebox files and unprojected tif:
        #            dep_path=path_to_dem + '\\' + breach_filename_dep
        #            tas_path=path_to_dem + '\\' + breach_filename_dep[:-3]+'tas'
        #            os.remove(dep_path)
        #            os.remove(tas_path)
        #            os.remove(breach_filepath_tif)

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
        #                logger.info("ERROR: return value={}".format(ret))

        # Rust version...
        #            name = "ConvertRasterFormat"
        #            args = ['--input='+dem_filename, '-o='+dem_filename_tif]
        #            ret = run_rust_whiteboxtool(name, args, str_whitebox_dir, str_whitebox_exe, path_to_dem, callback)

        if run_wg:
            create_wg_from_streamlines(str_nhdhires_path, p, str_danglepts_path)

        if run_taudem:

            # Testing...
            #            mpipath = 'r' +  '"' + mpipath + '"'
            #            fel = 'r' +  '"' + fel + '"'

            #            logger.info('mpipath: ' + mpipath)
            #            logger.info('fel: ' + fel_breach)
            #            logger.info(' ')

            #            # ==============  << 1. Pit Filling with TauDEM >> ================
            #            cmd = 'mpiexec' + ' -n ' + inputProc + ' PitRemove -z ' + '"' + breach_filepath_tif + '"' + ' -fel ' + '"' + fel_pitremove + '"'
            #
            #            # Submit command to operating system
            #            logger.info('Running TauDEM PitRemove...')
            #            os.system(cmd)
            #
            #            # Capture the contents of shell command and logger.info it to the arcgis dialog box
            #            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            #
            #            # Get some feedback from the process to logger.info out...
            #            message = "\n"
            #            for line in process.stdout.readlines():
            #                line = line.decode()
            #                if isinstance(line, bytes):	   # true in Python 3
            #                    line = line.decode()
            #                message = message + line
            #            logger.info(message)
            #
            #            # ==============  << 2. D8 FDR with TauDEM >> ================       YES
            #            cmd = '"' + mpipath + '"' + ' -n ' + inputProc + ' ' + d8flowdir + ' -fel ' + '"' + str_dem_path + '"' + ' -p ' + '"' + p + '"' + \
            #                  ' -sd8 ' + '"' + sd8 + '"'
            #
            #            cmd = 'mpiexec' + ' -n ' + inputProc + ' d8flowdir -fel ' + '"' + str_dem_path + '"' + ' -p ' + '"' + p + '"' + \
            #                  ' -sd8 ' + '"' + sd8 + '"'
            #
            #            # Submit command to operating system
            #            logger.info('Running TauDEM D8 Flow Direction...')
            #            os.system(cmd)
            #
            #            # Capture the contents of shell command and logger.info it to the arcgis dialog box
            #            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            #
            #            # Get some feedback from the process to logger.info out...
            #            message = "\n"
            #            for line in process.stdout.readlines():
            #                line = line.decode()
            #                if isinstance(line, bytes):	   # true in Python 3
            #                    line = line.decode()
            #                message = message + line
            #            logger.info(message)
            #
            #    #        # ============= << 3.a AD8 with weight grid >> ================        YES
            #    #        cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_wg + '"'  + ' -wg ' + '"' + wtgr + '"'  + ' -nc '
            #            cmd = 'mpiexec' + ' -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_wg + '"'  + ' -wg ' + '"' + str_danglepts_path + '"'  + ' -nc '
            #
            #            # Submit command to operating system
            #            logger.info('Running TauDEM D8 FAC (with weight grid)...')
            #            os.system(cmd)
            #            # Capture the contents of shell command and logger.info it to the arcgis dialog box
            #            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            #
            #            message = "\n"
            #            for line in process.stdout.readlines():
            #                if isinstance(line, bytes):	    # true in Python 3
            #                    line = line.decode()
            #                message = message + line
            #            logger.info(message)

            # ============= << 3.b AD8 no weight grid >> ================
            #        cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_no_wg + '"'  +  ' -nc '
            cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8_no_wg + '"' + ' -nc '

            # Submit command to operating system
            logger.info('Running TauDEM D8 FAC (no weights)...')
            os.system(cmd)
            # Capture the contents of shell command and logger.info it to the arcgis dialog box
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            message = "\n"
            for line in process.stdout.readlines():
                if isinstance(line, bytes):  # true in Python 3
                    line = line.decode()
                message = message + line
            logger.info(message)

            # ============= << 4 StreamReachandWatershed with TauDEM >> ================      
            cmd = 'mpiexec -n ' + inputProc + ' StreamNet -fel ' + '"' + breach_filepath_tif + '"' + ' -p ' + '"' + p + '"' + \
                  ' -ad8 ' + '"' + ad8_no_wg + '"' + ' -src ' + '"' + ad8_wg + '"' + ' -ord ' + '"' + ord_g + '"' + ' -tree ' + \
                  '"' + tree + '"' + ' -coord ' + '"' + coord + '"' + ' -net ' + '"' + net + '"' + ' -w ' + '"' + w + \
                  '"'
            # Submit command to operating system
            logger.info('Running TauDEM Stream Reach and Watershed...')
            os.system(cmd)

            # Capture the contents of shell command and logger.info it to the arcgis dialog box
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            message = "\n"
            for line in process.stdout.readlines():
                if isinstance(line, bytes):  # true in Python 3
                    line = line.decode()
                message = message + line
            logger.info(message)

            # Let's get rid of some output that we are not currently using...
            try:
                #            os.remove(w)
                os.remove(coord)
                os.remove(tree)
                os.remove(ord_g)
            except:
                logger.info('Warning: Problem removing files!')
                pass

    #            # ============= << 5. Dinf with TauDEM >> =============        YES
    #            logger.info('Running TauDEM Dinfinity...')
    #            cmd = 'mpiexec -n ' + inputProc + ' DinfFlowDir -fel ' + '"' + breach_filepath_tif + '"' + ' -ang ' + '"' + ang + '"' + \
    #                  ' -slp ' + '"' + slp + '"'
    #
    #            # Submit command to operating system
    #            os.system(cmd)
    #
    #            # Capture the contents of shell command and logger.info it to the arcgis dialog box
    #            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    #
    #            # Get some feedback from the process to logger.info out...
    #            message = "\n"
    #            for line in process.stdout.readlines():
    #                line = line.decode()
    #        #            if isinstance(line, bytes):	   # true in Python 3
    #        #                line = line.decode()
    #                message = message + line
    #            logger.info(message)
    #
    #            # ============= << 6. DinfDistanceDown (HAND) with TauDEM >> ============= YES
    #            distmeth = 'v'
    #            statmeth = 'ave'
    #
    #            # Use original DEM here...
    #            logger.info('Running TauDEM Dinf Distance Down...') # Use Breached or Raw DEM here?? Currently using Raw
    #            cmd = 'mpiexec -n ' + inputProc + ' DinfDistDown -fel ' + '"' + str_dem_path_tif + '"' + ' -ang ' + '"' + ang + '"' + \
    #                  ' -src ' + '"' + ad8_wg + '"' + ' -dd ' + '"' + dd + '"' + ' -m ' + statmeth + ' ' + distmeth
    #
    #            # Submit command to operating system
    #            os.system(cmd)
    #
    #            # Get some feedback from the process to logger.info out...
    #            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    #
    #            message = "\n"
    #            for line in process.stdout.readlines():
    #                line = line.decode()
    #        #            if isinstance(line, bytes):	   # true in Python 3
    #        #                line = line.decode()
    #                message = message + line
    #
    #            logger.info(message)

    except:
        logger.info("Unexpected error:", sys.exc_info()[0])
        raise

    return


def interpolate(arr_in, ind_val):
    if ind_val == np.ceil(ind_val):
        out_val = arr_in[int(np.ceil(ind_val))]
    else:
        out_val = arr_in[np.int(ind_val)] + (ind_val - np.int(ind_val)) * (arr_in[int(np.ceil(ind_val))] - arr_in[
            np.int(ind_val)])  # it will always be divided by 1

    return out_val


# Count the number of features in a vector file...    
def get_feature_count(str_shp_path):
    with fiona.open(str(str_shp_path), 'r') as features:
        i_count = len(features)

    return i_count


# ===============================================================================
#  Reach characteristics from HAND
# =============================================================================== 
def reach_characteristics_hand(str_sheds_path, str_hand_path, str_slp_path, logger):
    '''
    Calculate the reach geometry metrics necessary for synthetic
    rating curve derivation from HAND grids.
    
    Returns: Dataframe.  Channel properties (m) and Q in cms
    '''
    # Define the vertical slice array:
    interval = 0.3048  # vertical step height NOTE:  Make sure you're considering units here
    rng = 25  # maximum height NOTE:  What should this be?? (25 is the NFIE setting)
    arr_slices = np.arange(interval, rng, interval)

    roughness = 0.05  # Manning's n. Hardcode for now (units of meters)

    #    arr_inds = np.arange(0, len(arr_slices)-1)

    lst_props = []

    # Open the HAND layer:
    with rasterio.open(str(str_hand_path)) as ds_hand:

        #        hand_mask = ds_hand.read_masks(1)

        #        nodata = ds_hand.nodata
        #        out_meta = ds_hand.meta.copy()
        #        res = ds_hand.res[0]  # CAREFUL:  If WGS, resolution is in degrees
        res = ds_hand.res[0]  # meters
        #        arr_fim = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])
        #        arr_fim[:,:] = out_meta['nodata']

        # Open the slope layer:
        with rasterio.open(str(str_slp_path)) as ds_slp:

            # Open the catchment polygon layer:
            with fiona.open(np.str(str_sheds_path), 'r') as sheds:

                # Loop over catchments:
                for shed in sheds:

                    # Get the linkno:
                    linkno = shed['properties']['FEATUREID']

                    if linkno != 4509202: continue

                    logger.info(linkno)

                    # Mask the HAND grid for each catchment polygon:
                    w_hand, out_transform = rasterio.mask.mask(ds_hand, [shed['geometry']], crop=True, filled=False,
                                                               invert=False)

                    # ALSO mask the slp grid:
                    w_slp, out_tranform = rasterio.mask.mask(ds_slp, [shed['geometry']], crop=True, filled=False,
                                                             invert=False)

                    # Get reach length here from the attribute table?
                    length = shed['properties']['LENGTHKM'] * 1000  # Make sure units are correct (should be meters)

                    # Get reach slope:
                    slope = shed['properties']['SLOPE']

                    # Vertical slice loop:
                    for i_step in arr_slices:

                        try:

                            i_hand_w = w_hand - i_step
                            #                            i_hand_w = i_hand_w[i_hand_w<0]
                            #                            i_hand_w = w_hand[(w_hand<=i_step) & (w_hand>=0)]
                            arr_inun = i_hand_w[i_hand_w < 0]
                            i_num_pixels = arr_inun.size

                            test = i_hand_w <= 0
                            cell_height = i_hand_w[test]
                            volume = 10 * 10 * cell_height.sum()

                            # Create a boolean array to use for selecting the slope pixels for this step:
                            #                            i_hand_bool = ((w_hand<=i_step)&(w_hand>=0))

                            # Slope grid values for this HAND interval:                            
                            i_slp_w = w_slp[i_hand_w < 0]

                            # Surface area of inundated zone:
                            area_surf = i_num_pixels * (res ** 2)

                            # Channel bed area of inundated zone:
                            area_bed = ((res ** 2) * (1 + i_slp_w ** 2) ** 0.5).sum()

                            # Volume of this slice:
                            volume = ((res ** 2) * arr_inun).sum() * -1.  # cubic meters correct?
                            # volume = area_surf*i_interval 

                            # Cross-sectional area:
                            area_xn = volume / length

                            # Wetted perimeter:
                            wetted_p = area_bed / length

                            # Top width:
                            top_width = area_surf / length

                            # Hydraulic radius:
                            hyd_radius = area_xn / wetted_p

                            # Calculate Q from Manning's equation:
                            Q = (1 / roughness) * area_xn * (hyd_radius ** 0.667) * (slope ** 0.5)

                            # Add to list:
                            lst_props.append((linkno, i_step, area_surf, area_bed, volume, area_xn, wetted_p, top_width,
                                              hyd_radius, Q))

                        #                        except Exception as e:
                        #                            logger.info('Exception!: {}'.format(str(e)))
                        #                            sys.exit()
                        except:
                            pass  # end of loop

    # Create final output dataframe:
    df_props = pd.DataFrame(lst_props, columns=['linkno', 'height', 'area_surf', 'area_bed', 'volume',
                                                'area_xn', 'wetted_p', 'top_width', 'hyd_radius', 'Q'])
    # Write it to csv:
    df_props.to_csv(r"D:\hand\nfie\020700\difficult_run\dr_hydroprops_test.csv")

    return


def hand_analysis_chsegs(str_hand_path, str_chanmet_segs, str_src_path, str_fp_path, str_dem_path, logger):
    """
    Calculation of floodplain metrics by analyzing HAND using 2D cross-sections and
    differentiating between channel pixels and floodplain pixels.

    :param str_hand_path: Path to the HAND grid .tif
    :param str_chanmet_segs: Path to the output of the channel_width_from_bank_pixels() func
    :param str_src_path: Path to the .tif file where stream reaches correspond to linkno values
    :param str_fp_path: Path to the floodplain grid .tif
    :param logger: Logger instance for messaging
    :return: Writes out additional attributes to the file in str_chanmet_segs

    """
    # Open the stream network segments layer with channel metrics:
    gdf_segs = gpd.read_file(str_chanmet_segs)

    lst_linkno = []
    # Channel metrics:
    lst_bnk_ht = []
    lst_chn_wid = []
    lst_chn_shp = []
    lst_geom = []  # for testing/creating a new output file
    # FP metrics:
    lst_fpmin = []
    lst_fpmax = []
    lst_fpstd = []
    lst_fprug = []
    lst_fpwid = []
    lst_fprange = []
    lst_fpmin_e = []
    lst_fpmax_e = []
    lst_fpstd_e = []
    lst_fprange_e = []

    # Open the hand layer...
    with rasterio.open(str(str_hand_path)) as ds_hand:

        out_meta = ds_hand.meta.copy()
        arr_chn = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])
        arr_chn[:, :] = out_meta['nodata']

        # Access the src layer for excluding channel pixels:
        with rasterio.open(str(str_src_path)) as ds_src:

            res = ds_hand.res[0]

            # Access the floodplain grid:
            with rasterio.open(str(str_fp_path)) as ds_fp:

                with rasterio.open(str(str_dem_path)) as ds_dem:

                    # Loop over each segment:
                    for tpl in gdf_segs.itertuples():

                        try:
                            # if tpl.linkno != 3343:
                            #     continue

                            logger.info(f'\t{tpl.Index}')

                            # Get Xn length based on stream order:
                            p_xnlength, p_fitlength = get_xn_length_by_order(tpl.order, False)

                            x, y = zip(*mapping(tpl.geometry)['coordinates'])

                            # Get the segment midpoint: NOTE-can also use this to identify the channel blob?
                            midpt_x = x[int(len(x) / 2)]  # Also-this isn't actually midpt by distance
                            midpt_y = y[int(len(y) / 2)]

                            # Build a 1D cross-section from the end points:
                            lst_xy = build_xns(y, x, midpt_x, midpt_y, p_xnlength)

                            try:
                                # Turn the cross-section into a linestring:
                                fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
                            except Exception as e:
                                logger.info(f'Error converting Xn endpts to LineString: {str(e)}')
                                pass

                            # Buffer the cross section to form a 2D rectangle:
                            # about half of the line segment straight line distance --> AFFECTS LABELLED ARRAY!
                            buff_len = tpl.dist_sl / 2.5
                            geom_fpls_buff = fp_ls.buffer(buff_len, cap_style=2)
                            xn_buff = mapping(geom_fpls_buff)

                            # Mask the hand grid for each feature...
                            w_hand, trans_hand = rasterio.mask.mask(ds_hand, [xn_buff], crop=True)
                            w_hand = w_hand[0]
                            w_hand[w_hand == ds_hand.nodata] = -9999.

                            # Set up vertical intervals to slice using 2D cross-section horizontal plane:
                            w_min = w_hand[w_hand > -9999.].min()
                            w_max = w_hand.max()
                            arr_slices = np.linspace(w_min, w_max, 50)

                            # Also mask the src layer (1 time) to get the indices of the raster stream line:
                            w_src, trans_src = rasterio.mask.mask(ds_src, [xn_buff], crop=True)
                            w_src = w_src[0]

                            # Also mask the floodplain grid:
                            w_fp, trans_fp = rasterio.mask.mask(ds_fp, [xn_buff], crop=True)
                            w_fp = w_fp[0]

                            # Aaaand, mask the dem:
                            w_dem, trans_dem = rasterio.mask.mask(ds_dem, [xn_buff], crop=True)
                            w_dem = w_dem[0]

                            # NOTE: You need to make sure you're looking at channel pixels here only and not other stuff
                            # Channel pixel indices:
                            src_inds = np.where(w_src == tpl.linkno)
                            # Convert to a set to keep only unique values:
                            src_inds = set(zip(src_inds[0], src_inds[1]))

                            # To produce labeled array:
                            s = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

                            # Consider the blobs of hand pixels in each slice:
                            lst_count = []
                            lst_width = []
                            # set_inds=set([])
                            lst_inds = []
                            lst_height = []

                            # Begin vertical slicing using 2D cross-section horizontal plane:
                            w_inds = np.indices(w_hand.shape)
                            for i, i_step in enumerate(arr_slices[1:]):  # skip the first entry

                                # Create a binary array where within step height threshold:
                                w_step = w_hand.copy()
                                w_step[(w_step < i_step) & (w_step > -9999.)] = 1
                                w_step[w_step != 1] = 0

                                # scipy labeled array:
                                labeled_arr, num_feats = label(w_step, structure=s)

                                # You need to loop over num_feats here and do the test:
                                for feat in np.arange(0, num_feats):  # does this produce the correct number for feat?

                                    # Get the window indices of each feature:
                                    inds = set(zip(w_inds[0][labeled_arr == feat + 1], w_inds[1][labeled_arr == feat + 1]))

                                    # if they share indices, consider this blob connected to the channel
                                    if len(src_inds.intersection(inds)) > 0:
                                        lst_count.append(len(inds))
                                        lst_width.append(len(inds) * (res ** 2) / tpl.dist_sl)
                                        lst_height.append(i_step)
                                        lst_inds.append(inds)

                            # End slices here
                            df_steps = pd.DataFrame(
                                {'count': lst_count, 'height': lst_height, 'width': lst_width, 'inds': lst_inds})

                            if len(df_steps.index) < 3:
                                logger.info('Too few slices!')
                                lst_bnk_ht.append(-9999.)
                                lst_chn_wid.append(-9999.)
                                lst_chn_shp.append(-9999.)
                                lst_linkno.append(tpl.linkno)
                                lst_geom.append(tpl.geometry)
                                # FP metrics:
                                lst_fpmax.append(-9999.)
                                lst_fpmin.append(-9999.)
                                lst_fpstd.append(-9999.)
                                lst_fprug.append(-9999.)
                                lst_fpwid.append(-9999.)
                                lst_fprange.append(-9999.)
                                lst_fpmin_e.append(-9999.)
                                lst_fpmax_e.append(-9999.)
                                lst_fpstd_e.append(-9999.)
                                lst_fprange_e.append(-9999.)
                                continue

                            df_steps['dy'] = df_steps.height.diff()
                            df_steps['dx'] = df_steps.width.diff()
                            df_steps['delta_width'] = df_steps.dx / df_steps.dy
                            indx = df_steps.delta_width.iloc[1:].idxmax() - 1

                            chn_wid = df_steps.width.iloc[indx]
                            bnk_ht = df_steps.height.iloc[indx]
                            chn_shp = np.arctan(bnk_ht / chn_wid)  # sort of like entrenchment?

                            # Separate the FP and channel pixels using bnk_ht
                            # Channel pixels only:
                            for i_set in df_steps.inds.iloc[0:indx+1].tolist():
                                src_inds.update(i_set)

                            # NEED A TUPLE OF 1D ARRAYS FOR ARRAY INDEXING
                            lst1, lst2 = zip(*list(src_inds))

                            # the hand grid window values at indices corresponding to channel pixels only?
                            # w_chn = w_hand[lst1, lst2]  # TODO: Calculate additional chan metrics from this?

                            # Get the FP pixels without the channel pixels:
                            mask = np.ones_like(w_hand, dtype=bool)
                            mask[lst1, lst2] = False

                            # Relative elevation (HAND):
                            try:
                                w_fp = w_fp[mask]
                                w_fp = w_fp[w_fp != ds_fp.nodata]  # also remove nodata vals

                                if w_fp.size == 0:
                                    logger.info('No FP!')
                                    # There's nothing we can do here related to FP:
                                    lst_fpmax.append(-9999.)
                                    lst_fpmin.append(-9999.)
                                    lst_fpstd.append(-9999.)
                                    lst_fprug.append(-9999.)
                                    lst_fpwid.append(-9999.)
                                    lst_fprange.append(-9999.)
                                    lst_fpmin_e.append(-9999.)
                                    lst_fpmax_e.append(-9999.)
                                    lst_fpstd_e.append(-9999.)
                                    lst_fprange_e.append(-9999.)
                                    # Channel metrics:
                                    lst_bnk_ht.append(bnk_ht)
                                    lst_chn_wid.append(chn_wid)
                                    lst_chn_shp.append(chn_shp)
                                    lst_linkno.append(tpl.linkno)
                                    lst_geom.append(tpl.geometry)
                                    continue

                                # FP width:
                                num_pixels = w_fp.size
                                area_pixels = num_pixels * (ds_fp.res[0] ** 2)  # get grid resolution
                                # Calculate width by stretching it along the length of the 2D Xn:
                                fp_width = area_pixels / (buff_len * 2)

                                # FP roughness (Planar area vs. actual area):
                                fp_rug = rugosity(w_fp, ds_fp.res[0], logger)  # returns -9999. if error

                                # Depth range:
                                fp_max = w_fp.max()
                                fp_min = w_fp.min()
                                fp_std = w_fp.std()
                                fp_range = fp_max - fp_min

                            except Exception as e:
                                logger.error(f'Error calculating relative elevation FP metrics: {str(e)}')
                                fp_max = -9999.
                                fp_min = -9999.
                                fp_std = -9999.
                                fp_range = -9999.
                                fp_rug = -9999.
                                fp_width = -9999.

                            try:
                                # Absolute elevation (DEM):
                                w_dem = w_dem[mask]
                                w_dem = w_dem[w_dem != ds_dem.nodata]  # also remove nodata vals
                                # Elevation range
                                fp_max_e = w_dem.max()
                                fp_min_e = w_dem.min()
                                fp_std_e = w_dem.std()
                                fp_range_e = fp_max_e - fp_min_e

                            except Exception as e:
                                logger.error(f'Error calculating absolute elevation FP metrics: {str(e)}')
                                fp_min_e = -9999.
                                fp_max_e = -9999.
                                fp_std_e = -9999.
                                fp_range_e = -9999.

                            # <editor-fold desc="TODO: Write out channel pixels to tif">
                            # # TODO: Write out the channel pixels to a tif?
                            # shp = np.shape(w_hand)
                            # # window bounds in x-y space (west, south, east, north)
                            # bounds = rasterio.transform.array_bounds(shp[0], shp[1], trans_hand)
                            # # Rasterio v1.x upper left row and column of window?
                            # col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[3])
                            #
                            # row_min = np.int(row_min)
                            # col_min = np.int(col_min)
                            # row_max = np.int(row_min + shp[0])
                            # col_max = np.int(col_min + shp[1])
                            #
                            # arr_w = np.empty([row_max - row_min, col_max - col_min], dtype=out_meta['dtype'])
                            # arr_w[:, :] = arr_chn[row_min:row_max, col_min:col_max]
                            #
                            # inds_lt = np.where(arr_chn[row_min:row_max, col_min:col_max] < w)
                            # arr_w[inds_lt] = w[inds_lt]
                            #
                            # # assign the FIM window for this catchment to the total array
                            # arr_chn[row_min:row_max, col_min:col_max] = arr_w
                            # </editor-fold>

                            # Save metrics to lists
                            # Relative elevation:
                            lst_fpmax.append(fp_max)
                            lst_fpmin.append(fp_min)
                            lst_fpstd.append(fp_std)
                            lst_fprug.append(fp_rug)
                            lst_fpwid.append(fp_width)
                            lst_fprange.append(fp_range)
                            # Absolute elevation:
                            lst_fpmin_e.append(fp_min_e)
                            lst_fpmax_e.append(fp_max_e)
                            lst_fpstd_e.append(fp_std_e)
                            lst_fprange_e.append(fp_range_e)
                            # Channel metrics:
                            lst_bnk_ht.append(bnk_ht)
                            lst_chn_wid.append(chn_wid)
                            lst_chn_shp.append(chn_shp)
                            lst_linkno.append(tpl.linkno)
                            lst_geom.append(tpl.geometry)

                            # logger.info('hey')
                        except Exception as e:
                            logger.info(f'Error with segment {tpl.Index}; skipping. {str(e)}')
                            # sys.exit()
                            # FP metrics:
                            lst_fpmax.append(-9999.)
                            lst_fpmin.append(-9999.)
                            lst_fpstd.append(-9999.)
                            lst_fprug.append(-9999.)
                            lst_fpwid.append(-9999.)
                            lst_fprange.append(-9999.)
                            lst_fpmin_e.append(-9999.)
                            lst_fpmax_e.append(-9999.)
                            lst_fpstd_e.append(-9999.)
                            lst_fprange_e.append(-9999.)
                            # Channel metrics:
                            lst_bnk_ht.append(-9999.)
                            lst_chn_wid.append(-9999.)
                            lst_chn_shp.append(-9999.)
                            lst_linkno.append(tpl.linkno)
                            lst_geom.append(tpl.geometry)
                            continue

        # Re-save the channel metrics shapefile with FP metrics added:
        gdf_segs['bnk_ht3'] = lst_bnk_ht
        gdf_segs['chn_shp3'] = lst_chn_shp
        gdf_segs['chn_wid3'] = lst_chn_wid
        gdf_segs['fp3_min'] = lst_fpmin
        gdf_segs['fp3_max'] = lst_fpmax
        gdf_segs['fp3_std'] = lst_fpstd
        gdf_segs['fp3_wid'] = lst_fpwid
        gdf_segs['fp3_rug'] = lst_fprug
        gdf_segs['fp3_rng'] = lst_fprange
        gdf_segs['fp3_min_e'] = lst_fpmin_e
        gdf_segs['fp3_max_e'] = lst_fpmax_e
        gdf_segs['fp3_std_e'] = lst_fpstd_e
        gdf_segs['fp3_rng_e'] = lst_fprange_e

        path_to_file, filename = os.path.split(str_chanmet_segs)
        gdf_segs.to_file(f'E:\\DRB_Files\\drb_chan_fp_metrics_2019.03.27\\{filename}')

        # TESTING:
        # gdf_test = gpd.GeoDataFrame()
        # gdf_test.crs = gdf_segs.crs
        # gdf_test['geometry'] = lst_geom
        # gdf_test['bnk_ht2'] = lst_bnk_ht
        # gdf_test['chn_shp'] = lst_chn_shp
        # gdf_test['chn_wid'] = lst_chn_wid
        # gdf_test['linkno'] = lst_linkno
        # gdf_test['fp_min'] = lst_fpmin
        # gdf_test['fp_max'] = lst_fpmax
        # gdf_test['fp_std'] = lst_fpstd
        # gdf_test['fp_wid'] = lst_fpwid
        # gdf_test['fp_rug'] = lst_fprug
        # gdf_test['fp_rng'] = lst_fprange
        # gdf_test.to_file(str_chanmet_segs[:-4] + '_TEST_HAND_ANALYSIS.shp')


def fp_metrics_chsegs(str_fim_path, str_chwid, str_chanmet_segs, logger):
    """
    Calculate floodplain metrics from 2D cross-sections
    :param str_fim_path: path to the flood inundation raster
    :param str_chwid: name of the field in the channel segment layer containing pre-calculated
                      channel width values
    :param str_chanmet_segs: path to the segmented streamline layer containing pre-calculated channel metrics
                             (output form channel width from bank pixel method)
    :param logger:
    :return: Resaves the str_chanmet_segs file with additonal attributes

    NOTE: stream order a required field in str_chanmet_segs ('order')
    """
    # Open the stream network segments layer with channel metrics:
    gdf_segs = gpd.read_file(str_chanmet_segs)

    lst_fpwid = []
    lst_fprng = []
    lst_geom = []
    lst_min = []
    lst_max = []
    lst_std = []
    lst_rug = []

    # Keep only valid geometries:
    gdf_segs = gdf_segs[gdf_segs.geometry.is_valid]

    # Open the floodplain layer...
    with rasterio.open(str(str_fim_path)) as ds_fim:
        # Loop over each segment:
        for tpl in gdf_segs.itertuples():

            try:
                # if tpl.Index==931: # FOR TESTING
                #      logger.info('pause')

                # Get Xn length based on stream order:
                p_xnlength, p_fitlength = get_xn_length_by_order(tpl.order, True)

                x, y = zip(*mapping(tpl.geometry)['coordinates'])

                # Get the segment midpoint:                                                    
                midpt_x = x[int(len(x) / 2)]
                midpt_y = y[int(len(y) / 2)]

                # Build a 1D cross-section from the end points:
                lst_xy = build_xns(y, x, midpt_x, midpt_y, p_xnlength)

                try:
                    # Turn the cross-section into a linestring:
                    fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
                except:
                    logger.info('Error converting Xn endpts to LineString')
                    pass

                # Buffer the cross section to form a 2D rectangle:        
                buff_len = tpl.dist_sl / 1.85  # about half of the line segment straight line distance
                geom_fpls_buff = fp_ls.buffer(buff_len, cap_style=2)
                xn_buff = mapping(geom_fpls_buff)

                # Mask the fp for each feature:
                w_fim, trans_fim = rasterio.mask.mask(ds_fim, [xn_buff], crop=True)
                w_fim = w_fim[0]
                w_fim = w_fim[w_fim != ds_fim.nodata]

                if w_fim[w_fim > 0].size == 0:
                    lst_fpwid.append(-9999.)
                    lst_fprng.append(-9999.)
                    lst_geom.append(-9999.)
                    lst_min.append(-9999.)
                    lst_max.append(-9999.)
                    lst_std.append(-9999.)
                    lst_rug.append(-9999.)
                    continue

                # OR, Get indices of FIM pixels and use those for the DEM
                # inds_fim=np.where(w_fim!=ds_fim.nodata)

                fp_fim = w_fim[w_fim > 0]  # Assumes is along zero height pixels
                fp_min = fp_fim.min()
                fp_max = fp_fim.max()
                fp_std = fp_fim.std()

                # << Related to mapping the floodplain based on HAND height >>
                # Count the number of pixels in the buffered Xn...
                num_pixels = w_fim.size

                # Calculate area of FP pixels...
                area_pixels = num_pixels * (ds_fim.res[0] ** 2)  # get grid resolution

                # Calculate width by stretching it along the length of the 2D Xn...
                fp_width = area_pixels / (buff_len * 2)
                #    fp_width=0 # For testing purposes

                # Elevation range using HAND heights:
                try:
                    fp_range = fp_max - fp_min
                except:
                    fp_range = 0
                    pass

                # Subtract channel width from fp width...
                fp_width = fp_width - getattr(tpl, str_chwid)
                # If negative, just set it to zero:
                if fp_width < 0.: fp_width = 0

                # Try calculating roughness (Planar area vs. actual area):
                fp_rug = rugosity(w_fim, ds_fim.res[0], logger)  # returns -9999. if error

                lst_min.append(fp_min)
                lst_max.append(fp_max)
                lst_std.append(fp_std)
                lst_fpwid.append(fp_width)
                lst_fprng.append(fp_range)
                lst_rug.append(fp_rug)
                lst_geom.append(tpl.geometry)
            #                logger.info('hey')
            except Exception as e:
                logger.info(f'Error with segment {tpl.Index}: {str(e)}')
                lst_fpwid.append(-9999.)
                lst_fprng.append(-9999.)
                lst_geom.append(-9999.)
                lst_min.append(-9999.)
                lst_max.append(-9999.)
                lst_std.append(-9999.)
                lst_rug.append(-9999.)
                #                sys.exit() # TEST TEST TEST
                continue

    # Re-save the channel metrics shapefile with FP metrics added:         
    #    gdf_out=gpd.GeoDataFrame()
    #    gdf_out.crs=gdf_segs.crs
    #    gdf_out.geometry=gdf_segs.geometry
    gdf_segs['fp_width_2d'] = lst_fpwid
    gdf_segs['fp_range_2d'] = lst_fprng
    gdf_segs['fp_min_2d'] = lst_min
    gdf_segs['fp_max_2d'] = lst_max
    gdf_segs['fp_std_2d'] = lst_std
    gdf_segs['fp_rug_2d'] = lst_rug
    gdf_segs.to_file(str_chanmet_segs)  # [:-4]+'_TEST.shp')


#            gdf=gpd.GeoDataFrame()
#            gdf['geometry']=lst_geom
#            gdf['fp_width']=lst_fpwid
#            gdf['fp_range']=lst_fprng            
#            gdf.crs=gdf_segs.crs
#            gdf['buff']=1
#            gdf.to_file(r"E:\test_buffs.shp")

# ===============================================================================
#  Delineate a FIM from the HAND grid using depth at each polygon (eg, catchment)
#
#    1. Read in catchment polygons (assume these have attributes necessary for regression?)
#    2. Calculate a HAND height (h) for each polygon based on some attribute(s)
#    3. Delineate FIM per catchment
#
# =============================================================================== 
def fim_hand_poly(str_hand_path, str_sheds_path, str_reachid, logger):
    # Open the HAND layer...
    with rasterio.open(str(str_hand_path)) as ds_hand:

        out_meta = ds_hand.meta.copy()
        arr_fim = np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])
        arr_fim[:, :] = out_meta['nodata']

        lst_h = []
        lst_linkno = []
        lst_prov = []

        # Open the catchment polygon layer...
        with fiona.open(np.str(str_sheds_path), 'r') as sheds:

            for shed in sheds:

                # Get the linkno...
                linkno = shed['properties']['gridcode']

                # Get the Province...
                prov = shed['properties']['PROVINCE']

                # Get the Drainage Area in km^2...
                da_km2 = shed['properties']['DSContArea'] / 1000000

                if (prov == 'COASTAL PLAIN' and da_km2 >= 3 and da_km2 <= 3000):
                    h = 1.65
                elif (prov == 'PIEDMONT' and da_km2 >= 3 and da_km2 <= 3000):
                    h = (np.log10(da_km2) * 0.471 + 0.523) ** 2
                elif (prov == 'VALLEY AND RIDGE' and da_km2 >= 3 and da_km2 <= 3000):
                    h = (np.log10(da_km2) * 0.471 + 0.375) ** 2
                elif (prov == 'APPALACHIAN PLATEAUS' and da_km2 >= 3 and da_km2 <= 3000):
                    h = (np.log10(da_km2) * 0.471 + 0.041) ** 2
                else:
                    lst_h.append(-9999)
                    lst_linkno.append(linkno)
                    lst_prov.append(prov)
                    continue  # skip this catchment

                lst_h.append(h)
                lst_linkno.append(linkno)
                lst_prov.append(prov)

                #                buff = mapping(shed)

                try:
                    # Mask the bankpts file for each feature...
                    w, out_transform = rasterio.mask.mask(ds_hand, [shed['geometry']], crop=True)

                    #                w[(w<=h) & (w>=0.)] # Get HAND depth below h
                    w[(w > h)] = out_meta['nodata']  # Assign NoData to everywhere else
                    #                w[(w<0.)] = out_meta['nodata'] # Assign NoData to everywhere else

                    # Now write out the FIM for this shed...
                    w = w[0]
                    shp = np.shape(w)

                    bounds = rasterio.transform.array_bounds(shp[0], shp[1],
                                                             out_transform)  # window bounds in x-y space (west, south, east, north)

                    col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[
                        3])  # RAsterio v1.x upper left row and column of window?
                    #                    col_min, row_min = ~ds_hand.affine * (bounds[0], bounds[3])

                    row_min = np.int(row_min)
                    col_min = np.int(col_min)
                    row_max = np.int(row_min + shp[0])
                    col_max = np.int(col_min + shp[1])

                    arr_w = np.empty([row_max - row_min, col_max - col_min], dtype=out_meta['dtype'])
                    arr_w[:, :] = arr_fim[row_min:row_max, col_min:col_max]
                    #
                    inds_lt = np.where(arr_fim[row_min:row_max, col_min:col_max] < w)
                    arr_w[inds_lt] = w[inds_lt]

                    arr_fim[row_min:row_max,
                    col_min:col_max] = arr_w  # assign the FIM window for this catchment to the total array
                except:
                    logger.info('WARNING:  Problem masking HAND grid using catchment Linkno:  ' + str(linkno))

                    # Write out the final FIM grid...
    #    logger.info('Writing final FIM .tif file...')
    str_fim_path = str_hand_path[:-4] + '_3sqkm_fim.tif'
    out_meta.update(compress='lzw')
    with rasterio.open(str_fim_path, 'w', tiled=True, blockxsize=512, blockysize=512, **out_meta) as dest:
        dest.write(arr_fim, indexes=1)

        # Write HAND heights to csv...
    str_csv_path = str_hand_path[:-4] + '_3sqkm_fim_h.csv'
    df_h = pd.DataFrame({str_reachid: lst_linkno, 'prov': lst_prov, 'h': lst_h})
    df_h.to_csv(str_csv_path)

    return


# ===============================================================================
#  Analyze DEM in vertical slices using an individual 2-D cross section
# ===============================================================================        
def get_indx(srs):
    first_diff = srs.diff()
    first_diff[first_diff < 0] = 0  # assign 0 to fall
    first_diff[first_diff > 0] = 1  # assign 1 to rise
    sh_first_diff = first_diff.shift(periods=-1)
    sub_diff = sh_first_diff - first_diff
    idx = sub_diff[sub_diff == 1].idxmin()
    return idx


def analyze_hand_2Dxn(w_hnd, tpl_xy, parm_ivert, dist_sl, res, logger):
    try:
        w_min = w_hnd[w_hnd > -9999.].min()
        w_max = w_hnd.max()
        # OR for scale-independence:
        arr_slices = np.linspace(w_min, w_max, 40)

        lst_count = []
        lst_width = []
        lst_inds = []
        lst_height = []

        # OR channel inds are where w_hnd==0:
        chan_inds = np.where(w_hnd == 0)
        chan_inds = set(zip(chan_inds[0], chan_inds[1]))

        # NOTE: You need to make sure you're looking at channel pixels here only and not other stuff
        from scipy.ndimage import label
        s = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

        w_inds = np.indices(w_hnd.shape)  # .T[:,:,[1, 0]]
        for i, i_step in enumerate(arr_slices[1:]):  # skip the first entry

            # Create a binary array where within step height threshold:
            w_step = w_hnd.copy()
            w_step[(w_step < i_step) & (w_step > -9999.)] = 1
            w_step[w_step != 1] = 0

            # Start with the first blob (assume there is only one and it is the channel)
            # Go to the next step; find indices of all blobs
            # If any of these indices are within (or adjacent to?) the previous blob, append them
            # Continue

            # scipy labeled array:
            labeled_arr, num_feats = label(w_step, structure=s)

            # You need to loop over num_feats here and do the test:          
            for feat in np.arange(1, num_feats):

                # Get the window indices of each feature:
                inds = set(zip(w_inds[0][labeled_arr == feat], w_inds[1][labeled_arr == feat]))

                if (len(chan_inds & inds) > 0):  # consider it the channel

                    #                    set_inds.update(inds)
                    lst_count.append(len(inds))
                    lst_width.append(len(inds) * (res ** 2) / dist_sl)
                    lst_height.append(i_step)
                    lst_inds.append(i)

        df_steps = pd.DataFrame({'count': lst_count, 'height': lst_height, 'width': lst_width})

        if len(df_steps.index) < 3:
            return -9999., -9999., -9999.

        df_steps['dy'] = df_steps.height.diff()
        df_steps['dx'] = df_steps.width.diff()
        df_steps['delta_width'] = df_steps.dx / df_steps.dy
        indx = df_steps.delta_width.iloc[1:].idxmax() - 1

        chn_wid = df_steps.width.iloc[indx]
        bnk_ht = df_steps.height.iloc[indx]
        bnk_ang = np.arctan(bnk_ht / chn_wid)  # bank angle

    #        if len(df_steps.index)<3:
    #            return -9999., -9999., -9999., -9999., -9999., -9999.

    #        # NOTE:  The scale at which you analyze this should depend on stream order? (by scale I mean vertical distance of slope/curvature measurement)
    #        # Let the vertical step vary based on min/max HAND height (and potentially width variance or something?) to automatically
    #        # account for variations in scale.  Use np.gradient (central difference) to the power of 2 (?) to consider curvature
    #        dx=arr_slices[2]-arr_slices[1]
    #        df_steps['slp']=np.gradient(df_steps.width, dx) # change in width per interval of depth
    #        df_steps['curv']=abs(np.gradient(df_steps.slp, dx))
    #
    #        df_steps['d1']=df_steps.width.diff()
    #        df_steps['d2']=df_steps.width.diff(periods=2)
    #        df_steps['d3']=df_steps.width.diff(periods=3)
    #        df_steps['d4']=df_steps.width.diff(periods=4)
    #        df_steps['d5']=df_steps.width.diff(periods=5)
    #
    #
    #        # Inflection points:
    #        idx1=get_indx(df_steps.d1)
    #        idx2=get_indx(df_steps.d2)
    #        idx3=get_indx(df_steps.d3)
    #        idx4=get_indx(df_steps.d4)
    #        idx5=get_indx(df_steps.d5)
    #        lst_inds=[idx1] #,idx3,idx4,idx5]
    #
    #        chn_wid=df_steps.width.iloc[lst_inds].mean()
    #        bnk_ht=df_steps.height.iloc[lst_inds].mean()

    #        first_diff = df_steps['d1'].diff()
    #        first_diff[first_diff<0]=0 # assign 0 to fall
    #        first_diff[first_diff>0]=1 # assign 1 to rise
    #        sh_first_diff = first_diff.shift(periods=-1)
    #        sub_diff = sh_first_diff-first_diff
    #        idx1=sub_diff[sub_diff==2].idxmin()

    #        idx_bnk=df_steps[1:].curv.idxmax() # point of maximum curvature, skipping the first element
    #        bnk_ht=df_steps.height.iloc[idx_bnk] # bank height
    #        chn_wid = df_steps.width.iloc[idx_bnk] # width at bank height

    #        # Plotting for tests:
    #        ax=df_steps.plot(x='width',y='height', marker='.', logx=False, logy=False)
    #        ax.set(xlabel='Width', ylabel='Height')
    #        sys.exit()
    # IDEAS:
    # 1. Get the coefficients of the 2D polynomial that fits the curve?:
    #    a,b,c=np.polyfit(df_steps.height, df_steps.width, 2)    
    # 2. Split the curve in two parts equally between w_hnd.min and max, then take several ratios:
    # Find the break in slope of the curve (by every fourth point?):
    #        w1, w2, w3 = jenkspy.jenks_breaks(df_steps.width, nb_class=2)
    #        h1=df_steps.height[df_steps.width==w1].values[0]
    #        h2=df_steps.height[df_steps.width==w2].values[0]
    #        h3=df_steps.height[df_steps.width==w3].values[0]

    #        # Channel vs. FP slope:
    #        slope_chan=chn_wid/bnk_ht
    #        slope_fp=(df_steps.width.iloc[-1]-chn_wid)/(w_max-bnk_ht)
    #        slp_ratio=slope_chan/slope_fp
    #
    #        # rugosity upper/lower:  (actual area/planar area)
    #        arr2=w_hnd[w_hnd<=bnk_ht]
    #        arr3=w_hnd[w_hnd<=w_max]
    #
    #        rugosity2=rugosity(arr2, res)
    #        rugosity3=rugosity(arr3, res)
    #
    #        if rugosity2 and rugosity3:
    #            rug_ratio=rugosity2/rugosity3
    #        else:
    #            rug_ratio=-9999.
    #            rugosity2=-9999.
    #            rugosity3=-9999.

    except Exception as e:
        logger.info(f'Error in analyze_2D_xns:  {str(e)}')
        sys.exit()
    #        slope_chan=-9999.
    #        slope_fp=-9999.
    #        slp_ratio=-9999.
    #        rugosity2=-9999.
    #        rugosity3=-9999.
    #        rug_ratio=-9999.

    #    return slope_chan, slope_fp, slp_ratio, rugosity2, rugosity3, rug_ratio, bnk_ht, chn_wid, bnk_ang
    return bnk_ht, chn_wid, bnk_ang


# ===============================================================================
#  Calculates channel  sinuosity using parallel offset buffering
# ===============================================================================    
def chan_wid_bnk_pixel_worker(i_step, max_buff, lst_buff, ds_bankpixels, args):
    i_linkno, df_linkno = args

    i_linkno = int(i_linkno)

    # << Analysis by reach segments >>
    # Set up index array to split up df_linkno into segments (these dictate the reach segment length)...
    # NOTE:  Reach might not be long enough to break up
    arr_ind = np.arange(i_step, len(df_linkno.index) + 1,
                        i_step)  # NOTE: Change the step for resolution?
    lst_dfsegs = np.split(df_linkno, arr_ind)

    for i_seg, df_seg in enumerate(lst_dfsegs):  # looping over each reach segment

        order = df_seg.order.max()

        try:
            order = int(order)
        except:
            order = 1

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
            dist_sl = np.sqrt((arr_x[0] - arr_x[-1]) ** 2 + (arr_y[0] - arr_y[-1]) ** 2)
        except:
            print('Error calculated straight line distance')
            dist_sl = -9999.

        dist = ls.length
        sinuosity = dist / dist_sl  # ratio of sinuous length to straight line length

        lst_tally = []

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
                out_left = np.array([0])

            try:
                out_rt, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_rt)], crop=True)
            except:
                print('Right offset error')
                out_rt = np.array([0])

            num_pixels_left = len(out_left[out_left > 0.])
            num_pixels_rt = len(out_rt[out_rt > 0.])

            # You want the number of pixels gained by each interval...                    
            tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
            lst_tally.append(tpl_out)
            df_tally = pd.DataFrame(lst_tally, columns=['linkno', 'buffer', 'interval_left', 'interval_rt'])

        # Calculate weighted average                     
        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
        weighted_avg_left = 0
        weighted_avg_rt = 0
        n_top = 2

        try:
            for tpl in df_tally.nlargest(n_top, 'interval_left').iloc[0:2].itertuples():
                weighted_avg_left += tpl.buffer * (np.float(tpl.interval_left) / np.float(
                    df_tally.nlargest(n_top, 'interval_left').iloc[0:2].sum().interval_left))
        except Exception as e:
            weighted_avg_left = max_buff
            print('Left width set to max. Exception: {} \n'.format(e))

        try:
            for tpl in df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].itertuples():
                weighted_avg_rt += tpl.buffer * (np.float(tpl.interval_rt) / np.float(
                    df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].sum().interval_rt))
        except Exception as e:
            weighted_avg_rt = max_buff
            print('Right width set to max. Exception: {} \n'.format(e))

    return weighted_avg_left + weighted_avg_rt, weighted_avg_left, weighted_avg_rt, sinuosity


def channel_width_from_bank_pixels(df_coords, str_streamlines_path, str_bankpixels_path, str_reachid, i_step, max_buff,
                                   str_chanmet_segs, logger):
    logger.info('Channel width from bank pixels -- segmented reaches...')

    j = 0
    #    progBar = self.progressBar
    #    progBar.setVisible(True)

    gp_coords = df_coords.groupby('linkno')

    # Schema for the output properties file:
    schema_output = {'geometry': 'LineString',
                     'properties': {'linkno': 'int', 'ch_wid_total': 'float', 'ch_wid_1': 'float', 'ch_wid_2': 'float',
                                    'dist_sl': 'float', 'dist': 'float', 'sinuosity': 'float', 'order': 'int'}}

    # Access the bank pixel layer...
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:  # open with share=False for multithreading

        # Successive buffer-mask operations to count bank pixels at certain intervals
        lst_buff = range(int(ds_bankpixels.res[0]), max_buff, int(ds_bankpixels.res[0]))

        #        # TEST Multiprocessing
        #        func = partial(chan_wid_bnk_pixel_worker, i_step, max_buff, lst_buff, ds_bankpixels) # Can't send an open dataset reader?
        #        pool = mp.Pool(processes=2)
        #        lst_out = pool.map(func, gp_coords)
        #        pool.close()
        #        pool.join()
        #        for item in lst_out:
        #            print('pause')
        #        # END TEST MP

        # Access the streamlines layer...
        with fiona.open(str(str_streamlines_path), 'r') as streamlines:  # WHY IS THIS BEING OPENED??

            #            progBar.setRange(0, len(streamlines))

            # Get the crs...
            streamlines_crs = streamlines.crs

            # Open another file to write the output props:
            with fiona.open(str_chanmet_segs, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:

                for i_linkno, df_linkno in gp_coords:

                    #                    progBar.setValue(j)
                    j += 1

                    i_linkno = int(i_linkno)

                    #                    if i_linkno != 1368: continue

                    logger.info('linkno:  {}'.format(i_linkno))

                    # << Analysis by reach segments >>
                    # Set up index array to split up df_linkno into segments (these dictate the reach segment length)...
                    # NOTE:  Reach might not be long enough to break up
                    arr_ind = np.arange(i_step, len(df_linkno.index) + 1,
                                        i_step)  # NOTE: Change the step for resolution?
                    lst_dfsegs = np.split(df_linkno, arr_ind)

                    for i_seg, df_seg in enumerate(lst_dfsegs):  # looping over each reach segment

                        order = df_seg.order.max()

                        try:
                            order = int(order)
                        except:
                            order = 1

                        arr_x = df_seg.x.values
                        arr_y = df_seg.y.values

                        try:
                            # Create a line segment from endpts in df_seg...
                            ls = LineString(zip(arr_x, arr_y))
                        except:
                            logger.error('Cannot create a LineString using these points, skipping')
                            continue

                        try:
                            # Calculate straight line distance...
                            dist_sl = np.sqrt((arr_x[0] - arr_x[-1]) ** 2 + (arr_y[0] - arr_y[-1]) ** 2)
                        except:
                            logger.warning('Error calculated straight line distance')
                            dist_sl = -9999.

                        dist = ls.length
                        sinuosity = dist / dist_sl  # ratio of sinuous length to straight line length

                        lst_tally = []

                        for buff_dist in lst_buff:

                            try:
                                # Watch out for potential geometry errors here...
                                ls_offset_left = ls.parallel_offset(buff_dist, 'left')
                                ls_offset_rt = ls.parallel_offset(buff_dist, 'right')
                            except:
                                logger.warning('Error performing offset buffer')

                            # Buffer errors can result from complicated line geometry... 
                            try:
                                out_left, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_left)],
                                                                             crop=True)
                            except:
                                logger.warning('Left offset error')
                                out_left = np.array([0])

                            try:
                                out_rt, out_transform = rasterio.mask.mask(ds_bankpixels, [mapping(ls_offset_rt)],
                                                                           crop=True)
                            except:
                                logger.warning('Right offset error')
                                out_rt = np.array([0])

                            num_pixels_left = len(out_left[out_left > 0.])
                            num_pixels_rt = len(out_rt[out_rt > 0.])

                            # You want the number of pixels gained by each interval...                    
                            tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
                            lst_tally.append(tpl_out)
                            df_tally = pd.DataFrame(lst_tally,
                                                    columns=['linkno', 'buffer', 'interval_left', 'interval_rt'])

                        # Calculate weighted average                     
                        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
                        weighted_avg_left = 0
                        weighted_avg_rt = 0
                        n_top = 2

                        try:
                            for tpl in df_tally.nlargest(n_top, 'interval_left').iloc[0:2].itertuples():
                                weighted_avg_left += tpl.buffer * (np.float(tpl.interval_left) / np.float(
                                    df_tally.nlargest(n_top, 'interval_left').iloc[0:2].sum().interval_left))
                        except Exception as e:
                            weighted_avg_left = max_buff
                            logger.warning('Left width set to max. Exception: {} \n'.format(e))

                        try:
                            for tpl in df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].itertuples():
                                weighted_avg_rt += tpl.buffer * (np.float(tpl.interval_rt) / np.float(
                                    df_tally.nlargest(n_top, 'interval_rt').iloc[0:2].sum().interval_rt))
                        except Exception as e:
                            weighted_avg_rt = max_buff
                            logger.warning('Right width set to max. Exception: {} \n'.format(e))

                        # Write to the output shapefile here...
                        output.write({'properties': {'linkno': i_linkno,
                                                     'ch_wid_total': weighted_avg_left + weighted_avg_rt,
                                                     'ch_wid_1': weighted_avg_left, 'ch_wid_2': weighted_avg_rt,
                                                     'dist_sl': dist_sl, 'dist': dist, 'sinuosity': sinuosity,
                                                     'order': int(order)}, 'geometry': mapping(ls)})

                    #                    if j > 50: break
    return


# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size, use_wavelet_method,
                                     logger):
    logger.info('Bank pixels from curvature windows...')

    # Convert df_coords x-y to row-col via DEM affine
    # Loop over center row-col pairs accessing the window
    # Now loop over the linknos to get access grid by window...

    # << PARAMETERS >>
    cell_size = int(cell_size)

    # 3 m...
    w_height = 20  # number of rows
    w_width = 20  # number of columns
    buff = 3  # number of cells
    curve_thresh = 0.30  # good for 3m DEM?

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

    j = 0
    #    progBar = self.progressBar
    #    progBar.setVisible(True)

    try:

        sigma = 1.0  # select the scale sigma=1 --> NOTE:  FIND THIS ON THE FLY? (within each window?)
        g2x1, g2y1 = gauss_kern(sigma)

        with rasterio.open(str_dem_path) as ds_dem:

            # Transform to pixel space
            df_coords['col'], df_coords['row'] = ~ds_dem.transform * (df_coords['x'], df_coords['y'])
            df_coords[['row', 'col']] = df_coords[['row', 'col']].astype(np.int32)
            df_coords.drop_duplicates(['col', 'row'], inplace=True)  # rounding to integer
            #            total_len = len(df_coords.index)

            out_meta = ds_dem.meta.copy()
            out_meta['dtype'] = rasterio.uint8  # no need for float32 for bankpixels to save size of output
            out_meta['compress'] = 'lzw'

            arr_bankpts = np.zeros([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])

            #        progBar.setRange(0, len(df_coords.index))

            for tpl_row in df_coords.itertuples():

                #                if tpl_row.order != 6:
                #                    continue
                #
                if tpl_row.order == 5:
                    w_height = 40  # number of rows
                    w_width = 40  # number of columns
                if tpl_row.order >= 6:
                    w_height = 80
                    w_width = 80

                #                if tpl_row.linkno != 1368: continue

                #            progBar.setValue(j)
                j += 1

                #                logger.info('{} | {} -- {}'.format(tpl_row.linkno, j, total_len))

                row_min = np.int(tpl_row.row - np.int(w_height / 2))
                row_max = np.int(tpl_row.row + np.int(w_height / 2))
                col_min = np.int(tpl_row.col - np.int(w_width / 2))
                col_max = np.int(tpl_row.col + np.int(w_width / 2))

                # Now get the DEM specified by this window as a numpy array...
                w = ds_dem.read(1, window=((row_min, row_max), (col_min, col_max)))

                # Then extract the internal part of the window that contains the rotated window??
                w[w > 9999999.0] = 0.0  # NoData values may have been corrupted by preprocessing?
                w[w < -9999999.0] = 0.0

                if np.size(w) > 9:  # make sure a window of appropriate size was returned from the DEM

                    if use_wavelet_method:
                        # === Wavelet Curvature from Chandana ===
                        gradfx1 = signal.convolve2d(w, g2x1, boundary='symm', mode='same')
                        gradfy1 = signal.convolve2d(w, g2y1, boundary='symm', mode='same')

                        w_curve = gradfx1 + gradfy1

                        # Pick out bankpts...                    
                        w_curve[w_curve < np.max(w_curve) * curve_thresh] = 0.

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
                            w_curve = (Zx ** 2 + 1) * Zyy - 2 * Zx * Zy * Zxy + (Zy ** 2 + 1) * Zxx
                            w_curve = -w_curve / (2 * (Zx ** 2 + Zy ** 2 + 1) ** (1.5))
                        except:
                            logger.info('Error calculating Curvature in window...skipping')
                            continue

                        w_curve[w_curve < np.max(w_curve) * curve_thresh] = 0.
                        # =======================================

                    w_curve[w_curve < -99999999.] = 0.
                    w_curve[w_curve > 99999999.] = 0.

                    w_curve[w_curve > 0.] = 1.

                    # Note:  This assumes that the w_curve window is the specified size, which is not always the case for edge reaches...
                    arr_bankpts[row_min + buff:row_max - buff, col_min + buff:col_max - buff] = w_curve[
                                                                                                buff:w_height - buff,
                                                                                                buff:w_width - buff]

                    out_meta['nodata'] = 0.

            logger.info('Writing bank pixels .tif...')
            with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
                dest.write(arr_bankpts.astype(rasterio.uint8), indexes=1)
    #                dest.write(arr_bankpts, indexes=1)

    except Exception as e:
        logger.info('\r\nError in bankpixels_from_curvature_window. Exception: {} \n'.format(e))

    #        df_dist = pd.DataFrame(lst_dist)

    return


# ===============================================================================
#  Calculate angle from vertical of left and right banks
# ===============================================================================
def find_bank_angles(tpl_bfpts, lst_total_slices, xn_len, xn_elev_n, parm_ivert, xn_ptdistance, logger):
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
                x1 = lf_bottom_ind - 1
                x2 = lf_bottom_ind
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                lf_bottombank_ind = interp_bank(x1, x2, y1, y2, yuk)

                if abs(lf_bottombank_ind - tpl_bfpts[1]) > 0:
                    lf_angle = atan((abs(lf_bottombank_ind - tpl_bfpts[1])) / (
                                (tpl_bfpts[3] - parm_ivert) * xn_ptdistance)) * 57.29578  # convert radians to degrees
                else:
                    lf_angle = 0

            # RIGHT BANK: Interpolate to find left position along bank...
            rt_bottom_ind = lst_total_slices[1][-1]

            # Make sure we're within bounds here
            if rt_bottom_ind == 0 or rt_bottom_ind == xn_len:
                rt_angle = 0
            else:
                x1 = rt_bottom_ind
                x2 = rt_bottom_ind + 1
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                rt_bottombank_ind = interp_bank(x1, x2, y1, y2, yuk)

                if abs(rt_bottombank_ind - tpl_bfpts[2]) > 0:
                    rt_angle = atan((abs(rt_bottombank_ind - tpl_bfpts[2])) / (
                                (tpl_bfpts[3] - parm_ivert) * xn_ptdistance)) * 57.29578  # convert radians to degrees
                else:
                    rt_angle = 0

        else:  # if there's only one slice, just set it to 0? or -9999?
            #                            lf_angle = -9999
            #                            rt_angle = -9999
            # Use bottom slice for bank angle estimate...
            lf_bottom_ind = lst_total_slices[0][0]
            rt_bottom_ind = lst_total_slices[0][-1]

            if abs(lf_bottom_ind - tpl_bfpts[1]) > 0:
                lf_angle = atan((abs(lf_bottom_ind - tpl_bfpts[1]) * xn_ptdistance) / tpl_bfpts[
                    3]) * 57.29578  # convert radians to degrees
            else:
                lf_angle = 0

            if abs(rt_bottom_ind - tpl_bfpts[2]) > 0:
                rt_angle = atan((abs(rt_bottom_ind - tpl_bfpts[2]) * xn_ptdistance) / tpl_bfpts[
                    3]) * 57.29578  # convert radians to degrees
            else:
                rt_angle = 0

        # NOTE: For now, just set any resulting negative values to -9999, until we figure out what's going on (27Mar2015, SJL)
        if lf_angle < 0:
            lf_angle = -9999.0

        if rt_angle < 0:
            rt_angle = -9999.0

        tpl_angles = (lf_angle, rt_angle)

    except Exception as e:
        logger.info('\r\nError in find_bank_angles. Exception: {} \n'.format(e))

    return tpl_angles


# ===============================================================================
# Search Xn outward to the right to find the first point greater than the left bank
# ===============================================================================
def search_right_gt(xnelev, prev_ind, lf_bf):
    # Search outward from end of previous slice...
    for i in range(prev_ind + 1, len(xnelev), 1):

        if xnelev[i] > lf_bf:
            bank_ind = i
            break
        else:  # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind


# ===============================================================================
# Search Xn outward to the left to find the first point greater than the right bank
# ===============================================================================
def search_left_gt(xnelev, prev_ind, rt_bf):
    # Search outward from end of previous slice...
    for i in range(prev_ind, 0, -1):

        if xnelev[i] > rt_bf:
            bank_ind = i
            break
        else:  # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind


# ===============================================================================
# Interpolate to find positions of right/left bank
# ===============================================================================
def interp_bank(x1, x2, y1, y2, y_uk):
    x_uk = (((x2 - x1) * (y_uk - y1)) / (y2 - y1)) + x1;

    return x_uk


# ===============================================================================
# Search for banks via slope break and vertical slices
# ===============================================================================
def find_bank_ratio_method(lst_total, ratio_threshold, xnelev_zero, slp_thresh, logger):
    """
        Compares the length of the last gtzero slice (num of indices) vs. the previous slice

        Inputs:     lst_total   - a list of 1D array slice index values
                    ratio_threshold
                    xnelev_zero - the Xn elevation profile normalized to zero

        Output:     tpl_bfpts   - a tuple of bankfull points (left_index, right_index, height)
    """
    tpl_bfpts = ()  # output tuple
    num_slices = len(lst_total) - 1  # total number of slices, each at a height of param_vertstep
    xn_len = len(xnelev_zero) - 1  # length of Xn

    try:
        if num_slices > 2 and len(lst_total[num_slices - 1]) > 2:
            top_area = len(lst_total[num_slices])  # should really include point distance but this will cancel out?
            below_area = len(lst_total[num_slices - 1])

            # Check the ratio...
            this_ratio = float(top_area) / float(below_area)

            if (this_ratio > ratio_threshold):  # USE THIS TO DRIVE THE BANK BREAK DETERMINATION INSTEAD

                # Find end indices of this and of previous slice...
                prev_lf_ind = lst_total[num_slices - 1][0]
                prev_rt_ind = lst_total[num_slices - 1][-1]

                this_lf_ind = lst_total[num_slices][0]
                this_rt_ind = lst_total[num_slices][-1]

                # Bottom left and right for searching...
                bottom_lf = lst_total[0][0]
                bottom_rt = lst_total[0][-1]

                # First derivative is slope...
                lf_arr = np.array(xnelev_zero[this_lf_ind:prev_lf_ind + 1])
                rt_arr = np.array(xnelev_zero[prev_rt_ind - 1:this_rt_ind])
                firstdiff_left = np.diff(lf_arr)
                firstdiff_right = np.diff(rt_arr)

                # Set both indices to negative 1 initially...
                rt_bank_ind = -1
                lf_bank_ind = -1

                # Look for the first occurrence of a very small slope value in both directions...
                # slp_thresh = 0.03 # ? a parameter
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
                    if rt_bank_ind > 0 and lf_bank_ind < 0:  # only the right index exists

                        # Interpolate to find left bankfull...
                        bf_height = xnelev_zero[rt_bank_ind]
                        lf_x2 = search_left_gt(xnelev_zero, bottom_rt, bf_height)

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if (lf_y1 == lf_y2):
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(lf_x1, lf_x2, lf_y1, lf_y2, bf_height)

                            tpl_bfpts = (lfbf_ind, rt_bank_ind, bf_height)

                    elif lf_bank_ind > 0 and rt_bank_ind < 0:  # only the left index exists

                        # Interpolate to find right bank index...
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(xnelev_zero, bottom_lf, bf_height)

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if (rt_y1 == rt_y2):
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(rt_x1, rt_x2, rt_y1, rt_y2, bf_height)

                            tpl_bfpts = (lf_bank_ind, rtbf_ind, bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[rt_bank_ind] < xnelev_zero[
                        lf_bank_ind]:  # right is smaller than left

                        # Interpolate to find left bankfull...
                        bf_height = xnelev_zero[rt_bank_ind]
                        lf_x2 = search_left_gt(xnelev_zero, bottom_rt, bf_height)  # search all the way across?

                        # lf_x2 = search_left_lt(xnelev_zero, lf_bank_ind, bf_height) # find the index that's just smaller than bank height on left side FASTER?

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if (lf_y1 == lf_y2):
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(lf_x1, lf_x2, lf_y1, lf_y2, bf_height)

                            tpl_bfpts = (lfbf_ind, rt_bank_ind, bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[lf_bank_ind] < xnelev_zero[
                        rt_bank_ind]:  # left is smaller than right

                        # Interpolate to find right bank index...
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(xnelev_zero, bottom_lf,
                                                bf_height)  # Searches all the way across channel

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if (rt_y1 == rt_y2):
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(rt_x1, rt_x2, rt_y1, rt_y2, bf_height)

                            tpl_bfpts = (lf_bank_ind, rtbf_ind, bf_height)

                    elif rt_bank_ind > 0 and lf_bank_ind > 0 and xnelev_zero[lf_bank_ind] == xnelev_zero[
                        rt_bank_ind]:  # they're exactly equal
                        # logger.info 'they are the same!'
                        bf_height = xnelev_zero[lf_bank_ind]
                        tpl_bfpts = (lf_bank_ind, rt_bank_ind, bf_height)

    except Exception as e:
        logger.info('\r\nError in find_bank_ratio_method. Exception: {} \n'.format(e))

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
def analyze_xnelev(df_xn_elev, param_ivert, xn_ptdist, param_ratiothreshold, param_slpthreshold, nodata_val, logger):
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
    # USE NUMPY ARRAYS VS. LISTS??
    lst_bfmetrics = []  # list of tuples to contain output

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

            #            logger.info('this_linkno: {} | index: {}'.format(this_linkno, tpl_row.Index))

            #            if tpl_row.Index == 2254:
            #                logger.info('pause')

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
                gtzero_indices = np.nonzero((this_slice - thisxn_norm) > 0)[
                    0]  # Zero index get the first element of the returned tuple

                # NOTE:  DO THIS IN TERMS OF MAP COORDINATES HERE??

                # Use funtion to check if contiguous...
                if np.size(gtzero_indices) == 0:  # the first loop only

                    # get the index of the zero value...
                    lst_total_cnt.append(np.where(thisxn_norm == 0)[0])
                    prev_val = lst_total_cnt[0][0]

                elif is_contiguous(gtzero_indices):

                    # Yes, it is contiguous
                    # Save it to the total count...
                    lst_total_cnt.append(gtzero_indices)
                    prev_val = gtzero_indices[0]  # just need any value from the contiguous array

                else:
                    # No, it's not contiguous
                    # Find the contiguous part of the slice...
                    tpl_parts = np.array_split(gtzero_indices, np.where(np.diff(gtzero_indices) != 1)[
                        0] + 1)  # splits the contiguous elements into separate tuple elements

                    # Find the one that contains an element of the previous slice?...
                    #                if prev_val in [this_arr for this_arr in tpl_parts]: # use a list comprehension here?
                    for this_arr in tpl_parts:
                        if prev_val in this_arr[:]:
                            lst_total_cnt.append(this_arr)
                            prev_val = this_arr[0]
                            break

                # NOTE:  DO THIS IN TERMS OF MAP COORDINATES HERE??
                tpl_bankfullpts = find_bank_ratio_method(lst_total_cnt, param_ratiothreshold, thisxn_norm,
                                                         param_slpthreshold)

                if tpl_bankfullpts:  # test to see if there's anything in there?? THIS SEEMS LIKE KIND OF A HACK?

                    if tpl_bankfullpts[0] - tpl_bankfullpts[
                        1] == 0:  # BF points are so close that rounding to int gives identical values
                        # logger.info('pause')
                        break

                        # Add Xn number to the output...
                    # xn_elev_norm = tpl_thisxn[2] - np.min(tpl_thisxn[2]) # normalized elevation profile
                    xn_length = len(thisxn_norm)

                    # tpl_bankfullpts = (xn_cnt,) + tpl_bankfullpts + (lst_total_cnt,) + (this_linkno,) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),)

                    # Bank points tuple...
                    #                    tpl_bankfullpts = (tpl_thisxn[3],) + tpl_bankfullpts # add local ID
                    tpl_bankfullpts = (tpl_row.Index,) + tpl_bankfullpts

                    # Find bank angles...
                    tpl_bankangles = find_bank_angles(tpl_bankfullpts, lst_total_cnt, xn_length, thisxn_norm,
                                                      param_ivert, xn_ptdist)

                    # Estimate bankfull area...
                    # (Bank height - xn_elev_norm[i])*xn_ptdist
                    # Round up on left, round down on right
                    bf_area = 0
                    lst_bf_rng = range(int(ceil(tpl_bankfullpts[1])), int(tpl_bankfullpts[2]) + 1, 1)

                    for i in lst_bf_rng:
                        bf_area += (tpl_bankfullpts[3] - thisxn_norm[i]) * xn_ptdist

                        # Channel width...
                    ch_width = (tpl_bankfullpts[2] - tpl_bankfullpts[1]) * xn_ptdist

                    # Overbank ratio...
                    #                    if (tpl_bankfullpts[2]-tpl_bankfullpts[1]) <= 0:
                    #                        overbank_ratio = -9999.0
                    #                    else:
                    #                        overbank_ratio = len(lst_total_cnt[-1])/(tpl_bankfullpts[2]-tpl_bankfullpts[1])

                    try:
                        overbank_ratio = len(lst_total_cnt[-1]) / (tpl_bankfullpts[2] - tpl_bankfullpts[1])
                    except:
                        overbank_ratio = -9999.0

                    if bf_area == 0:
                        #                        logger.info('pause')
                        bf_area = -9999.0
                        total_arearatio = -9999.0
                    else:
                        # Also try area under entire Xn length relative to BF area...
                        total_xn_area = sum(thisxn_norm * xn_ptdist)
                        try:
                            total_arearatio = (total_xn_area - bf_area) / bf_area
                        except:
                            total_arearatio = -9999.0

                    tpl_metrics = tpl_bankfullpts + (lst_total_cnt,) + (this_linkno,) + (
                    tpl_bankfullpts[3] + np.min(tpl_row.elev),) + tpl_bankangles + (bf_area,) + (ch_width,) + (
                                  overbank_ratio,) + (total_arearatio,)

                    #                    # Output metrics tuple...       local ID      bank angles        bank height                        bank elevation                    xn area
                    #                    tpl_metrics = (this_linkno,) +  (xn_cnt,) + tpl_bankangles + (tpl_bankfullpts[3],) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),) + (bf_area,)

                    #                    # Output bank points tuple...
                    #                    tpl_bankpts = tpl_bankfullpts + (this_linkno,) + (tpl_bankfullpts[2] + np.min(tpl_thisxn[2]),)

                    lst_bfmetrics.append(tpl_metrics)

                    break  # no need to keep slicing here, unless we want to try for FP analysis
    #            else:
    #                logger.info 'no bank!'

    #        progDialog.close()

    except Exception as e:
        logger.info('\r\nError in analyze_xn_elev. Exception: {} \n'.format(e))
        pass

    return lst_bfmetrics


# ====================================================================================
#  Calculate channel metrics based on the bankpoint slope-threshold method at each Xn,
#   writing the bank points to a shapefile 
# ====================================================================================
def chanmetrics_bankpts(df_xn_elev, str_xnsPath, str_demPath, str_bankptsPath, parm_ivert, XnPtDist, parm_ratiothresh,
                        parm_slpthresh, logger):
    logger.info('Channel metrics from bank points...')

    # << BEGIN LOOP >>
    # Do the rest by looping in strides, rather than all at once, to conserve memory...(possibly using multiprocessing)
    xn_count = get_feature_count(str_xnsPath)

    # Striding...
    arr_strides = np.linspace(0, xn_count, xn_count / 100)
    arr_strides = np.delete(arr_strides, 0)

    # NOTE:  Use MultiProcessing here over arr_strides??

    #    progBar = self.progressBar
    #    progBar.setVisible(True)
    #    progBar.setRange(0, arr_strides.size)

    # Now loop over the linknos to get access grid by window...
    with rasterio.open(str_demPath) as ds_dem:

        dem_crs = ds_dem.crs
        nodata_val = ds_dem.nodata

        # Define the schema for the output bank points shapefile...
        schema = {'geometry': 'Point',
                  'properties': {'xn_num': 'int', 'linkno': 'int', 'bank_hght': 'float', 'bank_elev': 'float',
                                 'bnk_ang_1': 'float', 'bnk_ang_2': 'float', 'bf_area': 'float', 'chan_width': 'float',
                                 'obank_rat': 'float',
                                 'area_ratio': 'float'}}

        with fiona.open(str_bankptsPath, 'w', driver='ESRI Shapefile', crs=dem_crs, schema=schema) as bankpts:

            #            k=0
            j = 0
            for indx in arr_strides:

                #                progBar.setValue(k)
                #                k+=1

                #                logger.info('\tIndex {} - {}/{}'.format(j,int(indx),xn_count))
                df_xn_elev_n = df_xn_elev.iloc[j:int(indx)]
                j = int(indx) + 1

                #                logger.info('\tIndex 143613 - 143712/{}'.format(xn_count))
                #                df_xn_elev_n = df_xn_elev.iloc[143613:143712]

                # << INTERPOLATE XNs >>  
                df_bank_metrics = pd.DataFrame(
                    analyze_xnelev(df_xn_elev_n, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh, nodata_val),
                    columns=['xn_no', 'left_ind', 'right_ind', 'bank_height', 'slices', 'linkno', 'bank_elev',
                             'lf_bank_ang', 'rt_bank_ang', 'bankful_area', 'chan_width', 'overbank_ratio',
                             'area_ratio'])

                df_bank_metrics.set_index('xn_no', inplace=True)

                df_map = pd.merge(df_xn_elev, df_bank_metrics, left_index=True, right_index=True)

                lst_lfbank_row = []
                lst_lfbank_col = []
                lst_rtbank_row = []
                lst_rtbank_col = []

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
                #                logger.info('Writing to bank points shapefile...')
                for tpl_row in df_map.itertuples():
                    #                    progDialog.setValue(k)
                    #                    k+=1

                    tpl_left = (tpl_row.lfbank_x, tpl_row.lfbank_y)
                    tpl_right = (tpl_row.rtbank_x, tpl_row.rtbank_y)

                    # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                    lf_pt = {'type': 'Point', 'coordinates': tpl_left}
                    rt_pt = {'type': 'Point', 'coordinates': tpl_right}

                    prop_lf = {'xn_num': int(tpl_row.Index), 'linkno': int(tpl_row.linkno_x),
                               'bank_hght': tpl_row.bank_height, 'bank_elev': tpl_row.bank_elev,
                               'bnk_ang_1': tpl_row.lf_bank_ang, 'bf_area': tpl_row.bankful_area, 'bnk_ang_2': -9999.,
                               'chan_width': tpl_row.chan_width, 'obank_rat': tpl_row.overbank_ratio,
                               'area_ratio': tpl_row.area_ratio}

                    prop_rt = {'xn_num': int(tpl_row.Index), 'linkno': int(tpl_row.linkno_x),
                               'bank_hght': tpl_row.bank_height, 'bank_elev': tpl_row.bank_elev,
                               'bnk_ang_2': tpl_row.rt_bank_ang, 'bf_area': tpl_row.bankful_area, 'bnk_ang_1': -9999.,
                               'chan_width': tpl_row.chan_width, 'obank_rat': tpl_row.overbank_ratio,
                               'area_ratio': tpl_row.area_ratio}

                    bankpts.write({'geometry': lf_pt, 'properties': prop_lf})
                    bankpts.write({'geometry': rt_pt, 'properties': prop_rt})


#                sys.exit() # for testing

# ==================================================================================
#      Floodplain Xn analysis              
# ==================================================================================                    
def read_fp_xns_shp_and_get_1d_fp_metrics(str_xns_path, str_fp_path, str_dem_path, logger):
    '''
    1. Read Xn file with geopandas, groupby linkno
    2. Linkno extent window like below using rasterio
    3. For each Xn x-y pair interpolate additional points along the length with shapely
    4. Convert to array space
    5. Sample DEM and fp grids as numpy arrays
    6. Calculate metrics
    '''
    # Depth:
    lst_min_d = []
    lst_max_d = []
    lst_rng_d = []
    lst_mean_d = []
    lst_std_d = []
    lst_sum_d = []
    # Elevation:
    lst_min_e = []
    lst_max_e = []
    lst_rng_e = []
    lst_mean_e = []
    lst_std_e = []
    lst_sum_e = []

    lst_id = []
    lst_width = []
    lst_index = []

    # Read xn file:
    logger.info('Reading Xn file...')
    gdf_xns = gpd.read_file(str_xns_path)
    # Groupby linkno:
    gp_xns = gdf_xns.groupby('linkno')

    # Access the floodplain and DEM grids:
    with rasterio.open(str(str_dem_path)) as ds_dem:
        with rasterio.open(str(str_fp_path)) as ds_fp:
            # Loop over the linkno groups:
            for linkno, gdf in gp_xns:

                #                if linkno!=3343:continue # for testing

                # Loop over the Xns along this linkno:
                for i, tpl in enumerate(gdf.itertuples()):
                    try:
                        # Xn ID and index for saving:                    
                        lst_id.append(i)
                        lst_index.append(tpl.Index)
                    except Exception as e:
                        print(f'Error with itertuple: {str(e)}')

                    try:
                        # Mask the floodplain grid with each Xn:
                        w_fp, w_trans = rasterio.mask.mask(ds_fp, [mapping(tpl.geometry)], crop=True)
                        w_fp = w_fp[0]
                        w_fp = w_fp[w_fp != ds_fp.nodata]  # ignore nodata vals

                        num_pixels = w_fp.size  # number of fp pixels along the xn
                        tot_width = num_pixels * ds_fp.res[0]  # num pixels times cell resolution
                        lst_width.append(tot_width)
                        #                net_width=tot_width-ch_wid # where to get ch_wid for a specific cross-section?
                    except:
                        lst_width.append(-9999.)

                    try:
                        # Relative elevation (depth) metrics:
                        min_depth = w_fp.min()
                        max_depth = w_fp.max()
                        rng_depth = max_depth - min_depth
                        mean_depth = w_fp.mean()
                        std_depth = w_fp.std()
                        sum_depth = w_fp.sum()
                        lst_min_d.append(min_depth)
                        lst_max_d.append(max_depth)
                        lst_rng_d.append(rng_depth)
                        lst_mean_d.append(mean_depth)
                        lst_std_d.append(std_depth)
                        lst_sum_d.append(sum_depth)
                    except:
                        lst_min_d.append(-9999.)
                        lst_max_d.append(-9999.)
                        lst_rng_d.append(-9999.)
                        lst_mean_d.append(-9999.)
                        lst_std_d.append(-9999.)
                        lst_sum_d.append(-9999.)

                    try:
                        # Also mask the DEM to get absolute elevation metrics:
                        w_dem, w_trans = rasterio.mask.mask(ds_dem, [mapping(tpl.geometry)], crop=True)
                        w_dem = w_dem[0]
                        w_dem = w_dem[w_dem != ds_dem.nodata]

                        # Absolute elevation metrics:
                        min_elev = w_dem.min()
                        max_elev = w_dem.max()
                        rng_elev = max_elev - min_elev
                        mean_elev = w_dem.mean()
                        std_elev = w_dem.std()
                        sum_elev = w_dem.sum()
                        lst_min_e.append(min_elev)
                        lst_max_e.append(max_elev)
                        lst_rng_e.append(rng_elev)
                        lst_mean_e.append(mean_elev)
                        lst_std_e.append(std_elev)
                        lst_sum_e.append(sum_elev)
                    except:
                        lst_min_e.append(-9999.)
                        lst_max_e.append(-9999.)
                        lst_rng_e.append(-9999.)
                        lst_mean_e.append(-9999.)
                        lst_std_e.append(-9999.)
                        lst_sum_e.append(-9999.)

                        # Initialize fields:
            gdf_xns['xn_id_1dfp'] = -9999.
            gdf_xns['totwid_1dfp'] = -9999.
            # Depth:
            gdf_xns['mindep_1dfp'] = -9999.
            gdf_xns['maxdep_1dfp'] = -9999.
            gdf_xns['rngdep_1dfp'] = -9999.
            gdf_xns['meandep_1dfp'] = -9999.
            gdf_xns['stddep_1dfp'] = -9999.
            gdf_xns['sumdep_1dfp'] = -9999.
            # Elevation:
            gdf_xns['minele_1dfp'] = -9999.
            gdf_xns['maxele_1dfp'] = -9999.
            gdf_xns['rngele_1dfp'] = -9999.
            gdf_xns['meanele_1dfp'] = -9999.
            gdf_xns['stdele_1dfp'] = -9999.
            gdf_xns['sumele_1dfp'] = -9999.

            # Add new values:
            gdf_xns.loc[lst_index, 'xn_id_1dfp'] = lst_id
            gdf_xns.loc[lst_index, 'totwid_1dfp'] = lst_width
            # Depth:
            gdf_xns.loc[lst_index, 'mindep_1dfp'] = lst_min_d
            gdf_xns.loc[lst_index, 'maxdep_1dfp'] = lst_max_d
            gdf_xns.loc[lst_index, 'rngdep_1dfp'] = lst_rng_d
            gdf_xns.loc[lst_index, 'meandep_1dfp'] = lst_mean_d
            gdf_xns.loc[lst_index, 'stddep_1dfp'] = lst_std_d
            gdf_xns.loc[lst_index, 'sumdep_1dfp'] = lst_sum_d
            # Elevation:
            gdf_xns.loc[lst_index, 'minele_1dfp'] = lst_min_e
            gdf_xns.loc[lst_index, 'maxele_1dfp'] = lst_max_e
            gdf_xns.loc[lst_index, 'rngele_1dfp'] = lst_rng_e
            gdf_xns.loc[lst_index, 'meanele_1dfp'] = lst_mean_e
            gdf_xns.loc[lst_index, 'stdele_1dfp'] = lst_std_e
            gdf_xns.loc[lst_index, 'sumele_1dfp'] = lst_sum_e

            # Save it again:
            gdf_xns.to_file(str_xns_path)

    '''
    # Get total FP width...
    #fpelev_len = len(lst_thisxn_fpelev)
    this_fpwidth = len(lst_thisxn_fpelev)*xnpt_dist # NOTE: need to subtract out channel width --> WHY CAN'T I MODIFY THIS LINE??

    this_fpwidth_net = this_fpwidth - ch_wid

    if this_fpwidth > 0:
        # Calculate elevation metrics...
        this_fp_min = min(lst_thisxn_fpelev)
        this_fp_max = max(lst_thisxn_fpelev)
        this_fp_range = this_fp_max - this_fp_min
        this_fp_mean = np.mean(lst_thisxn_fpelev)
        this_fp_std = np.std(lst_thisxn_fpelev)
        this_fp_sum = np.sum(lst_thisxn_fpelev)
    else:
        this_fp_min = -9999
        this_fp_max = -9999
        this_fp_range = -9999
        this_fp_mean = -9999
        this_fp_std = -9999
        this_fp_sum = -9999
    '''


# ===================================================================================
#  Build the Xns for all reaches and write to shapefile
# ===================================================================================
def write_xns_shp(df_coords, streamlines_crs, str_xns_path, bool_isvalley, p_xngap, logger):
    """
        Builds Xns from x-y pairs representing shapely interpolations along a reach

        Input: a list of tuples (row, col, linkno) for a reach
        
        Output: list of tuples of lists describing the Xn's along a reach (row, col)
        
    """
    j = 0
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

    lst_xnrowcols = []  # the final output, a list of tuples of XY coordinate pairs for all Xn's for this reach

    XnCntr = 0
    #    m_init = 0

    gp_coords = df_coords.groupby('linkno')

    # Create the Xn shapefile for writing...
    #    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int', 'endpt1_x':'float', 'endpt1_y':'float', 'endpt2_x':'float', 'endpt2_y':'float'}}
    test_schema = {'geometry': 'LineString', 'properties': {'linkno': 'int', 'strmord': 'int'}}

    logger.info('Building and Writing Cross Section File...')
    with fiona.open(str_xns_path, 'w', driver='ESRI Shapefile', crs=streamlines_crs, schema=test_schema) as chan_xns:

        for i_linkno, df_linkno in gp_coords:

            i_linkno = int(i_linkno)
            i_order = int(df_linkno.order.iloc[0])

            #            if i_linkno != 177:
            #                continue

            #            progBar.setValue(j)
            j += 1

            #                # ======================================
            #                QtCore.QCoreApplication.processEvents()
            #                if progdialog.wasCanceled():
            #                    break
            #
            #                progdialog.setValue(j)
            #                # ======================================

            # NOTE:  Define Xn length (p_xnlength) -- and other parameters? -- relative to stream order
            # Settings for stream channel cross-sections:
            p_xnlength, p_fitlength = get_xn_length_by_order(i_order, bool_isvalley)

            reach_len = len(df_linkno['x'])

            if reach_len <= p_xngap:
                #                logger.info('Less than!')
                continue  # skip it for now

            # Loop along the reach at the specified intervals...(Xn loop)
            for i in range(p_xngap, reach_len - p_xngap, p_xngap):

                lstThisSegmentRows = []
                lstThisSegmentCols = []

                if p_fitlength > i or i + p_fitlength >= reach_len:  # if i + paramFitLength > reach_len
                    fitLength = p_xngap
                else:
                    fitLength = p_fitlength

                lstThisSegmentRows.append(df_linkno['y'].iloc[i + fitLength])
                lstThisSegmentRows.append(df_linkno['y'].iloc[i - fitLength])
                lstThisSegmentCols.append(df_linkno['x'].iloc[i + fitLength])
                lstThisSegmentCols.append(df_linkno['x'].iloc[i - fitLength])

                midPtRow = df_linkno['y'].iloc[i]
                midPtCol = df_linkno['x'].iloc[i]

                # Send it the endpts of what you to draw a perpendicular line to...
                lst_xy = build_xns(lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow,
                                   p_xnlength)  # returns a list of two endpoints

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
                line = {'type': 'LineString', 'coordinates': lst_xy}
                #                prop = {'linkno': i_linkno, 'endpt1_x':lst_xy[0][0], 'endpt1_y':lst_xy[0][1], 'endpt2_x':lst_xy[1][0], 'endpt2_y':lst_xy[1][1]}
                prop = {'linkno': i_linkno, 'strmord': i_order}
                chan_xns.write({'geometry': line, 'properties': prop})

            #                if XnCntr > 10:
    #                    break

    return lst_xnrowcols


# ===================================================================================
#  Build Xn's based on vector features
# ===================================================================================
def get_stream_coords_from_features(str_streams_filepath, cell_size, str_reachid, str_orderid, logger):
    #        try:
    #    lst_x=[]
    #    lst_y=[]
    #    lst_linkno=[]
    #    lst_order=[]
    lst_df_final = []

    p_interp_spacing = int(
        cell_size)  # 3 # larger numbers would simulate a more smoothed reach | NOTE: Hardcode this = grid resolution?
    j = 0  # prog bar

    # Open the streamlines shapefile...
    with fiona.open(str(str_streams_filepath), 'r') as streamlines:

        # Get the crs...
        streamlines_crs = streamlines.crs
        #        str_proj4 = crs.to_string(streamlines.crs)

        #            progBar.setRange(0,len(streamlines))

        #        # ==================================
        #        progdialog = QtGui.QProgressDialog("Getting stream coords from features...",
        #        "Cancel", 0, len(streamlines))
        #        progdialog.setWindowTitle("FACET")
        #        progdialog.setWindowModality(QtCore.Qt.WindowModal)
        #        progdialog.resize(350, 110)
        #        progdialog.show()
        #        # ==================================
        tot = len(streamlines)
        for line in streamlines:

            #           # ======================================
            #           QtCore.QCoreApplication.processEvents()
            #           if progdialog.wasCanceled():
            #               break
            #
            #           progdialog.setValue(j)
            #           # ======================================

            j += 1
            #               self.emit(QtCore.SIGNAL("update(int)"), int(100*len(streamlines)/j))

            line_shply = LineString(line['geometry']['coordinates'])

            length = line_shply.length  # units depend on crs

            if length > 9:  # NOTE: This value is dependent on CRS!!

                i_linkno = line['properties'][str_reachid]
                i_order = line['properties'][str_orderid]

                #               logger.info(f'{i_linkno}; {j}|{tot}')
                #               print(f'{i_linkno}; {j}|{tot}')

                #               if i_linkno == 3343:
                #                   print('he')
                #               else:
                #                   continue

                # Smoothing reaches via Shapely...
                if i_order <= 3:
                    line_shply = line_shply.simplify(5.0, preserve_topology=False)
                elif i_order == 4:
                    line_shply = line_shply.simplify(10.0, preserve_topology=False)
                elif i_order == 5:
                    line_shply = line_shply.simplify(20.0, preserve_topology=False)
                elif i_order >= 6:
                    line_shply = line_shply.simplify(30.0, preserve_topology=False)

                length = line_shply.length

                int_pts = np.arange(0, length, p_interp_spacing)  # p_interp_spacing in projection units?

                lst_x = []
                lst_y = []
                lst_linkno = []
                lst_order = []
                for i in int_pts:  # lambda here instead?
                    i_pt = np.array(line_shply.interpolate(i))
                    lst_x.append(i_pt[0])
                    lst_y.append(i_pt[1])
                    lst_linkno.append(i_linkno)
                    lst_order.append(i_order)

                df_coords = pd.DataFrame({'x': lst_x, 'y': lst_y, 'linkno': lst_linkno, 'order': lst_order})
                #           if j == 100:
                #               logger.info('pause')

                df_coords.drop_duplicates(subset=['x', 'y'],
                                          inplace=True)  # potential duplicates due to interpolation (?)
                lst_df_final.append(df_coords)

        df_final = pd.concat(lst_df_final)
    #        except Exception as e:
    #            logger.info('Error!: {}'.format(e.strerror))
    #            QtGui.QMessageBox.information(self, 'Error!', '{}'.format(e.strerror))

    return df_final, streamlines_crs  # A list of lists
