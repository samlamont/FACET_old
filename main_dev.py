# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:11:00 2016

@author: sam.lamont
"""

from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtWidgets import QMainWindow, QApplication

#from PyQt4 import QtCore, QtGui, uic
#from PyQt4 import QMainWindow, QApplication # Py

#import time
import glob
import timeit
#import sys
import numpy as np
#from numpy import array
#from scipy import ndimage
import os
import rasterio
#from rasterio.plot import show
#from rasterio import features
#from affine import Affine # Georeferenced raster datasets use affine transformations to map from image coordinates to world coordinates
#from shapely.geometry import LineString
import fiona
#from fiona.crs import to_string
#from math import isinf, sqrt #, hypot, modf, atan, ceil
#import matplotlib.pyplot as plt
import pandas as pd

#from functools import partial
#import multiprocessing

import funcs_v2 
#from funcs_v2 import workerfuncs

# OR from funcs_v2 import <<ClassName>>

#from shapely import speedups
#
#if speedups.available:
#    speedups.enable() # enable performance enhancements written in C

#plt.ion() 

#class ResultObj(QtCore.QObject):
#    def __init__(self, val):
#        self.val = val 
#
#class GenericThread(QtCore.QThread):
#    finished = QtCore.pyqtSignal(object)
#  
#    def __init__(self, function, *args):
#        QtCore.QThread.__init__(self)
#        self.function = function
#        self.args = args
##        self.kwargs = kwargs
#        self.finished.connect(function)
#     
#    def __del__(self):
#        self.wait()
#     
#    def run(self):
##        self.function(*self.args) #,**self.kwargs)
#        self.finished.emit(ResultObj(self.function(*self.args)))
##        return
  
# ===================================================================================
#                                       MAIN
# ===================================================================================
class facet(QtGui.QMainWindow): # 2.7, PyQt4
#class facet(QtGui.QDialog):
#class facet(QtGui.QTabWidget):
#class MainForm(QtGui.QDialog, FORM_CLASS):
    
    def __init__(self, parent=None):   
        
        super(facet, self).__init__(parent)
        
        # ================= UI ================
        uic.loadUi('facet.ui', self)  
              
        # << Browse buttons >>
        # Pre-processing...
        self.browse_dem_pre.clicked.connect(self.get_dem_file)                      # Raw DEM for pre-processing
        self.browse_nhd_pre.clicked.connect(self.get_streams_file)                  # Input NHD for weight grid     
        self.browse_hand_pre.clicked.connect(self.set_hand_file)                    # Output HAND grid 
        self.browse_streamlines_pre.clicked.connect(self.set_streamlines_pre_file)  # Output delineated streamlines
        self.browse_taudem_pre.clicked.connect(self.get_taudem_pre_path)            # Path to directory containing TauDEM executables
        self.browse_whitebox_pre.clicked.connect(self.get_whitebox_pre_exe)         # Path to the Whitebox exe
        self.browse_mpi_pre.clicked.connect(self.get_mpi_pre_exe)                   # Path to the MPI exe (for TauDEM parallelization)
        
        # Original Xn Method...
        self.browseDEM_Xns.clicked.connect(self.get_dem_file)               # DEM - Xns
        self.browseStreams_Xns.clicked.connect(self.get_streams_file)       # Streams - Xns
        self.browseXns.clicked.connect(self.get_or_set_xns_file)            # Xns
        self.browseBankPts.clicked.connect(self.set_bankpts_file)           # BankPts
        
        # Mean Curvature Method...
        self.browseDEM_Curve.clicked.connect(self.get_dem_file)             # DEM - curvature
        self.browseStreams_Curve.clicked.connect(self.get_streams_file)     # Streams - curvature
        self.browseBankPixels.clicked.connect(self.get_or_set_bankpixels)   # Bank Pixels   
        
        # Wavelet Curvature Method...
        
        
        # Floodplain analysis...
        
               
#        self.connect( self, QtCore.SIGNAL("add(QString)"), self.run_main )
#        self.genThread = GenericThread(self.run_main)        
        
        # << Button box >>      
#        self.buttonBox.accepted.connect(self.genThread.start) 
        self.buttonBox.accepted.connect(self.run_main)        # run_main in a thread?
        self.buttonBox.rejected.connect(self.close)
                     
        # << Defaults >>
        self.textFitLength.setText('3')
        self.textXnSpacing.setText('3')
        self.textXnLength.setText('30')
        self.textVertIncrement.setText('0.2')
        self.textRatioThresh.setText('1.5')
        self.textSlpThresh.setText('0.03')        
        
        # << Progress Bar >>
        self.progressBar.setVisible(False)
#        self.progressBar.setTextVisible(False)        
        
        # << Load Previous Paths from Settings >>
        self.settings = QtCore.QSettings('taudem_path', 'facet')        
        self.settings = QtCore.QSettings('mpi_path', 'facet')
        self.settings = QtCore.QSettings('whitebox_path', 'facet')
        
        td_fname = self.settings.value('taudem_path', type=str)        
        self.txt_taudem_pre.setText(td_fname)
        
        mpi_fname = self.settings.value('mpi_path', type=str)        
        self.txt_mpi_pre.setText(mpi_fname)
        
        wb_fname = self.settings.value('whitebox_path', type=str)        
        self.txt_whitebox_pre.setText(wb_fname)        
        
        
        # =====================================
        
    # =========================================    
    #              << FILE I/O >>    
    # =========================================        
    def get_dem_file(self):
        fname_dem = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
        
        if self.tabWidget.currentIndex() == 0:      # Pre-processing        
            self.txt_dem_pre.setText(fname_dem)        
        
        if self.tabWidget.currentIndex() == 1:      # Channel metrics via Xn's        
            self.textDEM_xns.setText(fname_dem)

        if self.tabWidget.currentIndex() == 2:      # Channel metrics via Curvature    
            self.textDEM_curvature.setText(fname_dem)
        
        # << Get DEM and set units, get cell size >>
        with rasterio.open(str(fname_dem)) as ds_dem:
            
            self.dem_crs = ds_dem.crs      
            self.dem_affine = ds_dem.affine            
#            self.dem_affine = ds_dem.transform # NOTE: Confused about affine vs. transform here
            self.cell_size = ds_dem.affine[0]
            self.nodata_val = ds_dem.nodata
            
#            print('DEM affine: {}'.format(ds_dem.affine))
            print('DEM crs: {}'.format(ds_dem.crs))
#            print('ds_dem.crs.wkt: {}'.format(ds_dem.crs.wkt))
            
            try:
#                self.units = self.dem_crs['units']  
                self.lbl_fit_units.setText(self.dem_crs['units'])
                self.lbl_xngap_units.setText(self.dem_crs['units'])
                self.lbl_xnlength_units.setText(self.dem_crs['units'])
                self.lbl_ivert_units.setText(self.dem_crs['units'])             
            except:
                self.lbl_fit_units.setText('unknown')
                self.lbl_xngap_units.setText('unknown')
                self.lbl_xnlength_units.setText('unknown')
                self.lbl_ivert_units.setText('unknown')                             
            
    def get_streams_file(self):
        fname_streams = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Shapefiles (*.shp)')

        if self.tabWidget.currentIndex() == 0:        
            self.txt_nhd_pre.setText(fname_streams) 
            
#            # Also load the reach ID combo box here...
#            with fiona.open(str(fname_streams), 'r') as streamlines:
#                lst_fieldnames = streamlines.schema['properties'].keys()      
#                self.streamlines_crs = streamlines.crs
#                print('Streamlines crs: {}'.format(streamlines.crs))
#            if lst_fieldnames:
#                for name in lst_fieldnames:
#                    self.comboBoxXns.addItem(str(name)) 
                    
        if self.tabWidget.currentIndex() == 1:        
            self.textStreams_Xns.setText(fname_streams) 
            
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.comboBoxXns.addItem(str(name))            

        if self.tabWidget.currentIndex() == 2:        
            self.textStreams_curvature.setText(fname_streams) 
        
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.comboBoxCurvature.addItem(str(name))

    def get_taudem_pre_path(self):
        fname = QtGui.QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.txt_taudem_pre.setText(fname) 
        self.settings.setValue('taudem_path', fname) # Save to settings
        
    def get_mpi_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.exe*)')
        self.txt_mpi_pre.setText(fname)  
        self.settings.setValue('mpi_path', fname) # Save to settings
    
    def get_whitebox_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.*)')
        self.txt_whitebox_pre.setText(fname) 
        self.settings.setValue('whitebox_path', fname) # Save to settings

    def set_hand_file(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Grids (*.*)')
        self.txt_hand_pre.setText(fname) 
        
    def set_streamlines_pre_file(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
        self.txt_streamlines_pre.setText(fname)        

    def set_bankpts_file(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
        self.textBankPts.setText(fname)  
        
#    def get_bankpts_file(self):
#        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
#        self.bankptsPath.setText(fname)  
        
    def get_or_set_xns_file(self):
        if self.chkCreateXns.isChecked():
            fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
            self.textXns.setText(fname)
        else:
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
            self.textXns.setText(fname)  
        
    def get_or_set_bankpixels(self):
        if self.chkBankPixels.isChecked():
            fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Grids (*.*)')
            self.textBankPixels.setText(fname)   
        else:
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
            self.textBankPixels.setText(fname)             
    
    # =========================================    
    #              << RUN MAIN >>    
    # =========================================     
    def run_main(self):        

        # Gray out the tabs...   
#        self.tabWidget.setEnabled(False) 
        self.hide()                             
        start_time_0 = timeit.default_timer()
        
        str_orderid='Order_' # TESTING
       
#         Check which tab is active??
#        i_tab = self.tab_widget.currentIndex() # (0: bankpts; 1: bankpixels; 2: FP-HAND)
#        print('selected tab:  {}'.format(i_tab))
              
        try:    
            
            # << Pre-processing >>
            if self.tabWidget.currentIndex() == 0: # Run pre-processing functions
            
                print('About to preprocess')
                funcs_v2.preprocess_dem(str(self.txt_dem_pre.text()), str(self.txt_nhd_pre.text()), str(self.txt_mpi_pre.text()), str(self.txt_taudem_pre.text()), str(self.txt_whitebox_pre.text()), self.chk_wg_pre.isChecked(), self.chk_whitebox_pre.isChecked(), self.chk_taudem_pre.isChecked())            
            
            # << Xn-based analysis >>
            if self.tabWidget.currentIndex() == 1:
                
                # Check coordinate systems (DEM vs. Streamlines)...
                if dict(self.dem_crs) != dict(self.streamlines_crs):
                    QtGui.QMessageBox.critical(self, 'Warning!', 'Coordinate Systems of DEM and Streams Layer May Not Match')                            
    
                # Build reach coords and get crs from a pre-existing streamline shapefile...
                funcs = funcs_v2.workerfuncs()            
                    
                # << Generate Cross Sections >>
                if self.chkCreateXns.isChecked():
                    
    #                str_streams=self.textStreams_Xns.text()
    #                cell_size = self.cell_size
    #                str_streamid = str(self.comboBoxXns.currentText())
    #                
    #                self.lblProgBar.setText('Building stream coordinates...') 
    #                self.genThread = GenericThread(funcs_v2.get_stream_coords_from_features, str_streams, cell_size, str_streamid, str_orderid)
    #                
    #                self.connect(self.genThread, QtCore.SIGNAL("finished()"), self.done)
    #
    #                self.genThread.start()
                    
                    # Need to wait here
                    
    #                self.lblProgBar.setText('Building stream coordinates...')               
    
                    
    #                QtCore.QObject.connect(func,QtCore.SIGNAL("update(int)"), self.updateUI)                
                    df_coords, streamlines_crs = funcs_v2.workerfuncs.get_stream_coords_from_features(funcs, self.textStreams_Xns.text(), self.cell_size, str(self.comboBoxXns.currentText()), str_orderid)
        
    #                self.lblProgBar.setText('Generating cross sections...')
    #                # Build and write the Xn shapefile...
                    funcs_v2.workerfuncs.write_xns_shp(funcs, df_coords, streamlines_crs, str(self.textXns.text()), False, int(self.textXnSpacing.text()), int(self.textFitLength.text()), float(self.textXnLength.text()))     
    #            
    #                self.genThread.terminate() # stop thread
                    
                # << Calculate Channel Metrics -- Xns Method >>
                if self.chkMetricsXns.isChecked():
        
    #                self.lblProgBar.setText('Extracting elevation along cross sections...')
                    # Read the cross section shapefile and interplote elevation along each one using the specified DEM
                    df_xn_elev = funcs_v2.workerfuncs.read_xns_shp_and_get_dem_window(funcs, self.textXns.text(), self.textDEM_xns.text())
                    
                    XnPtDist = self.cell_size
                    
    #                self.lblProgBar.setText('Calculating channel metrics...')
                    # Calculate channel metrics and write bank point shapefile...# parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh
                    funcs_v2.workerfuncs.chanmetrics_bankpts(funcs, df_xn_elev, str(self.textXns.text()), str(self.textDEM_xns.text()), str(self.textBankPts.text()), float(self.textVertIncrement.text()), XnPtDist, float(self.textRatioThresh.text()), float(self.textSlpThresh.text()))
            
            # << Mean Curvature analysis >>
            if self.tabWidget.currentIndex() == 2:            
            
                # << Create bank pixel layer (.tif) >>
                if self.chkBankPixels.isChecked():
                     
                    self.lblProgBar.setText('Reading streamline coordinates...') 
                    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(self.textStreams_curvature.text(), self.cell_size, str(self.comboBoxCurvature.currentText()), str_orderid, self)
                    
                    self.lblProgBar.setText('Creating bank pixel layer...') 
                    funcs_v2.bankpixels_from_streamline_window(df_coords, str(self.textDEM_curvature.text()), str(self.textBankPixels.text()), self)
                    
                # << Calculate Channel Metrics -- Bank Pixels Method >>
                if self.chkPixelMetrics.isChecked():
                    
                    self.lblProgBar.setText('Calculating channel width from bank pixel layer!')
                    funcs_v2.channel_width_bankpixels(self.textStreams_curvature.text(), self.textBankPixels.text(), self.comboBoxCurvature .currentText(), self.cell_size, self)
                
    #        # << Calculate Floodplain >>
    #        if self.ck_fp.isChecked():
    #            print('Create floodplain here!')
                
    #        # << Run All Functions >>
    #        if self.ck_runall_bankpts.isChecked():
    #            print('Run all is currently disabled')
    #
    #        if self.ck_runall_pixels.isChecked():
    #            print('Run all is currently disabled')
                
            print('Done with main!')
#            self.lblProgBar.setText('Done!                                Elasped Time:  {0:.2f} mins'.format(float(timeit.default_timer() - start_time_0))/60.)
            
        except BaseException as e:
            QtGui.QMessageBox.critical(self, 'Error!', '{}'.format(str(e)))
                
        self.tabWidget.setEnabled(True) 
#        self.close()
        
#        self.emit( QtCore.SIGNAL('add(QString)'), 'run_main complete')
        
        # Show message when complete...
        QtGui.QMessageBox.information(self, 'Code Complete', 'Successful Run! \n\n\t   {0:.2f} secs'.format(timeit.default_timer() - start_time_0))       

#    def done(self):
#        QtGui.QMessageBox.information(self, "Done!", "Done whatevah!")

 
#    def updateUI(self, val):
#        self.progressBar.setValue(val)   
        
if __name__ == '__main__':    
    
###    # ============================ << FOR DEBUGGING >> ==================================
    print('\n<<< Start >>>\r\n')
    start_time_0 = timeit.default_timer() 
    
#    funcs = funcs_v2.workerfuncs()
##    
    ## << PATHS >>
    ## DEM:
    str_dem_path = r"D:\hand\nfie\020700\020700dd.tif"
    
    ## Slope:
    str_slp_path = r"D:\drb\02040205\02040205_breach_dem_sd8.tif"
    
    ## Vector streamlines:
    str_net_in_path = r"D:\hand\nfie\020700\020700-flows.shp"
      
    ## Bank pixels:   
#    str_bankpixels_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpixels_PO.tif'
    
    ## Bank points:   
#    str_bankpts_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpts_TEST.shp'
 
    ## Floodplain:
#    str_floodplain_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM_breach_hand_slice2.3.tif'    
    
    ## Cross sections:    
    str_xns_path = r"D:\hand\nfie\020700\020700_xns_test.shp"

    ## HAND:
    str_hand_path = r"D:\drb\02040205\02040205_breach_hand.tif"
    
    ## Catchments:
    str_sheds_path = r"D:\drb\02040205\02040205_w_diss_physio.shp"
    
    ## Channel start points:
#    str_startptgrid_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMnet_UNIQUE_ID.shp'
    
    ## Openness: THIS NEEDS WORK -- Table this til later if you have time (July5)
#    str_pos_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\dr_pos_raw.tif'    
    
    ## Flow direction: For discerning between right/left bank?
#    str_fdr_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach_p.tif"    
        
    str_reachid='LINKNO'
#    str_reachid='ARCID'
#    str_reachid='COMID'
    str_orderid='strmOrder'
#    cell_size=3
    
#        self.p_fitlength.setText('3')
#        self.p_xngap.setText('3')
#        self.p_xnlength.setText('30')
#        self.p_ivert.setText('0.2')
#        self.p_ratiothresh.setText('1.5')
#        self.p_slpthresh.setText('0.03')    
    
    # << PARAMETERS >>
    parm_ivert = 0.2 # 0.2 default
    XnPtDist = 3 # 3 is default
    parm_ratiothresh = 1.5 # 1.5 default
    parm_slpthresh = 0.03 # 0.03 default
    p_buffxnlen = 30 # meters (if UTM) ??
    p_xngap = 3 # 3 default
    
    # Width from curvature via buffering...
    use_wavelet_curvature_method = False
    i_step = 100 # length of reach segments for measuring width from bank pixels (and others?)
    max_buff = 30 # maximum buffer length for measuring width based on pixels
    
    # By stream order:
    # - Line smoothing parameters
    # - Curvature window size
    # - Xn length
    # - Above Xn parameters    

    #=============================================================================================== 
    #                             BEGIN BULK PROCESSING LOOP
    #===============================================================================================    
    
    ## << FOR BULK PROCESSING >>
    lst_paths = glob.glob('/home/sam/nwm/bulk_processing/*')
    lst_paths.sort() # for testing
    
    for i, path in enumerate(lst_paths):
        
        if i <= 3: continue
        print('Processing:  ' + path)
        
        start_time_i = timeit.default_timer()
    
        try:
            str_dem_path = glob.glob(path + '/*dem*.tif')[0]
            str_hand_path = glob.glob(path + '/*hand*.tif')[0]
            str_net_path = glob.glob(path + '/*net*.shp')[0]    
            str_sheds_path = glob.glob(path + '/*w_diss_physio*.shp')[0]
        except:
            print('WARNING:  There is an error in the paths!')
            pass # depending on what's being run, it might not matter if a file doesn't exist
        
        path_to_dem, dem_filename = os.path.split(str_dem_path)
        csv_filename = dem_filename[:-8] + '.csv'
        str_csv_path = path_to_dem + '/' + csv_filename
        
        # Output layers...
        str_xns_path = path_to_dem + '/' + dem_filename[:-8] + '_xns.shp'
        str_bankpts_path = path_to_dem + '/' + dem_filename[:-8] + '_bankpts.shp'
        str_bankpixels_path = path_to_dem + '/' + dem_filename[:-8] + '_bankpixels.tif'
        
        # << GET CELL SIZE >>
        cell_size = int(funcs_v2.get_cell_size(str_dem_path)) # range functions need int?        

        # << BUILD STREAMLINES COORDINATES >>
#        df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_path, cell_size, str_reachid, str_orderid) # YES!        
#        df_coords.to_csv(str_csv_path)
        df_coords = pd.read_csv(str_csv_path, )    
        streamlines_crs = {'init': u'epsg:26918'} # NAD83, UTM18N     

        # ============================= << CROSS SECTION ANALYSES >> =====================================
#        # << CREATE Xn SHAPEFILES >>
#        funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(p_xngap), int(3), float(30))     
#
#        # << BUILD STREAMLINES COORDINATES >>
#        df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_path, cell_size, str_reachid, str_orderid) # YES!        
#        df_coords.to_csv(str_csv_path)
##        df_coords = pd.read_csv(str_csv_path, )    
##        streamlines_crs = {'init': u'epsg:26918'} # NAD83, UTM18N     
#
#        # ============================= << CROSS SECTION ANALYSES >> =====================================
##        # << CREATE Xn SHAPEFILES >>
##        funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(p_xngap), int(3), float(30))     
##
##        # << INTERPOLATE ELEVATION ALONG Xns >>
##        df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path)
##        
##        # Calculate channel metrics and write bank point shapefile...# NOTE:  Use raw DEM here??        
##        funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
#        
#        # ========================== << BANK PIXELS AND WIDTH FROM CURVATURE >> ====================================
#        funcs_v2.bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size, use_wavelet_curvature_method) # YES!        
#
#        funcs_v2.channel_width_from_bank_pixels(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size, i_step, max_buff)        
#        # Calculate channel metrics and write bank point shapefile...# NOTE:  Use raw DEM here??        
#        funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
        
        # ========================== << BANK PIXELS AND WIDTH FROM CURVATURE >> ====================================
#        funcs_v2.bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size, use_wavelet_curvature_method) # YES!        
#
#        funcs_v2.channel_width_from_bank_pixels(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size, i_step, max_buff)        
        
        # ============================= << DELINEATE FIM >> =====================================
        funcs_v2.fim_hand_poly(str_hand_path, str_sheds_path, str_reachid)
#        
#        # ============================= << DELINEATE FIM >> =====================================
##        funcs_v2.fim_hand_poly(str_hand_path, str_sheds_path, str_reachid)
##        
##        break # for testing
#        
#        # ==================== CHANNEL WIDTH, FLOODPLAIN WIDTH, HAND ANALYSIS ALL IN ONE ===========
##        funcs_v2.channel_and_fp_width_bankpixels_segments_po_2Dfpxns(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size, p_buffxnlen, str_hand_path, parm_ivert)    
#        
#        print('\nRun time for {}:  {}\r\n'.format(path, timeit.default_timer() - start_time_i))
#    #=============================================================================================== 
#    #                               END BULK PROCESSING LOOP
#    #===============================================================================================
    # << GET CELL SIZE >>
    cell_size = int(funcs_v2.get_cell_size(str_dem_path)) # range functions need int?
  
##    # << BUILD STREAMLINES COORDINATES >>
##    # Build reach coords and get crs from a pre-existing streamline shapefile...
    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_in_path, cell_size, str_reachid, str_orderid) # YES!
    
    df_coords.to_csv(r"D:\hand\nfie\020700\df_coords_020700.csv") # just for testing
    
#    print('NOTE:  Reading pre-calculated csv file...')
#    df_coords = pd.read_csv('df_coords_DifficultRun.csv')
#    df_coords = pd.read_csv('df_coords_Chillisquaque.csv', )
#    df_coords = pd.read_csv('df_coords_020802.csv', )    
#    streamlines_crs = {'init': u'epsg:26918'} # NAD83, UTM18N
    
    
    # << UNIFORM REACH POINTS AND CATCHMENTS >>    
#    str_diss_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_net_diss.shp'
#    funcs_v2.dissolve_line_features(str_net_path, str_diss_path)
    
#    str_pts_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_net_pts.shp'
#    funcs_v2.points_along_line_features(str_diss_path, str_pts_path)
    
#    str_d8fdr_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach_p.tif"
#    funcs_v2.taudem_gagewatershed(str_pts_path, str_d8fdr_path)
    
    # << FIM FROM HAND AND CATCHMENTS >>
#    funcs_v2.fim_hand_poly(str_hand_path, str_sheds_path) # NOTE:  Will need to know which regression eqn to use?
#   
#    # << BANK PIXELS FROM CURVATURE >>
#    funcs_v2.bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size) # YES!
#    funcs_v2.bankpixels_from_openness_window(df_coords, str_pos_path, str_bankpixels_path) 
#    funcs_v2.bankpixels_from_openness_window_buffer_all(df_coords, str_dem_path, str_net_path, str_pos_path, str_neg_path) 
    
#     << FLOODPLAIN WIDTH >> 
#    buff_dist = 40
#    funcs_v2.floodplain_width_2D_xns(str_xns_path, str_floodplain_path, buff_dist)
#    funcs_v2.floodplain_width_fppixels_segments_po(df_coords, str_net_in_path, str_floodplain_path, str_reachid, cell_size)
#    funcs_v2.floodplain_width_reach_buffers_po(funcs, str_net_path, str_fp_path, str_reachid, cell_size)
    
    # << CHANNEL WIDTH, FLOODPLAIN WIDTH, HAND ANALYSIS ALL IN ONE >>
#    funcs_v2.channel_and_fp_width_bankpixels_segments_po_2Dfpxns(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size, p_buffxnlen, str_hand_path, parm_ivert)    
    
    # << CREATE Xn SHAPEFILES >>
#    # Channel...
#    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(3), int(3), float(30))     
#    # FP...
##    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), True, int(30), int(30), float(100))  # For FP width testing
#
##    # << INTERPOLATE ELEVATIN ALONG Xns >>
#    df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path)
#    
##    print('Writing df_xn_elev to .csv for testing...')
##    df_xn_elev.to_csv(columns=['index','linkno','elev','xn_row','xn_col']) 
##    df_xn_elev2 = pd.read_csv('df_xn_elev.csv') #, dtype={'linko':np.int,'elev':np.float,'xn_row':np.float,'xn_col':np.float})
#    
#    # Loop over df_xn_elev here to determine the flow direction/slope from one Xn midpoint to the next?
# 
#    # Calculate channel metrics and write bank point shapefile...
#    print('Calculating channel metrics from bank points...')
#    funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
#     
     # << DEM PRE-PROCESSING using TauDEM and GoSpatial >>              
    # (1) Clip original streamlines layer (NHD hi-res 4 digit HUC to DEM of interest)...     
    # Build the output streamlines file name...
#    path_to_dem, dem_filename = os.path.split(str_dem_path)
#    str_output_nhdhires_path = path_to_dem + '\\' + dem_filename[:-4]+'_nhdhires.shp' 
    
#    funcs_v2.clip_features(str_net_in_path, str_output_nhdhires_path, str_dem_path)     
    
###    # (2) Do all Whitebox and TauDEM functions including HAND based on output r'"C:\Program Files\TauDEM\TauDEM5Exe\D8FlowDir.exe"'
#    str_mpi_path=r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe'
#    str_taudem_dir=r'C:\Program Files\TauDEM\TauDEM5Exe' #\D8FlowDir.exe"'
#    str_whitebox_path= r"C:\gospatial\go-spatial_win_amd64.exe" # Go version
####    str_whitebox_path= r'C:\Terrain_and_Bathymetry\Whitebox\WhiteboxTools\whitebox_tools.exe'   # Rust version
###    
#    run_whitebox = False
#    run_wg = False
#    run_taudem = True
####    
#    funcs_v2.preprocess_dem(str_dem_path, str_net_in_path, str_mpi_path, str_taudem_dir, str_whitebox_path, run_whitebox, run_wg, run_taudem)    

    print('\n<<< End >>>\r\n')
    print('Total Run Time:  {}'.format(timeit.default_timer() - start_time_0))
    
        # << BEGIN TEST MULTIPROCESSING >>  This did not help!
#        # Do the rest by looping in strides, rather than all at once, to conserve memory...(possibly using multiprocessing)
#        xn_count = funcs_v2.get_feature_count(str_xns_path)
#        
#        # Striding...
#        arr_strides = np.linspace(0, xn_count, xn_count/100)
#        arr_strides = np.delete(arr_strides,0)   
#        
#        lst_arr_strides = np.array_split(arr_strides, 10)
#
#        func = partial(funcs_v2.chanmetrics_bankpts, df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
#        
#        pool = multiprocessing.Pool(processes=10)
#        
#        lst_final = pool.map(func, lst_arr_strides)                
        # << END TEST MULTIPROCESSING >>    
#     =====================================================================================

#    # ===== Run GUI ==============================
#    app = QtGui.QApplication(sys.argv)
#
#    mainform = facet()
#    
#    # Show it in the center of the window...
#    fr_geo=mainform.frameGeometry()
#    cp = QtGui.QDesktopWidget().availableGeometry().center()
#    fr_geo.moveCenter(cp)
#    mainform.move(fr_geo.topLeft())    
#    
#    mainform.show()
#
#    sys.exit(app.exec_())
##    # ============================================