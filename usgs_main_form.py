# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:11:00 2016

@author: sam.lamont
"""
#import time
from glob import glob
import timeit
from PyQt4 import QtCore, QtGui, uic
import sys
#import numpy as np
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
class facet(QtGui.QMainWindow):
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
        
    def get_mpi_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.exe*)')
        self.txt_mpi_pre.setText(fname)         
    
    def get_whitebox_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.*)')
        self.txt_whitebox_pre.setText(fname)  

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
                funcs_v2.preprocess_dem(str(self.txt_dem_pre.text()), str(self.txt_nhd_pre), str(self.txt_mpi_pre), str(self.txt_taudem_pre), str(self.txt_whitebox_pre), self.chk_wg_pre.isChecked(), self.chk_whitebox_pre.isChecked(), self.chk_taudem_pre.isChecked())            
            
            # << Xn-based analysis >>
            if self.tabWidget.currentIndex() == 1:
                
                # Check coordinate systems (DEM vs. Streamlines)...
                if dict(self.dem_crs) <> dict(self.streamlines_crs):
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
    
####    # ============================ << FOR DEBUGGING >> ==================================
#    print('\n<<< Start >>>\r\n')
#    start_time_0 = timeit.default_timer() 
#    
##    funcs = funcs_v2.workerfuncs()
###    
#    ## << DEM (raw or fel?) >>
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_dem.tif'
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700dem_UTM18N_drclip_3mresample.tif'
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM_02050205_USGS_AEA.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMfel_dem_post.tif'
##    str_dem_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_raw_dem_clip_utm18.tif'
##    str_dem_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18.tif'
##    str_dem_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach.tif'
##    str_dem_path = r'D:\fernando\HAND\020802\020802_utm.tif' # 10 m
##    str_dem_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM.tif'
##    str_dem_path = r'D:\CFN_data\DEM_Files\020503010802_MahantangoCreek\DEM.tif'    
##    str_dem_path = r'D:\CFN_data\DEM_Files\020503050802_QuittapahillaCreek\DEM_utm18.tif'
##    str_dem_path = r'D:\CFN_data\DEM_Files\020600050203_ChoptankRiver\01_02_03_04_utm18_breach.tif'
##    str_dem_path = r'D:\Terrain_and_Bathymetry\USGS\DRB_2016\gis\drb_dems\02040101\02040101.dep'
##    str_dem_path = r'D:\Terrain_and_Bathymetry\OWP\grids\121003_SanAntonio\121003.dep'
##    str_dem_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\CFN_data\DEM_Files\020600020405_MorganCreek\facet\morgan_utm18.tif"
#    str_dem_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\CFN_data\DEM_Files\020600030902_GwynnsFalls\facet\DEM_utm18.tif"
#    
#    ## << VECTOR STREAMLINES (net) >>  
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/Network04.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net.shp'
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700_flowsUTM18N_drsheds_clip_diss.shp'
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_StreamNetwork_USGS_AEA_2.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMnet1000.shp'
##    str_net_in_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\nhdhires_chillisquaque.shp'
##    str_net_in_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM_nhdhires.shp'
##    str_net_in_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM_breach_net.shp'
##    str_site_nhdhires_path = r'D:\CFN_data\nhd_hires\nhdhires_choptank0206_utm18.shp'
##    str_site_nhdhires_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\nhd_hires\dr_nhd_hires_utm18.shp'
##    str_net_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach_net.shp'
##    str_net_path = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\headwaters\upperreaches.shp'
##    str_net_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\dr3m_net_raw.shp'
##    str_net_in_path = r'D:\fernando\HAND\020802\020802-flows-utm.shp'
##    str_net_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\breach\dr3m_raw_net.shp'
##    str_net_path = r'D:\CFN_data\DEM_Files\020600050203_ChoptankRiver\01_02_03_04_utm18_breach_net.shp'
#    str_net_in_path = r"D:\Terrain_and_Bathymetry\USGS\CBP_analysis\CFN_data\DEM_Files\020600030902_GwynnsFalls\facet\gwynnsfall_nhdhires_utm18.shp"
##    str_net_in_path = r'D:\Terrain_and_Bathymetry\USGS\DRB_2016\gis\drb_dems\02040101\02040101_nhdhires_utm18.shp'
##    str_net_in_path = r'D:\Terrain_and_Bathymetry\OWP\grids\121003_SanAntonio\121003\121003-flows.shp'
#      
#    ## << BANK PIXELS >>   
###    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/bank_pixels_02050205_H.tif'
##    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpixels.tif'
##    str_bankpixels_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpixels_PO.tif'
##    str_bankpixels_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\bankpixels_raw_dem.tif'
##    str_bankpixels_path =r'D:\fernando\HAND\020802\0202802_utm_bankpixels.tif'
##    str_bankpixels_path =r'D:\CFN_data\DEM_Files\020600050203_ChoptankRiver\01_02_03_04_utm18_breach_bankpixels.tif'
#    
#    ## << BANK POINTS >>    
##    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpoints.shp'
###    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_bankpoints_USGS_AEA.shp'
##    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_bankpts.shp'
##    str_bankpts_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpts_TEST.shp'
##    str_bankpts_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\headwaters\bankpts.shp'
# 
#    ## << FLOODPLAIN >>
##    str_floodplain_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/fp_pixels.tif'
##    str_floodplain_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM_breach_hand_slice2.3.tif'    
#    
#    ## << CROSS SECTIONS >>    
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_xns_USGS_AEA_2.shp'
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net_xnptdist1.shp'
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_xns.shp'
##    str_xns_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\chan_xns_TEST.shp'
##    str_xns_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\fp_xns.shp'
##    str_xns_path = r'C:\USGS_TerrainBathymetry_Project\CBP_analysis\DifficultRun\headwaters\xns.shp'
#
#    ## << HAND (dd) >>    
##    str_hand_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_hand.tif'
##    str_hand_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMhand.tif'
##    str_hand_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\raw\dr3m_raw_dem_clip_utm18_breach_hand.tif'
##    str_hand_path = r'D:\CFN_data\DEM_Files\020600050203_ChoptankRiver\01_02_03_04_utm18_breach_hand.tif'
#    
#    ## << SHEDS >>
##    str_sheds_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/dr3m_sheds_diss2.shp'   
#    
#    ## << DANGLE POINTS >>
##    str_startptgrid_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMnet_UNIQUE_ID.shp'
##    str_startptgrid_path= r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\breach\dr_nhd_hires_dangles.tif'
##    str_danglepts_path=r'D:\fernando\HAND\020802\sam_test\020802_dangles.tif'
#    
#    ## << OPENNESS >>  THIS NEEDS WORK -- Table this til later if you have time (July5)
##    str_pos_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\dr_pos_raw.tif'    
##    str_pos_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\openness\Openness\diffphix_pos' 
##    str_neg_path = r'D:\Terrain_and_Bathymetry\USGS\CBP_analysis\DifficultRun\facet_tests\dr_neg_raw.tif'    
#    
#    str_reachid='LINKNO'
##    str_reachid='ARCID'
##    str_reachid='COMID'
#    str_orderid='strmOrder'
##    cell_size=3
#    
##        self.p_fitlength.setText('3')
##        self.p_xngap.setText('3')
##        self.p_xnlength.setText('30')
##        self.p_ivert.setText('0.2')
##        self.p_ratiothresh.setText('1.5')
##        self.p_slpthresh.setText('0.03')    
#    
#    # << PARAMETERS >>
#    parm_ivert = 0.2
#    XnPtDist = 3
#    parm_ratiothresh = 1.5
#    parm_slpthresh = 0.03
#    p_fpxnlen = 100 # meters (if UTM) ??
#    
#    # << Get cell size directly from grid >>
##    cell_size = int(funcs_v2.get_cell_size(str_dem_path)) # range functions need int?
#  
###    # << Build Strealine Coordinates >>
###    # Build reach coords and get crs from a pre-existing streamline shapefile...
##    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_path, cell_size, str_reachid, str_orderid)
##    df_coords.to_csv('df_coords_Chillisquaque.csv') # just for testing
##    df_coords.to_csv('df_coords_DifficultRun.csv') # just for testing
##    df_coords.to_csv('df_coords_020802.csv') # just for testing
#    
##    print('NOTE:  Reading pre-calculated csv file...')
##    df_coords = pd.read_csv('df_coords_DifficultRun.csv')
##    df_coords = pd.read_csv('df_coords_Chillisquaque.csv', )
##    df_coords = pd.read_csv('df_coords_020802.csv', )    
##    streamlines_crs = {'init': u'epsg:26918'}
##   
##    # << Find bank pixels using moving window along streamline >>
##    funcs_v2.bankpixels_from_curvature_window(df_coords, str_dem_path, str_bankpixels_path, cell_size) 
##    funcs_v2.bankpixels_from_openness_window(df_coords, str_pos_path, str_bankpixels_path) 
##    funcs_v2.bankpixels_from_openness_window_buffer_all(df_coords, str_dem_path, str_net_path, str_pos_path, str_neg_path) 
#    
###    # << Channel Width via Bank Pixels >>    
##    funcs_v2.channel_width_bankpixels(str_net_path, str_bankpixels_path, str_reachid, cell_size)    
##    funcs_v2.channel_width_bankpixels_segments(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size)
##    funcs_v2.channel_and_fp_width_bankpixels_segments_po(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size, p_fpxnlen, str_hand_path)
#    
##     << FLOODPLAIN WIDTH >> 
##    buff_dist = 40
##    funcs_v2.floodplain_width_2D_xns(str_xns_path, str_floodplain_path, buff_dist)
##    funcs_v2.floodplain_width_fppixels_segments_po(df_coords, str_net_in_path, str_floodplain_path, str_reachid, cell_size)
##    funcs_v2.floodplain_width_reach_buffers_po(funcs, str_net_path, str_fp_path, str_reachid, cell_size)
#
##      << Find bank pixels using moving window along streamline >>
##    funcs_v2.fp_from_streamline_window(df_coords, str_dem_path, str_fp_path)
##       
##      << Analyze hand using moving window along streamline >>  WAY overestimates (this is more of a terrace finder)
##    funcs_v2.workerfuncs.analyze_hand_from_streamline_window(funcs, df_coords, str_hand_path)   ### USE THIS ONE
##    
##      << Analyze the DEM using 2D (buffered) cross-sections >>  ## OTHER TESTs
##    funcs_v2.analyze_hand_2D_xns(str_xns_path, str_hand_path, parm_ivert) # just use DEM here?   
##    funcs_v2.analyze_hand_reach_buffers(str_net_path, str_hand_path, str_reachid) 
##    funcs_v2.analyze_hand_sheds(str_sheds_path, str_hand_path, parm_ivert)
#    
#    # << CREATE Xn SHAPEFILES >>
#    # Channel...
##    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(3), int(3), float(30))     
#    # FP...
##    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), True, int(30), int(30), float(100))  # For FP width testing
#
##    # << INTERPOLATE XNs >>
##    df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path)
#    
##    print('Writing df_xn_elev to .csv for testing...')
##    df_xn_elev.to_csv(columns=['index','linkno','elev','xn_row','xn_col']) 
##    df_xn_elev2 = pd.read_csv('df_xn_elev.csv') #, dtype={'linko':np.int,'elev':np.float,'xn_row':np.float,'xn_col':np.float})
# 
#    # Calculate channel metrics and write bank point shapefile...
##    print('Calculating channel metrics from bank points...')
##    funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
# 
##    # << Run TauDEM HAND/Dinf Distance Down tool >>    
##    str_directory = r'D:\CFN_data\DEM_Files\*'
##    lst_subdirs = glob(str_directory+'/')
##    
##    for subdir in lst_subdirs:
##        
##        str_fel_path = subdir + 'DEMfel.tif'
##        str_src_path = subdir + 'DEMsrc.tif'
##        str_ang_path = subdir + 'DEMang.tif'
##        str_hand_path = subdir + 'DEMdd.tif'
##        str_slp_path = subdir + 'DEMslp.tif'
##        
##        print('Processing:  {}'.format(subdir))
##        
##        funcs_v2.workerfuncs.taudem_dinfdd(funcs, str_fel_path, str_src_path, str_ang_path, str_hand_path, str_slp_path)
#     
#     # << DEM Pre-processing using TauDEM and GoSpatial >>              
#    # (1) Clip original streamlines layer (NHD hi-res 4 digit HUC to DEM of interest)...     
#    # Build the output streamlines file name...
##    path_to_dem, dem_filename = os.path.split(str_dem_path)
##    str_output_nhdhires_path = path_to_dem + '\\' + dem_filename[:-4]+'_nhdhires.shp' 
#    
##    funcs_v2.clip_features(str_net_in_path, str_output_nhdhires_path, str_dem_path)     
#    
##    # (2) Do all Whitebox and TauDEM functions including HAND based on output r'"C:\Program Files\TauDEM\TauDEM5Exe\D8FlowDir.exe"'
#    str_mpi_path=r'"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"'
#    str_taudem_dir=r'"C:\Program Files\TauDEM\TauDEM5Exe' #\D8FlowDir.exe"'
#    str_whitebox_path= r'C:\Terrain_and_Bathymetry\Whitebox\GoSpatial\go-spatial_win_amd64.exe'
#    
#    
#    # def preprocess_dem(str_dem_path, str_streamlines_path, str_mpi_path, str_taudem_path, str_whitebox_path, run_wg, run_whitebox, run_taudem):
#    
#    funcs_v2.preprocess_dem(str_dem_path, str_net_in_path, str_mpi_path, str_taudem_dir, str_whitebox_path, True, True, True)
#    
#
#    
#    print('\n<<< End >>>\r\n')
#    print('Run time:  {}'.format(timeit.default_timer() - start_time_0))
##     =====================================================================================

    # ===== Run GUI ==============================
    app = QtGui.QApplication(sys.argv)

    mainform = facet()
    
    # Show it in the center of the window...
    fr_geo=mainform.frameGeometry()
    cp = QtGui.QDesktopWidget().availableGeometry().center()
    fr_geo.moveCenter(cp)
    mainform.move(fr_geo.topLeft())    
    
    mainform.show()

    sys.exit(app.exec_())
#    # ============================================