# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:11:00 2016

@author: sam.lamont
"""

import timeit
from PyQt4 import QtGui, uic
import sys
#import numpy as np
#from numpy import array
#from scipy import ndimage
#import os
import rasterio
#from rasterio.plot import show
#from rasterio import features
#from affine import Affine # Georeferenced raster datasets use affine transformations to map from image coordinates to world coordinates
#from shapely.geometry import LineString
import fiona
#from fiona.crs import to_string
#from math import isinf, sqrt #, hypot, modf, atan, ceil
#import matplotlib.pyplot as plt
#import pandas as pd

import funcs_v2

#from shapely import speedups
#
#if speedups.available:
#    speedups.enable() # enable performance enhancements written in C

#plt.ion()  

# ===================================================================================
#                                       MAIN
# ===================================================================================
class facet(QtGui.QDialog):
#class facet(QtGui.QTabWidget):
#class MainForm(QtGui.QDialog, FORM_CLASS):
    
    def __init__(self, parent=None):        
        super(facet, self).__init__(parent)
        
        # ================= UI ================
        uic.loadUi('facet.ui', self)  
        
        # << Browse buttons >>
        self.browseDEM_Xns.clicked.connect(self.get_dem_file)             # DEM - Xns
        self.browseStreams_Xns.clicked.connect(self.get_streams_file)     # Streams - Xns
        self.browseXns.clicked.connect(self.get_or_set_xns_file)          # Xns
        self.browseBankPts.clicked.connect(self.set_bankpts_file)         # BankPts
        self.browseDEM_Curve.clicked.connect(self.get_dem_file)           # DEM - curvature
        self.browseStreams_Curve.clicked.connect(self.get_streams_file)   # Streams - curvature
        self.browseBankPixels.clicked.connect(self.get_or_set_bankpixels) # Bank Pixels      
        
        # << Button box >>
        self.buttonBox.accepted.connect(self.run_main)        
        self.buttonBox.rejected.connect(self.close)
#        self.pushButton.clicked.connect(self.run_main)
        
        # << Defaults >>
        self.textFitLength.setText('3')
        self.textXnSpacing.setText('3')
        self.textXnLength.setText('30')
        self.textVertIncrement.setText('0.2')
        self.textRatioThresh.setText('1.5')
        self.textSlpThresh.setText('0.03')        
        
        # << Progress Bar >>
#        self.progressBar.setVisible(False)
#        self.progressBar.setTextVisible(False)
        
        # =====================================
        
    # =========================================    
    #              << FILE I/O >>    
    # =========================================        
    def get_dem_file(self):
        fname_dem = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
        
        if self.tabWidget.currentIndex() == 0:        
            self.textDEM_xns.setText(fname_dem)

        if self.tabWidget.currentIndex() == 1:        
            self.textDEM_curvature.setText(fname_dem)
            
#        if self.tab_widget.currentIndex() == 2:        
#            self.demPath_hand.setText(fname_dem)            
        
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
            self.textStreams_Xns.setText(fname_streams) 
            
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.comboBoxXns.addItem(str(name))            

        if self.tabWidget.currentIndex() == 1:        
            self.textStreams_curvature.setText(fname_streams) 
        
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.comboBoxCurvature.addItem(str(name))

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
        if self.chkPixelMetrics.isChecked():
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
            self.textBankPixels.setText(fname)   
        else:
            fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Grids (*.*)')
            self.textBankPixels.setText(fname)             
    
    # =========================================    
    #              << RUN MAIN >>    
    # =========================================
    def run_main(self):        
        
#        self.progDialog = QtGui.QProgressDialog("", "Cancel", 0, len(self.rename_list), self)
               
        start_time_0 = timeit.default_timer()
#        self.close() # keep open to display progress bar?
#        self.show()
        
#         Check which tab is active??
#        i_tab = self.tab_widget.currentIndex() # (0: bankpts; 1: bankpixels; 2: FP-HAND)
#        print('selected tab:  {}'.format(i_tab))
        
        # Check coordinate systems (DEM vs. Streamlines)...
        if dict(self.dem_crs) <> dict(self.streamlines_crs):
            QtGui.QMessageBox.critical(self, 'Warning!', 'Coordinate Systems of DEM and Streams Layer May Not Match')
#            sys.exit()
            
        str_orderid='Order_'
            
        # << Generate Cross Sections >>
        if self.chkCreateXns.isChecked():
            
#            print('here')
                   
            # Build reach coords and get crs from a pre-existing streamline shapefile...
            df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(self.textStreams_Xns.text(), self.cell_size, str(self.comboBoxXns.currentText()), str_orderid, self)

            # Build and write the Xn shapefile...
            funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(self.textXns.text()), False, int(self.textXnSpacing.text()), int(self.textFitLength.text()), float(self.textXnLength.text()), self)     
        
        # << Calculate Channel Metrics -- Xns Method >>
        if self.chkMetricsXns.isChecked():

            # Read the cross section shapefile and interplote elevation along each one using the specified DEM
            df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(self.textXns.text(), self.textDEM_xns.text(), self)
            
            XnPtDist = self.cell_size
            
            # Calculate channel metrics and write bank point shapefile...# parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh
            funcs_v2.chanmetrics_bankpts(df_xn_elev, str(self.textXns.text()), str(self.textDEM_xns.text()), str(self.textBankPts.text()), float(self.textVertIncrement.text()), XnPtDist, float(self.textRatioThresh.text()), float(self.textSlpThresh.text()), self)
                   
        # << Create bank pixel layer (.tif) >>
        if self.chkBankPixels.isChecked():
             
            print('Reading streamline coordinates...') 
            df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(self.textStreams_curvature.text(), self.cell_size, str(self.comboBoxCurvature.currentText()))
            
            print('Creating bank pixel layer...') 
            funcs_v2.bankpixels_from_streamline_window(df_coords, str(self.textDEM_curvature.text()), str(self.textBankPixels.text()))
            
        # << Calculate Channel Metrics -- Bank Pixels Method >>
        if self.chkPixelMetrics.isChecked():
            print('Calculating channel width from bank pixel layer!')
            # def channel_width_bankpixels(str_streamlines_path, str_bankpixels_path, str_reachid):                        
            funcs_v2.channel_width_bankpixels(self.textStreams_curvature.text(), self.textBankPixels.text(), self.comboBoxCurvature .currentText(), self.cell_size)
        
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
                
        # Show message when complete...
        # '\tTotal time:  ' + str(timeit.default_timer() - start_time_0
        QtGui.QMessageBox.information(self, 'Code Complete', 'Successful Run! \n\n\t   {0:.2f} secs'.format(timeit.default_timer() - start_time_0))
                 
if __name__ == '__main__':    
    
###    # ============================ << FOR DEBUGGING >> ==================================
#    print('\n<<< Start >>>\r\n')
#    start_time_0 = timeit.default_timer() 
###    
#    ## << DEM >>
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_dem.tif'
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700dem_UTM18N_drclip_3mresample.tif'
###    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM_02050205_USGS_AEA.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMfel_dem_post.tif'
#    str_dem_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMfel.tif'
#    
#    ## << STREAMLINES >>  
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/Network04.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net.shp'
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700_flowsUTM18N_drsheds_clip_diss.shp'
###    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_StreamNetwork_USGS_AEA_2.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMnet1000.shp'
#    str_net_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEMnet_UNIQUE_ID.shp'
#
#    ## << BANK PIXELS >>   
###    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/bank_pixels_02050205_H.tif'
##    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpixels.tif'
#    str_bankpixels_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpixels.tif'
#    
#    ## << BANK POINTS >>    
##    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpoints.shp'
###    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_bankpoints_USGS_AEA.shp'
##    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_bankpts.shp'
#    str_bankpts_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\bankpts.shp'
# 
#    ## << FLOODPLAIN >>
#    str_fp_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/fp_pixels.tif'
#    
#    ## << CROSS SECTIONS >>    
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_xns_USGS_AEA_2.shp'
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net_xnptdist1.shp'
##    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_xns.shp'
#    str_xns_path = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\chan_xns.shp'
#
#    ## << HAND >>    
#    str_hand_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_hand.tif'
#    
#    ## << SHEDS >>
#    str_sheds_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/dr3m_sheds_diss2.shp'    
#    
#    str_reachid='LINKNO'
##    str_reachid='ARCID'
#    str_orderid='Order_'
#    cell_size=3
#    
##        self.p_fitlength.setText('3')
##        self.p_xngap.setText('3')
##        self.p_xnlength.setText('30')
##        self.p_ivert.setText('0.2')
##        self.p_ratiothresh.setText('1.5')
##        self.p_slpthresh.setText('0.03')    
#    
#    parm_ivert = 0.2
#    XnPtDist = 3
#    parm_ratiothresh = 1.5
#    parm_slpthresh = 0.03
#  
###    # << Build Strealine Coordinates >>
###    # Build reach coords and get crs from a pre-existing streamline shapefile...
##    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_path, cell_size, str_reachid, str_orderid)
##    df_coords.to_csv('df_coords_Chillisquaque.csv') # just for testing
#    
#    print('Reading pre-calculated csv file...')
#    df_coords = pd.read_csv('df_coords_Chillisquaque.csv', )    
##   
##        << Find bank pixels using moving window along streamline >>
##    funcs_v2.bankpixels_from_streamline_window(df_coords, str_dem_path, str_bankpixels_path) 
#    
###    # << Channel Width via Bank Points >>    
###    buff_dist = 10.0
##    funcs_v2.channel_width_bankpixels(str_net_path, str_bankpixels_path, str_reachid, cell_size)    
#    funcs_v2.channel_width_bankpixels_segments(df_coords, str_net_path, str_bankpixels_path, str_reachid, cell_size)
##
#
##
##      << Find bank pixels using moving window along streamline >>
##    funcs_v2.fp_from_streamline_window(df_coords, str_dem_path, str_fp_path)
##       
##      << Analyze hand using moving window along streamline >>
##    funcs_v2.analyze_hand_from_streamline_window(df_coords, str_hand_path)   ### USE THIS ONE
##    
##      << Analyze the DEM using 2D (buffered) cross-sections >>  ## OTHER TESTs
##    funcs_v2.analyze_hand_2D_xns(str_xns_path, str_hand_path, parm_ivert) # just use DEM here?   
##    funcs_v2.analyze_hand_reach_buffers(str_net_path, str_hand_path, str_reachid) 
##    funcs_v2.analyze_hand_sheds(str_sheds_path, str_hand_path, parm_ivert)
#    
#    # Build and write the Xn shapefile...
##    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(3), int(3), float(30))     
#
##    # << INTERPOLATE XNs >>
##    df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path)
#    
##    print('Writing df_xn_elev to .csv for testing...')
##    df_xn_elev.to_csv(columns=['index','linkno','elev','xn_row','xn_col']) 
##    df_xn_elev2 = pd.read_csv('df_xn_elev.csv') #, dtype={'linko':np.int,'elev':np.float,'xn_row':np.float,'xn_col':np.float})
# 
##    # Calculate channel metrics and write bank point shapefile...
##    print('Calculating channel metrics from bank points...')
##    funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
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
    # ============================================