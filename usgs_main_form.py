# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:11:00 2016

@author: sam.lamont
"""

import timeit
from PyQt4 import QtGui, QtCore
import sys
import numpy as np
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
import pandas as pd

import funcs_v2

#from shapely import speedups
#
#if speedups.available:
#    speedups.enable() # enable performance enhancements written in C

#plt.ion()  
# ===================================================================================
#                                       MAIN
# ===================================================================================
class MainForm(QtGui.QWidget):
    
    def __init__(self, parent=None):
        super(MainForm, self).__init__(parent)
        
        # << TABS >>
        # Establish the main tab widget...
        self.tab_widget = QtGui.QTabWidget()
#        self.tab_widget.setStyleSheet("QTabBar::tab:disabled {"+\
#                        "width: 300px;"+\
#                        "color: transparent;"+\
#                        "background: transparent;}")
#        
#        self.tab_widget.blockSignals(True)
#        self.tab_widget.currentChanged.connect(self.tab_click)
#        self.tab_widget.currentChanged.connect(self.ck_chanmetrics_click)

        # Create each tab...        
        self.tab_bankpts = QtGui.QWidget()
        self.tab_bankpixels = QtGui.QWidget()
        self.tab_floodplain = QtGui.QWidget()
        
        # ========================================================================
        # ========================== BANK POINTS LAYOUT ==========================
        # ========================================================================
        # << Define text labels, line edits, buttons >>
        self.ck_createXns = QtGui.QCheckBox('Create Xns')
        self.ck_chanmetrics_bankpts = QtGui.QCheckBox('Calculate Channel Metrics')
        self.ck_runall_bankpts = QtGui.QCheckBox('All')        
        
        # DEM...
        self.lbl_dem = QtGui.QLabel('DEM')
        self.demPath_bankpts = QtGui.QLineEdit()
        self.btn_dem = QtGui.QPushButton('Browse')
        self.btn_dem.clicked.connect(self.get_dem_file)
#        self.btn_dem.setFixedSize(120, 25)
        
        # Streams...
        self.lbl_streams = QtGui.QLabel('Streamlines')
        self.streamsPath_bankpts = QtGui.QLineEdit()
        self.btn_streams = QtGui.QPushButton('Browse')
        self.btn_streams.clicked.connect(self.get_streams_file)
        
        # Reach ID combo box...
        self.lbl_reachid = QtGui.QLabel('Reach ID')
        self.reachid_bankpts = QtGui.QComboBox()

        # Xns...
        self.lbl_xns = QtGui.QLabel('Cross Sections')
        self.xnsPath = QtGui.QLineEdit()
        self.btn_xns = QtGui.QPushButton('Browse')
        self.btn_xns.clicked.connect(self.get_or_set_xns_file) 
        
        # Bank Points...
        self.lbl_bankpts = QtGui.QLabel('Bank Points')
        self.bankptsPath = QtGui.QLineEdit()
        self.btn_bankpts = QtGui.QPushButton('Browse')
        self.btn_bankpts.clicked.connect(self.set_bankpts_file)
        
        # Run...        
        self.btn_run = QtGui.QPushButton('Run')
        self.btn_run.clicked.connect(self.run_main) 
        self.btn_run.setFixedSize(80, 25)

        # Cancel...        
        self.btn_cancel = QtGui.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(QtCore.QCoreApplication.instance().quit) 
        self.btn_cancel.setFixedSize(80, 25)
               
        # Labels...
        self.lbl_p_fitlength = QtGui.QLabel('Fit Length')
        self.lbl_p_xngap = QtGui.QLabel('Cross Section Spacing')
        self.lbl_p_xnlength = QtGui.QLabel('Cross Section Length')
        self.lbl_p_ivert = QtGui.QLabel('Vertical Increment')
        self.lbl_p_ratiothresh = QtGui.QLabel('Ratio Threshold')
        self.lbl_p_slpthresh = QtGui.QLabel('Slope Threshold')
        
        self.lbl_fit_units = QtGui.QLabel(' ')
        self.lbl_xngap_units = QtGui.QLabel(' ')
        self.lbl_xnlength_units = QtGui.QLabel(' ')
        self.lbl_ivert_units = QtGui.QLabel(' ')
        self.lbl_ratiothresh_units = QtGui.QLabel(' ')
        self.lbl_slpthresh_units = QtGui.QLabel(' ')
        
        # Right justify labels...        
        self.lbl_reachid.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_fitlength.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_xngap.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_xnlength.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_streams.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_xns.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_fitlength.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_dem.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_ivert.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_ratiothresh.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_p_slpthresh.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_bankpts.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        
        # Parameter value input text boxes...
        self.p_fitlength = QtGui.QLineEdit()
        self.p_xngap = QtGui.QLineEdit()
        self.p_xnlength = QtGui.QLineEdit()    
        self.p_ivert = QtGui.QLineEdit()
        self.p_ratiothresh = QtGui.QLineEdit()
        self.p_slpthresh = QtGui.QLineEdit()
        
        self.p_fitlength.setText('3')
        self.p_xngap.setText('3')
        self.p_xnlength.setText('30')
        self.p_ivert.setText('0.2')
        self.p_ratiothresh.setText('1.5')
        self.p_slpthresh.setText('0.03')
        
        # << Form Layout >>
        gridLayout = QtGui.QGridLayout(self.tab_bankpts)
#        gridLayout.setSpacing(1)
        gridLayout.setColumnStretch(0, 1)
        gridLayout.setColumnStretch(1, 1)        
        gridLayout.setColumnStretch(2, 1)
        gridLayout.setColumnStretch(3, 1) 
#        gridLayout.setColumnStretch(4, 1)
        
        # Groupbox with check boxes for selecting program control...
        gpBox = QtGui.QGroupBox('Program Controls')
#        gpBox.setFlat(True)
        
        vBox = QtGui.QVBoxLayout()
#        self.ck_createXns = QtGui.QCheckBox('Create Xns')
#        self.ck_chanmetrics_bankpts = QtGui.QCheckBox('Calculate Channel Metrics')
#        self.ck_runall_bankpts = QtGui.QCheckBox('All')
                
        vBox.addWidget(self.ck_createXns)
        vBox.addWidget(self.ck_chanmetrics_bankpts)
        vBox.addWidget(self.ck_runall_bankpts)
#        vBox.addStretch(1)
        
        gpBox.setLayout(vBox)
        
        row_gpbox = 0
        row_dem = 1
        row_streams = 2
        row_reachid = 3
        row_fitlength = 4
        row_xngap = 5
        row_xnlen = 6
        row_ivert = 7
        row_ratthresh = 8
        row_slpthresh = 9
        row_xnlyr = 10
        row_bankpts = 11
        row_runcnl = 12
                
        gridLayout.addWidget(gpBox, row_gpbox, 1, 1, 1)
        
        # Input DEM...
        gridLayout.addWidget(self.lbl_dem, row_dem, 0, 1, 1)
        gridLayout.addWidget(self.demPath_bankpts, row_dem, 1, 1, 2)
        gridLayout.addWidget(self.btn_dem, row_dem, 3, 1, 1)
                  
        # Input Streamlines...
        gridLayout.addWidget(self.lbl_streams, row_streams, 0, 1, 1)
        gridLayout.addWidget(self.streamsPath_bankpts, row_streams, 1, 1, 2)
        gridLayout.addWidget(self.btn_streams, row_streams, 3, 1, 1)
        
        # Reach ID...
        gridLayout.addWidget(self.lbl_reachid, row_reachid, 0, 1, 1)
        gridLayout.addWidget(self.reachid_bankpts, row_reachid, 1, 1, 1)
        
        # << Input Parameters >>
        # Fit length...
        gridLayout.addWidget(self.lbl_p_fitlength, row_fitlength, 0, 1, 1)
        gridLayout.addWidget(self.p_fitlength, row_fitlength, 1, 1, 1) 
        gridLayout.addWidget(self.lbl_fit_units, row_fitlength, 2, 1, 1)
        # Xn gap...
        gridLayout.addWidget(self.lbl_p_xngap, row_xngap, 0, 1, 1)
        gridLayout.addWidget(self.p_xngap, row_xngap, 1, 1, 1) 
        gridLayout.addWidget(self.lbl_xngap_units, row_xngap, 2, 1, 1) 
        # Xn length...
        gridLayout.addWidget(self.lbl_p_xnlength, row_xnlen, 0, 1, 1)
        gridLayout.addWidget(self.p_xnlength, row_xnlen, 1, 1, 1)     
        gridLayout.addWidget(self.lbl_xnlength_units, row_xnlen, 2, 1, 1)
        # IVert...
        gridLayout.addWidget(self.lbl_p_ivert, row_ivert, 0, 1, 1)
        gridLayout.addWidget(self.p_ivert, row_ivert, 1, 1, 1)
        gridLayout.addWidget(self.lbl_ivert_units, row_ivert, 2, 1, 1)
        # Ratio Thresh...
        gridLayout.addWidget(self.lbl_p_ratiothresh, row_ratthresh, 0, 1, 1)
        gridLayout.addWidget(self.p_ratiothresh, row_ratthresh, 1, 1, 1)
        gridLayout.addWidget(self.lbl_ratiothresh_units, row_ratthresh, 2, 1, 1)
        # Slope Thresh...
        gridLayout.addWidget(self.lbl_p_slpthresh, row_slpthresh, 0, 1, 1)
        gridLayout.addWidget(self.p_slpthresh, row_slpthresh, 1, 1, 1)
        gridLayout.addWidget(self.lbl_slpthresh_units, row_slpthresh, 2, 1, 1)
        
        # Cross Sections...
        gridLayout.addWidget(self.lbl_xns, row_xnlyr, 0, 1, 1)
        gridLayout.addWidget(self.xnsPath, row_xnlyr, 1, 1, 2)
        gridLayout.addWidget(self.btn_xns, row_xnlyr, 3, 1, 1)
        
        # Bank Points...
        gridLayout.addWidget(self.lbl_bankpts, row_bankpts, 0, 1, 1)
        gridLayout.addWidget(self.bankptsPath, row_bankpts, 1, 1, 2)
        gridLayout.addWidget(self.btn_bankpts, row_bankpts, 3, 1, 1)
               
        # Run and Cancel...
        gridLayout.addWidget(self.btn_run, row_runcnl, 1, 1, 1)
        gridLayout.addWidget(self.btn_cancel, row_runcnl, 2, 1, 1)        
        
        
        # ========================================================================
        # =========================== BANK PIXELS LAYOUT =========================
        # ========================================================================
        # << Define text labels, line edits, buttons >>
        int_pathwidth = 300

        self.ck_createPixels = QtGui.QCheckBox('Create Bank Pixel Layer')
        self.ck_chanmetrics_pixels = QtGui.QCheckBox('Calculate Channel Metrics')
        self.ck_runall_pixels = QtGui.QCheckBox('All')
        
        # DEM...
        self.lbl_dem = QtGui.QLabel('DEM')
        self.demPath_bankpixels = QtGui.QLineEdit()
        self.demPath_bankpixels.setFixedSize(int_pathwidth, 20)
        self.btn_dem = QtGui.QPushButton('Browse')
        self.btn_dem.clicked.connect(self.get_dem_file)
        self.btn_dem.setFixedSize(80, 25)
        
        # Streams...
        self.lbl_streams = QtGui.QLabel('Streamlines')
        self.streamsPath_bankpixels = QtGui.QLineEdit()
        self.streamsPath_bankpixels.setFixedSize(int_pathwidth, 20)
        self.btn_streams = QtGui.QPushButton('Browse')
        self.btn_streams.clicked.connect(self.get_streams_file)
        self.btn_streams.setFixedSize(80, 25)
        
        # Reach ID combo box...
        self.lbl_reachid = QtGui.QLabel('Reach ID')
        self.reachid_bankpixels = QtGui.QComboBox()    
        
        # Bank Pixels layer...
        self.lbl_bankpixels = QtGui.QLabel('Bank Pixels')
        self.bankpixelsPath = QtGui.QLineEdit()
        self.bankpixelsPath.setFixedSize(int_pathwidth, 20)
        self.btn_bankpixels = QtGui.QPushButton('Browse')
        self.btn_bankpixels.clicked.connect(self.get_or_set_bankpixels)
        
        # Run...        
        self.btn_run = QtGui.QPushButton('Run')
        self.btn_run.clicked.connect(self.run_main) 
        self.btn_run.setFixedSize(80, 25)

        # Cancel...        
        self.btn_cancel = QtGui.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(QtCore.QCoreApplication.instance().quit) 
        self.btn_cancel.setFixedSize(80, 25)
        
        # Right justify labels...        
        self.lbl_reachid.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_streams.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_dem.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl_bankpixels.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter) 
        
        # Spacer item...
        spacerItem = QtGui.QSpacerItem(int_pathwidth, 275, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        
        # << Main VBox >>
        gridLayout = QtGui.QGridLayout(self.tab_bankpixels)
        gridLayout.setColumnStretch(0, 1)
        gridLayout.setColumnStretch(1, 1)        
        gridLayout.setColumnStretch(2, 1)
        gridLayout.setColumnStretch(3, 1) 
       
        # Groupbox with check boxes for selecting program control...
        gpBox = QtGui.QGroupBox('Program Controls')
        
#        gpBox.setMaximumHeight(150)
#        gpBox.setMinimumWidth(275)
       
        vBox = QtGui.QVBoxLayout()
               
        vBox.addWidget(self.ck_createPixels)
        vBox.addWidget(self.ck_chanmetrics_pixels)
        vBox.addWidget(self.ck_runall_pixels)
#        vBox.addSpacing(100)
        
        gpBox.setLayout(vBox)
        
        row_gpbox = 0
        row_spacer = 1
        row_dem = 2
        row_streams = 3
        row_reachid = 4
        row_bankpixels = 5
        row_runcnl = 6
                
        # Group box...                
        gridLayout.addWidget(gpBox, row_gpbox, 1, 1, 1)
        
        gridLayout.addItem(spacerItem, row_spacer, 1, 1, 2)

                
        # Input DEM...
        gridLayout.addWidget(self.lbl_dem, row_dem, 0, 1, 1)
        gridLayout.addWidget(self.demPath_bankpixels, row_dem, 1, 1, 2)
        gridLayout.addWidget(self.btn_dem, row_dem, 3, 1, 1)
                  
        # Input Streamlines...
        gridLayout.addWidget(self.lbl_streams, row_streams, 0, 1, 1)
        gridLayout.addWidget(self.streamsPath_bankpixels, row_streams, 1, 1, 2)
        gridLayout.addWidget(self.btn_streams, row_streams, 3, 1, 1)
        
        # Reach ID...
        gridLayout.addWidget(self.lbl_reachid, row_reachid, 0, 1, 1)
        gridLayout.addWidget(self.reachid_bankpixels, row_reachid, 1, 1, 1)
           
        # Bank Points...
        gridLayout.addWidget(self.lbl_bankpixels, row_bankpixels, 0, 1, 1)
        gridLayout.addWidget(self.bankpixelsPath, row_bankpixels, 1, 1, 2)
        gridLayout.addWidget(self.btn_bankpixels, row_bankpixels, 3, 1, 1)
               
        # Run and Cancel...
        gridLayout.addWidget(self.btn_run, row_runcnl, 1, 1, 1)
        gridLayout.addWidget(self.btn_cancel, row_runcnl, 2, 1, 1)
        # ========================================================================
        # ======================= END LAYOUTS ====================================
        # ========================================================================
        
        self.ck_runall_bankpts.setEnabled(False)
        self.ck_runall_pixels.setEnabled(False)
        
        # Format the layout of each tab...
#        self.tab_bankpts_layout()
#        self.tab_bankpixels_layout()       
        
        # Connect the tabs to the widget...
        self.tab_widget.addTab(self.tab_bankpts, 'Banks-Xn\'s')
        self.tab_widget.addTab(self.tab_bankpixels, 'Banks-Gridded')
        self.tab_widget.addTab(self.tab_floodplain, 'Floodplain-HAND')

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.tab_widget)
        
#        self.tab_widget.blockSignals(False)

        self.setLayout(mainLayout)
        self.setWindowTitle('Bank Detection and FP Analysis')
#        self.tab_widget.pos()
#        self.tab_widget.setFixedSize(550, 500)
        self.setFixedSize(570, 550)
        self.show()

#        if self.ck_chanmetrics_bankpts.isChecked():
#            self.btn_xns.clicked.connect(self.get_xns_file)
#        else:
#            self.btn_xns.clicked.connect(self.set_xns_file)
#            
#        if self.ck_chanmetrics_pixels.isChecked():
#            self.btn_bankpixels.clicked.connect(self.get_bankpixels_file)
#        else:
#            self.btn_bankpixels.clicked.connect(self.set_bankpixels_file)            
#        
#        if self.tab_widget.currentIndex()==1:
#            self.tab_bankpixels_layout()
#            print('hey!')
            
#    def tab_click(self):
#        self.tab_bankpts.layout
##        if self.tab_widget.currentIndex() == 0:
##            self.tab_bankpts_layout()
#            
#        if self.tab_widget.currentIndex() == 1:
#            self.tab_bankpixels_layout()
            
#        print('hey!')                     
    # =========================================    
    #              << FILE I/O >>    
    # =========================================        
    def get_dem_file(self):
        fname_dem = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
        
        if self.tab_widget.currentIndex() == 0:        
            self.demPath_bankpts.setText(fname_dem)

        if self.tab_widget.currentIndex() == 1:        
            self.demPath_bankpixels.setText(fname_dem)
            
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

        if self.tab_widget.currentIndex() == 0:        
            self.streamsPath_bankpts.setText(fname_streams) 
            
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
    #            print(streamlines.schema['properties'].keys())
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.reachid_bankpts.addItem(str(name))            

        if self.tab_widget.currentIndex() == 1:        
            self.streamsPath_bankpixels.setText(fname_streams) 
        
            # Also load the reach ID combo box here...
            with fiona.open(str(fname_streams), 'r') as streamlines:
    #            print(streamlines.schema['properties'].keys())
                lst_fieldnames = streamlines.schema['properties'].keys()      
                self.streamlines_crs = streamlines.crs
                print('Streamlines crs: {}'.format(streamlines.crs))
            if lst_fieldnames:
                for name in lst_fieldnames:
                    self.reachid_bankpixels.addItem(str(name))

    def set_bankpts_file(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
        self.bankptsPath.setText(fname)  
        
#    def get_bankpts_file(self):
#        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
#        self.bankptsPath.setText(fname)  
        
    def get_or_set_xns_file(self):
        if self.ck_createXns.isChecked():
            fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
            self.xnsPath.setText(fname)
        else:
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Shapefiles (*.shp)')
            self.xnsPath.setText(fname)  
        
    def get_or_set_bankpixels(self):
        if self.ck_chanmetrics_pixels.isChecked():
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')
            self.bankpixelsPath.setText(fname)   
        else:
            fname = QtGui.QFileDialog.getSaveFileName(self, 'Browse', '', 'Grids (*.*)')
            self.bankpixelsPath.setText(fname)             
    
    # =========================================    
    #              << RUN MAIN >>    
    # =========================================
    def run_main(self):
               
        start_time_0 = timeit.default_timer()
        self.close()
        
#         Check which tab is active??
#        i_tab = self.tab_widget.currentIndex() # (0: bankpts; 1: bankpixels; 2: FP-HAND)
#        print('selected tab:  {}'.format(i_tab))
        
        # Check coordinate systems (DEM vs. Streamlines)...
        if dict(self.dem_crs) <> dict(self.streamlines_crs):
            QtGui.QMessageBox.critical(self, 'Warning!', 'Coordinate Systems of DEM and Streams Layer May Not Match')
#            sys.exit()
            
        # << Generate Cross Sections >>
        if self.ck_createXns.isChecked():
                   
            # Build reach coords and get crs from a pre-existing streamline shapefile...
            df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(self.streamsPath_bankpts.text(), self.cell_size, str(self.reachid_bankpts.currentText()))
            
            # Build and write the Xn shapefile...
            funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(self.xnsPath.text()), False, int(self.p_xngap.text()), int(self.p_fitlength.text()), float(self.p_xnlength.text()))     
        
        # << Calculate Channel Metrics -- Bank Points Method >>
        if self.ck_chanmetrics_bankpts.isChecked():

#            print('here!')

            # Read the cross section shapefile and interplote elevation along each one using the specified DEM
            df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(self.xnsPath.text(), self.demPath_bankpts.text())
            
            XnPtDist = self.cell_size
            
            # Calculate channel metrics and write bank point shapefile...# parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh
            funcs_v2.chanmetrics_bankpts(df_xn_elev, str(self.xnsPath.text()), str(self.demPath_bankpts.text()), str(self.bankptsPath.text()), float(self.p_ivert.text()), XnPtDist, float(self.p_ratiothresh.text()), float(self.p_slpthresh.text()))
                   
        # << Create bank pixel layer (.tif) >>
        if self.ck_createPixels.isChecked():
             
            print('Reading streamline coordinates...') 
            df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(self.streamsPath_bankpixels.text(), self.cell_size, str(self.reachid_bankpixels.currentText()))
            
            print('Creating bank pixel layer...') 
            funcs_v2.bankpixels_from_streamline_window(df_coords, str(self.demPath_bankpixels.text()), str(self.bankpixelsPath.text()))
            
        # << Calculate Channel Metrics -- Bank Pixels Method >>
        if self.ck_chanmetrics_pixels.isChecked():
            print('Calculating channel width from bank pixel layer!')
            # def channel_width_bankpixels(str_streamlines_path, str_bankpixels_path, str_reachid):                        
            funcs_v2.channel_width_bankpixels(self.streamsPath_bankpixels.text(), self.bankpixelsPath.text(), self.reachid_bankpixels.currentText(), self.cell_size)
        
#        # << Calculate Floodplain >>
#        if self.ck_fp.isChecked():
#            print('Create floodplain here!')
            
        # << Run All Functions >>
        if self.ck_runall_bankpts.isChecked():
            print('Run all is currently disabled')

        if self.ck_runall_pixels.isChecked():
            print('Run all is currently disabled')
            
        print('Done with main!')
                
        # Show message when complete...
        # '\tTotal time:  ' + str(timeit.default_timer() - start_time_0
        QtGui.QMessageBox.information(self, 'Code Complete', 'Successful Run! \n\n\t   {0:.2f} secs'.format(timeit.default_timer() - start_time_0))
                 
if __name__ == '__main__':    
    
##    # ============================ << FOR DEBUGGING >> ==================================
##    print('\n<<< Start >>>\r\n')
##    start_time_0 = timeit.default_timer() 
##    
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM.tif'
    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_dem.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700dem_UTM18N_drclip_3mresample.tif'
##    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEM_02050205_USGS_AEA.tif'
#    str_dem_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMfel_dem_post.tif'
##    
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/Network04.shp'
    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/10m/020700_flowsUTM18N_drsheds_clip_diss.shp'
##    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_StreamNetwork_USGS_AEA_2.shp'
#    str_net_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/DEMnet1000.shp'
##    
##    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/bank_pixels_02050205_H.tif'
    str_bankpixels_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpixels.tif'
###    
    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_bankpoints.shp'
##    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_bankpoints_USGS_AEA.shp'
#    str_bankpts_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_bankpts.shp'
##   
    str_fp_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/fp_pixels.tif'
##     
#    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/02050205_xns_USGS_AEA_2.shp'
    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_net_xnptdist1.shp'
#    str_xns_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/020700_xns.shp'
##    
    str_hand_path = '/home/sam.lamont/USGSChannelFPMetrics/drb/test_gis/dr3m_hand.tif'
    
##    
    str_sheds_path = '/home/sam.lamont/USGSChannelFPMetrics/testing/dr3m_sheds_diss2.shp'    
    
#    str_reachid='LINKNO'
    str_reachid='ARCID'
    cell_size=3
    
#        self.p_fitlength.setText('3')
#        self.p_xngap.setText('3')
#        self.p_xnlength.setText('30')
#        self.p_ivert.setText('0.2')
#        self.p_ratiothresh.setText('1.5')
#        self.p_slpthresh.setText('0.03')    
    
    parm_ivert = 0.2
    XnPtDist = 3
    parm_ratiothresh = 1.5
    parm_slpthresh = 0.03

##    # << Channel Width via Bank Points >>    
###    buff_dist = 10.0
    funcs_v2.channel_width_bankpixels(str_net_path, str_bankpixels_path, str_reachid, cell_size)
##    
##    # << Build Strealine Coordinates >>
##    # Build reach coords and get crs from a pre-existing streamline shapefile...
#    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(str_net_path, cell_size, str_reachid)
#    df_coords.to_csv('df_coords_DR.csv') # just for testing
#   
#        << Find bank pixels using moving window along streamline >>
#    funcs_v2.bankpixels_from_streamline_window(df_coords, str_dem_path, str_bankpixels_path)  
#
#    print('Reading pre-calculated csv file...')
#    df_coords = pd.read_csv('df_coords_DR.csv', )
#
#      << Find bank pixels using moving window along streamline >>
#    funcs_v2.fp_from_streamline_window(df_coords, str_dem_path, str_fp_path)
#       
#      << Analyze hand using moving window along streamline >>
#    funcs_v2.analyze_hand_from_streamline_window(df_coords, str_hand_path)   ### USE THIS ONE
#    
#      << Analyze the DEM using 2D (buffered) cross-sections >>  ## OTHER TESTs
#    funcs_v2.analyze_hand_2D_xns(str_xns_path, str_hand_path, parm_ivert) # just use DEM here?   
#    funcs_v2.analyze_hand_reach_buffers(str_net_path, str_hand_path, str_reachid) 
#    funcs_v2.analyze_hand_sheds(str_sheds_path, str_hand_path, parm_ivert)
    
    # Build and write the Xn shapefile...
#    funcs_v2.write_xns_shp(df_coords, streamlines_crs, str(str_xns_path), False, int(3), int(3), float(150))     

#    # << INTERPOLATE XNs >>
#    df_xn_elev = funcs_v2.read_xns_shp_and_get_dem_window(str_xns_path, str_dem_path)
    
#    print('Writing df_xn_elev to .csv for testing...')
#    df_xn_elev.to_csv(columns=['index','linkno','elev','xn_row','xn_col']) 
#    df_xn_elev2 = pd.read_csv('df_xn_elev.csv') #, dtype={'linko':np.int,'elev':np.float,'xn_row':np.float,'xn_col':np.float})
 
#    # Calculate channel metrics and write bank point shapefile...
#    print('Calculating channel metrics from bank points...')
#    funcs_v2.chanmetrics_bankpts(df_xn_elev, str_xns_path, str_dem_path, str_bankpts_path, parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh)
#     =====================================================================================

#    # ===== Run GUI ==============================
#    app = QtGui.QApplication(sys.argv)
#
#    mainform = MainForm()
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
#    # ============================================