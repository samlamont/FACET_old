# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:11:00 2016

@author: sam.lamont
"""

from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtWidgets import QMainWindow, QApplication, QDesktopWidget, QFileDialog
import glob
import timeit
import numpy as np
import os
import rasterio
import fiona
import pandas as pd
import sys
import funcs_v2


class Facet(QMainWindow):

    def __init__(self, parent=None):

        super(Facet, self).__init__(parent)

        uic.loadUi('facet.ui', self)

        # << Browse buttons >>
        # Pre-Processing tab:
        self.browse_dem_pre.clicked.connect(self.get_dem_file)  # Raw DEM for pre-processing
        self.browse_nhd_pre.clicked.connect(self.get_streams_file)  # Input NHD for weight grid
        self.browse_hand_pre.clicked.connect(self.set_hand_file)  # Output HAND grid
        self.browse_streamlines_pre.clicked.connect(self.set_streamlines_pre_file)  # Output delineated streamlines
        self.browse_taudem_pre.clicked.connect(
            self.get_taudem_pre_path)  # Path to directory containing TauDEM executables
        self.browse_whitebox_pre.clicked.connect(self.get_whitebox_pre_exe)  # Path to the Whitebox exe
        self.browse_mpi_pre.clicked.connect(self.get_mpi_pre_exe)  # Path to the MPI exe (for TauDEM parallelization)

        # Channel Metrics - Xns tab:
        self.browseDEM_Xns.clicked.connect(self.get_dem_file)  # DEM - Xns
        self.browseStreams_Xns.clicked.connect(self.get_streams_file)  # Streams - Xns
        self.browseXns.clicked.connect(self.get_or_set_xns_file)  # Xns
        self.browseBankPts.clicked.connect(self.set_bankpts_file)  # BankPts

        # Channel Metrics - Curvature tab:
        self.browseDEM_Curve.clicked.connect(self.get_dem_file)  # DEM - curvature
        self.browseStreams_Curve.clicked.connect(self.get_streams_file)  # Streams - curvature
        self.browseBankPixels.clicked.connect(self.get_or_set_bankpixels)  # Bank Pixels

        # Floodplain Analysis tab:
        # TODO

        # << Button box >>      
        #        self.buttonBox.accepted.connect(self.genThread.start)
        self.buttonBox.accepted.connect(self.run_main)  # run_main in a thread?
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

    # =========================================    
    #              << FILE I/O >>    
    # =========================================        
    def get_dem_file(self):

        fname_dem = QFileDialog.getOpenFileName(self, 'Browse', '', 'Grids (*.*)')

        # fname_dem = 'blah'

        if self.tabWidget.currentIndex() == 0:  # Pre-processing
            self.txt_dem_pre.setText(fname_dem)

        if self.tabWidget.currentIndex() == 1:  # Channel metrics via Xn's
            self.textDEM_xns.setText(fname_dem)

        if self.tabWidget.currentIndex() == 2:  # Channel metrics via Curvature
            self.textDEM_curvature.setText(fname_dem)

        # # << Get DEM and set units, get cell size >>
        # with rasterio.open(str(fname_dem)) as ds_dem:
        #
        #     self.dem_crs = ds_dem.crs
        #     self.dem_affine = ds_dem.affine
        #     #            self.dem_affine = ds_dem.transform # NOTE: Confused about affine vs. transform here
        #     self.cell_size = ds_dem.affine[0]
        #     self.nodata_val = ds_dem.nodata
        #
        #     #            print('DEM affine: {}'.format(ds_dem.affine))
        #     print('DEM crs: {}'.format(ds_dem.crs))
        #     #            print('ds_dem.crs.wkt: {}'.format(ds_dem.crs.wkt))
        #
        #     try:
        #         # self.units = self.dem_crs['units']
        #         self.lbl_fit_units.setText(self.dem_crs['units'])
        #         self.lbl_xngap_units.setText(self.dem_crs['units'])
        #         self.lbl_xnlength_units.setText(self.dem_crs['units'])
        #         self.lbl_ivert_units.setText(self.dem_crs['units'])
        #     except:
        #         self.lbl_fit_units.setText('unknown')
        #         self.lbl_xngap_units.setText('unknown')
        #         self.lbl_xnlength_units.setText('unknown')
        #         self.lbl_ivert_units.setText('unknown')

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
        self.settings.setValue('taudem_path', fname)  # Save to settings

    def get_mpi_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.exe*)')
        self.txt_mpi_pre.setText(fname)
        self.settings.setValue('mpi_path', fname)  # Save to settings

    def get_whitebox_pre_exe(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Browse', '', 'Executable (*.*)')
        self.txt_whitebox_pre.setText(fname)
        self.settings.setValue('whitebox_path', fname)  # Save to settings

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

        str_orderid = 'Order_'  # TESTING

        #         Check which tab is active??
        #        i_tab = self.tab_widget.currentIndex() # (0: bankpts; 1: bankpixels; 2: FP-HAND)
        #        print('selected tab:  {}'.format(i_tab))

        try:

            # << Pre-processing >>
            if self.tabWidget.currentIndex() == 0:  # Run pre-processing functions

                print('About to preprocess')
                funcs_v2.preprocess_dem(str(self.txt_dem_pre.text()), str(self.txt_nhd_pre.text()),
                                        str(self.txt_mpi_pre.text()), str(self.txt_taudem_pre.text()),
                                        str(self.txt_whitebox_pre.text()), self.chk_wg_pre.isChecked(),
                                        self.chk_whitebox_pre.isChecked(), self.chk_taudem_pre.isChecked())

                # << Xn-based analysis >>
            if self.tabWidget.currentIndex() == 1:

                # Check coordinate systems (DEM vs. Streamlines)...
                if dict(self.dem_crs) != dict(self.streamlines_crs):
                    QtGui.QMessageBox.critical(self, 'Warning!',
                                               'Coordinate Systems of DEM and Streams Layer May Not Match')

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
                    df_coords, streamlines_crs = funcs_v2.workerfuncs.get_stream_coords_from_features(funcs,
                                                                                                      self.textStreams_Xns.text(),
                                                                                                      self.cell_size,
                                                                                                      str(
                                                                                                          self.comboBoxXns.currentText()),
                                                                                                      str_orderid)

                    #                self.lblProgBar.setText('Generating cross sections...')
                    #                # Build and write the Xn shapefile...
                    funcs_v2.workerfuncs.write_xns_shp(funcs, df_coords, streamlines_crs, str(self.textXns.text()),
                                                       False, int(self.textXnSpacing.text()),
                                                       int(self.textFitLength.text()), float(self.textXnLength.text()))
                    #
                #                self.genThread.terminate() # stop thread

                # << Calculate Channel Metrics -- Xns Method >>
                if self.chkMetricsXns.isChecked():
                    #                self.lblProgBar.setText('Extracting elevation along cross sections...')
                    # Read the cross section shapefile and interplote elevation along each one using the specified DEM
                    df_xn_elev = funcs_v2.workerfuncs.read_xns_shp_and_get_dem_window(funcs, self.textXns.text(),
                                                                                      self.textDEM_xns.text())

                    XnPtDist = self.cell_size

                    #                self.lblProgBar.setText('Calculating channel metrics...')
                    # Calculate channel metrics and write bank point shapefile...# parm_ivert, XnPtDist, parm_ratiothresh, parm_slpthresh
                    funcs_v2.workerfuncs.chanmetrics_bankpts(funcs, df_xn_elev, str(self.textXns.text()),
                                                             str(self.textDEM_xns.text()), str(self.textBankPts.text()),
                                                             float(self.textVertIncrement.text()), XnPtDist,
                                                             float(self.textRatioThresh.text()),
                                                             float(self.textSlpThresh.text()))

            # << Mean Curvature analysis >>
            if self.tabWidget.currentIndex() == 2:

                # << Create bank pixel layer (.tif) >>
                if self.chkBankPixels.isChecked():
                    self.lblProgBar.setText('Reading streamline coordinates...')
                    df_coords, streamlines_crs = funcs_v2.get_stream_coords_from_features(
                        self.textStreams_curvature.text(), self.cell_size, str(self.comboBoxCurvature.currentText()),
                        str_orderid, self)

                    self.lblProgBar.setText('Creating bank pixel layer...')
                    funcs_v2.bankpixels_from_streamline_window(df_coords, str(self.textDEM_curvature.text()),
                                                               str(self.textBankPixels.text()), self)

                # << Calculate Channel Metrics -- Bank Pixels Method >>
                if self.chkPixelMetrics.isChecked():
                    self.lblProgBar.setText('Calculating channel width from bank pixel layer!')
                    funcs_v2.channel_width_bankpixels(self.textStreams_curvature.text(), self.textBankPixels.text(),
                                                      self.comboBoxCurvature.currentText(), self.cell_size, self)

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
        #            self.lblProgBar.setText('Done!
        #            Elasped Time:  {0:.2f} mins'.format(float(timeit.default_timer() - start_time_0))/60.)

        except BaseException as e:
            QtGui.QMessageBox.critical(self, 'Error!', '{}'.format(str(e)))

        self.tabWidget.setEnabled(True)
        #        self.close()

        #        self.emit( QtCore.SIGNAL('add(QString)'), 'run_main complete')

        # Show message when complete...
        QtGui.QMessageBox.information(self, 'Code Complete', 'Successful Run! \n\n\t   {0:.2f} secs'.format(
            timeit.default_timer() - start_time_0))

    #    def done(self):


#        QtGui.QMessageBox.information(self, "Done!", "Done whatevah!")


#    def updateUI(self, val):
#        self.progressBar.setValue(val)   

if __name__ == '__main__':
    app = QApplication(sys.argv)

    mainform = Facet()

    # Show it in the center of the window...
    fr_geo = mainform.frameGeometry()
    cp = QDesktopWidget().availableGeometry().center()
    fr_geo.moveCenter(cp)
    mainform.move(fr_geo.topLeft())

    mainform.show()

    sys.exit(app.exec_())
