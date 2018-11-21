# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:16:17 2018

@author: sam.lamont
"""

# =================================================================================
#  Calculate channel width based on bank pixels and stream line parallel offsets, 
#  potentially subdividing using Xn's
#  NOTE: Make this a generic metric calculator by segment?  (ie., sinuosity, etc)
# =================================================================================
def channel_and_fp_2Dxn_analysis(df_coords, str_streamlines_path, str_bankpixels_path, str_hand_path, str_fim_path, str_reachid, cell_size, i_step, max_buff, p_fpxnlen):
    
    print('Channel and floodplain metrics from 2D cross sections along reach segments...')    

    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 
    
    lst_geom=[] # TESTING.  For saving 2D Xn buffer polygons

    # Successive buffer-mask operations to count bank pixels at certain intervals
    lst_buff=range(cell_size,max_buff,cell_size)    
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_ch_fp_width.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')
    
    ## Schema for the output properties file:
    schema_output = {'geometry': 'LineString', 'properties': {'linkno':'int','ch_wid_total':'float', 'ch_wid_1':'float', 'ch_wid_2':'float',
                                                              'dist_sl':'float', 'dist':'float', 'sinuosity':'float','fp_width':'float','fp_range':'float'}}                                  

    ## Access the hand grid:
    with rasterio.open(str_hand_path, 'r') as ds_hand:                        
    
        ## Access the FIM:
        with rasterio.open(str_fim_path, 'r') as ds_fim:
    
            # Access the bank pixel layer:
            with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:    
                
                # Access the streamlines layer...
                with fiona.open(str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
              
        #            progBar.setRange(0, len(streamlines)) 
                    
                    # Get the crs...
                    streamlines_crs = streamlines.crs                
                                
                    # Open another file to write the output props:
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
                                         
                                ## Call 2D Xn channel and FP analysis here:
                                parm_ivert=0.1
                                fp_width, fp_range = metrics_from_2D_xns(ds_fim, ds_hand, df_seg, p_fpxnlen, dist_sl, parm_ivert)
                            
#                                lst_geom.append(xn_buff)
                                
                                # Write to the output shapefile here...
                                output.write({'properties':{'linkno':i_linkno,'ch_wid_total': weighted_avg_left+weighted_avg_rt,'ch_wid_1': weighted_avg_left,'ch_wid_2': weighted_avg_rt, 'dist_sl':dist_sl,
                                                            'dist':dist,'sinuosity': sinuosity,'fp_width':fp_width,'fp_range':fp_range.astype(np.float64)},'geometry':mapping(ls)})                                                    
                                                              
                                
                            break
                                
                                    
#                    if j > 50: break
#    print('Saving test 2D Xn buffer layer...')                           
#    gdf_buff=gpd.GeoDataFrame()
#    gdf_buff['geometry']=lst_geom
#    gdf_buff.crs = from_epsg(26918)
#    gdf_buff.to_file(r"D:\facet\dr_working_data\dr_working_data\dr3m_2d_buffs_test.shp" )
                            
#    ## Save output metrics in segmented streamline file:
#    print('Saving channel and floodplain metric stream segment layer...')                           
#    gdf_buff=gpd.GeoDataFrame()
#    gdf_buff['geometry']=lst_geom
#    gdf_buff.crs = from_epsg(26918)
#    gdf_buff.to_file(r"D:\facet\dr_working_data\dr_working_data\dr3m_2d_buffs_test.shp" )                            
    
    return

# ===============================================================================
#  Metrics using polygons from buffering cross-sections
# ===============================================================================
def metrics_from_2D_xns(ds_fim, ds_dem, tpl_seg, p_fpxnlen):
    try:
        if tpl_seg.geometry:
            x,y=zip(*mapping(tpl_seg.geometry)['coordinates'])
        
            ## Get the segment midpoint:                                                    
            midpt_x = x[int(len(x)/2)]
            midpt_y = y[int(len(y)/2)]                              
            
            # Build a cross-section from the end points:
            lst_xy = build_xns(y, x, midpt_x, midpt_y, p_fpxnlen)
            
            try:            
                # Turn the cross-section into a linestring:
                fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
            except:
                print('Error converting Xn endpts to LineString')
                pass
                                            
            ## Buffer the cross section:        
            buff_len=tpl_seg.dist_sl/1.85 # about half of the line segment straight line distance
            geom_fpls_buff = fp_ls.buffer(buff_len, cap_style=2)
            xn_buff = mapping(geom_fpls_buff)                    
            
            # Mask the fp for each feature...
            w_fim, trans_fim = rasterio.mask.mask(ds_fim, [xn_buff], crop=True)
            w_fim=w_fim[0]
            w_fim=w_fim[w_fim!=ds_fim.nodata]
            
            # << Related to mapping the floodplain based on HAND height >>
            # Count the number of pixels in the buffered Xn...
            num_pixels = w_fim.size
             
            # Calculate area of FP pixels...
            area_pixels = num_pixels*(ds_fim.res[0]**2) # get grid resolution               
            
            # Calculate width by stretching it along the length of the 2D Xn...
            fp_width = area_pixels/(buff_len*2)   
        #    fp_width=0 # For testing purposes
            
            ## TO DO:  Get other properties by analyzing values in w_fim (depth) 
            try:
                fp_range=w_fim.max()-w_fim.min()
            except:
                fp_range=0
                pass
                
        #    # Subtract channel width from fp width...
            fp_width = fp_width - tpl_seg.ch_wid_tot
            
            if fp_width<0.: fp_width = 0 # don't have negatives              
    
        # << CALL HAND VERTICAL SLICE ANALYSIS HERE >>                                                       
    #    bank_height_hand, chan_width_hand, bank_ang_hand = analyze_hand_poly(ds_hand, xn_buff, p_fpxnlen, parm_ivert)
    except Exception as e:
        print(f'Error in metrics_from_2D_xns: {str(e)}')

    return fp_width, fp_range, geom_fpls_buff

# ===============================================================================
#  Analyze DEM in vertical slices using an individual polygon
# =============================================================================== 
def analyze_hand_poly(ds_hand, xn_buff, p_fpxnlen, parm_ivert):

    ## Mask hand grid using the 2D Xn buffer polygon:
    w_hnd, trans_hnd = rasterio.mask.mask(ds_hand, [xn_buff], crop=True)
    w_hnd=w_hnd[0]    
                    
    i_rng=10
    arr_slices = np.arange(parm_ivert, i_rng, parm_ivert)    
    
    lst_count=[]
    lst_width=[]

    # List comprehension here instead??
    for i_step in arr_slices: 

        num_pixels = w_hnd[(w_hnd<=i_step) & (w_hnd>=0.)].size           
        lst_count.append(num_pixels) # number of pixels greater than or equal to zero and less than the height interval
            
        # Calculate area of FP pixels:
        area_pixels = num_pixels*(ds_hand.res[0]**2)
        
        # Calculate width by stretching it along the length of the 2D Xn:
        lst_width.append(area_pixels/p_fpxnlen)                 

    df_steps = pd.DataFrame({'count':lst_count, 'height':arr_slices, 'width':lst_width})
    
#    df_steps.plot(x='width',y='height', marker='.')  

    ## IDEA:  You basically need to fit a straight line to segments along the curve and test the slope of that line          
    
    # Slope of width...
    df_steps['width_diff'] = df_steps['width'].diff() # Use .diff(2) or 3 to average?
    
    # Slope of slope of count...
    df_steps['width_diff_2nd'] = df_steps['width_diff'].diff()  # Max val for width_diff_2nd is the bank?
    
    # Find the top three maximum diff_2nd values and select the one with the lowest height?
    df_top3 = df_steps.nlargest(3, columns='width_diff_2nd')
    
    bank_height = df_steps['height'].iloc[df_top3['height'].idxmin()]
    chan_width = df_steps['width'].iloc[df_top3['height'].idxmin()]
    bank_ang = np.arctan(bank_height/chan_width)
                    
    return bank_height, chan_width, bank_ang

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
    
# ===============================================================================
#  Calculate FP width using 2D Xn's (buffered Xn lines)
#  Requires existing cross section file
# =============================================================================== 
def floodplain_width_2D_xns(str_xns_path, str_grid_path, buff_dist):
    
    ###### TESTING
    str_streamspath = r'D:\CFN_data\DEM_Files\020502061102_ChillisquaqueRiver\DEM_breach_net.shp'
    ######
           
    # Open the dem...
    with rasterio.open(str(str_grid_path)) as ds_grid:
        
        # Open the xn layer...
        with fiona.open(np.str(str_xns_path), 'r') as xns: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)

            # Create the filename and path for the output file... # NOT OS INDEPENDENT??
            head, tail = os.path.split(str_xns_path)    
            str_outname = tail[:-4] + '_streams_fpwidth.shp'
            str_outpath = head + '/' + str_outname 
       
            # Get the crs...
            schema = xns.schema.copy()
            schema['properties']['fp_width'] = 'float'
           
            lst_lyr=[]
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema, xns.crs) as output:  
                with fiona.open(str_streamspath, 'r') as streams:
                    
                    for line in streams:                
                        lst_lyr.append(shape(line['geometry']))                
            
                    lyr = MultiLineString(lst_lyr)
    
                    for line in xns:    
        
        #                test_id=777
        #                if line['properties']['linkno'] <> test_id:
        #                    continue
                        
                        print(line['properties']['linkno'])
                        
                        # Buffer each feature...
                        geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
                                            
                        geom_buff = geom.buffer(buff_dist,cap_style=2)
                        buff = mapping(geom_buff)           
                                    
                        # Mask the bankpts file for each feature...
                        out_image, out_transform = rasterio.tools.mask.mask(ds_grid, [buff], crop=True)
                        
                        # Count the number of pixels in the buffered Xn...
                        num_pixels = len(out_image[out_image==1])
                        
        #                w[w<0]=9999. # handle no data vals?
                        
                        # Calculate area of FP pixels...
                        area_pixels = num_pixels*(ds_grid.res[0]**2) # get grid resolution               
                        
                        # Calculate width by stretching it along the length of the 2D Xn...
                        fp_width = area_pixels/(buff_dist*2)
                        
                        # Intersect buffer polygon with streamlines...
                        line_seg = geom_buff.intersection(lyr)
                        
                        output.write({'properties':{'linkno':0,'fp_width':fp_width}, 'geometry':mapping(line_seg)})
                                                
#                        line['properties']['fp_width'] = fp_width
                        # 'properties':{'width_total': 1, 'width_left': 0, 'width_right':1}
#                        for ls in line_seg:
#                            output.write({'properties':{'linkno':0,'fp_width':fp_width}, 'geometry':mapping(ls)})

        # Remove the original FP Xn file?
        os.remove(str_xns_path)
                    
    return


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
def floodplain_width_fppixels_segments_po(df_coords, str_streamlines_path, str_fpgrid_path, str_reachid, cell_size):
    
    print('Floodplain width from bank pixels -- segmented...')
    
#    lst_output=[]
    
#    j=0   
#    progBar = self.progressBar
#    progBar.setVisible(True) 
    
    # Create the filename and path for the output file... # NOT OS INDEPENDENT??
    head, tail = os.path.split(str_streamlines_path)    
    str_outname = tail[:-4] + '_fpwidth_po_avg.shp'
    str_outpath = head + '/' + str_outname 
    
    gp_coords = df_coords.groupby('linkno')

    schema_output = {'geometry': 'LineString', 'properties': {'width_total':'float', 'width_left':'float', 'width_right':'float'}}
    
    # Open the bankpts layer...
    with rasterio.open(str(str_fpgrid_path)) as ds_fppixels:    
        
        # Open the streamlines layer...
        with fiona.open(np.str(str_streamlines_path), 'r') as streamlines: # NOTE: For some reason you have to explicitly convert the variable to a string (is it unicode?)
       
#            progBar.setRange(0, len(streamlines)) 
            
            # Get the crs...
            streamlines_crs = streamlines.crs                
            
            # Open another file to write the width...
            with fiona.open(str_outpath, 'w', 'ESRI Shapefile', schema_output, streamlines_crs) as output:
                
                for i_linkno, df_linkno in gp_coords:
                    
#                    progBar.setValue(j)
#                    j+=1                    
            
                    i_linkno = int(i_linkno)
                    
                    if i_linkno <> 605:
                        continue
                    
                    print('linkno:  {}'.format(i_linkno))
                    
                    # testing...
                    arr_x = df_linkno.x.values
                    arr_y = df_linkno.y.values    
                    
                    # Create a line segment from endpts in df_linkno...
                    ls = LineString(zip(arr_x, arr_y))                    
                    
                    # Successive buffer-mask operations to count bank pixels at certain intervals
                    lst_left=[]
                    lst_rt=[]                    
                    lst_buff=range(cell_size,60,cell_size)
                    for buff_dist in lst_buff:                                                
                        
                        ls_offset_left = ls.parallel_offset(buff_dist, 'left', resolution=16, join_style=2, mitre_limit=1)
                        ls_offset_rt = ls.parallel_offset(buff_dist, 'right', resolution=16, join_style=2, mitre_limit=1)
                        
#                        try:
#                            arr_os_left = asarray(ls_offset_left)
#                            arr_os_rt = asarray(ls_offset_rt)
#                        except:
#                            continue
#                        
#                        # Testing...
#                        left_pts = MultiPoint(arr_os_left)
#                        out_left, out_transform = rasterio.tools.mask.mask(ds_fppixels, [mapping(left_pts)], crop=True)
                        
                        # Now segment the offsets and do the masking?
#                        i_step=30
#                        arr_ind_left = np.arange(0, len(arr_os_left), int(len(arr_os_left)/i_step)) # NOTE: Change the step for longer/shorter segments                        
#                        arr_ind_rt = np.arange(0, len(arr_os_rt), int(len(arr_os_rt)/i_step))
                        
                        ###############################################         
                        # Left...                                    
                        int_pts = np.arange(0, ls_offset_left.length, 10)                        
                        left_pts = MultiPoint([ls_offset_left.interpolate(i_dist) for i_dist in int_pts])  # list of points at distance intervals
                        left_splitted = split(ls_offset_left, left_pts)
                        
                        # Right...
                        int_pts = np.arange(0, ls_offset_rt.length, 10)                        
                        rt_pts = MultiPoint([ls_offset_rt.interpolate(i_dist) for i_dist in int_pts])  # list of points at distance intervals
                        rt_splitted = split(ls_offset_rt, rt_pts)                        
                        ##############################################
                        
#                        # Get points at specified indices...
#                        left_coords = arr_os_left[arr_ind_left]
#                        rt_coords = arr_os_rt[arr_ind_rt]
#                        
#                        left_pts = MultiPoint(left_coords)
#                        rt_pts = MultiPoint(rt_coords)
#                                                
#                        left_splitted = split(ls_offset_left, left_pts)
#                        rt_splitted = split(ls_offset_rt, rt_pts)
                                                
#                        out_left, out_transform = rasterio.tools.mask.mask(ds_fppixels, [mapping(left_splitted)], crop=True)                         
#                        out_rt, out_transform = rasterio.tools.mask.mask(ds_fppixels, [mapping(rt_splitted)], crop=True) 
                        
                        for i, line_left in enumerate(left_splitted):
                            out_left, out_transform = rasterio.tools.mask.mask(ds_fppixels, [mapping(line_left)], crop=True) 
                            output.write({'properties':{'width_total': 1, 'width_left': 1, 'width_right':0}, 'geometry':mapping(line_left)})
                            
                            num_pixels_left = len(out_left[out_left==1])
                            
                            tpl_left = (buff_dist, i, num_pixels_left)
                            lst_left.append(tpl_left)
                            
                        for i, line_rt in enumerate(rt_splitted):
                            out_rt, out_transform = rasterio.tools.mask.mask(ds_fppixels, [mapping(line_rt)], crop=True) 
                            output.write({'properties':{'width_total': 1, 'width_left': 0, 'width_right':1}, 'geometry':mapping(line_rt)})      
                            
                            num_pixels_rt = len(out_rt[out_rt==1])
                            
                            tpl_rt = (buff_dist, i, num_pixels_rt)
                            lst_rt.append(tpl_rt)
                            
                    break
                            
#                    # Now create one big dataframe...
#                    df_left = pd.DataFrame(lst_left, columns=['buff_dist','step','count'])
#                    df_rt = pd.DataFrame(lst_rt, columns=['buff_dist','step','count'])
#                    
#                    df_test = df_left[df_left['step']==0]
                    
                    print('pause')
        
                                    
#                            tpl_out = i_linkno, buff_dist, num_pixels_left, num_pixels_rt
#                            lst_tally.append(tpl_out)                    
#                            df_tally = pd.DataFrame(lst_tally, columns=['linkno','buffer','interval_left','interval_rt'])
                   
        
                        # Avoid division by zero                     
#                        if df_tally.interval.sum() == 0:
#                            continue                    
                        
#                        # Calculate weighted average                     
#                        # Only iterate over the top 3 or 2 (n_top) since distance is favored...    
#                        try:                        
#                            weighted_avg_left=0 
#                            weighted_avg_rt=0
#                        
#                            avg_left = df_tally[df_tally.interval_left>0].nsmallest(3, 'interval_left').interval_left.mean()                            
#                            fp_left = np.interp(avg_left, df_tally.buffer, df_tally.interval_left)
#                            
#                            avg_rt = df_tally[df_tally.interval_rt>0].nsmallest(3, 'interval_rt').interval_rt.mean()                            
#                            fp_rt = np.interp(avg_rt, df_tally.buffer, df_tally.interval_rt)                            
#                            
##                            weighted_avg_left = np.float(df_tally[df_tally.interval_left>0].nsmallest(1, 'interval_left').buffer.max())
##                            weighted_avg_rt = np.float(df_tally[df_tally.interval_rt>0].nsmallest(1, 'interval_rt').buffer.max())
#                            
#                            n_top=10
#                            
#                            for tpl in df_tally[df_tally.interval_left>0].nsmallest(n_top, 'interval_left').itertuples():
#                                weighted_avg_left += tpl.buffer*(np.float(tpl.interval_left)/np.float(df_tally[df_tally.interval_left>0].nsmallest(n_top, 'interval_left').sum().interval_left))
#
#                            for tpl in df_tally[df_tally.interval_rt>0].nsmallest(n_top, 'interval_rt').itertuples():
#                                weighted_avg_rt += tpl.buffer*(np.float(tpl.interval_rt)/np.float(df_tally[df_tally.interval_rt>0].nsmallest(n_top, 'interval_rt').sum().interval_rt))
#                                
##                                weighted_avg=weighted_avg*2   # Multiply buffer by 2 to get width
#                        except:
##                            weighted_avg=-9999.
#                            continue
                        
                        # Write to an output file here...
#                        output.write({'properties':{'width_total': fp_left+fp_rt, 'width_left': fp_left, 'width_right': fp_rt}, 'geometry':mapping(ls)})
#                        output.write({'properties':{'width_total': weighted_avg_left+weighted_avg_rt, 'width_left': weighted_avg_left, 'width_right': weighted_avg_rt}, 'geometry':mapping(ls)})
                                                    
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
    
#def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
#    """Kernel Density Estimation with Scipy"""
#    # Note that scipy weights its bandwidth by the covariance of the
#    # input data.  To make the results comparable to the other methods,
#    # we divide the bandwidth by the sample standard deviation here.
#    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
#    return kde.evaluate(x_grid)   
#    
#def gauss_func(x, *p):
#    A, mu, sigma = p
#    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
# ===============================================================================
#  Analyze DEM in vertical slices successive buffers stream reaches
# ===============================================================================    
def analyze_hand_reach_buffers(str_net_path, str_hand_path, str_reachid):    
    
#    cell_size = 3
    
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
    
    str_reachid='gridcode'
           
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
#  Calculates openness from DEM window
# ===============================================================================  
def openness_from_elev_window(arr_elev_norm, arr_mask_dist, cell_size, searchR, ind_diag_nw, ind_diag_ne, ind_diag_sw, ind_diag_se):
    
    arr_mask_dist[50,50]=1 # TEST TEST
    
    arr_ang = np.arctan(arr_elev_norm/arr_mask_dist) # NOTE:  Division by zero here
    
    # Find maximum and minimum in every direction...OR construct an array of each direction and use np.amin/amax??
    df_dirs = pd.DataFrame([arr_ang[ind_diag_nw],
                            arr_ang[:searchR, searchR],
                            arr_ang[ind_diag_ne],
                            arr_ang[searchR,searchR+1:],
                            arr_ang[ind_diag_se],
                            arr_ang[searchR+1:,searchR],
                            arr_ang[ind_diag_sw],
                            arr_ang[searchR,:searchR]])

#                    nw_max = np.max(arr_ang[ind_diag_nw])
#                    n_max = np.max(arr_ang[:w_height/2, w_width/2])
#                    ne_max = np.max(arr_ang[ind_diag_ne])
#                    e_max = np.max(arr_ang[w_height/2,w_width/2:])
#                    se_max = np.max(arr_ang[ind_diag_se])
#                    s_max = np.max(arr_ang[w_height/2+1:,w_width/2])
#                    sw_max = np.max(arr_ang[ind_diag_sw])
#                    w_max = np.max(arr_ang[w_height/2,:w_width/2])
    
    # Get max and mins along columns...
    srs_max = df_dirs.max(axis=1)
#    srs_min = df_dirs.min(axis=1)
    
    # Convert to zenith and nadir (add/subtract from 90)
    srs_zenith = 90.0 - srs_max
#    srs_nadir = 90 + srs_min
                        
    # Openness is average of each (positive: zenith; negative: nadir)
    po = np.mean(srs_zenith)   
    
    return po
    
# ===============================================================================
#  Calculates openness from DEM window (alternate method)
# ===============================================================================  
def openness_from_elev_window_numpy(arr_elev_norm, arr_mask_dist, cell_size, searchR, ind_diag_nw, ind_diag_ne, ind_diag_sw, ind_diag_se):
    
    arr_mask_dist[50,50]=1 # TEST TEST
    
    arr_ang = np.arctan(arr_elev_norm/arr_mask_dist) # NOTE:  Division by zero here
       
    # NOTE:  Convert these to degrees...
    degcon = 180 / math.pi
    
    nw_max = degcon*np.max(arr_ang[ind_diag_nw])
    n_max = degcon*np.max(arr_ang[:searchR, searchR])
    ne_max = degcon*np.max(arr_ang[ind_diag_ne])
    e_max = degcon*np.max(arr_ang[searchR,searchR+1:])
    se_max = degcon*np.max(arr_ang[ind_diag_se])
    s_max = degcon*np.max(arr_ang[searchR+1:,searchR])
    sw_max = degcon*np.max(arr_ang[ind_diag_sw])
    w_max = degcon*np.max(arr_ang[searchR,:searchR])
    
    nw_min = degcon*np.min(arr_ang[ind_diag_nw])
    n_min = degcon*np.min(arr_ang[:searchR, searchR])
    ne_min = degcon*np.min(arr_ang[ind_diag_ne])
    e_min = degcon*np.min(arr_ang[searchR,searchR+1:])
    se_min = degcon*np.min(arr_ang[ind_diag_se])
    s_min = degcon*np.min(arr_ang[searchR+1:,searchR])
    sw_min = degcon*np.min(arr_ang[ind_diag_sw])
    w_min = degcon*np.min(arr_ang[searchR,:searchR])    
    
    # Get max and mins along columns...
#    srs_max = df_dirs.max(axis=1)
#    srs_min = df_dirs.min(axis=1)
    
    # Convert to zenith and nadir (add/subtract from 90)
#    srs_zenith = 90.0 - srs_max
#    srs_nadir = 90 + srs_min
                        
    # Openness is average of each (positive: zenith; negative: nadir)
#    po = np.mean(srs_zenith)   
    po_pos = np.mean([90. - nw_max,90. - n_max,90. - ne_max,90. - e_max,90. - se_max,90. - s_max,90. - sw_max,90. - w_max])
    po_neg = np.mean([90. + nw_min,90. + n_min,90. + ne_min,90. + e_min,90. + se_min,90. + s_min,90. + sw_min,90. + w_min])
    
    return po_pos, po_neg    

# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_openness_window_buffer_all(df_coords, str_dem_path, str_net_path, str_pos_path, str_neg_path):
    
    cell_size=3    
    searchRadius = 150 # m
    buff_dist=30    
    
    searchR = searchRadius/cell_size
    
    w_height=searchR*2 # number of rows.       This is actually half the window size
    w_width=searchR*2  # number of columns.    This is actually half the window size
    
    # Construct a distance array based on searchR and cell size...
    arr_dist_xy = cell_size*np.array(range(0,searchR+1,1))      
     
    arr_dist_ang = np.sqrt(cell_size**2 + cell_size**2)*np.array(range(0,searchR+1,1))  
    arr_mask_dist = np.ones((searchR*2+1, searchR*2+1)) # Add 1 only if it's an even number?
    
    arr_mask_dist[:,searchR] = np.hstack((arr_dist_xy[::-1], arr_dist_xy[1:]))
    arr_mask_dist[searchR,:] = np.hstack((arr_dist_xy[::-1], arr_dist_xy[1:]))
    
    # Get indices of cardinal directions for later...
    ind_diag = np.diag_indices_from(arr_mask_dist) # Will I ever use a non-square window?
    
    ind_diag_nw = (ind_diag[0][:searchR], ind_diag[1][:searchR])    
    ind_diag_se = (ind_diag[0][searchR+1:], ind_diag[1][searchR+1:])
    
    ind_diag_ne = (ind_diag[0][:searchR], np.abs(ind_diag_nw[1]-w_width))
    ind_diag_sw = (np.abs(ind_diag_nw[0]-w_height), ind_diag[1][:searchR])   
    
    arr_mask_dist[ind_diag] = np.hstack((arr_dist_ang[::-1], arr_dist_ang[1:]))
        
    # Flip it and fill the other diagonal...
    arr_mask_dist = np.fliplr(arr_mask_dist)
    arr_mask_dist[np.diag_indices_from(arr_mask_dist)] = np.hstack((arr_dist_ang[::-1], arr_dist_ang[1:]))      
    
    print('Bank pixels from openness windows, buffered')
    
    lst_lyr=[]
    
    
#    j=0
        
    # IDEA:  Buffer all streamlines first and extract DEM
    # Open the DEM...
    with rasterio.open(str(str_dem_path)) as ds_dem:    
        
        out_meta = ds_dem.meta.copy()  
        arr_pos_openness=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])          
        arr_neg_openness=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype']) 
    
        # Open the streamlines and build a multiline string for buffering...
        with fiona.open(np.str(str_net_path), 'r') as streamlines:
            
            for line in streamlines:                
                lst_lyr.append(shape(line['geometry']))                
        
            lyr = MultiLineString(lst_lyr)
            
            # Apply the buffer...
            print('Buffering stream lines...')
            lyr_buffer = lyr.buffer(buff_dist)
            
            buff = mapping(lyr_buffer)            
            
            # Clipping DEM...
            print('Clipping DEM...')
            dem_clip, out_transform = rasterio.tools.mask.mask(ds_dem, [buff], crop=False) 
            dem_clip = dem_clip[0]
            
            
#            out_meta['affine']=out_transform
#            out_meta['height']=np.shape(dem_clip)[0] # rows
#            out_meta['width']=np.shape(dem_clip)[1] # cols
            
#            plt.imshow(dem_clip[0])
            
#            # DO THIS PER REACH?
#            for line in streamlines:
#                
#                j += 1
#                print('j: {}'.format(j))
#                if j > 1: break
#            
#                # Clip the DEM using the buffered streams...            
#                geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
##                bounds = geom.bounds
#                geom_buff = geom.buffer(buff_dist,cap_style=2)
#                buff = mapping(geom_buff)                
#                                
#                dem_clip, out_transform = rasterio.tools.mask.mask(ds_dem, [buff], crop=True) 
#            
#        #        plt.imshow(dem_clip[0])
#            
#                # Transform to pixel space
#        #       col_max, row_max = ~ds_dem.affine * (df_coords['x1'], df_coords['y1'])
#                                        
#                print('Calculating Openness...')
#                
            shp=np.shape(dem_clip)
#                
#                bounds = rasterio.transform.array_bounds(shp[0],shp[1],out_transform) # window bounds in x-y space (west, south, east, north)
#                
#                # df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y']) 
#                col_min, row_min = ~ds_dem.affine * (bounds[0], bounds[3]) # upper left row and column of window?
              
            print('Calculating Openness...') # Gotta be a quicker way here
            for row in range(searchR, shp[0]-searchR):
                for col in range(searchR, shp[1]-searchR):
                    
                    if dem_clip[int(row),int(col)]<-9999.0: continue # skip nodata values
                    
                    row_min_w = np.int(row - searchR)
                    row_max_w = np.int(row + searchR)
                    col_min_w = np.int(col - searchR)
                    col_max_w = np.int(col + searchR)
                    
                    arr_elev = ds_dem.read(1, window=((row_min_w, row_max_w+1),(col_min_w, col_max_w+1)))                 
                    
                    # Get elevation difference...
                    arr_elev = arr_elev - dem_clip[row,col]
                    
                    po_pos, po_neg = openness_from_elev_window_numpy(arr_elev, arr_mask_dist, cell_size, searchR, ind_diag_nw, ind_diag_ne, ind_diag_sw, ind_diag_se)                    

                    arr_pos_openness[int(row),int(col)]=po_pos
                    arr_neg_openness[int(row),int(col)]=po_neg
                    
                    
         
        arr_pos_openness[arr_pos_openness<=0.] = out_meta['nodata'] 
        arr_neg_openness[arr_neg_openness<=0.] = out_meta['nodata'] 

        print('Writing Positive Openness .tif...')
        with rasterio.open(str_pos_path, "w", **out_meta) as dest:
            dest.write(arr_pos_openness, indexes=1)
            
        print('Writing Negative Openness .tif...')
        with rasterio.open(str_neg_path, "w", **out_meta) as dest:
            dest.write(arr_neg_openness, indexes=1)            

    
    return

# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_openness_window_buffer_per_reach(df_coords, str_dem_path, str_net_path, str_bankpixels_path):
    
    cell_size=3    
    searchRadius = 150 # m
    searchR = searchRadius/cell_size
    
    w_height=searchR*2 # number of rows.       This is actually half the window size
    w_width=searchR*2  # number of columns.    This is actually half the window size
    
    # Construct a distance array based on searchR and cell size...
    arr_dist_xy = cell_size*np.array(range(0,searchR+1,1))      
     
    arr_dist_ang = np.sqrt(cell_size**2 + cell_size**2)*np.array(range(0,searchR+1,1))  
    arr_mask_dist = np.ones((searchR*2+1, searchR*2+1)) # Add 1 only if it's an even number?
    
    arr_mask_dist[:,searchR] = np.hstack((arr_dist_xy[::-1], arr_dist_xy[1:]))
    arr_mask_dist[searchR,:] = np.hstack((arr_dist_xy[::-1], arr_dist_xy[1:]))
    
    # Get indices of cardinal directions for later...
    ind_diag = np.diag_indices_from(arr_mask_dist) # Will I ever use a non-square window?
    
    ind_diag_nw = (ind_diag[0][:searchR], ind_diag[1][:searchR])    
    ind_diag_se = (ind_diag[0][searchR+1:], ind_diag[1][searchR+1:])
    
    ind_diag_ne = (ind_diag[0][:searchR], np.abs(ind_diag_nw[1]-w_width))
    ind_diag_sw = (np.abs(ind_diag_nw[0]-w_height), ind_diag[1][:searchR])   
    
    arr_mask_dist[ind_diag] = np.hstack((arr_dist_ang[::-1], arr_dist_ang[1:]))
        
    # Flip it and fill the other diagonal...
    arr_mask_dist = np.fliplr(arr_mask_dist)
    arr_mask_dist[np.diag_indices_from(arr_mask_dist)] = np.hstack((arr_dist_ang[::-1], arr_dist_ang[1:]))      
    
    print('Bank pixels from openness windows, buffered...')
    
#    lst_lyr=[]
    buff_dist=10
    
    j=0
        
    # IDEA:  Buffer all streamlines first and extract DEM
    # Open the DEM...
    with rasterio.open(str(str_dem_path)) as ds_dem:    
        
        out_meta = ds_dem.meta.copy()  
        arr_openness=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])          
    
        # Open the streamlines and build a multiline string for buffering...
        with fiona.open(np.str(str_net_path), 'r') as streamlines:
            
#            for line in streamlines:                
#                lst_lyr.append(shape(line['geometry']))                
#        
#            lyr = MultiLineString(lst_lyr)
#            
#            # Apply the buffer...
#            lyr_buffer = lyr.buffer(buff_dist)
            
            # DO THIS PER REACH?
            for line in streamlines:
                
                j += 1
                print('j: {}'.format(j))
                if j > 1: break
            
                # Clip the DEM using the buffered streams...            
                geom = shape(line['geometry'])   # Convert to shapely geometry to operate on it
#                bounds = geom.bounds
                geom_buff = geom.buffer(buff_dist,cap_style=2)
                buff = mapping(geom_buff)                
                                
                dem_clip, out_transform = rasterio.tools.mask.mask(ds_dem, [buff], crop=True) 
            
        #        plt.imshow(dem_clip[0])
            
                # Transform to pixel space
        #       col_max, row_max = ~ds_dem.affine * (df_coords['x1'], df_coords['y1'])
                          
                dem_clip = dem_clip[0]
                
                print('Calculating Openness...')
                
                shp=np.shape(dem_clip)
                
                bounds = rasterio.transform.array_bounds(shp[0],shp[1],out_transform) # window bounds in x-y space (west, south, east, north)
                
                # df_coords['col'], df_coords['row'] = ~ds_dem.affine * (df_coords['x'], df_coords['y']) 
                col_min, row_min = ~ds_dem.affine * (bounds[0], bounds[3]) # upper left row and column of window?
            #    
                for row in range(shp[0]):
                    for col in range(shp[1]):
                        
                        if dem_clip[int(row),int(col)]<-9999.0: continue # skip nodata values
                        
                        row_min_w = np.int(row - searchR + row_min)
                        row_max_w = np.int(row + searchR + row_min)
                        col_min_w = np.int(col - searchR + col_min)
                        col_max_w = np.int(col + searchR + col_min)
                        
                        arr_elev = ds_dem.read(1, window=((row_min_w, row_max_w+1),(col_min_w, col_max_w+1)))                 
                        
                        # Get elevation difference...
                        arr_elev = arr_elev - dem_clip[row,col]
                        
                        po = openness_from_elev_window(arr_elev, arr_mask_dist, cell_size, searchR, ind_diag_nw, ind_diag_ne, ind_diag_sw, ind_diag_se)                    

                        arr_openness[int(row+row_min),int(col+col_min)]=po
                    

          
        arr_openness[arr_openness<=0.] = out_meta['nodata'] 
        print('Writing Openness .tif...')
        with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
            dest.write(arr_openness, indexes=1)

    
    return
    
# ===============================================================================
#  Searchable window with center pixel defined by get_stream_coords_from_features
# ===============================================================================  
def bankpixels_from_openness_window(df_coords, str_openness_path, str_bankpixels_path):
    
    print('Bank pixels from streamline windows of openness...')
    
    w_height=20 # number of rows
    w_width=20  # number of columns
    
    j=0   
    
    with rasterio.open(str_openness_path) as ds_open:
        
        # Transform to pixel space
        df_coords['col'], df_coords['row'] = ~ds_open.affine * (df_coords['x'], df_coords['y'])   
        
        df_coords[['row','col']] = df_coords[['row','col']].astype(np.int32)  
        df_coords.drop_duplicates(['col','row'], inplace=True) # rounding to integer
        
        out_meta = ds_open.meta.copy()  
        buff=3 # cell size?
#        cell_size=3
        
        arr_bankpts=np.empty([out_meta['height'], out_meta['width']], dtype=out_meta['dtype'])         
        
        for tpl_row in df_coords.itertuples():

            j+=1
            
#            print(tpl_row.linkno)
            
#            if tpl_row.linkno <> 12:
#                continue
            
#            print('pause')
            
            row_min = np.int(tpl_row.row - np.int(w_height/2))
            row_max = np.int(tpl_row.row + np.int(w_height/2))
            col_min = np.int(tpl_row.col - np.int(w_width/2))
            col_max = np.int(tpl_row.col + np.int(w_width/2))            
            
            # Now get the DEM specified by this window as a numpy array...
            w = ds_open.read(1, window=((row_min, row_max),(col_min, col_max))) 
            
            w = np.ma.masked_less(w, 0) 
            
            # Just pick out the low values (less than the Nth percentile)...
            w[w>np.percentile(w,20)] = out_meta['nodata']            
            
            
            w = np.ma.masked_less(w, 0) # for plotting?
            
#            plt.imshow(w)
#
#            sys.exit()
            
            if np.size(w) > 9: # make sure a window of appropriate size was returned from the DEM
#                            
#                
##                w_curve[w_curve<np.max(w_curve)*0.30] = out_meta['nodata']
#
#                
                arr_bankpts[row_min+buff:row_max-buff, col_min+buff:col_max-buff] = w[buff:w_height-buff, buff:w_width-buff]
#                
#
#
#
        arr_bankpts[arr_bankpts<=0.] = out_meta['nodata']            
        
        print('Writing bankpts .tif...')
        with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
            dest.write(arr_bankpts, indexes=1)
            
#        df_dist = pd.DataFrame(lst_dist)
    
    return    