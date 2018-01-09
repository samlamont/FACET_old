# FACET
[Floodplain and Channel Evaluation Toolkit]

FACET is a standalone Python tool that uses open source modules to map the floodplain extent and compute stream channel and floodplain geomorphic metrics such as channel width, streambank height, active floodplain width, and stream slope from DEMs. 

FACET currently contains several automated methods for extracting metrics including:

Cross Section Automation
- Cross sections automatically created parallel to reach at user specified spacing.
- Geomorphic metrics are calculated at each cross section and summarized by reach.

Identifying Channel Banks using Curvature
- Stream channel banks are delineated by calculating curvature (two methods exist) along the stream channel network.
- Curvature is calculated within a moving window that traverses the stream network, followed by the application of a threshold to identify grid cells representing the bank.
- Channel width for the reach is then computed using a buffering technique applied to each reach segments at user-specified length (e.g., 100m).

Identifying Floodplain Extent using Height Above Nearest Drainage (HAND)
- The active floodplain is delineated based on the HAND technique. The HAND grid is the elevation of every grid cell relative to the nearest stream channel cell it drains to. 
- Field data is used to calibrate the HAND elevation threshold that defines floodplain extent.
- Floodplain and channel metrics can be calculated using 1D and/or 2D cross sections of the HAND grid.
- Analysis of the HAND grid also provides a method for mapping the stream channel at bankfull stage.

