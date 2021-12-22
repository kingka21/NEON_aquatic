# NEON_aquatic
This repository uses NEON data to provide case studies that examine spatial patterns and temporal variability in stream water chemistry as well as predictions of the River Continuum Concept, including metabolism, carbon chemistry, and macroinvertebrate community composition.   


The code supports all analyses in the manuscript **Jennifer W. Edmonds, Katelyn B.S. King, Merrie Beth Neely, Robert Hensley, Keli J. Goodman, and Kaelin Cawley. Using large, open-datasets to understand spatial and temporal patterns in lotic ecosystems: NEON case studies. Ecosphere.**   


**The 'Code' folder includes the following:**  
01_download_grab_data that includes methods for downloading the stream surface water chemistry grab sample data
02_plot_surfacewater_averages that includes methods for plotting and investigating water chemistry average concentrations 
03_PCA_analysis that includes methods for the PCA and HCA analysis that are presented in case study 1
04_discharge_converstion_level1to2 that includes methods for converting the level 1 discharge measures to level 2 discharge 
05_grab_discharge_analysis that includes methods used for the concentration-discharge (C-Q) analysis using the grab discharge measures presented in case study 2 
06_continuous_discharge_analysis that includes methods used for the C-Q regression analysis using the continuous discharge measures presented in case study 2
07_metabolism_analysis that includes methods used for the stream metabolism analysis presented in case study 3 

**The 'Data' folder includes the following:**  
Discharge_Level2_QAQC.csv that includes the level 2 grab discharge measures 
neon_sites_latlon.csv that includes the lat and lon of NEON stream and river sites 
surface_water_grab_QCQC_ph.csv that includes the stream water chemistry grab data, including pH from the Domain lab 
