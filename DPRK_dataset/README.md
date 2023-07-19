This directory contains the following files:

(1) **DPRK_CC_times.txt**  

This file contains the differential time measurements provided in the supplementary data to Gibbons et al. (2017),
but converted to the same format as the time difference file provided for the Finland dataset in Gibbons et al. (2020).

(2) **DPRK_phaselist.txt**  
  
This file just contains a list with  
**station** **phase** **station_latitude** **station_longitude**  
for all phases in the DPRK_CC_times.txt file.  

(3) **calculate_DPRK_theoretical_slovecs.sh**  

A script which calculates the outgoing slowness vectors at the North Korean nuclear test site from the phases
provided in the file DPRK_phaselist.txt using the AK135 velocity model (Kennett et al., 1995).  

(4) **DPRK_ak135_slovecs_expected.txt**  

The anticipated outcome of running the calculate_DPRK_theoretical_slovecs.sh script.
This is in the correct format for use with the VIScERaL_GD.py program.  

(5) **Yang_etal_initial_locations.txt**  

The best fit location estimates from the study of Yang et al. (2021) for the first five DPRK nuclear tests.
Taken from Table 1 of Yang et al. (2021).  


**References**  

Gibbons, S.J., Pabian, F., Näsholm, S.P., Kværna, T., and Mykkeltveit, S. 2017.  
Accurate relative location estimates for the North Korean nuclear tests using empirical slowness corrections,  
Geophys J Int, 208, 101-117 https://doi.org/10.1093/gji/ggw379  

Gibbons, S.J., Kvaerna, T., Tiira, T., Kozlovskaya, E., 2020.  
A benchmark case study for seismic event relative location,  
Geophys J Int, 223, 1313-1326 https://doi.org/10.1093/gji/ggaa362  

Kennett, B.L.N. Engdahl, E.R. & Buland R., 1995.  
Constraints on seismic velocities in the Earth from travel times,  
Geophys J Int, 122, 108-124  

Yang, G. and Zhao, L. F. and Xie, X. B. and Zhang, L. and Yao, Z. X., 2021.  
High‐Precision Relocation With the Burial Depths of the North Korean Underground
Nuclear Explosions by Combining Pn and Pg Differential Traveltimes,  
Journal of Geophysical Research: Solid Earth, 126, DOI 10.1029/2020JB020745  
