# NESFCDrifterArchive

NEFSC Drifter Archive: comparison with FVCOM

Download data from the dataset:
http://www.nefsc.noaa.gov/drifter/
http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.html

generated ERDAP URL
http://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?id,time,latitude,longitude,depth,sea_water_temperature&distinct()&orderBy(%22id,time%22)
in browser download
drifters_3680_38e9_bb39.csv 
33.6MB 2017-01-27
all archive is ~30MB in 2015
drifters_archive.csv

s01_plt_drift.py - load all dataset, browse through drifters, plot tracks, velocities

test_dtr_v27.py - comparison of drifter tracks and FVCOM tracks based on dtr.py (derived from v25).

