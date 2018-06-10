README_IBIS_INPUT.txt
last update: 6/2/04 (John Bachan)


                  ###   ######    ###    #####
                   #    #     #    #    #     #
                   #    #     #    #    #
                   #    ######     #     #####
                   #    #     #    #          #
                   #    #     #    #    #     #
                  ###   ######    ###    #####

                  Integrated BIosphere Simulator

======================================================================


IBIS INPUT FILES
-----------------------
For IBIS Version 2.5 - 2.6 
http://www.sage.wisc.edu


GLOBAL FILES
------------
Thirteen input files are required to run IBIS in the standard global
mode. These files are not readily available in **netcdf** format. We have
them available for redistribution if you agree to properly reference 
the original author. 

Contact us for access to these input files.


NETCDF
------
NetCDF (network Common Data Form) is an interface for array-oriented
data access and a library that provides an implementation of the
interface. The netCDF library also defines a machine-independent format
for representing scientific data. Together, the interface, library, and
format support the creation, access, and sharing of scientific data. The
netCDF software was developed at the Unidata Program Center in Boulder,
Colorado. The freely available source can be obtained as a compressed
tar file or a zip file from Unidata or from other mirror sites. 
http://www.unidata.ucar.edu/packages/netcdf

Version 3.3 or better is needed.


RESOLUTION
----------
The raw data we use for the global model is at a resolution of 0.5¡x 0.5¡.
We have repackaged it to 1.0¡, 2.0¡ and 4.0¡ resolutions.

You can create your own input data at any resolution. 
Some IBIS users have created input data at one-meter resolution to 
describe a crop field. 
  

THIRTEEN INPUT FILES
--------------
cld.mon.nc
deltat.mon.nc
rh.mon.nc
temp.mon.nc
prec.mon.nc
soita.clay.nc
soita.sand.nc
surta.nc
trange.mon.nc
vegtype.nc
wetd.mon.nc
wspd.mon.nc


Name    Description                             Units
----    -----------                             -----
cld     monthly mean cloudiness                 %
deltat  *(see below)                            C 
prec    monthly mean precipitation rate         mm/day
rh      monthly mean relative humidity          %
sand    percentage of sand                      %
clay    percentage of clay                      %
surta   land mask                               1=land,0=ocean
temp    monthly mean temperature                C
topo    topography                              m
trange  monthly mean temperature range          C
vegtype initial vegetation types                indexd 
wetd    mean "wet" days per month               days
wspd    monthly mean wind speed at sig=0.995    m/s

* deltat = minimum temp ever recorded at that location minus
              avg temp of coldest month

Name    Source
----    ----------- 
cld     CRU
deltat  Oregon
prec    CRU
rh      CRU
sand    IGBP/CONUS
clay    IGBP/CONUS	
surta   adapted from ETOPO5
temp    CRU
topo    ETOPO5
trange  CRU
vegtype SAGE
wetd    CRU
wspd    CRU


CRU - GLOBAL CLIMATE DATASET
----------------------------
Institute:    Climate Research Unit - East Anglia, UK
Website:      http://ipcc-ddc.cru.uea.ac.uk/cru_data/examine/cru_climate.html
Contact:      David Viner <d.viner@uea.ac.uk>
Temporal:     1961-1990 mean climatology
Spatial:      0.5 degree 
Reference: 
New M., Hulme, M and Jones, P.D, 1999: 
       Representating twentieth century space-time climate 
       variability. Part 1: development of a 1961-90 mean 
       monthly terrestrial climatology. 
       J. Climate, 12, 829-856
SAGE Comment: 
       Anomaly files also available for 1901-1996 
       (anomaly = individual year monthly value minus climatology monthly value)
       http://www.cru.uea.ac.uk/link
       SAGE cleaned, flipped and netcdf, reformated anomaly files.


OREGON - deltaT - absolute min temp minus avg of coldest mon mean temp
--------------
Institute:    Dept. Geography, Univ. Oregon, Eugene OR
Contact:      Pat Bartlein <bartlein@oregon.uoregon.edu>
Source:       The data sources were on the World Weather Disc CD:
       Globe:  Worldwide Airfield Summaries (NCDC TD-9647), 
       variable lengths of record (really variable)
       US:  Climatology of the U.S. No. 20 (1951-80)
Temporal:  
Spatial:      0.5 degree
Reference:    Unpublished personal communication
SAGE Comment: 


IGBP - GLOBAL SOILS DATASET
---------------------------
Institute:    Global Soil Data Task of the International Geosphere-Biosphere Programme 
              Data and Information System (IGBP-DIS) Potsdam, Germany. 
Website:      http://www.daac.ornl.gov/SOILS/igbp.html
Contact:      webmaster@www.daac.ornl.gov 
Temporal:  
Spatial:      5 min
Reference:    Global Soil Data Task. 2000. Global Soil Data Products CD-ROM (IGBP-DIS). 
International Geosphere-Biosphere Programme - Data and Information Services. Available
on-line [http://www.daac.ornl.gov ] from Oak Ridge National Laboratory Distributed 
Active Archive Center, Oak Ridge, Tennessee, U.S.A. 
SAGE Comment: Sage has reformatted this data to soita.clay and soita.sand for use in IBIS.


CONUS - US SOILS DATASET
------------------------
Institute:    Penn. State. Univ.
Website:      http://www.essc.psu.edu/soil_info/index.cgi?soil_data&conus 
Contact:      -- 
Temporal:     --
Spatial:      30 arc-seconds  
Reference: 
Miller, D.A. and R.A. White, 1998: A Conterminous United States Multi-Layer
       Soil Characteristics Data Set for Regional Climate and Hydrology Modeling.
       Earth Interactions, 2. [Available on-line at http://EarthInteractions.org]
SAGE Comment: Should be used when running over the US.


ETOP05  - SURTA and TOPO
------------------------
Institute:    NOAA/NGDC
Website:      http://www.ngdc.noaa.gov/mgg/global/global.html
Contact:      Peter.W.Sloss@noaa.gov
Temporal:     --
Spatial:      5 min (reprocessed to 0.5)
Reference:    Various contributors
SAGE Comment: 


SAGE - VEGTYPE
------------------------
Institute:    SAGE/UW-MADISON
Website:      http://sage.sage.wisc.edu
Contact:      Navin Ramankutty <nramanku@facstaff.wisc.edu>
Temporal:     
Spatial:      5 min
Reference:    Ramankutty, N., and J.A. Foley (1999). 
       Estimating historical changes in global land cover: croplands from 1700 to 1992. 
       Global Biogeochemical Cycles 13(4), 997-1027.
SAGE Comment: Initial Vegetation typel as needed by IBIS.



CREATING YOUR OWN DATA
----------------------
Arcview, R, PyClimate all can be used to create your own datasets.

======================================================================


