# compile

Compile IBIS by typing '**make ibis**' in the directory which contains the ibis FORTRAN code.  However, before doing so, you may need to edit makefile so that the **F77_OPTIONS**, **INCLUDE_DIRS**, and **LD_OPTIONS_NETCDF** match those on your computer.  The executable will be named 'ibisnetcdf'.  You may move the executable to the directory in which your input/output directories exist.  If so, be sure to move ibis.infile and diag.infile there too.

IBIS was written and tested on SGI workstations.  We have also successfully run IBIS on a Sun, and on a PC running RedHat Linux.  To our knowledge, only two portions of the code may not be consistent with your platform. The first is the function '**date**'.  The MipsPro compiler on our SGI's returns an 8-character string containing the date.  Some Sun compilers return a 9-character string.  Space is reserved for up to 10 characters (look for the variable cdate in **io.f**).  If your computer's FORTRAN compiler does not have an intrinsic function for date, has a function of a different name, or returns something other than 10 or fewer characters, you'll have to edit io.f accordingly.  Look in wrestart, wdaily, wmonthly, and wyearly for places where you may need to change the code.  If nothing is available for you to use, you may eliminate the reference to the function 'date' and the variable cdate.  It is only used to create a history attribute in the output files.

The second portion of code which may not consistent with your platform is the command '**flush**' in **stats.f**.  This command is used to update the output in the files 'ibis.out.global' and 'ibis.out.vegtype'.  This prevents your computer from buffering output for these files.  If you do not have a version of flush in your compiler, you may comment out those lines of code.  You may not be able to read output for a particular year in ibis.out.global or ibis.out.vegtype until the program has completed a few extra years.

# input

## GLOBAL FILES
Thirteen input files are required to run IBIS in the standard global mode. These files are not readily available in **netcdf** format. We have them available for redistribution if you agree to properly reference  the original author. Contact us for access to these input files.

## NETCDF
NetCDF (network Common Data Form) is an interface for array-oriented data access and a library that provides an implementation of the interface. The netCDF library also defines a machine-independent format for representing scientific data. Together, the interface, library, and format support the creation, access, and sharing of scientific data. The netCDF software was developed at the Unidata Program Center in Boulder, Colorado. The freely available source can be obtained as a compressed tar file or a zip file from Unidata or from other mirror sites.  http://www.unidata.ucar.edu/packages/netcdf
Version 3.3 or better is needed.

## RESOLUTION
The raw data we use for the global model is at a resolution of 0.5�x 0.5�. We have repackaged it to 1.0�, 2.0� and 4.0� resolutions. You can create your own input data at any resolution.  Some IBIS users have created input data at one-meter resolution to  describe a crop field. 

## THIRTEEN INPUT FILES
8 met files and 5 site files
| Name | Filename | Description | Units | Source | Type |
| ---- | -------- | ----------- | ----- | ------ |-|
| cld     | cld.mon.nc    | monthly mean cloudiness                                                     | %              | CRU                 | Met |
| deltat  | deltat.mon.nc | minimum temp ever recorded at that location minus avg temp of coldest month | C              | Oregon              | Met |
| prec    | prec.mon.nc   | monthly mean precipitation rate                                             | mm/day         | CRU                 | Met |
| rh      | rh.mon.nc     | monthly mean relative humidity                                              | %              | CRU                 | Met |
| temp    | temp.mon.nc   | monthly mean temperature                                                    | C              | CRU                 | Met |
| trange  | trange.mon.nc | monthly mean temperature range                                              | C              | CRU                 | Met |
| wetd    | wetd.mon.nc   | mean "wet" days per month                                                   | days           | CRU                 | Met |
| wspd    | wspd.mon.nc   | monthly mean wind speed at sig=0.995                                        | m/s            | CRU                 | Met |
| sand    | soita.sand.nc | percentage of sand                                                          | %              | IGBP/CONUS          | Site |
| clay    | soita.clay.nc | percentage of clay                                                          | %              | IGBP/CONUS          | Site |
| vegtype | vegtype.nc    | initial vegetation types                                                    | indexd         | SAGE                | Site |
| surta   | surta.nc      | land mask                                                                   | 1=land,0=ocean | adapted from ETOPO5 | Site |
| topo    | topo.nc       | topography                                                                  | m              | ETOPO5              | Site |

### CRU - GLOBAL CLIMATE DATASET
|-|-|
|-|-|
|Institute:|    Climate Research Unit - East Anglia, UK|
|Website:|      http://ipcc-ddc.cru.uea.ac.uk/cru_data/examine/cru_climate.html|
|Contact:|      David Viner <d.viner@uea.ac.uk>|
|Temporal:|     1961-1990 mean climatology|
|Spatial:|      0.5 degree |
|Reference:| New M., Hulme, M and Jones, P.D, 1999:  Representating twentieth century space-time climate variability. Part 1: development of a 1961-9 mean monthly terrestrial climatology.  J. Climate, 12, 829-856|
|SAGE Comment: | Anomaly files also available for 1901-1996  (anomaly = individual year monthly value minus climatology monthly value) http://www.cru.uea.ac.uk/link SAGE cleaned, flipped and netcdf, reformated anomaly files.|


### OREGON - deltaT - absolute min temp minus avg of coldest mon mean temp
|-|-|
|-|-|
|Institute:|    Dept. Geography, Univ. Oregon, Eugene OR|
|Contact:|      Pat Bartlein <bartlein@oregon.uoregon.edu>|
|Source:|       The data sources were on the World Weather Disc CD: Globe:  Worldwide Airfield Summaries (NCDC TD-9647),  variable lengths of record (really variable) US:  Climatology of the U.S. No. 2(1951-80)|
|Temporal: | |
|Spatial: |     0.5 degree|
|Reference: |   Unpublished personal communication|
|SAGE Comment:| |


### IGBP - GLOBAL SOILS DATASET
|-|-|
|-|-|
| Institute:   | Global Soil Data Task of the International Geosphere-Biosphere Programme Data and Information System (IGBP-DIS) Potsdam, Germany. |
| Website:      |http://www.daac.ornl.gov/SOILS/igbp.html|
| Contact:      |webmaster@www.daac.ornl.gov |
| Temporal:  ||
| Spatial:  |    5 min|
| Reference: |   Global Soil Data Task. 2000. Global Soil Data Products CD-ROM (IGBP-DIS). International Geosphere-Biosphere Programme - Data and Information Services. Available on-line [http://www.daac.ornl.gov ] from Oak Ridge National Laboratory Distributed  Active Archive Center, Oak Ridge, Tennessee, U.S.A. |
| SAGE Comment:| Sage has reformatted this data to soita.clay and soita.sand for use in IBIS.|


### CONUS - US SOILS DATASET
|-|-|
|-|-|
|Institute:|    Penn. State. Univ.|
|Website:|      http://www.essc.psu.edu/soil_info/index.cgi?soil_data&conus |
|Contact:|      -- |
|Temporal:|     --|
|Spatial:|      30 arc-seconds  |
|Reference:| Miller, D.A. and R.A. White, 1998: A Conterminous United States Multi-Layer Soil Characteristics Data Set for Regional Climate and Hydrology Modeling. Earth Interactions, 2. [Available on-line at http://EarthInteractions.org]|
|SAGE Comment:| Should be used when running over the US.|


### ETOP05  - SURTA and TOPO
|-|-|
|-|-|
|Institute: |   NOAA/NGDC|
|Website: |     http://www.ngdc.noaa.gov/mgg/global/global.html|
|Contact: |     Peter.W.Sloss@noaa.gov|
|Temporal: |    --|
|Spatial: |     5 min (reprocessed to 0.5)|
|Reference: |   Various contributors|
|SAGE Comment:| |


### SAGE - VEGTYPE
|-|-|
|-|-|
|Institute:|    SAGE/UW-MADISON|
|Website:|      http://sage.sage.wisc.edu|
|Contact:|      Navin Ramankutty <nramanku@facstaff.wisc.edu>|
|Temporal:|     |
|Spatial:|      5 min|
|Reference:|    Ramankutty, N., and J.A. Foley (1999). Estimating historical changes in global land cover: croplands from 1700 to 1992.  Global Biogeochemical Cycles 13(4), 997-1027.|
|SAGE Comment: |Initial Vegetation typel as needed by IBIS.|

### CREATING YOUR OWN DATA
Arcview, R, PyClimate all can be used to create your own datasets.

## TEST DATA
We provide some nonsense test files at 10x10 degree resolution that you may use in order to test whether your compilation of ibis works or not.  Please download input10-nc.tar.gz and ungzip and untar it. 
 
You'll have the 13 files you need to run ibis (make sure your value for nanom in ibis.infile or ibismac.infile is large). Place the files in the subdirectory 'input'.  Don't expect to have meaningful results from using this data as ibis input - they are just some numbers we made up.  You should not, however, crash the program, core dump, or get NaN's (not-a-number) or Inf's (infinity) using this input.

You may also use the test files as templates to guide you in creating your own input data.  Due to the proprietary nature of some of the climate data we have, we cannot freely distribute all of our input files.  This test data, however, can be freely modified or given away to anyone.

## REAL DATA
We are able to provide real data at resolutions of 0.5, 1.0, 2.0 and 4.0.
Contact us to download the 13 input files.

Good Luck!

# run

1. Set up proper input/output directories in the location from which you will run IBIS:
     ```
     mkdir output
     mkdir output/yearly output/monthly output/daily
     mkdir input
     mkdir input/anom
     mkdir restart
     ```

2. Place your input files in input (and input/anom if you have monthly anomalies).  Input files are expected to be in netcdf format and be readable by subroutine **readvar** (ies-io.f).  Files should have the dimensions of longitude, latitude, an optional 3rd dimension ('level', usually), and time. See subroutines readit, rdanom, rdday, and inird for details.  Subroutine inird is special in that it reads an integer year from the units attribute of the time dimension. Time units for input, therefore are required to be in the form  'days since yyyy-12-31'.  Optionally, you may have a clock time after the yyyy-mm-dd portion. Since the yyyy portion of the string is converted to integer with a formatted read statement, it is essential that yyyy sits in the 12th through 15th place in the units string. Also, it is expected that the monthly anomalies begin in January of the year after yyyy (one day since yyyy-12-31) as stored in time's units attribute in all monthly anomaly files, and daily fields begin on January 1st of the year after yyyy as stored in time's units attribute in all daily anomaly files.

    *subroutine readit*:
    | input file name     | **when used** | variable name | 3rd dimension |
    | ------------------- | ------------- | ------------- | ------------- |
    | input/surta.nc      | always        | surta         | layer         |
    | input/topo.nc       | always        | topo          | level         |
    | input/vegtype.nc    | isimveg <= 1  | vegtype       | level         |
    | input/soita.sand.nc | always        | sandpct       | layer         |
    | input/soita.clay.nc | always        | claypct       | layer         |
    | input/deltat.nc     | always        | deltat        | level         |
    | input/wetd.mon.nc   | always        | wetd          | level         |
    | input/temp.mon.nc   | always        | temp          | level         |
    | input/trange.mon.nc | always        | trange        | level         |
    | input/prec.mon.nc   | always        | prec          | level         |
    | input/wspd.mon.nc   | always        | wspd          | level         |
    | input/cld.mon.nc    | always        | cld           | level         |
    | input/rh.mon.nc     | always        | rh            | level         |

    *subroutine rdanom*:
    | input file name            | **when used**     | variable name | 3rd dimension |
    | -------------------------- | ----------------- | ------------- | ------------- |
    | input/anom/temp.danom.nc   | monthly anomalies | temp          | level         |
    | input/anom/trange.danom.nc | monthly anomalies | trange        | level         |
    | input/anom/prec.danom.nc   | monthly anomalies | prec          | level         |
    | input/anom/cld.danom.nc    | monthly anomalies | cld           | level         |
    | input/anom/rhum.danom.nc   | monthly anomalies | rhum          | level         |
    | input/anom/wspd.danom.nc   | monthly anomalies | wspd          | level         |
    | input/anom/wetd.danom.nc   | monthly anomalies | wetd          | level         |

    *subroutine rdday*:
    | input file name       | **when used** | variable name | 3rd dimension |
    | --------------------- | ------------- | ------------- | ------------- |
    | input/prec.daily.nc   | daily means   | prec          | level         |
    | input/temp.daily.nc   | daily means   | temp          | level         |
    | input/trange.daily.nc | daily means   | trange        | level         |
    | input/cld.daily.nc    | daily means   | cld           | level         |
    | input/wspd.daily.nc   | daily means   | wspd          | level         |
    | input/sphum.daily.nc  | daily means   | sphum         | level         |

    *subroutine inird*:
    | input file name        | **when used**     | variable name | 3rd dimension |
    | ---------------------- | ----------------- | ------------- | ------------- |
    | input/prec.daily.nc    | daily means       | prec          | not used      |
    | input/anom/temp.mon.nc | monthly anomalies | temp          | not used      |

3. Edit **compar.h** so that the following parameters are valid for your particular input data:
    - **nlon** - integer number of grid cells in the east-west direction.
    - **nlat** - integer number of grid cells in the north-south direction.
    - **npoi** - integer number of land grid cells in the lat/lon box over which you will be running IBIS(see ibis.infile explanation of snorth, etc.).
    - **xres** - real number of degrees longitude separating the first and second grid points in the east-west direction.
    - **yres** - real number of degrees latitude separating the first and second grid points in the north-south direction.
 
4. Edit **ibis.infile** so that the following input values are to your liking (note: input values must be in the exact order mentioned below).  Unless specifically mentioned, do not change any of the values if you restart the model.
    - **irestart** - enter 0 if this is a new run, 1 if it is a continuation of a previous run.
    - **iyear0** - enter the calendar year for the very first run of the simulation.  Do not change the values of this year if you are restarting the run (irestart = 1).
    - **nrun** - the number of years you wish to run the model.  You should change this value if you are restarting the model. For example, if you initially wanted to run the model for 100 years, but had to restart during year 80, change nrun form 100 (initial value) to 21 (number of years left to complete).
    - **nanom** - number of the year for which to start reading **monthly anomalies**.  Set this value to something very large if you don't want to read in anomalies during the simulation (greater than iyear0+nrun-1).  If you will be reading in monthly anomalies, set this to the first anomaly year.  For example, if iyear0 = 1950 and you wish to start reading monthly anomalies in 1965, then set nanom equal to 1965.
    - **ndprecy** - number of the year for which to start reading daily means. Set this value to something very large if you don't want to read in daily means during the simulation (greater than iyear0+nrun-1).  If you will be reading in daily means, set this to the first year in which you want to use daily means.  For example, if iyear0 = 1950 and you wish to start reading daily fields in 1965, then set ndprecy equal to 1965.
    - **soilcspin** - enter 1 if you want to use the accelerated soil spinup, 0 if you don't.
    - **iyearout** - enter 1 if you want yearly output, 0 if not.
    - **imonthout** - enter 1 if you want monthly output, 0 if not.
    - **idailyout** - enter 1 if you want daily output, 0 if not.
    - **isimveg** - enter 0 for static vegetation (vegetation will not grow or die), 1 for dynamic vegetation, or 2 for dynamic vegetation starting from a cold start (no vegetation to start out with, just 'seeds').  Note that for isimveg equal to 0 or 1, you must have the file input/vegtype.nc to provide initial conditions for vegetation distribution.
    - **isimfire** - enter 0 for fixed fire conditions, 1 for dynamic fire.  We are not satisfied with the fire scheme, so we recommend that you don't use dynamic fire.
    - **isimco2** - enter 0 for fixed co2 concentrations (value in co2init), or 1 for ramped co2 concentrations.  If you choose 1, the co2 concentration will change each year according to the polynomial equation found in subroutine co2 (physiology.f), which approximates the change in co2 from pre-industrial times to present (this currently commented out, so you'll need to remove the c's).
    - **co2init** - initial value for co2 concentration (mol/mol) to be used in the simulation.
    - **o2init** -  initial value for o2 concentration (mol/mol) to be used in the simulation.
    - **dtime** - time step (seconds) to use for the simulation, must be an even divisor of 86400 (24 hours).
    - **idiag** - enter 0 for no diagnostic output, 1 through 10 for number of diagnostic files (see diag.infile explanation for details).
    - **snorth** - latitude of northern edge of a subset box over which you wish to run IBIS.  A subset is considered to be any rectangular (in lat/lon) region that is less than the full lat/lon grid contained in your input files.  The value of snorth will be used to calculate the index of the gridcell which contains this value. If the lat/lon boundaries of the subset box that you indicate contains fractions of gridboxes, the entire gridboxes are used.  Subsetting ONLY works for  rectangular grids - grids in which latitude is constant for all grid cells in each row and longitude is constant for all grid cells in each column.
    - **ssouth** - latitude of southern edge of a subset box over which you wish to run IBIS (see snorth).
    - **swest** - longitude of western edge of a subset box over which you wish to run IBIS (see snorth).
    - **seast** - longitude of eastern edge of a subset box over which you wish to run IBIS (see snorth).

5. If you wish, edit diag.infile. Diag.infile is used to print out selected variables for every nth time step.  For example, if you wish to know how the temperature of the top soil layer changes for every other time step in 5 separate grid cells, you would first of all enter 5 for idiag in ibis.infile (each file represents a single grid cell).

    Then select five locations and enter the latitude and longitude for each grid cell center in diag.infile's spaces for diaglat0 and diaglon0, diaglat1 and diaglon1, diaglat2 and diaglon2, diaglat3 and diaglon3, and diaglat4 and diaglon4.

    Then enter the beginning and end years for the desired time span in diag.infile's spaces for diagstart# and diagend# (where # means 0 through 4 in this example).

    Also enter the frequency for writing the diagnostic output in the space for nfreq# in ibis.infile.  For every time step, nfreq# = 1, every other means nfreq# = 2, every third means nfreq# = 3, etc.  Note that if you want to print out every nfreq steps, and that is not an even factor of the number of time steps per day, the diagnostics won't be evenly spaced between the last diagnostic write of the day and the first diagnostic write of the next day.

    Lastly, in diag.infile, edit the array at the bottom so that for each diagnostic file (column), there is a 1 in the row for the diagnostic variable(s) that you want.  In our example, you would change the 0's to 1's in columns 1 through 5 in the row for tsoi, nsoilay=1.  All other places in the array would contain zeroes.

    10 is the maximum number of diagnostic files (10 individual grid cells), and a maximum of 12 variables may be chosen for each grid cell.  If you want more than 12 variables for a single grid cell, you can request two diagnostic files which merely have the same coordinates, but different variables.  Diagnostic files may be written out over the same or different time periods, and the same or different frequency.

    Read the comments at the beginning of diag.infile for another explanation of how to use diag.infile.  Do not, however, remove those comments.  The values entered in diag.infile must be in exactly the right rows.

    Each time you run ibis with diagnostics output turned on, you may wish to remove any old diagnostic files.  Any existing diagnostic files will have output appended to them.

6. To run IBIS, type 'ibis' while in the directory which contains the ibis executable and the input/output directories.  To run ibis in the background and capture all screen output to a file, type 'ibis >& ibis.out &' instead (without the quotes).  All screen output, including system error messages, will be written to ibis.out.

1. You do not need to recompile if you do a restart.  Just edit ibis.infile.  *You must recompile if you want to run over a different number of land points (subsetting) or for a different resolution*.  In either case, you must change compar.h and therefore recompile.  You may not restart if you change the number of points, location of points, or resolution. If you run ibis a second time in the same working directory, but do not do a restart, all output files (except diagnostic files) will be overwritten.

8. *If you wish to change input or output files, read the comments at the end of **io.f** before making any changes*.  Those comments explain the syntax of read/write calls and give examples.

# parameter
## choice of 'wpudmax' parameter
In soil.f subroutine soilctrl, runoff is calculated and raing (rainfall reaching ground) gets apportioned before infiltration is calculated.  This causes a problem in which raing gets alternately assigned either mostly to puddle or mostly to runoff.  This oscillating puddle depth causes the infiltration and soil moisture in the top layer to go up and down each iteration.  A shallow puddle exacerbates the fluctuation.  You are encouraged to adjust wpudmax - the current 4.5mm is less than what a silt-loam soil can infiltrate in 1 hour.  Zobek and Onstad measured random roughness for a "smooth" soil surface = 6mm.  A very rough surface (chisel-plowed soil) is 25mm.  Your puddle depth should be larger than the definition for a "heavy" rain for your soil type.  A deeper puddle will help reduce the fluctuations in infiltration and surface soil moisture by allowing a reserve of water to exist on the surface from one time step to the next.

A heavy rain event is greater than or equal to 0.6 mm/hr for a clay soil 7 mm/hr for a silt loam soil 210 mm/hr for a sand soil You can calculate the definition of "heavy" rain for other soil types by converting zdpud to units of kg/m2 (equivalent to mm of water).

## choice of 'nsoilay' parameter
In this verion, IBIS is set up to run with a 4m soil depth. Although, we find that this gives reasonable results in our global simulation, please change this  to suit your specific region. This version is set up to have only 6 layers of soil:  of depths (variable hsoi), 0.10m,0.15m,0.25m,0.50m,1.0m,2.0m from top to bottom.  However, you can easily choose different depths by modifying 'nsoilay' and/or 'hsoi'.

# HOW TO READ/WRITE NETCDF FILES IN IBIS
Reading/writing files in ibis is done through subroutines in **ies-io.f**.  The most important concept to understand is the relationship between the locations of points in an n-dimensional array and the values of istart and icount.  In FORTRAN, values in an array are stored so that the first dimension varies the fastest.  In C, the last dimension varies the fastest.  If you were to use ncdump on a netcdf file whose variable 'mydata' has 4 dimensions (latitude, longitude, level, time), the variable would be shown as follows: mydata(time, level, latitude, longitude) because ncdump was written in C and reflects C's ordering of dimensions. In this example time varies the the slowest and longitude the fastest. If you were to define this variable in a FORTRAN program so that time again varies the slowest and longitude the fastest,（高维数组的存放次序可由二维数组类推，即最右边指标（相当于最外层循环变量）变动最慢，最左边的第一个指标变动最快） it would look like this:
    real mydata(nlons, nlats, nlevels, ntimes)
where nlons, nlats, nlevels, and ntimes are integer parameters of some specified size.  Looping through the data in the order of storage in memory would look like this:
    do l = 1, ntimes
      do k = 1, nlevels
        do j = 1, nlats
          do i = 1, nlons
            mydata(i,j,k,l) = ...
          enddo
        enddo
      enddo
    enddo
Since ies-io.f is FORTRAN code, keep in mind the FORTRAN representation will be flipped in comparison to what you see from ncdump.
The netcdf interface reads and writes data using two integer vectors, istart, and icount.  Istart refers to the starting point for reading/writing along each dimension.  If you have a 4-d variable as above and want to write starting at the first point, istart would be a vector of length 4 and have the values (1,1,1,1). Icount refers to how far along each dimension data will be read/written.  If you wish to write the entire array in the sample FORTRAN code above, icount would be a vector of length 4 and have the values (nlons,nlats,nlevels,ntimes).
Things get a little more complicated when you want to read/write only a portion of the variable.  In ibis we typically read/write a single time step, but possibly more than one level/pft and either an entire or partial lat/lon grid. Below are examples of istart and icount for various situations of reading/writing.
1) an entire lat/lon grid for only one (6th) time step and one (2nd) pft: istart is (1,1,2,6), icount is (nlons,nlats,1,1)
2) entire lat/lon arrays for all 9 pfts at one (6th) time step: istart is (1,1,1,6), icount is (nlons,nlats,npfts,1)
3) a single lat/lon point (at the ilon-th longitude and the ilat-th latitude) of a 3-d variable at all 100 time steps: istart is (ilon,ilat,1), icount is (1,1,100) Note that if istart and icount have been declared length 4, the 4th value in each is ignored when referencing a 3-d variable.
4) a subsection of a lat/lon grid, 20 points wide in longitude, 15 points high in latitude, starting at point (ilon,ilat), at one (18th) level and 12 times, beginning at time step itime: istart is (ilon,ilat,18,itime) icount is (20,15,1,12)

HOW TO ADD NEW CODE TO READ A FILE:
To read a file, use subroutine readvar in ies-io.f. This subroutine assumes that the variable being read has dimensions called longitude, latitude, possibly time, and possibly another dimension, which is named in the call.  Only the bare essentials are returned.

General call:
    call readvar(filen,varname,name3d,istart,icount,values,
   > alons,alats,vals3d,times,ierror)

INPUT
    filen - character*(*) - file name from which to read data
    varname - character*(*) - name of variable from which to read
    name3d - character*(*) - name of 3rd, nontime dimension variable.
     Ignored if varname variable is only 3-d.
    istart - integer(*) - starting points along each dimension of
     the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
     and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
    icount - integer(*) - number of points along each dimension to read,
     for example, to read in a single lat/lon grid from a 4-d variable,
     icount would be (nlon,nlat,1,1)
OUTPUT
    values - real(*) - returned real values of the designated hyperslab
    alons - real(*) - longitude vector for those points read in
    alats - real(*) - latitude vector for those points read in
    vals3d - real(*) - vector of values of the 3rd dimension (unchanged if
     variable is 2- or 3-d, or if 3rd dimension has character values).
    times - real(*) - time vector vector for those points read in (unchanged
     if variable is 2-d).
    ierror - integer - error code = 0 if no error, < 0 if error occurred

Example: Read in an entire lat/lon grid for one (9th) pft (dimension 
name is 'pft') and one (24th) time
    parameter (nlons = 360, nlats=180)
    real x(nlons,nlats), alons(nlons), alats(nlats), xjunk, time
    integer istart(4), icount(4), ierr
    istart(1) = 1
    istart(2) = 1
    istart(3) = 9
    istart(4) = 24
    icount(1) = nlons
    icount(2) = nlats
    icount(3) = 1
    icount(4) = 1
    call readvar('myfile.nc','myvar','pft',istart,icount,x,alons,alats,
   > xjunk,time,ierr)
    if (ierr .lt. 0) then
      print *, 'Error occurred in readvar'
    end if
Note that above, I used a junk variable for the values of the 3rd dimension, because in ibis, the pft dimension is a character dimension.  If the 3rd dimension type is character, readvar does not return the value.

Ibis-specific example: read in a variable which has 3rd dimension 'layer', using cdummy to store full array before extracting land points into the variable used directly in ibis.  The variable work is used to store returned info that won't be used later.
Previously declared
     istart(1) = 1        !begin at 1st longitude point
     istart(2) = 1        !begin at 1st latitude point
     istart(3) = 1        !begin at 1st layer
     istart(4) = 1        !begin at 1st time
     icount(1) = nlon     !read nlon points along longitude dimension
     icount(2) = nlat     !read nlat points along latitude dimension
Read in data in file input/soita.nc to ibis variable tex
     icount(3) = nsoilay  !read nsoilay points along layer dimension
     icount(4) = 1        !read one time step
     filen = 'input/soita.nc'
     aname = 'tex'
     call readvar(filen,aname,'layer',istart,icount,
    > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
     if (istat.lt.0) then
        write(*,9000)
        print *, 'while reading soita'
        stop 1
     end if
For each lat/lon grid, strip off ocean points. tex only has land points.
     do 12 j = 1, nsoilay
       call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), tex(1,j))
12   continue

HOW TO ADD NEW CODE TO WRITE TO A FILE:
Writing to a file involves several steps.  First, if the file does not yet exist, you must create the file using inifile or inifilec. Then, you must create one or more variables using inivar.  Then after all initializing is done, declare the end of the initialization phase using endini.
Finally, once the file and variable exist, you can write data into the variable using writevar.
Creating and writing to files in seperate steps may seem complicated at first, but this 4-step method allows much flexibility. Not only do you have the choice of creating a variable without a 3rd dimension, with a real 3rd dimension, or with a character 3rd dimension, you may also have more than one variable within the same file.
Note that when using ies-io.f routines if you have more than one variable in the file, and both have 4 dimensions, that 3rd dimension (level, pft, etc.) must be shared by both variables.  For example, you may have one variable whose 3rd dimension is pft and another variable which does not have a 3rd non-time dimension in the same file, but you may not have in the same file one variable whose 3rd dimension is pft and another variable whose 3rd dimension is level.

1) Initialize a file
If you wish a file to contain only latitude, longitude, and time dimensions, or want it to also have a 3rd real dimension, use inifile.  If you want to the file to have a 3rd character dimension, use inifilec.
General call for inifile:
     call inifile(idies,filen,title,source,history,nlon,alons,
    > nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,
    > calendar,ierror)

INPUT
    filen - character*(*) - name for new file
    title - character*(*) - title for file
    source - character*(*) - source for file
    history - character*(*) - date of creation, and possibly author
    nlon - integer - number of point along longitude direction
    alons - real(nlon) - values of longitudes
    nlat - integer - number of point along latitude direction
    alats - real(nlat) - values of latitudes
    name3rd - character*(*) - name of 3rd dimension - use '' if you
     don't want this dimension ('' is the empty string)
    long3rd - character*(*) - long name for 3rd dimension variable
     (ignored if nam3rd='')
    n3rd - integer - size of 3rd dimension (ignored if name3rd='')
    vals3rd - real(n3rd) - values along 3rd dimension (ignored if 
     name3rd='')
    pos3rd - character*(*) - if the 3rd dimension exists and refers to
     a measured level or height, use 'up' if values increase with distance
     above earth (as in height), or 'down' if values decrease with distance 
     (as in pressure) (ignored if name3rd=''). If 3rd dimension does not
     refer to a height or level, use ''.
    tunits - character*(*) - units for time, must be in the form of days
     since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
     month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00.
    calendar - character*(*) - type of calendar.  Choose from 'noleap',
     'gregorian','n ka BP', etc.  Use iescal (in ies.f) if orbital 
     parameters differ from modern.
OUTPUT
    idies - integer - id number of new file for later use
    ierror - integer - error code, 0 = no error, < 0 = an error occured

Example: Create a file which will hold fractional snow
cover (only 3-d) and snow thickness for each layer (4-d).
previously defined: 

    parameter (nlons = 360, nlats = 180, nsnolay = 6)
    real lonscale(nlons), latscale(nlats), slayers(nsnolay)
    real snowc(nlons,nlats), snowh(nlons,nlats,nsnolay)
    integer istat
    character*80 cdate, tunits
    alats = ...   ! create values for latitude somehow
    alons = ...   ! create values for longitude somehow
    slayers = ... ! create values for snow layers somehow
    cdate = 'created on 6/29/97'
    tunits = 'days since 1969-12-31'
Now initialize a file with a 3rd real variable
    filen = 'snow.nc'
    call inifile(idies,filen,
   > 'file for snow variables',
   > 'C Molling, program.f v1.01',cdate,nlons,lonscale,nlats,latscale,
   > 'snowlayer','snow layers top to bottom','',nsnolay,slayers,
   > 'down',tunits,'gregorian',istat)
    if (istat .lt. 0)then
      print *, 'Error in inifile'
    end if
Note above, the empty string ('') is used because the 3rd dimension, snowlayer, does not have any units.  The 'positive' attribute for the 3rd dimension is 'down' because snowlayer 1 is at the top and the last layer is at the bottom (distance above the center of the earth is decreasing as snowlayer increases).  The calendar is gregorian because there are leap years included.  If all years are 365 days, you should use 'noleap'.  Units of time are days dince a date of the form yyyy-mm-dd.  Time units should use this form to be compatible with GrADS.  Other units should be compatible with Udunits (see www.unidata.ucar.edu).  The returned variable idies will be used in the subroutine that initializes a variable.

General call for inifilec
    call inifilec(idies,filen,title,source,history,nlon,alons,
   > nlat,alats,name3rd,long3rd,units3rd,n3rd,len3rd,chars3rd,
   > tunits,calendar,ierror)

INPUT
    filen - character*(*) - name for new file
    title - character*(*) - title for file
    source - character*(*) - source for file
    history - character*(*) - date of creation, and possibly author
    nlon - integer - number of point along longitude direction
    alons - real(nlon) - values of longitudes
    nlat - integer - number of point along latitude direction
    alats - real(nlat) - values of latitudes
    name3rd - character*(*) - name of 3rd dimension
    long3rd - charqcter*(*) - long name for 3rd dimension variable
    n3rd - integer - size of 3rd dimension
    len3rd - integer length of chracter strings in vals3rd
    chars3rd - character*len3rd(n3rd) - values along 3rd dimension
    tunits - character*(*) - units for time, must be in the form of days
     since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
     month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
    calendar - charater*(*) - type of calendar.  Choose from 'noleap',
     'gregorian','n ka BP', etc.  Use iescal if orbital parameters
     differ from modern.
OUTPUT
    idies - integer - id number of new file for later use
    ierror - integer - error code, 0 = no error, < 0 = an error occured

Example: Create a file that has a character 3rd dimension.
Defined previously - most variables as in previous example, plus...
     parameter (npft = 9)
     character*80 pftdef(npft)
     pftdef(1) = 'boreal evergreens'
     pftdef(2) = ...
     etc...
Initialize file with inifilec
     filen = 'exist.nc'
     call inifilec(idies,filen,
    > 'annual existence of each plant functional type',
    > 'ibis wyearly',cdate,nlon,lonscale,nlat,latscale,
    > 'pft','plant fuctional type','',npft,80,pftdef,
    > tunits,'gregorian',istat)
The '80' above refers to the length of the character strings in pftdef.
        dimnames(3) = 'pft'
        call inivar(idies,'exist','existence for each pft',
    >    '',4,dimnames,OCEAN,istat)
End initialization phase
        call endini(idies,istat)

2) Initialize a variable and end initialization phase After you initialize a file, you need to initialize a variable. Initializing reserves space for data to be written into the file, and saves information about the data (such as units, descriptive name, and a missing value).  You can use inivar to initialize any variable of 1 or more dimensions, provided that those dimensions already exist in the file.  You may initialize more than one variable in a single file, as stated above.
If you have a special value that denotes 'missing', such as using 1.e+36 to denote ocean grid cells, use this value for valmissing.  If you use 0. for valmissing, it is ignored.  Pick a value for valmissing which is well outside the valid range for the data.  For example, -99. would be a fine missing value for pressure, but not a good value for topography.
The character array dimnames is used to store the names of the dimensions upon which the new variable will depend. For example, a 3-d variable would only depend on longitude, latitude, and time.  A 4-d variable would depend on those and also pft, or snowlater, or level, or some other dimension.  The dimension names must be in the same order of varying: fastest to slowest.  See the discussions for istart and icount above.  For example, the 4-d variable to go in the file snow.nc (above inifile example)  would have dimnames as follows
     dimnames(1) = 'longitude'
     dimnames(2) = 'latitude'
     dimnames(3) = 'snowlayer'
     dimnames(4) = 'time'

General call for inivar:
     call inivar(idies,varname,longname,units,ndims,dimnames,
    > valmissing,ierror)

INPUT
    idies - integer - id number of a new file from a previous call
     to inifile, inifilec, or iesopen
    varname - charcter*(*) - name for new variable
    longname - character*(*) - descriptive name for variable
    units - character*(*) - units for variable, compatible with udunits
    ndims - integer - number of dimensions for this variable
    dimnames - character*(*)(ndims) - name of each dimension, in order
    valmissing - real - missing value, ignored if valmissing=0.
OUTPUT
    ierror - integer - error code, 0 = no error, < 0 = error occurred

Example: Define a 4-d and a 3-d variable to go into the file snow.nc.
defined previously
       character*80 dimnames(4)
       real OCEAN
       dimnames(1) = 'longitude'
       dimnames(2) = 'latitude'
       OCEAN = 9.e+20
define 4-d variable
       dimnames(3) = 'snowlayer'
       dimnames(4) = 'time'
       call inivar(idies,'hsno','snow layer thickness',
    >   'meters',4,dimnames,OCEAN,istat)
define 3-d variable in same file
       dimnames(3) = 'time'
       call inivar(idies,'snowf','snow cover fraction',
    >   'fraction',3,dimnames,OCEAN,istat)
end initialization phase
        call endini(idies,istat)
Notice that you need to change the value of dimnames(3) for the 3-d
variable.  The value in dimnames(4) gets ignored.

3) End the initialization phase for the file When you are done initializing the file and the variables in the file, you must call endini.  This subroutine merely closes the file.  But by closing the file, this is a signal to the computer to write out all changes to this file, synchronizing the instructions in the program with what is written on the hard disk.  Since most computers buffer, that is save up, data until a there is a large chunk to write, it is essential that all buffered information for a netcdf file be written to disk before attempting to write data to the file.  The subroutine endini accomplishes this.  The subroutine endini should be called after initializing the last variable in the netcdf file (inivar) and before calling writevar.  See the example above.

4) Write to a variable
After you have completed the initialization phase by initializing the file and initializing one or more variables in the file, you can write data  into the space reserved for each variable with subroutine writevar. You need to supply istart and icount, just as when you read a variable.  It is perfectly legal to write values in only a portion of the space reserved for the variable.  If nothing is written to a particular grid point, it has a special fill value supplied by netcdf when the variable was initialized. Notice that for ease of use, I have designed writevar to also write in the values of the time steps you are writing as well as the weighting for that time step. Time weighting refers to the number of days in that time sample. For example, some monthly means will be a mean value over 31 days, while others will be over only 30, or 28, or 29 days.  I assumed that most persons will calculate data and write it out one time step at a time, so it is nice to write out the time values and weight as you go along.  Another time-related item written is the character label for the time step.  Since after a few years, it's hard to keep track of what month is represented by a large number of days, I invented the date variable.  Date is a 10-character long string in which you can put a time label.  For example the time value 4745 may not be informative, but the data label 'JAN1913   ' is.  I usually use 9 characters plus a null character at the end (char(0)).  Use strings like 'DJF001005 ' to denote time averages (Winter years 1 through 5).

General call
    call writevar(filen,varname,istart,icount,values,
   > times,tweights,dates,ierror)

INPUT
    filen - character*(*) - file name to which to write data
    varname - character*(*) - name of variable to which to write
    istart - integer(*) - starting points along each dimension of
     the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
     and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
    icount - integer(*) - number of points along each dimension to write,
     for example, to write a single lat/lon grid to a 4-d variable,
     icount would be (nlon,nlat,1,1).
    values - real(*) - real values of the designated hyperslab
    times - real(*) - time vector vector for those points written (ignored
     if variable is 2-d).
    tweights - real(*) - number of days in each time step in times (ignored
     if variable is 2-d).
    dates - character*10(*) - character labels for each time step in times,
     last character should be a null (dates ignored if variable 2-d).
OUTPUT
    ierror - integer - error code = 0 if no error, < 0 if error occurred

Example: write data to the 4-d variable initialized above
previously defined
    parameter (ndim3=nlons*nlats*nsnolay)
    real cdummy(ndim3), hsno(npts,nsnolay)
    real time, timewght
    integer istart(4), icount(4)
    character*11 cdate
    time = 730.
    timewght = 365.
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = 2
    icount(1) = nlons
    icount(2) = nlats
    icount(3) = nsnolay
    icount(4) = 1
    cdate = 'ANN1980  '//char(0)
put land-only data in hsno into land and sea lat/lon grid in workspace
    do 20 k = 1, nsnolay
       call vec2arr (hsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
20  continue
    call writevar('snow.nc','hsno',istart,icount,cdummy,time,
   > timewght,cdate,istat)
    if (istat .ne. 0) then
       write(*,*) 'ERROR writing hsno'
       stop 1
    end if

Some other thoughts:
It's a good idea to end all strings stored in netcdf file with a null character.  Some less-robust C programs may fail when encountering non-null-ended strings, because in C, all strings end with a null character.
Udunits provides several ways to express multiplication, powers, etc.  I like to use ^ to denote exponents (as in m^2 instead of just m2) because programs like Matlab see ^ as a LaTeX command to write the next character as a superscript.  Likewise I try to avoid using _, because this is the LaTeX command to make the next character a subscript.  So instead of using deg_C, I use degC, to prevent the C from becoming a subscript.
The format of the netcdf files written by ies-io.f subroutines is meant to be compatible with the COARDS and CSM conventions (see www.cgd.ucar.edu:80/csm/experiments/output/format.html). Theoretically, you should be able to use NCO (see www.cgd.ucar.edu:80/cms/nco/index.html) and NCL (see ngwww.ucar.edu/ngdoc/ng4.1alpha/ug/ncl/ncloview.html) on the files.
A few other packages that work with this format are GrADS (grads.iges.org/grads/head.html) and Ncview (meteora.ucsd.edu/~pierce/ncview_home_page.html).
The routines in ies-io.f are just the beginning.  You can put many more things in a netcdf file by using routines in ies.f or by using the low-level netcdf commands directly.  Use one of the pre-existing subroutines in ies-io.f as a template and go on from there.  You are free to change or distribute my code as long as you 1) acknowlege me, Christine C. Molling (cmolling@facstaff.wisc.edu) as the author of the original code and 2) do not sell it.
Of course if there is anything wrong with the code, or you somehow encounter damage to your reputation/wallet because of its use, you are forbidden to sue me or anyone else.
Good Luck!
