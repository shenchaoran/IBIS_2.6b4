VERSION = 2.6b4

COMFILES = com1d.h    comage.h   comatm.h   combcs.h   combgc.h   comhyd.h   \
           commac.h   compar.h   comsat.h   comsno.h   \
           comsoi.h   comsum.h   compft.h   comtex.h   comveg.h   comwork.h  \
           comdiag.h  implicit.h

# SGI                       ** make ibis **
#
# F77 = f77
# F77_OPTIONS = -O2 -col120    
#  WARNING - DO NOT USE -Ofast OPTION DURING COMPILING
#  Useful options for debugging on SGI:
#  F77_OPTIONS = -col120 -c -O0 -g -C -trapuv
#INCLUDE_DIRS  = -I/usr/local/include
#LD_OPTIONS_NETCDF  = -L/usr/local/lib -lnetcdf


# OS X ppc with g77         ** make ibis **
#
  F77 = gfortran
#   debug mode
# F77_OPTIONS =   -g -ffixed-line-length-132  
#  
  F77_OPTIONS = -C -ffixed-line-length-132 -funroll-loops
# -m32
# 
# with local include and netcdf files.
  INCLUDE_DIRS  = -I/usr/local/include
  LD_OPTIONS_NETCDF  = /usr/local/lib/libnetcdff.so
#  LD_OPTIONS_NETCDF  = /usr/lib/x86_64-linux-gnu/libnetcdf.so


# OS X ppc with absoft f77   ** make ibisabs **
#
# F77 = f77
#   debug mode
# F77_OPTIONS =   -g -W
#       
# F77_OPTIONS =  -O -N2 -W -s -N11 -N15 -f -N90  
#
# with local include and netcdf files.
# INCLUDE_DIRS  = -I./include
# LD_OPTIONS_NETCDF  = ./lib/libnetcdf.a



# Linux on a PIII using g77   ** make ibis **
#
# F77 = g77 -fdebug-kludge -mcpu=pentiumpro -ffixed-line-length-132
# F77_OPTIONS = -mcpu=pentiumpro -ffast-math -malign-double -O5 -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -ffast-math -malign-double -O5 -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -O1 -malign-double -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -O2 -malign-double -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -g -ffixed-line-length-132
#
#INCLUDE_DIRS  = -I/usr/local/include
#LD_OPTIONS_NETCDF  = -L/usr/local/lib -lnetcdf

ibis: biogeochem.o canopy.o climate.o diagnostics.o ies-io.o \
      initial.o io.o main.o math.o physiology.o radiation.o \
      readpars.o snow.o soil.o stats.o utilities.o vegetation.o weather.o \
      $(COMFILES)
        

ibisabs: biogeochem.o canopy.o climate.o diagnostics.o ies-io.o \
         initial.o io-absoft.o main.o math.o physiology.o radiation.o \
         readpars.o snow.o soil.o stats-abs.o utilities.o vegetation.o weather.o \
         $(COMFILES)


biogeochem.o:	biogeochem.f    $(COMFILES)
canopy.o:	canopy.f        $(COMFILES)
climate.o:	climate.f       $(COMFILES)
diagnostics.o:	diagnostics.f   $(COMFILES)
ies-io.o:       ies-io.f        $(COMFILES)
initial.o:	initial.f       $(COMFILES)
io.o:		io.f            $(COMFILES)
io-absoft.o:	io-absoft.f     $(COMFILES)
main.o:		main.f          $(COMFILES)
math.o:		math.f          $(COMFILES)
physiology.o:	physiology.f    $(COMFILES)
radiation.o:	radiation.f     $(COMFILES)
readpars.o:     readpars.f      $(COMFILES)
snow.o:		snow.f          $(COMFILES)
soil.o:		soil.f          $(COMFILES)
stats.o:	stats.f         $(COMFILES)
stats-abs.o:	stats-abs.f     $(COMFILES)
utilities.o:	utilities.f     $(COMFILES)
vegetation.o:	vegetation.f    $(COMFILES)
weather.o:	weather.f	$(COMFILES)

ibis:
	$(F77) biogeochem.o canopy.o climate.o diagnostics.o ies-io.o \
	initial.o io.o main.o math.o physiology.o radiation.o \
	readpars.o snow.o soil.o stats.o utilities.o vegetation.o weather.o \
	$(F77_OPTIONS) $(INCLUDE_DIRS) $(LD_OPTIONS_NETCDF) -o ibis

ibisabs: 
	$(F77) biogeochem.o canopy.o climate.o diagnostics.o ies-io.o \
	initial.o io-absoft.o main.o math.o physiology.o radiation.o \
	readpars.o snow.o soil.o stats-abs.o utilities.o vegetation.o weather.o \
	$(F77_OPTIONS) $(INCLUDE_DIRS) $(LD_OPTIONS_NETCDF) \
	-o ibisabs


.f.o:
	$(F77) $(F77_OPTIONS) $(INCLUDE_DIRS) -c -o $*.o $*.f

clean:
	rm -f *.o ibis ibisabs ibis_$(VERSION).tar

#------------------------------------------------------------------------
# To put all necessary files in one file ready for distribution,
# make sure the value for VERSION (at top of this file) is correct.
# then, type:
# make tardist
# you may then gzip or compress the file ibis_VERSION.tar
# If you make changes to the names of .f files (above), or add more .f
# files, you must change/add to the list below.

tardist:
	tar -cjf ibis_$(VERSION).tar.bz2 \
	README-1st.txt READMEnet.txt \
	README_IBIS_INPUT.txt READMEnotes.txt history.txt \
	diag.infile ibis.infile \
	makefile $(COMFILES) \
	biogeochem.f canopy.f climate.f diagnostics.f ies-io.f initial.f \
	io.f  io-absoft.f main.f math.f physiology.f radiation.f readpars.f snow.f \
	soil.f stats.f utilities.f vegetation.f weather.f 
  

