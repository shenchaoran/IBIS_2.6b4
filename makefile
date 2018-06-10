VERSION = 2.6b4

# SGI                       ** make ibis **
#
# COMPILER = f77
# F77_OPTIONS = -O2 -col120    
#  WARNING - DO NOT USE -Ofast OPTION DURING COMPILING
#  Useful options for debugging on SGI:
#  F77_OPTIONS = -col120 -c -O0 -g -C -trapuv
#INCLUDE_DIRS  = -I/usr/local/include
#LD_OPTIONS_NETCDF  = -L/usr/local/lib -lnetcdf


# OS X ppc with g77         ** make ibis **
#
  COMPILER = gfortran
#   debug mode
# F77_OPTIONS =   -g -ffixed-line-length-132  
#  
  F77_OPTIONS = -C -ffixed-line-length-132 -funroll-loops
# 
# with local include and netcdf files.
  INCLUDE_DIRS  = -I/usr/local/include
  LD_OPTIONS_NETCDF  = /usr/local/lib/libnetcdff.so
#  LD_OPTIONS_NETCDF  = /usr/lib/x86_64-linux-gnu/libnetcdf.so


# OS X ppc with absoft f77   ** make ibisabs **
#
# COMPILER = f77
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
# COMPILER = g77 -fdebug-kludge -mcpu=pentiumpro -ffixed-line-length-132
# F77_OPTIONS = -mcpu=pentiumpro -ffast-math -malign-double -O5 -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -ffast-math -malign-double -O5 -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -O1 -malign-double -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -mpentiumpro -O2 -malign-double -ffixed-line-length-132 -funroll-loops
# F77_OPTIONS = -g -ffixed-line-length-132
#
#INCLUDE_DIRS  = -I/usr/local/include
#LD_OPTIONS_NETCDF  = -L/usr/local/lib -lnetcdf

SPATH = ./src
OBJPATH = ./build
DPATH = ./debug

HEADERS = $(wildcard $(SPATH)/*.h)
SRC = $(wildcard $(SPATH)/*.f)
SRCNAME = $(notdir $(SRC))
OBJS = $(patsubst %.f, $(OBJPATH)/%.o, $(SRCNAME))

ibis: $(OBJS) $(HEADERS)
	@mkdir -p $(DPATH)
	$(COMPILER) $(OBJS) $(F77_OPTIONS) $(INCLUDE_DIRS) $(LD_OPTIONS_NETCDF) -o $(DPATH)/ibis

$(OBJS): $(OBJPATH)/%.o: $(SPATH)/%.f $(HEADERS)
	@mkdir -p $(OBJPATH)
	$(COMPILER) $(F77_OPTIONS) $(INCLUDE_DIRS) $(LD_OPTIONS_NETCDF) -c $< -o $@

all:
	@echo $(OBJS)

clean:
	rm -rf $(OBJPATH)/*.o $(DPATH)/ibis

#------------------------------------------------------------------------
# To put all necessary files in one file ready for distribution,
# make sure the value for VERSION (at top of this file) is correct.
# then, type:
# make tardist
# you may then gzip or compress the file ibis_VERSION.tar
# If you make changes to the names of .f files (above), or add more .f
# files, you must change/add to the list below.

# tardist:
# 	tar -cjf ibis_$(VERSION).tar.bz2 \
# 	README-1st.txt READMEnet.txt \
# 	README_IBIS_INPUT.txt READMEnotes.txt history.txt \
# 	diag.infile ibis.infile \
# 	makefile $(HEADERS) \
# 	biogeochem.f canopy.f climate.f diagnostics.f ies-io.f initial.f \
# 	io.f  io-absoft.f main.f math.f physiology.f radiation.f readpars.f snow.f \
# 	soil.f stats.f utilities.f vegetation.f weather.f 
  

