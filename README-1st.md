README-1st.txt
last update: 6/2/04 (John Bachan)


                  ###   ######    ###    #####
                   #    #     #    #    #     #
                   #    #     #    #    #
                   #    ######     #     #####
                   #    #     #    #          #
                   #    #     #    #    #     #
                  ###   ######    ###    #####

                  Integrated BIosphere Simulator

AGREEMENT
---------

In giving you the IBIS code, we ask you to agree to the following:

1.  In publications resulting from your use of the IBIS code please
acknowledge us and cite the appropriate papers:

Foley, J.A. et al., 1996:  An integrated biosphere model of land surface
processes, terrestrial carbon balance, and vegetation dynamics.  Global
Biogeochemical Cycles, Vol.10, No.4, pages 603-628.

Kucharik, C.J. et al., 2000: Testing the performance of a dynamic global 
ecosystem model: Water balance, carbon balance and vegetation structure. 
Global Biogeochemical Cycles 14(3), 795-825.

We would also appreciate a copy of your publication so we stay informed
about the uses of the model.  Send printed copies to

Jonathan Foley
Center for Sustainability and the Global Enviroment
University of Wisconsin-Madison
1710 University Avenue
Madison, WI 53726
USA

2.  In order to keep track of where the code is, who is using which
version, and what it is being used for, we request that you not distribute
the code to any other persons or groups.  Any requests you receive for the
code can be forwarded to Jonathan Foley at (608) 265-5144, or email
jfoley@wisc.edu.  See also the SAGE web site:
http://www.sage.wisc.edu

3.  In an effort to keep the model up to date and in the most useful
form, we would like feedback from you.  If you make any substantial
changes to the code please let us know.  Your change may be one that
should be incorporated in the next version of the code.  The code is
in constant development and feedback from you will help us improve the
code for your future use as well as the use of others.  If you
experience any problems with the model or have any questions please do
not hesitate to contact us.  If you think you have found a bug, please
contact Jonathan Foley at jfoley@wisc.edu.

4.  The standard disclaimer applies.

Thank you for your cooperation.


MODES
-----

This version of IBIS can be compiled and run in only one:
* NetCDF input/output on a Unix/Linux computer

If you plan to run IBIS on a Unix/Linux platform, but are not familiar
with netCDF, we recommend that you look in the documentation of your 
visualization software to see if it is supported.

We have successfully run IBIS in netCDF mode on SGI, Sun, Linux, and
Mac OS X. If you successfully port IBIS to any other
platform or popular i/o format, we'd really like to know about it.
Please let us know if you needed to do anything special to get it to
compile and run.

To learn details of compiling/running IBIS, read
* READMEnet.txt for Unix/Linux netCDF mode
The file 'makefile' shows the file dependencies.

The file 'READMEnotes.txt' has some notes about certain parameter 
choices for the model. Please read this note, and make parameter choices 
based on your specific situation.

TEST DATA
---------

We provide some nonsense test files at 10x10 degree resolution that
you may use in order to test whether your compilation of ibis works or
not.  Please download input10-nc.tar.gz and ungzip and untar it. 
 
You'll have the 13 files you need to run ibis (make sure
your value for nanom in ibis.infile or ibismac.infile is large).
Place the files in the subdirectory 'input'.  Don't expect to have
meaningful results from using this data as ibis input - they are just
some numbers we made up.  You should not, however, crash the program,
core dump, or get NaN's (not-a-number) or Inf's (infinity) using this
input.

You may also use the test files as templates to guide you in creating
your own input data.  Due to the proprietary nature of some of the
climate data we have, we cannot freely distribute all of our input
files.  This test data, however, can be freely modified or given away
to anyone.

REAL DATA
---------
We are able to provide real data at resolutions of 0.5, 1.0, 2.0 and 4.0.
Contact us to download the 13 input files.

Good Luck!
