* ies-io.f last update 4/1/98 Christine C. Molling, cmolling@facstaff.wisc.edu
* University of Wisconsin - Madison; Institute for Environmental Studies;
* Center for Climatic Research; Climate, People, and Environment Program
*
* You are free to change or distribute my code as long as you
* 1) acknowlege me, Christine C. Molling (cmolling@facstaff.wisc.edu)
* as the author of the original code and 
* 2) do not sell it.
* Of course if there is anything wrong with the code, or you somehow
* encounter damage to your reputation/wallet because of its use, you
* are forbidden to sue me or anyone else.
*
*     These subroutines can be used as a primary interface to the ies
* format netcdf files.  You can also use the slightly lower level functions
* found in ies.f or the low level netcdf commands (see the netcdf V3
* users guide).
*SUBROUTINE LIST:
* inifile - create a new file with dimensions and variables for
*  longitude, latitude, an optional 3rd real dimension, and time
*  Time must have units of days since a date!!!
* inifilec - create a new file with dimensions and variables for
*  longitude, latitude, a 3rd character-variable dimension, and time
*  Time must have units of days since a date!!!
* inivar - initialize a new variable
* endini - end initialization of file: sync file to disk, close file
* writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
* readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
*
*.....................................................................
* inifile - create a new file with dimensions and variables for
*  longitude, latitude, an optional 3rd real dimension, and time.
*  Time must have units of days since a date!!!
*.....................................................................
      subroutine inifile(idies,filen,title,source,history,nlon,alons,
     . nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,
     . calendar,ierror)

*INPUT
*     filen - character*(*) - name for new file
*     title - character*(*) - title for file
*     source - character*(*) - source for file
*     history - character*(*) - date of creation, and possibly author
*     nlon - integer - number of point along longitude direction
*     alons - real(nlon) - values of longitudes
*     nlat - integer - number of point along latitude direction
*     alats - real(nlat) - values of latitudes
*     name3rd - character*(*) - name of 3rd dimension - use '' if you
*      don't want this dimension
*     long3rd - charqcter*(*) - long name for 3rd dimension variable
*      (ignored if nam3rd='')
*     n3rd - integer - size of 3rd dimension (ignored if name3rd='')
*     vals3rd - real(n3rd) - values along 3rd dimension (ignored if 
*      name3rd='')
*     pos3rd - character*(*) - if the 3rd dimension exists and refers to
*      a measured level or height, use 'up' if values increase with distance
*      above earth (as in height), or 'down' if values decrease with distance 
*      (as in pressure) (ignored if name3rd=''). If 3rd dimension does not
*      refer to a height or level, use ''.
*     tunits - character*(*) - units for time, must be in the form of days
*      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
*      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
*     calendar - character*(*) - type of calendar.  Choose from 'noleap',
*      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
*      differ from modern.
*OUTPUT
*     idies - integer - id number of new file for later use
*     ierror - integer - error code, 0 = no error, < 0 = an error occured

      integer idies, nlon, nlat, n3rd, ierror
      character*(*) filen, title, source, history, name3rd
      character*(*) long3rd, units3rd, pos3rd, tunits, calendar
      real alons(nlon), alats(nlat), vals3rd(n3rd)

      include 'netcdf.inc'

      integer ierr, iddlon, idlon, iddlat, idlat, idd3rd, id3rd
      integer iddtime, idtime, iddlend, iddate, idtw, ll
      integer idims(2)
      character junk*180, cnull*1

      cnull = char(0)

      ierror = 0

*     Open file
*     ---------
      ierr = NF_CREATE(filen,NF_CLOBBER,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define global attributes
*     ------------------------
      ll = lenchr(title)
      junk = title(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll+1,
*     > //char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      ll = lenchr(source)
      junk = source(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll+1,
*     . source(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(history)
      junk = history(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll+1,
*     . history(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(calendar)
      junk = calendar(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll+1,
*     . calendar(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'conventions',9,
     . 'NCAR-CSM'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimensions
*     -----------------
      ierr = NF_DEF_DIM(idies,'longitude',nlon,iddlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'latitude',nlat,iddlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idd3rd = -999
      if (name3rd .ne. '') then
         ierr = NF_DEF_DIM(idies,name3rd,n3rd,idd3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_DEF_DIM(idies,'time',NF_UNLIMITED,iddtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'lengthd',10,iddlend)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimension variables
*     --------------------------
      ierr = NF_DEF_VAR(idies,'longitude',NF_FLOAT,1,iddlon,idlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'latitude',NF_FLOAT,1,iddlat,idlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ierr = NF_DEF_VAR(idies,name3rd,NF_FLOAT,1,idd3rd,id3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_DEF_VAR(idies,'time',NF_FLOAT,1,iddtime,idtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define variables time_weights and date
*     --------------------------------------
      ierr = NF_DEF_VAR(idies,'time_weights',NF_FLOAT,1,iddtime,idtw)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1)=iddlend
      idims(2)=iddtime
      ierr = NF_DEF_VAR(idies,'date',NF_CHAR,2,idims,iddate)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if


*     Attributes for dimension variables
*     ----------------------------------
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'long_name',10,
     . 'longitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'units',13,
     . 'degrees_east'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'long_name',9,
     . 'latitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'units',14,
     . 'degrees_north'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ll = lenchr(long3rd)
         junk = long3rd(1:ll)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
     >    ll+1,junk//char(0))
*         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
*     .    ll+1,long3rd(1:ll)//char(0))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ll = lenchr(units3rd)
         junk = units3rd(1:ll)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     >    ll+1,junk//char(0))
*         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
*     .    ll+1,units3rd(1:ll)//char(0))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'long_name',5,
     . 'time'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      junk = tunits(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll+1,
     > junk//char(0))
*      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll+1,
*     . tunits(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'long_name',29,
     . 'number of days per time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'units',5,
     . 'days'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'long_name',25,
     . 'label for each time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'units',1,
     . char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     End define mode, enter data mode
*     --------------------------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put dimension variables except for time
*     ---------------------------------------
      ierr = NF_PUT_VAR_REAL(idies,idlon,alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_REAL(idies,idlat,alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ierr = NF_PUT_VAR_REAL(idies,id3rd,vals3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

*     Don't close file
*     ----------------

      return
      end
*.....................................................................
* inifilec - create a new file with dimensions and variables for
*  longitude, latitude, a 3rd character-variable dimension, and time
*  Time must have units of days since a date!!!
*.....................................................................
      subroutine inifilec(idies,filen,title,source,history,nlon,alons,
     . nlat,alats,name3rd,long3rd,units3rd,n3rd,len3rd,chars3rd,
     . tunits,calendar,ierror)

*INPUT
*     filen - character*(*) - name for new file
*     title - character*(*) - title for file
*     source - character*(*) - source for file
*     history - character*(*) - date of creation, and possibly author
*     nlon - integer - number of point along longitude direction
*     alons - real(nlon) - values of longitudes
*     nlat - integer - number of point along latitude direction
*     alats - real(nlat) - values of latitudes
*     name3rd - character*(*) - name of 3rd dimension
*     long3rd - charqcter*(*) - long name for 3rd dimension variable
*     n3rd - integer - size of 3rd dimension
*     len3rd - integer length of chracter strings in vals3rd
*     chars3rd - character*len3rd(n3rd) - values along 3rd dimension
*     tunits - character*(*) - units for time, must be in the form of days
*      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
*      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
*     calendar - charater*(*) - type of calendar.  Choose from 'noleap',
*      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
*      differ from modern.
*OUTPUT
*     idies - integer - id number of new file for later use
*     ierror - integer - error code, 0 = no error, < 0 = an error occured

      integer idies, nlon, nlat, n3rd, len3rd, ierror
      character*(*) filen, title, source, history, name3rd
      character*(*) long3rd, units3rd, tunits, calendar
      real alons(nlon), alats(nlat)
      character chars3rd(len3rd,n3rd)

      include 'netcdf.inc'

      integer ierr, iddlon, idlon, iddlat, idlat, idd3rd, id3rd
      integer iddlen, iddtime, idtime, iddate, idtw, ll, idims(2)
      integer iddlend

      character junk*180, cnull*1

      cnull = char(0)

      ierror = 0

*     Open file
*     ---------
      ierr = NF_CREATE(filen,NF_CLOBBER,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define global attributes
*     ------------------------
      ll = lenchr(title)
      junk = title(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll+1,
*     . title(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(source)
      junk = source(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll+1,
     . junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll+1,
*     . source(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(history)
      junk = history(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll+1,
*     . history(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(calendar)
      junk = calendar(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll+1,
*     . calendar(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'conventions',9,
     . 'NCAR-CSM'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimensions
*     -----------------
      ierr = NF_DEF_DIM(idies,'longitude',nlon,iddlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'latitude',nlat,iddlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,name3rd,n3rd,idd3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'length',len3rd,iddlen)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'lengthd',10,iddlend)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      ierr = NF_DEF_DIM(idies,'time',NF_UNLIMITED,iddtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimension variables
*     --------------------------
      ierr = NF_DEF_VAR(idies,'longitude',NF_FLOAT,1,iddlon,idlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'latitude',NF_FLOAT,1,iddlat,idlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1) = iddlen
      idims(2) = idd3rd
      ierr = NF_DEF_VAR(idies,name3rd,NF_CHAR,2,idims,id3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'time',NF_FLOAT,1,iddtime,idtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define variables time_weights and date
*     --------------------------------------
      ierr = NF_DEF_VAR(idies,'time_weights',NF_FLOAT,1,iddtime,idtw)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1)=iddlend
      idims(2)=iddtime
      ierr = NF_DEF_VAR(idies,'date',NF_CHAR,2,idims,iddate)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Attributes for dimension variables
*     ----------------------------------
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'long_name',10,
     . 'longitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'units',13,
     . 'degrees_east'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'long_name',9,
     . 'latitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'units',14,
     . 'degrees_north'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ll = lenchr(long3rd)
         junk = long3rd(1:ll)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
     >    ll+1,junk//cnull)
*         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
*     .    ll+1,long3rd(1:ll)//char(0))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifilec'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ll = lenchr(units3rd)
         junk = units3rd(1:ll)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     >    ll+1,junk//cnull)
*         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
*     .    ll+1,units3rd(1:ll)//char(0))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifilec'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'long_name',5,
     . 'time'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      junk = tunits(1:ll)
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll+1,
     > junk//cnull)
*      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll+1,
*     . tunits(1:ll)//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'long_name',29,
     . 'number of days per time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'units',5,
     . 'days'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'long_name',25,
     . 'label for each time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'units',1,
     . char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     End define mode, enter data mode
*     --------------------------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put dimension variables except for time
*     ---------------------------------------
      ierr = NF_PUT_VAR_REAL(idies,idlon,alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_REAL(idies,idlat,alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_TEXT(idies,id3rd,chars3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Don't close file
*     ----------------

      return
      end
*.....................................................................
* inivar - initialize a new variable
*.....................................................................
      subroutine inivar(idies,varname,longname,units,ndims,dimnames,
     . valmissing,ierror)

*INPUT
*     idies - integer - id number of a new file from a previous call
*      to inifile, inifilec, or iesopen
*     varname - charcter*(*) - name for new variable
*     longname - character*(*) - descriptive name for variable
*     units - character*(*) - units for variable, compatible with udunits
*     ndims - integer - number of dimensions for this variable
*     dimnames - character*(*)(ndims) - name of each dimension, in order
*     valmissing - real - missing value, ignored if valmissing=0.
*OUTPUT
*     ierror - integer - error code, 0 = no error, < 0 = error occurred

      integer idies, ndims, ierror
      character*(*) varname, longname, units, dimnames(ndims)
      real valmissing

      include 'netcdf.inc'

      integer ierr, idvar, iddims(4)
      character junk*180, cnull*1

      cnull = char(0)

      ierror = 0

*     Put into redefine mode, but don't fail if already in define mode
*     ----------------------------------------------------------------
      ierr = NF_REDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in inivar'
         print *, NF_STRERROR(ierr)
      end if

*     Find id's of dimensions
*     -----------------------
      do i = 1, ndims
         ierr = NF_INQ_DIMID(idies,dimnames(i),iddims(i))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inivar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end do

*     Define variable
*     ---------------
      ierr = NF_DEF_VAR(idies,varname,NF_FLOAT,ndims,iddims,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define attributes
*     -----------------
      ierr = NF_PUT_ATT_TEXT(idies,idvar,'long_name',
     . lenchr(longname),longname)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      ierr = NF_PUT_ATT_TEXT(idies,idvar,'units',lenchr(units),
     . units)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (valmissing .ne. 0.) then
         ierr = NF_PUT_ATT_REAL(idies,idvar,'missing_value',NF_FLOAT,
     .    1,valmissing)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inivar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

*     Exit define mode
*     ----------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
      end

*.....................................................................
* endini - end initialization of file: sync file to disk, close file
*.....................................................................
      subroutine endini(idies, ierror)

*INPUT
*     idies - integer - id of netcdf file (generated by inifile or inifilec)
*
*OUTPUT
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

      integer idies, ierror

      include 'netcdf.inc'

      integer ierr

*     Close file, sync to disk
*     ------------------------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
      end

*.....................................................................
* writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
*.....................................................................
      subroutine writevar(filen,varname,istart,icount,values,
     . times,tweights,dates,ierror)

*INPUT
*     filen - character*(*) - file name to which to write data
*     varname - character*(*) - name of variable to which to write
*     istart - integer(*) - starting points along each dimension of
*      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
*      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
*     icount - integer(*) - number of points along each dimension to write,
*      for example, to write a single lat/lon grid to a 4-d variable,
*      icount would be (nlon,nlat,1,1).
*     values - real(*) - real values of the designated hyperslab
*     times - real(*) - time vector vector for those points written (ignored
*      if variable is 2-d).
*     tweights - real(*) - number of days in each time step in times (ignored
*      if variable is 2-d).
*     dates - character*10(*) - character labels for each time step in times,
*      last character should be a null (dates ignored if variable 2-d).
*OUTPUT
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

      character*(*) filen, varname
      character*10 dates(*)
      integer istart(*), icount(*), ierror
      real values(*), times(*), tweights(*)


      include 'netcdf.inc'

      integer idies, idvar, ndims, ierr, itype, idtime, idtw, iddate
      integer ist(2), ict(2)
      character junk*180, cnull*1

      cnull = char(0)

      ierror = 0

*     Open file, but don't fail if already open
*     -----------------------------------------
      ierr = NF_OPEN(filen,NF_WRITE,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in writevar'
         print *, NF_STRERROR(ierr)
      end if

*     Get id of variable
*     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Inquire number of dimensions
*     ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put variable
*     ------------
      ierr = NF_PUT_VARA_REAL(idies,idvar,istart,icount,values)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put time value(s), time weights, dates if necesary
*     --------------------------------------------------
      if (ndims .gt. 2) then

         ierr = NF_INQ_VARID(idies,'time',idtime)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_PUT_VARA_REAL(idies,idtime,istart(ndims),
     .    icount(ndims),times)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_INQ_VARID(idies,'time_weights',idtw)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_PUT_VARA_REAL(idies,idtw,istart(ndims),
     .    icount(ndims),tweights)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_INQ_VARID(idies,'date',iddate)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ist(1)=1
         ist(2)=istart(ndims)
         ict(1)=10
         ict(2)=icount(ndims)
         ierr = NF_PUT_VARA_TEXT(idies,iddate,ist,
     .    ict,dates)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if

      end if

*     Close file
*     ----------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if


      return
      end
*.....................................................................
* readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
*.....................................................................

      subroutine readvar(filen,varname,name3d,istart,icount,values,
     . alons,alats,vals3d,times,ierror)

*INPUT
*     filen - character*(*) - file name from which to read data
*     varname - character*(*) - name of variable from which to read
*     name3d - character*(*) - name of 3rd, nontime dimension variable.
*      Ignored if varname variable is only 3-d.
*     istart - integer(*) - starting points along each dimension of
*      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
*      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
*     icount - integer(*) - number of points along each dimension to read,
*      for example, to read in a single lat/lon grid from a 4-d variable,
*      icount would be (nlon,nlat,1,1)
*OUTPUT
*     values - real(*) - returned real values of the designated hyperslab
*     alons - real(*) - longitude vector for those points read in
*     alats - real(*) - latitude vector for those points read in
*     vals3d - real(*) - vector of values of the 3rd dimension (unchanged if
*      variable is 2- or 3-d, or if 3rd dimension has character values).
*     times - real(*) - time vector vector for those points read in (unchanged
*      if variable is 2-d).
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

      character*(*) filen, varname, name3d
      integer istart(*), icount(*), ierror
      real values(*), alons(*), alats(*), vals3d(*), times(*)


      include 'netcdf.inc'


      integer idies, idvar, ndims, idlon, idlat, id3d, ierr, itype
      character junk*180, cnull*1

      cnull = char(0)

      ierror = 0

*     Open file, but don't fail if already open
*     -----------------------------------------
      ierr = NF_OPEN(filen,NF_NOWRITE,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in readvar'
         print *, NF_STRERROR(ierr)
      end if

*     Get id of variable
*     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, idies
         print *, varname
c TODO         
         print *, idvar
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Inquire number of dimensions
*     ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Get variable
*     ------------
      ierr = NF_GET_VARA_REAL(idies,idvar,istart,icount,values)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Get values of dimension variables
*     ---------------------------------
      ierr = NF_INQ_VARID(idies,'longitude',idlon)
      ierr = NF_GET_VARA_REAL(idies,idlon,istart(1),
     . icount(1),alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_INQ_VARID(idies,'latitude',idlat)
      ierr = NF_GET_VARA_REAL(idies,idlat,istart(2),
     . icount(2),alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (ndims .ge. 3) then
         ierr = NF_INQ_VARID(idies,'time',idtime)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_GET_VARA_REAL(idies,idtime,istart(ndims),
     .    icount(ndims),times)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

      if (ndims .eq. 4) then
         ierr = NF_INQ_VARID(idies,name3d,id3d)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_INQ_VARTYPE(idies,id3d,itype)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         if (itype .ne. NF_CHAR) ierr = NF_GET_VARA_REAL(idies,id3d,
     .    istart(3),icount(3),vals3d)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

*     Close file
*     ----------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
      end

