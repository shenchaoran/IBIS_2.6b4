c
c #####   #####    ##    #####   #####     ##    #####    ####   
c #    #  #       #  #   #    #  #    #   #  #   #    #  #      
c #    #  #      #    #  #    #  #    #  #    #  #    #   ####       
c #####   #####  ######  #    #  #####   ######  #####        #  
c #   #   #      #    #  #    #  #       #    #  #   #   #    # 
c #    #  #####  #    #  #####   #       #    #  #    #   ####  
c

c----------------------------------------------------------
      subroutine RD_PARAM
c----------------------------------------------------------
c
c  Read various parameters from ibis.params
c 
      use implicit

      use compar    ! npoi, npft
      use comage    ! npftu
      use comveg    ! woodnorm
      use comsoi    ! wpudmax
      use compft    ! PFT parameters incl physiological "constants"
      use comtex    ! Soil texture-related parameters
      use combgc    ! Soil biogeochemistry parameters
      use argvs
c  
c Local variables
c      
      integer*4 parm_unit,  ! file unit assignment for input
     >          j,          ! PFT number (in range 1 to npftu)
     >          npft2,      ! number of PFTs reported in params.veg
     >          npftu2,     ! number of upper canopy PFTs reported in params.veg
     >          nsoil2,     ! number of soil texture classes reported in params.soi
     >          nveg        ! number of vegetation classes reported in params.veg

      real*4    dummyvar    ! use this to read in integers that might be pretending
                            ! to be reals. Also use it as a filler when reading 
                            ! non-existent variables in subroutine ReadItems

      character*100 parm_file
      
      parameter (parm_unit = 9)


c ******************************************************************************
c open the parameter file 'params.can' for input; read in canopy parameters...

      parm_file = 'params.can'
      parm_file = params_can_
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      call readitem(parm_unit, parm_file, tau15)
      call readitem(parm_unit, parm_file, kc15)
      call readitem(parm_unit, parm_file, ko15)
      call readitem(parm_unit, parm_file, cimax)
      call readitem(parm_unit, parm_file, dummyvar)

      npft2 = nint(dummyvar)
      if (npft2 .ne. npft) then
        write (*, 9003) parm_file, npft2, npft
        goto 9006 ! In the circumstances this seems the best thing to do! 
      end if
9003  format ('RD_PARAM Warning: Number of PFTs in ', A10, ' is: ',
     >        I2, ' number in compar.h is: ', I2)  

      call readitem(parm_unit, parm_file, dummyvar)

      npftu2 = nint(dummyvar)
      if (npftu2 .ne. npftu) then
        write (*,9004) parm_file, npftu2, npftu
        goto 9006 ! In the circumstances this seems the best thing to do! 
      end if
9004  format ('RD_PARAM Warning: Number of upper canopy (tree) PFTs ',
     >        'in ', A10, ' is: ', 
     >        I2, ' number in comage.h is: ', I2)  

* Standard photosynthesis parameters for C3 and C4 physiological pathways
      call readitem(parm_unit, parm_file, alpha3)
      call readitem(parm_unit, parm_file, theta3)
      call readitem(parm_unit, parm_file, beta3)
      call readitem(parm_unit, parm_file, alpha4)
      call readitem(parm_unit, parm_file, theta4)
      call readitem(parm_unit, parm_file, beta4)

* Physiological parameters for broadleaf trees
      call readitems(parm_unit, parm_file, 4,
     >         gammaub, coefmub, coefbub, gsubmin, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

* Physiological parameters for conifer trees
      call readitems(parm_unit, parm_file, 4,
     >         gammauc, coefmuc, coefbuc, gsucmin, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

* Physiological parameters for shrubs
      call readitems(parm_unit, parm_file, 4,
     >         gammals, coefmls, coefbls, gslsmin, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

* Physiological parameters for C4 grasses
      call readitems(parm_unit, parm_file, 4,
     >         gammal4, coefml4, coefbl4, gsl4min, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

* Physiological parameters for C3 grasses
      call readitems(parm_unit, parm_file, 4,
     >         gammal3, coefml3, coefbl3, gsl3min, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)


*** DTP 2001/06/05: Note that I have included a new variable here: vmax_pft
*                   This reads in values of vmax assigned initially to each
*                   PFT. These values are then transferred to vmaxub, vmaxuc,
*                   vmaxls, vmaxl4, vmaxl3, in the modified PHYSIOLOGY.F 
*                   module (subroutine STOMATA).

      do j = 1, npft 
        call readitems(parm_unit, parm_file, 8, 
     >           vmax_pft(j), specla(j), 
     >           tauleaf(j), tauroot(j), tauwood0(j), 
     >           aleaf(j), aroot(j), awood(j),
     >           dummyvar, dummyvar) 
      end do

      call readitem(parm_unit, parm_file, woodnorm)

      do j = 1, nband
        call readitems(parm_unit, parm_file, 2, 
     >           rhoveg(j,1), rhoveg(j,2),
     >         dummyvar, dummyvar, dummyvar, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar)
      end do

      do j = 1, nband
        call readitems(parm_unit, parm_file, 2, 
     >           tauveg(j,1), tauveg(j,2),
     >         dummyvar, dummyvar, dummyvar, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar)
      end do

      call readitem(parm_unit, parm_file, chifuz)
      call readitem(parm_unit, parm_file, chiflz)

      call readitems(parm_unit, parm_file, 2, dleaf(1), dleaf(2),
     >         dummyvar, dummyvar, dummyvar, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar)

      call readitems(parm_unit, parm_file, 2, dstem(1), dstem(2),
     >         dummyvar, dummyvar, dummyvar, dummyvar,
     >         dummyvar, dummyvar, dummyvar, dummyvar)

      call readitem(parm_unit, parm_file, alaimu)
      call readitem(parm_unit, parm_file, alaiml)

      call readitem(parm_unit, parm_file, cleaf)
      call readitem(parm_unit, parm_file, cstem)
      call readitem(parm_unit, parm_file, cgrass)

      call readitem(parm_unit, parm_file, chs)
      call readitem(parm_unit, parm_file, chu)
      call readitem(parm_unit, parm_file, chl)

      call readitem(parm_unit, parm_file, wliqumax)
      call readitem(parm_unit, parm_file, wliqsmax)
      call readitem(parm_unit, parm_file, wliqlmax)

      call readitem(parm_unit, parm_file, wsnoumax)
      call readitem(parm_unit, parm_file, wsnosmax)
      call readitem(parm_unit, parm_file, wsnolmax)

      call readitem(parm_unit, parm_file, tdripu)
      call readitem(parm_unit, parm_file, tdrips)
      call readitem(parm_unit, parm_file, tdripl)

      call readitem(parm_unit, parm_file, tblowu)
      call readitem(parm_unit, parm_file, tblows)
      call readitem(parm_unit, parm_file, tblowl)

      close (parm_unit)
      

c ******************************************************************************
c open the parameter file 'params.veg' for input; read in vegetation PFT parameters...

      parm_file = 'params.veg'
      parm_file = params_veg_
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      do j = 1, npft
        call readitems(parm_unit, parm_file, 4, 
     >           TminL(j), TminU(j), Twarm(j), GDD(j),
     >           dummyvar, dummyvar, dummyvar,
     >           dummyvar, dummyvar, dummyvar)
      end do 

      call readitem(parm_unit, parm_file, dummyvar)
      nveg = nint(dummyvar)

      do j = 1, nveg
        call readitems(parm_unit, parm_file, 4,
     >           plai_init(1,j), plai_init(2,j), 
     >           plai_init(3,j), plai_init(4,j),
     >           dummyvar, dummyvar, dummyvar,
     >           dummyvar, dummyvar, dummyvar)
      end do

      call readitem(parm_unit, parm_file, plaiupper)
      call readitem(parm_unit, parm_file, plailower)
      call readitem(parm_unit, parm_file, xminlai)
      call readitem(parm_unit, parm_file, sapfrac_init)
      call readitem(parm_unit, parm_file, beta1)
      call readitem(parm_unit, parm_file, beta2)

      close (parm_unit)

      
c ******************************************************************************
c open the parameter file 'params.soi' for input; read in soil parameters...

      parm_file = 'params.soi'
      parm_file = params_soi_
      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

      do j = 1, nsoilay
        call readitem(parm_unit, parm_file, hsoi(j))
      end do

      call readitem(parm_unit, parm_file, bperm)
      call readitem(parm_unit, parm_file, wpudmax)
      call readitem(parm_unit, parm_file, zwpmax)

      call readitem(parm_unit, parm_file, dummyvar)

      nsoil2 = nint(dummyvar)

      if (nsoil2 .ne. ndat) then ! Is this needed????
        write (*, 9031) parm_file, nsoil2
        write (*, 9032) ndat
        goto 9006 ! In the circumstances this seems the best thing to do! 
      end if

9031  format (' RD_PARAM Warning: Number of soil types in ', 
     >        A10, ' is: ', I2)
9032  format (' Number of soil types in comtex.h is: ', I2)  

      do j = 1, ndat
        call readitems(parm_unit, parm_file, 9,
     >           texdat(1,j), texdat(2,j), texdat(3,j),
     >           porosdat(j), sfielddat(j), swiltdat(j), 
     >           bexdat(j), suctiondat(j), hydrauldat(j),
     >           dummyvar)
      end do

      call readitem(parm_unit, parm_file, lig_frac)
      call readitem(parm_unit, parm_file, fbsom)
      call readitem(parm_unit, parm_file, effac)

      call readitems(parm_unit, parm_file, 8,
     >         cnr(1), cnr(2), cnr(3), cnr(4),
     >         cnr(5), cnr(6), cnr(7), cnr(8),
     >         dummyvar, dummyvar)

      call readitems(parm_unit, parm_file, 5,
     >         fmax, rconst, cnleaf, cnroot, cnwood,
     >         dummyvar, dummyvar, dummyvar, dummyvar, dummyvar)

      call readitems(parm_unit, parm_file, 9,
     >         klm, kls, kll, krm, krs, krl, kwm, kws, kwl,
     >         dummyvar)

      call readitems(parm_unit, parm_file, 7,
     >         kbn, kbp, knb, kns, kpb, kps, ksb,
     >         dummyvar, dummyvar, dummyvar)

      call readitems(parm_unit, parm_file, 9,
     >         ylm, yrm, ywm, yls, yrs, yws, yll, yrl, ywl,
     >         dummyvar)

      call readitems(parm_unit, parm_file, 7,
     >         ybn, ybp, yps, yns, ysb, ypb, ynb,
     >         dummyvar, dummyvar, dummyvar)

      close (parm_unit)
      

c ******************************************************************************
c open the parameter file 'params.hyd' for input; read in vegetation PFT parameters...

*      parm_file = 'params.hyd' ! What is this file supposed to contain???
*      parm_file = params_hyd_
*      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)

*      close (parm_unit) 

c ******************************************************************************

      write (*,*) 'RD_PARAM: All data read in from parameter files ',
     >            'successfully.'

      return ! subroutine RD_PARAM


9000  write (*,2003) parm_file, parm_unit
2003  format ('RD_PARAM: Unexpected EOF encountered in file ',
     >        A12, ' on unit ', I2)
      stop

9001  write (*,2001) parm_file, parm_unit 
2001  format ('RD_PARAM: Error opening parameter file ', A12, 
     >        ' on unit ', I2)      
      stop

9002  write (*,2002) parm_file, parm_unit 
2002  format ('RD_PARAM: Error reading parameter file ', A12, 
     >        ' on unit ', I2)      
      stop

9006  write (*,2006) parm_file, parm_unit
2006  format ('RD_PARAM: Data inconsistency in parameter file ',
     >        A12, ' on unit ', I2) 
      close (parm_unit)
      stop
      end ! subroutine RD_PARAM

c*******************************************************************************

      subroutine ReadItem (funit, fname, item)

c Simple routine to read in data in free format, with a built-in error handler 
c designed to locate and skip over comment strings buried in the data stream.

      use implicit
     
      integer*4    funit,
     >             iocode     ! dummy variable required by system for error handling?
      real*4       item
      character*100 fname

 101  read (funit, *, end=999, err=911, iostat=iocode) item
      return

 911  call CommentHandler (funit, fname)
      goto 101

 999  write (*,2) fname, funit
   2  format ('RD_PARAM: Unexpected EOF encountered in file ',
     >        A20, ' on unit ', I2)
      stop

      end ! subroutine ReadItem

c*******************************************************************************

      subroutine ReadItems (funit, fname, N, item1, item2, item3, item4,
     >                      item5, item6, item7, item8, item9, item10)

c Hokey routine to read in up to N data items in free format, with a built-in 
c error-handler designed to locate and skip over comment strings buried in the
c data stream. Item1 to Item10 are holders for a sequence of legitimate data
c values.  N specifies how many of these are to be used. Therefore N must be 
c less than or equal to 10 (unless you increase the number of items). Note
c that N = 1 is allowed, but it is easier to use subroutine ReadItem for 
c reading single values. 

      use implicit
     
      integer*4    funit,
     >             iocode,   ! dummy variable required by system for error handling?
     >             N

      real*4       item1, item2, item3, item4, item5,
     >             item6, item7, item8, item9, item10
      character*100 fname

* Check to see that N is within acceptable range. Stop if not.

      if ((N.gt.10).or.(N.lt.0)) then
        write (*,3)
   3    format ('READITEMS: Invalid number of items specified.')
        stop
      end if

 100  goto (101,102,103,104,105,106,107,108,109,110) N
  
 101  read (funit, *, end=999, err=911, iostat=iocode) 
     >      item1
      return

 102  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2
      return

 103  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3
      return

 104  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4
      return

 105  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5
      return

 106  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5, item6
      return

 107  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5, item6,
     >      item7
      return

 108  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5, item6,
     >      item7, item8
      return

 109  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5, item6,
     >      item7, item8, item9
      return

 110  read (funit, *, end=999, err=911, iostat=iocode)
     >      item1, item2, item3, item4, item5, item6,
     >      item7, item8, item9, item10
      return

 911  call CommentHandler (funit, fname)
      goto 100

 999  write (*,2) fname, funit
   2  format ('RD_PARAM: Unexpected EOF encountered in file ',
     >        A20, ' on unit ', I2)
      stop

      end ! subroutine ReadItems

c*******************************************************************************

      subroutine CommentHandler (parm_unit, parm_file)
      
c #, *, c, C, !, // and /* denote comments in the parameter file. Therefore the 
c input reading routine must ignore these symbols and any text to their right. 
C This subroutine checks for the occurrence of these symbols, whenever RD_PARAM
c detects an error reading numeric data. If the error was caused by text starting
c with one of these symbols, the text is skipped over and the routine returns
c so that the next line of data can be read in.

c If an error is detected, the line containing the offending characters is 
c echoed to the screen and execution stops. But if someone were so-inclined,
c they could write some extra code to analyse the problem and possibly recover
c from it.

      use implicit

      integer*4     parm_unit,
     >              iocode    ! dummy variable required by system for error handling?
      integer*2     is_comment
      character*100  parm_file ! name of file from which input is being read in
      character*255 inputline ! arbitrary string. Input lines cannot exceed 255 chars!

      backspace(parm_unit)    ! back up to start of field that caused this error
 
      read (parm_unit, 1002, end=800, err=801, iostat=iocode) inputline
 1002 format (A255)

      is_comment = 0          ! must initialize this correctly

      if ((index (inputline, 'C') .eq. 1) .or. 
     >    (index (inputline, 'c') .eq. 1) .or.
     >    (index (inputline, '*') .eq. 1) .or.
     >    (index (inputline, '#') .eq. 1) .or.
     >    (index (inputline, '!') .eq. 1) .or.
     >    (index (inputline, '/') .eq. 1)) is_comment = 1


      if (is_comment .NE. 1) goto 801 ! we have a problem....
 
      return  ! But if we get to here, all's well and we can move on.

  800 write (*, 1000) parm_file, inputline ! For now, just echo the line and stop
 1000 format ('CommentHandler: Unexpected EOF encountered in ',
     >        A12, '\nLast input was:', A255)

      stop

  801 write (*, 1001) parm_file, inputline ! For now, just echo the line and stop.
 1001 format ('CommentHandler: Error encountered reading data ',
     >        'from ', A12, '\nLast input was:', A255)

      stop

      end ! subroutine CommentHandler

 


