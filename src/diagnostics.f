c diagnostics.f  last update 12/18/00 CCM
c
c This code contains subroutines relevant to diagnostics.  Code
c fragments were removed from main.f and io.f and collected here.
c
c---------------------------------------------------------------------
      subroutine inidiag(idiag)
c---------------------------------------------------------------------
c CALLED BY main
c CALLS diaginit
c This subroutine reads in diagnostics info from diag.infile and calls
c the subroutine diaginit.
c
c INPUT
      integer idiag          ! number of diagnostic files
c COMMON
      use comdiag
      use argvs
c INTERNAL
      integer lun, i, ivar
c
      real    diaglat(nfiles),       ! latitude of diagnostic point
     >        diaglon(nfiles)        ! longitude of diagnostic point
c
        lun = 12
        open (lun, status='old', file=diag_infile_)
c
        do 100 i = 1, 26
          read (lun,*)
 100    continue
c
        do 110 i = 1, nfiles
          read (lun,*) diaglat(i)
          read (lun,*) diaglon(i)
          read (lun,*) diagstart(i)
          read (lun,*) diagend(i)
          read (lun,*) nfreq(i)
          read (lun,*)
 110    continue
c
        do 120 i = 1, 8
          read (lun,*)
 120    continue
c
        do 130 ivar = 1, nvars
          read (lun,*) (ldiag(ivar,ifile),ifile=1,nfiles)
 130    continue
c
        close (lun)
c
         call diaginit (idiag, diaglat, diaglon)
c
         return
         end
c
c ---------------------------------------------------------------------
      subroutine diaginit(idiag, diaglat, diaglon)
c ---------------------------------------------------------------------
c CALLED BY inidiag
c CALLS nothing
c Initialize diagnostic output files
c
      use compar
      use comwork
      use comdiag
      use argvs
c
c local variables
c
      character*100 filen          ! name of diagnostic file
      character*50 varname(nvars) ! variable names
c
      character*50 header(16,nfiles)  ! column headers for output file
c
      integer idiag,              ! number of diganostic files
     >        ifile,              ! counter
     >        ivar,               ! counter
     >        kolumn,             ! column in output file
     >        lat(nfiles),
     >        lon(nfiles),
     >        countvar
c
      real    alat0,
     >        alon0
c
      real    diaglat(nfiles),       ! latitude of diagnostic point
     >        diaglon(nfiles)        ! longitude of diagnostic point
c
      data varname /'adco2mic','adco2root','adco2soi','adneetot',
     >              'adtsoic','adwsoic','amco2root','amco2soi',
     >              'ancl3','ancl4','ancub','ancuc','asurd1',
     >              'asurd2','asuri1','asuri2','ayco2soi',
     >              'biomass1','biomass2','biomass3','biomass4',
     >              'biomass5','biomass6','biomass7','biomass8',
     >              'biomass9','biomass10','biomass11','biomass12',
     >              'cloud','coszen','fi','fira','firb','frac1','frac2',
     >              'frac3','frac4','frac5','frac6','frac7','frac8',
     >              'frac9','frac10','frac11','frac12','gadjust',
     >              'gdrain','ginvap','grunof','gsuvap','gtrans',
     >              'gtransl','gtransu','hsno1','hsno2','hsno3',
     >              'hsnotop','plai1','plai2','plai3','plai4','plai5',
     >              'plai6','plai7','plai8','plai9','plai10','plai11',
     >              'plai12','precip','psurf','qa','raina','snowa',
     >              'solad1','solad2','solai1','solai2','ta','td','tl',
     >              'totalit','totcmic','totcsoi','totrlit','trng',
     >              'ts','tsno1','tsno2','tsno3','tsoi1','tsoi2',
     >              'tsoi3','tsoi4','tsoi5','tsoi6','tu',
     >              'wet','wipud','wisoi1','wisoi2','wisoi3','wisoi4',
     >              'wisoi5','wisoi6','wliql',
     >              'wliqs','wliqu',
     >              'wpud','wsnol','wsnos','wsnou','wsoi1','wsoi2',
     >              'wsoi3','wsoi4','wsoi5','wsoi6' /
c
c Initialize diagnostic output files
c
      do 10 ifile = 1,nfiles
        kolumn = 4
        countvar = 0
        do 20 ivar = 1,nvars
          if ( ldiag(ivar,ifile) .eq. 1 ) then
            kolumn = kolumn + 1
            header(kolumn,ifile) = varname(ivar)
            countvar = countvar + 1
          end if
 20     continue
          if (countvar .lt. 12) then
            do 25 i = countvar+5,16
              header(i,ifile) = 'empty'
 25         continue
          end if
 10   continue
c    
      do 30 i = 1, idiag
c       filen = 'ibis.diag'
        if(i==1) then
            filen = out_diag_0_
        elseif(i==2) then
            filen = out_diag_1_
        elseif(i==3) then
            filen = out_diag_2_
        elseif(i==4) then
            filen = out_diag_3_
        elseif(i==5) then
            filen = out_diag_4_
        elseif(i==6) then
            filen = out_diag_5_
        elseif(i==7) then
            filen = out_diag_6_
        elseif(i==8) then
            filen = out_diag_7_
        elseif(i==9) then
            filen = out_diag_8_
        elseif(i==10) then
            filen = out_diag_9_
        endif
c       write(filen(10:10),'(i1)')i-1
        open (i+77, status='unknown', file=filen, access='append')
        write (i+77,3000) diaglat(i), diaglon(i)
        write (i+77,3100) diagstart(i), diagend(i)
        write (i+77,3200) (header(kolumn,i),kolumn=5,16)
c
 3000 format ('%Cell_latitude ', f8.2,1x,' Cell_longitude ',f8.2)
 3100 format ('%Begin_year ',i8,' End_year ',i8)
 3200 format ('%Year',2x,'Month',4x,'Day',4x,'Step',3x,12(a12,1x))
c
 30   continue
c
c Determine the grid cells corresponding to the chosen coordinates
c
      alat0 = latscale(1)
      alon0 = lonscale(1)
c
      do 40 i = 1, idiag
        lat(i) = 1 + int( (alat0 - diaglat(i))/yres + 0.5 )
        lat(i) = min(lat(i),nlatsub)
        lat(i) = max(lat(i),1)
        lon(i) = 1 + int( (diaglon(i) - alon0)/xres + 0.5 )
        lon(i) = min(lon(i),nlonsub)
        lon(i) = max(lon(i),1)
        ndiagpt(i) = -999.
        do 50 n = 1, npoi
          if (latindex(n) .eq. lat(i) .and. 
     >        lonindex(n) .eq. lon(i)) then
             ndiagpt(i) = n
          end if
 50    continue
       if (ndiagpt(i) .lt. 1) then
          print *, 'Warning: diag point ',i,':',diaglon(i),diaglat(i),
     >     ' is not a land point.  Ignoring'
       end if
 40   continue
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine wdiag(i, iyear, imonth, iday, istep)
c ---------------------------------------------------------------------
c CALLED by main
c CALLS nothing
c
c This subroutine writes diagnostic output to the "diag" files
c
      use compar
      use com1d
      use comveg
      use comhyd
      use comatm
      use comsno
      use comsoi
      use comsum
      use comdiag
c
      real      diag(nvars)              ! diagnostic variables
c
      real      rowdata(16)              ! rows of data for output
c
      integer   i,         ! # diagnostic file requested
     >          ivar,      ! index
     >          kolumn,    ! columns for output
     >          countvar   !
c
      diag(1) = adco2mic(ndiagpt(i))
      diag(2) = adco2root(ndiagpt(i))
      diag(3) = adco2soi(ndiagpt(i))
      diag(4) = adneetot(ndiagpt(i))
      diag(5) = adtsoic(ndiagpt(i))
      diag(6) = adwsoic(ndiagpt(i))
      diag(7) = amco2root(ndiagpt(i))
      diag(8) = amco2soi(ndiagpt(i))
      diag(9) = ancl3(ndiagpt(i))
      diag(10) = ancl4(ndiagpt(i))
      diag(11) = ancub(ndiagpt(i))
      diag(12) = ancuc(ndiagpt(i))
      diag(13) = asurd(ndiagpt(i),1)
      diag(14) = asurd(ndiagpt(i),2)
      diag(15) = asuri(ndiagpt(i),1)
      diag(16) = asuri(ndiagpt(i),2)
      diag(17) = ayco2soi(ndiagpt(i))
      diag(18) = biomass(ndiagpt(i),1)
      diag(19) = biomass(ndiagpt(i),2)
      diag(20) = biomass(ndiagpt(i),3)
      diag(21) = biomass(ndiagpt(i),4)
      diag(22) = biomass(ndiagpt(i),5)
      diag(23) = biomass(ndiagpt(i),6)
      diag(24) = biomass(ndiagpt(i),7)
      diag(25) = biomass(ndiagpt(i),8)
      diag(26) = biomass(ndiagpt(i),9)
      diag(27) = biomass(ndiagpt(i),10)
      diag(28) = biomass(ndiagpt(i),11)
      diag(29) = biomass(ndiagpt(i),12)
      diag(30) = cloud(ndiagpt(i))
      diag(31) = coszen(ndiagpt(i))
      diag(32) = fi(ndiagpt(i))
      diag(33) = fira(ndiagpt(i))
      diag(34) = firb(ndiagpt(i))
      diag(35) = frac(ndiagpt(i),1)
      diag(36) = frac(ndiagpt(i),2)
      diag(37) = frac(ndiagpt(i),3)
      diag(38) = frac(ndiagpt(i),4)
      diag(39) = frac(ndiagpt(i),5)
      diag(40) = frac(ndiagpt(i),6)
      diag(41) = frac(ndiagpt(i),7)
      diag(42) = frac(ndiagpt(i),8)
      diag(43) = frac(ndiagpt(i),9)
      diag(44) = frac(ndiagpt(i),10)
      diag(45) = frac(ndiagpt(i),11)
      diag(46) = frac(ndiagpt(i),12)
      diag(47) = gadjust(ndiagpt(i))
      diag(48) = gdrain(ndiagpt(i))
      diag(49) = ginvap(ndiagpt(i))
      diag(50) = grunof(ndiagpt(i))
      diag(51) = gsuvap(ndiagpt(i))
      diag(52) = gtrans(ndiagpt(i))
      diag(53) = gtransl(ndiagpt(i))
      diag(54) = gtransu(ndiagpt(i))
      diag(55) = hsno(ndiagpt(i),1)
      diag(56) = hsno(ndiagpt(i),2)
      diag(57) = hsno(ndiagpt(i),3)
      diag(58) = hsnotop
      diag(59) = plai(ndiagpt(i),1)
      diag(60) = plai(ndiagpt(i),2)
      diag(61) = plai(ndiagpt(i),3)
      diag(62) = plai(ndiagpt(i),4)
      diag(63) = plai(ndiagpt(i),5)
      diag(64) = plai(ndiagpt(i),6)
      diag(65) = plai(ndiagpt(i),7)
      diag(66) = plai(ndiagpt(i),8)
      diag(67) = plai(ndiagpt(i),9)
      diag(68) = plai(ndiagpt(i),10)
      diag(69) = plai(ndiagpt(i),11)
      diag(70) = plai(ndiagpt(i),12)
      diag(71) = precip(ndiagpt(i))
      diag(72) = psurf(ndiagpt(i))
      diag(73) = qa(ndiagpt(i))
      diag(74) = raina(ndiagpt(i))
      diag(75) = snowa(ndiagpt(i))
      diag(76) = solad(ndiagpt(i),1)
      diag(77) = solad(ndiagpt(i),2)
      diag(78) = solai(ndiagpt(i),1)
      diag(79) = solai(ndiagpt(i),2)
      diag(80) = ta(ndiagpt(i))
      diag(81) = td(ndiagpt(i))
      diag(82) = tl(ndiagpt(i))
      diag(83) = totalit(ndiagpt(i))
      diag(84) = totcmic(ndiagpt(i))
      diag(85) = totcsoi(ndiagpt(i))
      diag(86) = totrlit(ndiagpt(i))
      diag(87) = tmax(ndiagpt(i))
      diag(88) = ts(ndiagpt(i))
      diag(89) = tsno(ndiagpt(i),1)
      diag(90) = tsno(ndiagpt(i),2)
      diag(91) = tsno(ndiagpt(i),3)
      diag(92) = tsoi(ndiagpt(i),1)
      diag(93) = tsoi(ndiagpt(i),2)
      diag(94) = tsoi(ndiagpt(i),3)
      diag(95) = tsoi(ndiagpt(i),4)
      diag(96) = tsoi(ndiagpt(i),5)
      diag(97) = tsoi(ndiagpt(i),6)
      diag(98) = tu(ndiagpt(i))
      diag(99) = iwet(ndiagpt(i))
      diag(100) = wipud(ndiagpt(i))
      diag(101) = wisoi(ndiagpt(i),1)
      diag(102) = wisoi(ndiagpt(i),2)
      diag(103) = wisoi(ndiagpt(i),3)
      diag(104) = wisoi(ndiagpt(i),4)
      diag(105) = wisoi(ndiagpt(i),5)
      diag(106) = wisoi(ndiagpt(i),6)
      diag(107) = wliql(ndiagpt(i))
      diag(108) = wliqs(ndiagpt(i))
      diag(109) = wliqu(ndiagpt(i))
      diag(110) = wpud(ndiagpt(i))
      diag(111) = wsnol(ndiagpt(i))
      diag(112) = wsnos(ndiagpt(i))
      diag(113) = wsnou(ndiagpt(i))
      diag(114) = wsoi(ndiagpt(i),1)
      diag(115) = wsoi(ndiagpt(i),2)
      diag(116) = wsoi(ndiagpt(i),3)
      diag(117) = wsoi(ndiagpt(i),4)
      diag(118) = wsoi(ndiagpt(i),5)
      diag(119) = wsoi(ndiagpt(i),6)
c
      rowdata(1) = float(iyear)
      rowdata(2) = float(imonth)
      rowdata(3) = float(iday)
      rowdata(4) = float(istep)
c
      kolumn = 4
      countvar = 0
c
      do 30 ivar = 1, nvars
        if( ldiag(ivar,i) .eq. 1) then
           countvar = countvar + 1
           kolumn = kolumn + 1
           rowdata(kolumn) = diag(ivar)
        end if
 30   continue
c
      if (countvar .lt. 12) then
        do 40 n = countvar + 5, 16
          rowdata(n) = -999.
 40     continue
      end if
c
      write(i+77, 3000) rowdata
 3000 format(1x,f5.0,2x,f3.0,3x,f7.0,1x,f5.0,1x,12(e12.3,1x))
c
      return
      end
c
