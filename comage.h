c
c ---------
c comage  header file  last update 04/22/00 (DTP)
c ---------

c ------------------------------------------------------------------------------
c Parameters for age-class management engine (ACME) 
c ------------------------------------------------------------------------------

      integer*4 maxage,            ! range of cohort ages to be simulated
     >          ageint,            ! age-class width (or interval) in years 
*     >          interval,          ! age-class reporting interval (years; mod (maxage, interval) = 0)
     >          ncohort,           ! number of age cohorts = maxage/ageint
     >          ac_var,            ! number of age-class dependent variables for
                                   ! which data will be tracked from year to year     
     >          npftu,             ! number of upper canopy pfts (used in modified 
                                   ! version of subroutine DYNAVEG) 
     >          ac_size            ! size (in REAL*4) of single age class data set

      integer*4 ac_recl,           ! size (in BYTES) of complete single grid point
                                   ! of age-class data  (Could be BIG!!!!)
     >          sum_recl           ! size (in BYTES) of single ageclass summary data set                                   
       real*4   fu_acmin           ! minimum cover fraction occupied by an upper
     >                             ! canopy age cohort (= 1/ncohort)                                       


*** N.B. maxage must be divisible by ageint with no remainder!!!

      parameter (maxage = 100)                 ! Maximum stand age allowed in years 
      parameter (ageint = 1)                   ! Width of each age-class in years  
      parameter (ncohort = maxage/ageint)      ! Max number of age cohorts allowed 
                                               ! oldest cohort includes stands >= maxage 

      parameter (ac_var = 3)                   ! cbiol_ac, cbior_ac, cbiow_ac
      parameter (npftu = 8)                    ! number of upper canopy PFTs                                   
      
      parameter(ac_size = npftu * (ncohort + 1))    ! size of single data arrays used to
                                                    ! handle age-dependent data - in REAL*4.
      parameter(ac_recl = ((ncohort + 1) + 
     >                     (ac_size * ac_var)) * 4) ! size of complete single file record in bytes
      
      parameter(sum_recl = ac_size * 4)             ! size of single ageclass summary data set in bytes                                     
      real*4 f_Aveg (0:ncohort),             ! land area fraction occupied by upper canopy 
     >       cbiow_ac (npftu, 0:ncohort),    ! biomass components in each age class and PFT 
     >       cbiol_ac (npftu, 0:ncohort),
     >       cbior_ac (npftu, 0:ncohort)
     
      real*4 cbiow_init(npftu),              ! These arrays contain copies of the initial
     >       cbiol_init(npftu),              ! values of c pools for wood, foliage and roots
     >       cbior_init(npftu)               ! assigned to each PFT at the start of the run

c ------------------------------------------------------------------------------
c File names and unit assignments for ACME binary and text files
c ------------------------------------------------------------------------------

      integer*4 ac_unit,    ! binary file to store age-class data
     >          pt_unit     ! text file to store age-class printouts
     
      character*12 ac_file, ! name of binary file to store age-class data  
     >             pt_file, ! name of text file to store age-class printouts
     >           acsum_file ! name of text file to store age-class summary        
     
      PARAMETER (ac_unit = 51, pt_unit = 52) ! Should check these don't conflict  
                                             ! with other I/O unit numbers....
                                            
      PARAMETER (ac_file = 'ageclass.dat', 
     >           pt_file = 'ageclass.txt',
     >           acsum_file = 'ac_summ.dat')
      
      common /comage1/ fu_acmin, cbiow_init, cbiol_init, cbior_init
      common /comage2/ f_Aveg, cbiow_ac, cbiol_ac, cbior_ac
                                            
c ------------------------------------------------------------------------------
c Variables needed to do summaries across landclasses.
c ------------------------------------------------------------------------------

      integer*4 n_lcs,         ! The max number of landclasses
     >          lc_unit        ! use a single unit number repeatedly.... 

      parameter (n_lcs = 100,  ! land classes can be set from 1 to 100
     >           lc_unit = 44)
    
      integer*4 used(n_lcs)    ! flag which lc values are 
                               ! in use in the current run

c ------------------------------------------------------------------------------
   
      common /comage3/ used 
