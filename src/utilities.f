c
c #    #   #####     #    #          #     #####     #    ######   ####
c #    #     #       #    #          #       #       #    #       #
c #    #     #       #    #          #       #       #    #####    ####
c #    #     #       #    #          #       #       #    #            #
c #    #     #       #    #          #       #       #    #       #    #
c  ####      #       #    ######     #       #       #    ######   ####
c
c ---------------------------------------------------------------------
      subroutine scopy (nt, arr, brr)
c ---------------------------------------------------------------------
c
c copies array arr to brr,for 1st nt words of arr
c
      include 'implicit.h'
c
c Arguments
c
      integer nt     
      real arr(nt),    ! input
     >     brr(nt)     ! output
c
c Local variables
c
      integer ia
c
      do 100 ia = 1, nt
        brr(ia) = arr(ia)
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine const (arr, nar, value)
c ---------------------------------------------------------------------
c
c sets all elements of real vector arr to value
c
      include 'implicit.h'
c
c Arguments
c
      integer nar
c     
      real value
      real arr(nar)
c
c Local variables
c
      integer j
c
      do 100 j = 1, nar
        arr(j) = value
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine endrun
c ---------------------------------------------------------------------
c
c stops gracefully
c
      stop
      end
c
c 
c ---------------------------------------------------------------------
      real function cvmgt (x,y,l)
c ---------------------------------------------------------------------
c
c chooses between two things.  Used in canopy.f
c
      include 'implicit.h'
c
      logical l
      real x, y
c
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine arr2vec(array, vect)
c ---------------------------------------------------------------------
c
c extracts land points from array and puts them into a land-only vector
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
c
c Arguments
c
      real array(nlonsub,nlatsub), vect(npoi)
c
c Local variables
c
      integer j,i,npts
c
      npts = 0
c
      do 10 j = 1, nlatsub
        do 20 i = 1, nlonsub
          if (lmask(i,j) .eq. 1) then
            npts = npts + 1
            vect(npts) = array(i,j)
          end if
 20     continue
 10   continue
c
      if (npts .ne. npoi) then
        write (*,*) 'ERROR in arr2vec'
        write (*,*) 'npts not equal to npoi'
        write (*,*) 'npts = ', npts, ' npoi = ', npoi
        stop 1
      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vec2arr(vect, array)
c ---------------------------------------------------------------------
c
c puts vector of land-only points back into array
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comwork.h'
c
c Arguments
c
      real array(nlonsub,nlatsub), vect(npoi)
c
c Local variables
c
      integer j, i, ii, jj
c
      do 10 j = 1, nlatsub
        do 20 i = 1, nlonsub
          array(i,j) = OCEAN
 20     continue
 10   continue
c
      do 30 i = 1, npoi
        ii = lonindex(i)
        jj = latindex(i)
        array(ii,jj) = vect(i)
 30   continue
c
      return
      end
c
c ---------------------------------------------------------------------
c lenchr - find index of last non-blank, non-null
c ---------------------------------------------------------------------
c
      integer function lenchr (ch)
c
c returns position of last non-blank,null character in ch,
c or 1 if ch is all blanks
c
      include 'implicit.h'
c
c Arguments
c
      character*(*) ch
c
      integer i
c
      do i = len(ch), 1, -1
        if (ch(i:i).ne.' '.and.ch(i:i).ne.char(0)) then
           lenchr = i
           return
        endif
      enddo
c
      lenchr = 1
c
      return
      end
c
