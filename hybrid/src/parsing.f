      module parsing

      implicit none

      contains

      subroutine parse( line, nn, lc, names, nv, values,
     .                  ni, integers, nr, reals )

c **********************************************************************
c Extracts the names and numbers in an input line.
c Written by J.M.Soler. August 1997.
c ********* Input ******************************************************
c character line*(*) : Input string.
c ********* Output *****************************************************
c integer   nn          : Number of names present in line
c integer   lc(0:*)     : Last character of each name in names string.
c                         Notice that first index is 0.
c character names*(*)   : Non-number strings present in line
c integer   nv          : Number of values present in line
c real*8    values(*)   : Numbers present in line (integer or real)
c integer   ni          : Number of integers present in line
c integer   integers(*) : Integer numbers present in line
c integer   nr          : Number of reals present in line
c real*8    reals(*)    : Real numbers present in line
c ********* Behaviour **************************************************
c Names are returned consecutively in string names. Name i may be
c   extracted as name = names(lc(i-1)+1:lc(i))
c Names and numbers in line may be separated by blanks (including tabs), 
c or commas. Other characters are considered part of names or numbers,
c not separators.
c Names are any non-number strings present in the line.
c Integers are numbers without a dot.
c Real numbers are i.e. 5.0, 2., -.25, +1.d6, 2.E-3
c Apostrophes enclosing strings are removed in output names.
c Double apostrophes inside a name are left unchanged.
c The size of arrays names, values, integers and reals must be long
c  enough to contain all instances present in line. No check is made.

      implicit          none
      integer, intent(out)           :: integers(:), lc(0:)
      integer, intent(out)           :: ni, nn, nr, nv
      character(len=*), intent(out)  :: names
      character(len=*), intent(in)   :: line
      double precision, intent(out)  :: reals(:), values(:)


c     Internal variables
      logical
     .  isinteger, isreal, opened
!
!     Automatic arrays
!
      logical, dimension(0:len(line)+1) ::
     .  isblank, isdigit, isdot, isexp, issign, isother

      character  c
      character fmtstr*10
      integer    i, i1, i2, ic, iw, ll, ndot, nexp, nsign
c
      names = ' '

c     Check line length
      ll = len(line)

c     Clasify line characters
      isblank(0) = .true.
      isdigit(0) = .false.
      isdot(0)   = .false.
      isexp(0)   = .false.
      issign(0)  = .false.
      isother(0) = .false.
      do i = 1,ll
        c = line(i:i)
        ic = ichar(c)
        isblank(i) = (c.eq.' ' .or. c.eq.',' .or. ic.eq.9)
        isdigit(i) = (ic.ge.48 .and. ic.le.57)
        isdot(i)   = (c.eq.'.')
        isexp(i)   = (c.eq.'e' .or. c.eq.'E' .or.
     .                c.eq.'d' .or. c.eq.'D')
        issign(i)  = (c.eq.'+' .or. c.eq.'-')
        isother(i) = (.not.isblank(i) .and. .not.isdigit(i) .and.
     .                .not.isdot(i)   .and. .not.isexp(i)   .and.
     .                .not.issign(i))
      enddo
      isblank(ll+1) = .true.
      isdigit(ll+1) = .false.
      isdot(ll+1)   = .false.
      isexp(ll+1)   = .false.
      issign(ll+1)  = .false.
      isother(ll+1) = .false.

c     Initialize number of instances
      nn = 0
      nv = 0
      ni = 0
      nr = 0
      lc(0) = 0

c     Iterate on 'words'
      i2 = 0
      do iw = 1,ll

c       Locate first character of word
        do i1 = i2+1,ll
          if (.not.isblank(i1)) goto 10
        enddo
        goto 999
   10   continue

c       Locate last character of word, while checking if
c       a pair of enclosing apostrophes is 'open'
        opened = .false.
        do i2 = i1,ll
          if (line(i2:i2).eq.'''') opened = .not.opened
          if (.not.opened .and. isblank(i2+1)) goto 15
        enddo
   15   continue

c       Find if word is a number
        isinteger = .true.
        isreal    = .true.
        ndot  = 0
        nexp  = 0
        nsign = 0
        do i = i1,i2

c         Non-allowed characters for any number (label 30 is loop exit)
          if (isother(i)) then
            isinteger = .false.
            isreal    = .false.
            goto 30
          endif

c         Allowed for integer and real numbers (label 20 is loop cont)
          if (isdigit(i)) goto 20
          if (issign(i) .and. isblank(i-1) .and.
     .        .not.isblank(i+1)) goto 20

          isinteger = .false.

c         Allowed only for real number
          if (isdot(i) .and.
     .        (isblank(i-1) .or. isdigit(i-1) .or. issign(i-1)) .and.
     .        (isblank(i+1) .or. isdigit(i+1) .or. isexp(i+1))) goto 20
          if (isexp(i) .and. 
     .        (isdigit(i-1) .or. isdot(i-1)) .and.
     .        (isdigit(i+1) .or. issign(i+1))) goto 20
          if (issign(i) .and. isexp(i-1) .and. isdigit(i+1)) goto 20

          isreal = .false.
          goto 30

   20     continue
          if (isdot(i))  ndot  = ndot  + 1
          if (isexp(i))  nexp  = nexp  + 1
          if (issign(i)) nsign = nsign + 1
        enddo

c       Check that there are not too many dots and so on
        if (ndot.gt.1)  isreal = .false.
        if (nexp.gt.1)  isreal = .false.
        if (nsign.gt.2) isreal = .false.

   30   continue

c       Some printout for debugging
*       write(6,'(a,2i4,3x,a10,2l3)')
*    .    'parse: i1, i2, word, isinteger, isreal =',
*    .    i1, i2, line(i1:i2), isinteger, isreal

c       Add word to proper class
        if (isinteger) then
          ni = ni + 1
          nv = nv + 1
          write(fmtstr,9000) i2-i1+1
 9000     format('(i',i2.2,')')
          read(line(i1:i2),fmt=fmtstr) integers(ni)
          values(nv) = integers(ni)
        elseif (isreal) then
          nr = nr + 1
          nv = nv + 1
          write(fmtstr,9020) i2-i1+1
 9020     format('(g',i2.2,'.0)')
          read(line(i1:i2),fmt=fmtstr) reals(nr)
          values(nv) = reals(nr)
        else
          nn = nn + 1
c         Remove enclosing apostrophes if present
          if (line(i1:i1).eq.'''' .and. line(i2:i2).eq.'''') then
            lc(nn) = lc(nn-1) + i2 - i1 - 1
            names(lc(nn-1)+1:lc(nn)) = line(i1+1:i2-1)
          else
            lc(nn) = lc(nn-1) + i2 - i1 + 1
            names(lc(nn-1)+1:lc(nn)) = line(i1:i2)
          endif
        endif
      enddo

C Exit point
  999 continue

      end subroutine parse

      end module parsing






