program main
use sharedvars
implicit none
include "fftw3.f"
   !----------------------------------------------------------------------------
   ! Some compiler needs to declare iargc.
   !integer, external             :: iargc
   !----------------------------------------------------------------------------
   ! read command line options
   narg = iargc()
   if ( narg.lt.1 ) call helpinfo
   !
   ! set default values
   outfile = 'acf.dat'
   nacf    = -1
   npad    = -1
   fave    = .true.
   fread   = .false.
   tmcol   = 1
   !
   ncol    = 1
   CMax    = max(2,narg)
   allocate( icol(CMax) )
   lcmplx = .false.
   lcross = .false.
   !
   iarg = 1
   do while ( iarg.le.narg )
      call getarg(iarg, oneline)
      select case ( trim(oneline) )
      case ( "-o", "-O" )
         iarg = iarg + 1
         if (iarg.gt.narg) stop 'No enough command line option!'
         call getarg(iarg, outfile)
      case ( "-n", "-N" )
         iarg = iarg + 1
         if (iarg.gt.narg) stop 'No enough command line option!'
         call getarg(iarg, oneline)
         read(oneline,*,iostat=rderr) nacf
         if (rderr.ne.0.or.nacf.lt.1) nacf = -1
      case ( "-p", "-P" )
         iarg = iarg + 1
         if (iarg.gt.narg) stop 'No enough command line option!'
         call getarg(iarg, oneline)
         read(oneline,*,iostat=rderr) npad
         if (rderr.ne.0.or.npad.lt.0) npad = -1
      case ( "-t", "-T" )
         iarg = iarg + 1
         if (iarg.gt.narg) stop 'No enough command line option!'
         call getarg(iarg, oneline)
         read(oneline,*,iostat=rderr) tmcol
         if (rderr.ne.0.or.tmcol.lt.1) tmcol = 1
      case ( "-ave", "-AVE" )
         fave = .true.
      case ( "-nave", "-nAVE", "-NAVE" )
         fave = .false.
      case ( "-c", "-C", "-complex", "-cmplx" )
         lcmplx = .true.
      case ( "-cr", "-cross", "-CR" )
         lcross = .true.
      case ( "-h", "-H", "-help" )
         call helpinfo
      case default
         if (.not.fread) then
            infile = trim(oneline)
            fread = .true.
         else
            read(oneline,*,iostat=rderr) idum
            if (rderr.eq.0.and.idum.ge.1) then
               ncol = ncol + 1
               icol(ncol) = idum
            endif
         endif
      end select
      iarg = iarg + 1
   enddo
   if (ncol.le.1) then
      ncol = 2
      icol(2) = 2
   endif
   icol(1) = tmcol
   icmax = maxval(icol)
   !
   if ( lcross ) then
      if ( ncol.lt.3 ) then
         write(*,'(/,5x,"For cross-correlation, two and only two columns are needed!")')
         stop
      else
         ncol = 3
      endif
   endif
   !
   call CPU_TIME(tm1)
   !
   if (lcmplx) then
      icmax = icmax + 1
      if ( lcross ) then
         call cmplx_fft_2_ccf
      else
        call cmplx_fft_2_acf
      endif
   else
      if ( lcross ) then
         call real_fft_2_ccf
      else
         call real_fft_2_acf
      endif
   endif
   !
   call free_all
   call CPU_TIME(tm2)
   !
   write(*,'(/,10x,"Total CPU time used: ", F10.2)') tm2-tm1
   !
contains
   subroutine helpinfo
   implicit none
      write(*,'(/,"acf    To evaluate the autocorrelation of one column data from file.")')
      write(*,'(  "Usage: fftacf [options] infile [col1 col2 ...]")')
      write(*,'(  "Options:")')
      write(*,'(  "   -o outfile; default: acf.dat")')
      write(*,'(  "   -n nacf;  maximum lag, default: same as # of data")')
      write(*,'(  "   -t tcom;  column number for time; by default is 1")')
      write(*,'(  "   -ave; to subtract average from data; default: true")')
      write(*,'(  "   -nave; not to subtract average from data; default: false")')
      write(*,'(  "   -c; to specify that the input data are complex")')
      write(*,'(  "   -cross; cross-correlation is required; two columns only are needed.")')
      write(*,'(  "   col1 ...: column numbers to evaluate ACF; by default, column 2.",/)')
      stop
   end subroutine
end program
