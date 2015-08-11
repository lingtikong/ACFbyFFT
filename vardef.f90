module sharedvars
implicit none
   !----------------------------------------------------------------------------
   integer                       :: i, II, npts, nfft, nacf=-1, idum, npad = 0
   integer                       :: tmcol, ncol, CMax, NMax, icmax, NINC
   integer, allocatable          :: icol(:)
   double precision              :: r_nfft, ave
   !
   double precision, allocatable :: data_in(:,:), data_tmp(:,:), data_read(:)
   double precision, allocatable :: fft_in(:)
   double complex, allocatable   :: fft_out(:), fft_cmplx_in(:), fft_out2(:)
   !
   integer*8                     :: pFW, pBW
   logical                       :: fexst, fread, fave, lcmplx, lcross
   integer                       :: narg, iarg, ioerr, rderr
   character(len=512)            :: infile, outfile, fmtstr
   character(len=1024)           :: oneline
   real                          :: tm1, tm2
   !----------------------------------------------------------------------------
contains
   !----------------------------------------------------------------------------
   subroutine free_all
   implicit none
      if (allocated( icol      ) ) deallocate( icol )
      if (allocated( data_in   ) ) deallocate( data_in )
      if (allocated( data_tmp  ) ) deallocate( data_tmp )
      if (allocated( data_read ) ) deallocate( data_read )
      if (allocated( fft_in    ) ) deallocate( fft_in )
      if (allocated( fft_out   ) ) deallocate( fft_out )
      if (allocated( fft_out2  ) ) deallocate( fft_out2 )
      if (allocated( fft_cmplx_in ) ) deallocate( fft_cmplx_in )
   end subroutine
   !----------------------------------------------------------------------------
end module
