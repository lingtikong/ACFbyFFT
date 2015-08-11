! Program to evaluate the auto-correlation function of columns of complex date from file
! You should link it with fftw 3 librarys while compiling
subroutine cmplx_fft_2_acf
use sharedvars
implicit none
include "fftw3.f"
   !----------------------------------------------------------------------------
   integer :: nsize
   !----------------------------------------------------------------------------
   nsize = ncol + ncol - 1
   !
   ! read in data
   inquire(file=infile, exist=fexst)
   if ( .not.fexst ) stop 'Input file not found!'
   !
   write(*,'(/,5x,"To read data from file: ",A)') trim(infile)
   !
   NMax = 100000
   NINC = 20000
   allocate( data_in(nsize,NMax), data_read(icmax) )
   open(10, file= infile, action='read', iostat=ioerr)
   read(10, '(A)', iostat=ioerr) oneline
   npts = 0
   do while (ioerr.eq.0)
      read(oneline, *, iostat=rderr) data_read
      do while (rderr.eq.0)
         npts = npts + 1
         if (npts.gt.NMax) then
            if ( allocated(data_tmp) ) deallocate(data_tmp)
            allocate(data_tmp(nsize,NMax))
            data_tmp = data_in
            deallocate(data_in)
            NMax = NMax + NINC
            allocate( data_in(nsize,NMax) )
            data_in(:,1:NMax-NINC) = data_tmp
         endif
         data_in(1, npts) = data_read(icol(1))
         II = 2
         do i = 2, ncol
            data_in(II,  npts) = data_read(icol(i)  )
            data_in(II+1,npts) = data_read(icol(i)+1)
            II = II + 2
         enddo
         read(10, *, iostat=rderr) data_read
      enddo
      read(10, '(A)', iostat=ioerr) oneline
   enddo
   close(10)
   if ( allocated(data_tmp) ) deallocate( data_tmp )
   write(*,'(5x,"Total # of data read from ",A,":",I10)') trim(infile), npts
   if (npts.lt.1) stop
   !
   ! initialization
   if (nacf.lt.1) nacf = npts
   nfft = npts + nacf
   !
   if ( allocated( fft_cmplx_in )) deallocate( fft_cmplx_in  )
   if ( allocated( fft_out)) deallocate( fft_out )
   allocate( fft_cmplx_in(nfft), fft_out(nfft) )
   !
   ! subtract average value if required
   if ( fave ) then
      do II = 2, nsize
         ave           = sum(data_in(II,:))/dble(npts)
         data_in(II,:) = data_in(II,:) - ave
      enddo
   endif
   !
   call dfftw_plan_dft_1d(pFW, nfft, fft_cmplx_in, fft_out, FFTW_FORWARD,FFTW_ESTIMATE)
   call dfftw_plan_dft_1d(pBW, nfft, fft_out, fft_cmplx_in, FFTW_BACKWARD,FFTW_ESTIMATE)
   !
   r_nfft  = 1.D0 / sqrt(dble(nfft))
   !
   do II = 2, ncol
      do i =1, npts
         fft_cmplx_in(i) = complex(data_in(2*II-2,i), data_in(2*II-1,i) ) * r_nfft
      enddo
      fft_cmplx_in(npts+1:nfft) = complex(0.D0, 0.D0)
      !
      ! do fft to get auto-correlation
      call dfftw_execute_dft(pFW, fft_cmplx_in, fft_out)
      forall (i=1:nfft) fft_out(i) = conjg(fft_out(i))*fft_out(i)
      call dfftw_execute_dft(pBW, fft_out, fft_cmplx_in)
      !
      do i = 1, nacf
         data_in(2*II-2,i) = real(fft_cmplx_in(i))/dble(npts-i+1)
         data_in(2*II-1,i) = imag(fft_cmplx_in(i))/dble(npts-i+1)
      enddo
   enddo
   !
   call dfftw_destroy_plan(pFW)
   call dfftw_destroy_plan(pBW)
   !
   ! write out auto-correlation
   open(11, file=outfile, action='write')
   fmtstr = '("# time-lag",10x,??(4x,"acf_of_col_",I2.2,23x))'
   write(fmtstr(19:20),'(I2.2)') ncol-1
   write(11,fmtstr) icol(2:ncol)
   fmtstr = '(G19.10,??(1x,G19.10))'
   write(fmtstr(9:10),'(I2.2)') nsize-1
   data_in(1,1:nacf) = data_in(1,1:nacf) - data_in(1,1)
   do i = 1, nacf
      write(11, fmtstr) data_in(:,i)
   enddo
   close(11)
   write(*,'(/,5x,"Autocorrelation info written to file: ",A,/)') trim(outfile)
   !
end subroutine
