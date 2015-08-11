subroutine real_fft_2_ccf
use sharedvars
implicit none
include "fftw3.f"
   !----------------------------------------------------------------------------
   ! read in data
   inquire(file=infile, exist=fexst)
   if ( .not.fexst ) stop 'Input file not found!'
   !
   write(*,'(/,5x,"To read data from file: ",A)') trim(infile)
   !
   NMax = 100000
   NINC = 20000
   allocate( data_in(ncol,NMax), data_read(icmax) )
   open(10, file= infile, action='read', iostat=ioerr)
   read(10, '(A)', iostat=ioerr) oneline
   npts = 0
   do while (ioerr.eq.0)
      read(oneline, *, iostat=rderr) data_read
      do while (rderr.eq.0)
         npts = npts + 1
         if (npts.gt.NMax) then
            if ( allocated(data_tmp) ) deallocate(data_tmp)
            allocate(data_tmp(ncol,NMax))
            data_tmp = data_in
            deallocate(data_in)
            NMax = NMax + NINC
            allocate( data_in(ncol,NMax) )
            data_in(:,1:NMax-NINC) = data_tmp
         endif
         forall (i=1:ncol) data_in(i,npts) = data_read(icol(i))
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
   if ( allocated( fft_in   )) deallocate( fft_in   )
   if ( allocated( fft_out  )) deallocate( fft_out  )
   if ( allocated( fft_out2 )) deallocate( fft_out2 )
   allocate( fft_in(nfft), fft_out(nfft/2+1), fft_out2(nfft/2+1) )
   !
   ! subtract average value if required
   if ( fave ) then
      do II = 2, ncol
         ave           = sum(data_in(II,:))/dble(npts)
         data_in(II,:) = data_in(II,:) - ave
      enddo
   endif
   !
   call dfftw_plan_dft_r2c_1d(pFW, nfft, fft_in,  fft_out, FFTW_ESTIMATE)
   call dfftw_plan_dft_c2r_1d(pBW, nfft, fft_out, fft_in,  FFTW_ESTIMATE)
   !
   r_nfft  = 1.D0 / sqrt(dble(nfft))
   !
   fft_in(1:npts) = data_in(2,1:npts) * r_nfft
   fft_in(npts+1:nfft) = 0.D0
   call dfftw_execute(pFW)
   fft_out2 = fft_out
   !
   fft_in(1:npts) = data_in(3,1:npts) * r_nfft
   fft_in(npts+1:nfft) = 0.D0
   call dfftw_execute(pFW)
   !
   forall (i=1:nfft/2+1) fft_out(i) = conjg(fft_out2(i))*fft_out(i)
   call dfftw_execute(pBW)
   !
   forall (i=1:nacf) data_in(2,i) = fft_in(i)/dble(npts-i+1)
   !
   call dfftw_destroy_plan(pFW)
   call dfftw_destroy_plan(pBW)
   !
   ! write out cross-correlation
   open(11, file=outfile, action='write')
   fmtstr = '("# time-lag",10x,"ccf_of_col_",I2.2,"_and_",I2.2)'
   write(11,fmtstr) icol(2:3)
   fmtstr = '(G19.10,1x,G19.10)'
   data_in(1,1:nacf) = data_in(1,1:nacf) - data_in(1,1)
   do i = 1, nacf
      write(11, fmtstr) data_in(1:2,i)
   enddo
   close(11)
   write(*,'(/,5x,"Cross-correlation info written to file: ",A,/)') trim(outfile)
   !
end subroutine
