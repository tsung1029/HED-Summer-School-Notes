!-----------------------------------------------------------------------
!
      module wbpush2
!
! Fortran90 wrappers to 2-1/2d PIC library bpush2.f
! wdistr2h calculates initial particle co-ordinates with uniform plasma
!          maxwellian velocities
!          calls DISTR2H
! wgbpush2 non-relativistic push particles with leap-frog
!          calls GBPUSH23L
! wgrbpush2 relativistic push particles with leap-frog
!           calls GRBPUSH23L
! wgpost2 deposit charge
!         calls GPOST2L
! wgjpost2 deposit non-relativistic current
!          calls GJPOST2L
! wgrjpost2 deposit relativistic current
!           calls GRJPOST2L
! wdsortp2yl sort particles
!            calls DSORTP2YL
! wcguard2 copy guard cells for periodic 2-1/2d vector data
!          calls BGUARD2L
! wacguard1 add guard cells for periodic 2-1/2d vector data
!          calls ACGUARD2L
! waguard1 add guard cells for periodic 2d scalar data
!          calls AGUARD2L
! pois2_init calculates table needed by 2d poisson solver
!            calls POIS23
! wpois2 poisson solver for periodic 2d electric field
!        calls POIS23
! wcuperp2 calculates the transverse current in fourier space
!          calls CUPERP2
! wibpois2 solves 2d poisson's equation for unsmoothed magnetic
!          field
!          calls IBPOIS23
! wmaxwel2 solves 2d maxwell's equation for unsmoothed transverse
!          electric and magnetic fields using verlet algorithm
!          calls MAXWEL2
! wemfield2 adds and smooths or copies and smooths complex vector fields
!           in fourier space
!           calls EMFIELD2
! fft2_init calculates tables needed by 2d FFTs
!           calls WFFT2RINIT
! wfft2r wrapper function for scalar 2d real/complex FFT
!        calls WFFT2RX
! w3fft2r wrapper function for vector 2d real/complex FFT
!         calls WFFT2R3
! fprecision determines if default reals are actually doubles
! written by viktor k. decyk, ucla
! copyright 2021, regents of the university of california
! update: july 21, 2023
!
      use bpush2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wdistr2h(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,nx,ny)
! calculates initial particle co-ordinates with uniform plasma and
! maxwellian velocities
      implicit none
      integer, intent(in) :: npx, npy, nx, ny
      double precision, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      double precision, dimension(:,:), intent(inout) :: part
! local data
      integer :: ipbc = 1
      integer :: idimp, nop
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
! call low level procedure
      call DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgbpush2(part,fxy,bxy,qbm,dt,dtc,ek,nx,ny)
! push non-relativistic particles with leap-frog
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: qbm, dt, dtc
      double precision, intent(inout) :: ek
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:,:,:), intent(in) :: fxy, bxy
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fxy,2); nyv = size(fxy,3)
! call low level procedure
      call GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgrbpush2(part,fxy,bxy,qbm,dt,dtc,ci,ek,nx,ny)
! push relativistic particles with leap-frog
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: qbm, dt, dtc, ci
      double precision, intent(inout) :: ek
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:,:,:), intent(in) :: fxy, bxy
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fxy,2); nyv = size(fxy,3)
! call low level procedure
      call GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgpost2(part,q,qm)
! deposit charge
      implicit none
      double precision, intent(in) :: qm
      double precision, dimension(:,:), intent(in) :: part
      double precision, dimension(:,:), intent(inout) :: q
! local data
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(q,1); nyv = size(q,2)
! call low level procedure
      call GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgjpost2(part,cu,qm,dt,nx,ny)
! deposit non-relativistic current
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: qm, dt
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(cu,2); nyv = size(cu,3)
! call low level procedure
      call GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgrjpost2(part,cu,qm,dt,ci,nx,ny)
! deposit relativistic current
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: qm, dt, ci
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(cu,2); nyv = size(cu,3)
! call low level procedure
      call GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wdsortp2yl(parta,partb,npic)
! sort particles
      implicit none
      double precision, dimension(:,:), intent(in) :: parta
      double precision, dimension(:,:), intent(inout) :: partb
      integer, dimension(:), intent(inout) :: npic
! local data
      integer :: idimp, nop, ny1
! extract dimensions
      idimp = size(parta,1); nop = size(parta,2)
      ny1 = size(npic,1)
! call low level procedure
      call DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wcguard2(fxy,nx,ny)
! copy guard cells for periodic 2-1/2d vector data
      implicit none
      integer, intent(in) :: nx, ny
      double precision, dimension(:,:,:), intent(inout) :: fxy
! local data
      integer :: nxe, nye
! extract dimensions
      nxe = size(fxy,2); nye = size(fxy,3)
! call low level procedure
      call BGUARD2L(fxy,nx,ny,nxe,nye)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wacguard2(cu,nx,ny)
! add guard cells for periodic 2-1/2d vector data
      implicit none
      integer, intent(in) :: nx, ny
      double precision, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: nxe, nye
! extract dimensions
      nxe = size(cu,2); nye = size(cu,3)
! call low level procedure
      call ACGUARD2L(cu,nx,ny,nxe,nye)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine waguard2(q,nx,ny)
! add guard cells for periodic 2d scalar data
      implicit none
      integer, intent(in) :: nx, ny
      double precision, dimension(:,:), intent(inout) :: q
! local data
      integer :: nxe, nye
! extract dimensions
      nxe = size(q,1); nye = size(q,2)
! call low level procedure
      call AGUARD2L(q,nx,ny,nxe,nye)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pois2_init(ffc,ax,ay,affp,nx,ny)
! calculates table needed by 2d poisson solver
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: ax, ay, affp
      double complex, dimension(:,:), intent(inout) :: ffc
! local data
      integer :: isign = 0
      integer :: nxvh, nyv, nxhd, nyhd
      double precision :: we
      double precision, dimension(2,1) :: q
      double precision, dimension(1,1,1) :: fxy
! extract dimensions
      nxvh = size(q,1)/2; nyv = size(q,2)
      nxhd = size(ffc,1); nyhd = size(ffc,2)
! call low level procedure
      call POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wpois2(q,fxy,ffc,we,nx,ny)
! poisson solver for periodic 2d electric field
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(inout) :: we
      double precision, dimension(:,:), intent(in)  :: q
      double precision, dimension(:,:,:), intent(inout) :: fxy
      double complex, dimension(:,:), intent(inout) :: ffc
! local data
      integer :: isign = -1
      integer :: nxvh, nyv, nxhd, nyhd
      double precision :: ax, ay, affp
! extract dimensions
      nxvh = size(q,1)/2; nyv = size(q,2)
      nxhd = size(ffc,1); nyhd = size(ffc,2)
! call low level procedure
      call POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wcuperp2(cu,nx,ny)
! calculates the transverse current in fourier space
      implicit none
      integer, intent(in) :: nx, ny
      double precision, dimension(:,:,:), intent(inout) :: cu
! local data
      integer :: nxvh, nyv
! extract dimensions
      nxvh = size(cu,2)/2; nyv = size(cu,3)
! call low level procedure
      call CUPERP2(cu,nx,ny,nxvh,nyv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wibpois2(cu,bxy,ffc,ci,wm,nx,ny)
! solves 2d poisson's equation for unsmoothed magnetic field
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: ci
      double precision, intent(inout) :: wm
      double precision, dimension(:,:,:), intent(in) :: cu
      double complex, dimension(:,:,:), intent(inout) :: bxy
      double complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nxvh, nyv, nxhd, nyhd
! extract dimensions
      nxvh = size(cu,2)/2; nyv = size(cu,3)
      nxhd = size(ffc,1); nyhd = size(ffc,2)
! call low level procedure
      call IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmaxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny)
! solves 2d maxwell's equation for unsmoothed transverse electric and
! magnetic fields using verlet algorithm
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: ci, dt
      double precision, intent(inout) :: wf, wm
      double complex, dimension(:,:,:), intent(inout) :: exy, bxy
      double precision, dimension(:,:,:), intent(in)  :: cu
      double complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nxvh, nyv, nxhd, nyhd
! extract dimensions
      nxvh = size(cu,2)/2; nyv = size(cu,3)
      nxhd = size(ffc,1); nyhd = size(ffc,2)
! call low level procedure
      call MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wemfield2(fxy,exy,ffc,isign,nx,ny)
! adds and smooths or copies and smooths complex vector fields in
! fourier space
      implicit none
      integer, intent(in) :: isign, nx, ny
      double precision, dimension(:,:,:), intent(inout) :: fxy
      double complex, dimension(:,:,:), intent(in) :: exy
      double complex, dimension(:,:), intent(in) :: ffc
! local data
      integer :: nxvh, nyv, nxhd, nyhd
! extract dimensions
      nxvh = size(fxy,2)/2; nyv = size(fxy,3)
      nxhd = size(ffc,1); nyhd = size(ffc,2)
! call low level procedure
      call EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fft2_init(mixup,sct,indx,indy)
! calculates tables needed by 2d FFTs
      implicit none
      integer, intent(in) :: indx, indy
      integer, dimension(:), intent(inout) :: mixup
      double complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhyd, nxyhd
! extract dimensions
      nxhyd = size(mixup,1); nxyhd = size(sct)
! call low level procedure
      call WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wfft2r(f,isign,mixup,sct,indx,indy)
! wrapper function for scalar 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy
      double precision, dimension(:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      double complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxhd, nyd, nxhyd, nxyhd
! extract dimensions
      nxhd = size(f,1)/2; nyd = size(f,2)
      nxhyd = size(mixup,1); nxyhd = size(sct)
! call low level procedure
      call WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine w3fft2r(f,isign,mixup,sct,indx,indy)
! wrapper function for vector 2d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx, indy
      double precision, dimension(:,:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      double complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxhd, nyd, nxhyd, nxyhd
! extract dimensions
      nxhd = size(f,2)/2; nyd = size(f,3)
      nxhyd = size(mixup,1); nxyhd = size(sct)
! call low level procedure-
      call WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
      end subroutine
!
!-----------------------------------------------------------------------
      integer function fprecision()
! function determines if default reals are actually doubles
      implicit none
      real :: prec
! ndprec = (0,1) = (no,yes) = (normal,autodouble) precision used
      if (digits(prec) > 24) then
         fprecision = 1
      else
         fprecision = 0
      endif
      end function fprecision
!
      end module
