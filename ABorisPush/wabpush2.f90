!-----------------------------------------------------------------------
!
      module wabpush2
!
! Fortran90 wrappers to 2-1/2d PIC library abpush2.f
! contains additional analytic pushers
! wagbpush2 non-relativistic analytic push procedure with leap-frog
!           calls AGBPUSH23L
! wagrbpush2 relativistic analytic push procedure with leap-frog
!            calls AGRBPUSH23L
! weagrbpush2 exact analytic push procedure with leap-frog
!            calls EAGRBPUSH23L
! written by viktor k. decyk, ucla
! copyright 2021, regents of the university of california
! update: july 21, 2023
!
      use abpush2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wagbpush2(part,fxy,bxy,qbm,dt,dtc,ek,nx,ny)
! non-relativistic analytic boris push procedure with leap-frog
      implicit none
      integer, intent(in) :: nx, ny
      real, intent(in) :: qbm, dt, dtc
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: fxy, bxy
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fxy,2); nyv = size(fxy,3)
! call low level procedure
      call AGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wagrbpush2(part,fxy,bxy,qbm,dt,dtc,ci,ek,nx,ny)
! relativistic analytic boris push procedure with leap-frog
      implicit none
      integer, intent(in) :: nx, ny
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: fxy, bxy
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fxy,2); nyv = size(fxy,3)
! call low level procedure
      call AGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny,nxv&
     &,nyv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine weagrbpush2(part,fxy,bxy,qbm,dt,dtc,ci,ek,nx,ny)
! exact analytic push procedure with leap-frog
      implicit none
      integer, intent(in) :: nx, ny
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: fxy, bxy
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv, nyv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fxy,2); nyv = size(fxy,3)
! call low level procedure
      call EAGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny,  &
     &nxv,nyv,ipbc)
      end subroutine
!
      end module
