!-----------------------------------------------------------------------
! Interface file for abpush2.f
      module abpush2_h
      implicit none
!
      interface
         subroutine AGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx, &
     &ny,nxv,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine AGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine EAGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop& 
     &,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      end module
