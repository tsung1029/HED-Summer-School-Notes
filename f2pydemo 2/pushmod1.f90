!-----------------------------------------------------------------------
!
      module wpush1
      implicit none
!
      interface
         subroutine PUSH1ZF(part,dt,ek,idimp,nop,nx)
         implicit none
         integer, intent(in) ::  nop, idimp, nx
         real, intent(in) :: dt
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine wpush1zf(part,dt,ek,nx)
! push free-streaming particles
      implicit none
      integer, intent(in) :: nx
      real, intent(in) :: dt
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
! call low level procedure
      call PUSH1ZF(part,dt,ek,idimp,nop,nx)
      end subroutine
!
      end module
