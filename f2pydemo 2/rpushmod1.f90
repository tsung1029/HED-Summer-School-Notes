!-----------------------------------------------------------------------
!
      module wrpush1
      use rpush1
      implicit none
!
      type (part1d), save, private :: wpart
!
      contains  
!
      subroutine set_part1d(qm,qbm,idimp,nop)
! this subroutine set a private descriptor for 1d particles
! wpart = part1d descriptor of particle data
      implicit none
      real, intent(in) :: qm, qbm
      integer, intent(in) :: idimp, nop
! create descriptor
      wpart%qm = qm; wpart%qbm = qbm
      wpart%idimp = idimp; wpart%nop = nop
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine w1rpush1zf(qm,qbm,idimp,nop,part,dt,ci,ek,nx)
! push free-streaming relativistic particles, unpacked arguments
      integer, intent(in) :: idimp, nop, nx
      real, intent(in) :: qm, qbm, dt, ci
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
! local data
      type (part1d) :: this
! set descriptor
      this%qm = qm; this%qbm = qbm
      this%idimp = idimp; this%nop = nop
! call low level procedure
      call rpush1zf(this,part,dt,ci,ek,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine w2rpush1zf(part,dt,ci,ek,nx)
! push free-streaming relativistic particles, using wpart
      integer, intent(in) :: nx
      real, intent(in) :: dt, ci
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
! call low level procedure
      call rpush1zf(wpart,part,dt,ci,ek,nx)
      end subroutine
!
      end module
