module data_mod
  implicit none
  
  real(8), allocatable :: G(:,:)
  real(8), allocatable :: H(:,:)
  real(8), allocatable :: A1(:,:)
  real(8), allocatable :: B1(:,:)
  real(8), allocatable :: U(:,:)
  real(8), allocatable :: T(:,:)
  
  real(8), allocatable :: X(:), Y(:)
  integer, allocatable :: KODE(:), IFUNC(:)
  real(8), allocatable :: BC(:), FI(:)
  
  integer, allocatable :: NPOIF(:)
  real(8), allocatable :: FSCALE(:)
  real(8), allocatable :: FUNC(:,:,:)
  
  ! Global variables needed by multiple subroutines
  integer :: N_nodes, NE, NFUNC_val, INTIMP, NSTEP, mid_node_number
  real(8) :: AT, GE, XNU, RO, CP, CS

end module data_mod
