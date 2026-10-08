module solver_mod
  use data_mod
  implicit none
  
  integer, allocatable :: IPIV(:)

contains

  subroutine RSOLVER(A, B, N, KKK)
    ! Replace manual RSOLVER with LAPACK DGETRF and DGETRS
    real(8), intent(inout) :: A(:,:)
    real(8), intent(inout) :: B(:)
    integer, intent(in) :: N, KKK
    integer :: INFO
    
    if (.not. allocated(IPIV)) then
       allocate(IPIV(N))
    end if
    
    if (KKK == 1) then
       ! Factorize
       call DGETRF(N, N, A, size(A,1), IPIV, INFO)
       if (INFO /= 0) then
          write(*,*) 'SOLVER: DGETRF failed with INFO = ', INFO
          stop 'SOLVER: DGETRF failure'
       end if
    end if
    
    ! Solve
    call DGETRS('N', N, 1, A, size(A,1), IPIV, B, size(B), INFO)
    if (INFO /= 0) then
       write(*,*) 'SOLVER: DGETRS failed with INFO = ', INFO
       stop 'SOLVER: DGETRS failure'
    end if
  end subroutine RSOLVER

end module solver_mod
