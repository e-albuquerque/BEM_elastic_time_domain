module io_mod
  use data_mod
  implicit none
contains

  subroutine INPUT12(INP, IPR)
    integer, intent(in) :: INP, IPR
    character*80 TITLE
    integer :: I, J, NN, NG

    write(IPR,'('' '',79(''*''))')
    read(INP,'(A)') TITLE
    write(IPR,'(A)') TITLE

    read(INP,*) NE, NFUNC_val, INTIMP, NSTEP, AT, GE, XNU, RO
    if(INTIMP == 0) INTIMP = NSTEP + 1
    CS = sqrt(GE/RO)
    CP = CS*sqrt((2.d0-2.d0*XNU)/(1.d0-2.d0*XNU))
    
    write(IPR,100) NE,NFUNC_val,INTIMP,NSTEP,AT,CP,CS,GE,XNU,RO
100 FORMAT(//2X,'DATA'//2X, &
           'NUMBER OF BOUNDARY ELEMENTS = ',I3/2X, &
           'NUMBER OF TIME FUNCTIONS    = ',I3/2X, &
           'INTERVAL OF PRINTING        = ',I3/2X, &
           'NUMBER OF TIME STEPS        = ',I3/2X, &
           'TIME INTERVAL               = ',E14.7/2X, &
           'P-WAVE VELOCITY             = ',E14.7/2X, &
           'S-WAVE VELOCITY             = ',E14.7/2X, &
           'SHEAR MODULUS               = ',E14.7/2X, &
           'POISSON RATIO               = ',E14.7/2X, &
           'DENSITY                     = ',E14.7)

    N_nodes = 2 * NE
    NN = 2 * N_nodes
    NG = 6 * NE

    ! Allocate global arrays
    allocate(X(N_nodes), Y(N_nodes))
    allocate(G(NN, NG), H(NN, NN))
    G = 0.0d0; H = 0.0d0
    allocate(A1(NN, NN), B1(NN, NG))
    A1 = 0.0d0; B1 = 0.0d0
    allocate(U(NN, NSTEP), T(NG, NSTEP))
    U = 0.0d0; T = 0.0d0
    allocate(BC(NG), IFUNC(NG), KODE(NG), FI(NN))
    BC = 0.0d0; IFUNC = 0; KODE = 0; FI = 0.0d0
    allocate(NPOIF(NFUNC_val), FSCALE(NFUNC_val))
    allocate(FUNC(NFUNC_val, 1000, 2))
    FUNC = 0.0d0

    write(IPR,200)
200 FORMAT(//2X,'BOUNDARY NODE COORDINATES'///4X, &
           'NODE',10X,'X',18X,'Y'/)
    read(INP,*) (X(I), Y(I), I=1, N_nodes)
    do I = 1, N_nodes
       write(IPR,300) I, X(I), Y(I)
    end do
300 FORMAT(5X,I3,2(5X,E14.7))

    write(IPR,400)
400 FORMAT(//2X,'TIME FUNCTIONS')
    do I = 1, NFUNC_val
       read(INP,*) NPOIF(I), FSCALE(I)
       read(INP,*) (FUNC(I,J,1), FUNC(I,J,2), J=1, NPOIF(I))
    end do
    do I = 1, NFUNC_val
       write(IPR,500) I, NPOIF(I), FSCALE(I)
       do J = 1, NPOIF(I)
          write(IPR,600) J, FUNC(I,J,1), FUNC(I,J,2)
       end do
    end do
500 FORMAT(/2X,'FUNCTION NUMBER  = ',I3/2X, &
           'NUMBER OF POINTS = ',I3/2X, &
           'FUNCTION SCALE   = ',E14.7//4X, &
           'POINT',8X,'TIME',11X,'VALUE')
600 FORMAT(5X,I3,4X,E14.7,2X,E14.7)

    write(IPR,700)
    do I = 1, NE
       read(INP,*) (KODE(6*I-6+J), BC(6*I-6+J), IFUNC(6*I-6+J), J=1,6)
       write(IPR,800) I, (BC(6*I-6+J), KODE(6*I-6+J), IFUNC(6*I-6+J), J=1,6)
    end do
700 FORMAT(//2X,'BOUNDARY CONDITIONS'//6X, &
           '------FIRST NODE------',2X,'-----SECOND NODE------',2X, &
           '------THIRD NODE------'/1X,'ELE', &
           3(2X,'X DIR. C F  Y DIR. C F'))
800 FORMAT(1X,I3,3(F8.2,1X,I1,I2,F8.2,1X,I1,I2))

    read(INP,*) mid_node_number

  end subroutine INPUT12

  subroutine OUTPUT12(ISTEP, TIME_VAL, IPR)
    integer, intent(in) :: ISTEP, IPR
    real(8), intent(in) :: TIME_VAL
    integer :: I, J, JB, JA, KODEB, KODEA, NN, NG
    integer :: trac_x_idx, trac_y_idx
    integer, save :: i_trac = 0

    NN = 2 * N_nodes
    NG = 6 * NE

    do J = 1, NN
       JB = J - 2 + 2 * ((J + 1) / 4)
       if (JB < 1) JB = JB + NG
       JA = J + 2 * ((J - 1) / 4)
       if (JA > NG) JA = JA - 2
       KODEB = KODE(JB)
       KODEA = KODE(JA)

       if (KODEB == 1 .and. KODEA == 1) then
          U(J, ISTEP) = FI(J)
       else if (KODEB == 0 .and. KODEA == 0) then
          T(JB, ISTEP) = GE * FI(J)
          T(JA, ISTEP) = GE * FI(J)
       else if (KODEB == 0 .and. KODEA == 1) then
          T(JB, ISTEP) = GE * FI(J)
       else if (KODEB == 1 .and. KODEA == 0) then
          T(JA, ISTEP) = GE * FI(J)
       else
          stop ' ERROR 3 IN BOUNDARY CONDITIONS'
       end if
    end do

    ! Open trac file once
    if (i_trac == 0) then
       open(unit=20, file='disp_trac.dat', status='replace')
       i_trac = 1
    end if

    ! Save tracking data to disp_trac.dat for mid_node_number only
    I = mid_node_number
    if (I == 1) then
        trac_x_idx = 1; trac_y_idx = 2
    else if (mod(I, 2) == 0) then
        trac_x_idx = 3 * I - 3
        trac_y_idx = trac_x_idx + 1
    else
        trac_x_idx = 3 * I - 2
        trac_y_idx = trac_x_idx + 1
    end if
    
    write(20, '(5(ES14.6, 2X))') &
       TIME_VAL, U(2*I-1, ISTEP), U(2*I, ISTEP), T(trac_x_idx, ISTEP), T(trac_y_idx, ISTEP)

    if (mod(ISTEP, INTIMP) /= 0 .and. ISTEP /= 1) return

    write(IPR,100) TIME_VAL
    do I = 1, N_nodes
       write(IPR,200) X(I), Y(I), U(2*I-1, ISTEP), U(2*I, ISTEP)
    end do
    write(IPR,300)
    write(IPR,400) X(1), Y(1), T(6*NE-1, ISTEP), T(6*NE, ISTEP), T(1, ISTEP), T(2, ISTEP)
    write(IPR,400) X(2), Y(2), T(3, ISTEP), T(4, ISTEP), T(3, ISTEP), T(4, ISTEP)
    do I = 2, NE
       write(IPR,400) X(2*I-1), Y(2*I-1), T(6*I-7, ISTEP), T(6*I-6, ISTEP), T(6*I-5, ISTEP), T(6*I-4, ISTEP)
       write(IPR,400) X(2*I),   Y(2*I),   T(6*I-3, ISTEP), T(6*I-2, ISTEP), T(6*I-3, ISTEP), T(6*I-2, ISTEP)
    end do

100 FORMAT(' ',79('*')//' RESULTS AT TIME = ',E12.5// &
           '  BOUNDARY NODES'//10X,'X',17X,'Y',11X,'DISPLACEMENT X', &
           4X,'DISPLACEMENT Y'/)
200 FORMAT(4(4X,E14.7))
300 FORMAT(///28X,'TRACTION X',3X,'TRACTION Y',3X,'TRACTION X',3X, &
           'TRACTION Y'/7X,'X',12X,'Y',7X,'BEFORE NODE',2X, &
           'BEFORE NODE',2X,'AFTER NODE',3X,'AFTER NODE'/)
400 FORMAT(6(1X,E12.5))
  end subroutine OUTPUT12

end module io_mod
