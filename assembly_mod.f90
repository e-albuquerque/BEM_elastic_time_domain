module assembly_mod
  use data_mod
  implicit none
contains

  subroutine GHMAT12(ISTEP, IAU)
    integer, intent(in) :: ISTEP, IAU
    real(8) :: CPT, DIST, C11, C12, C22
    integer :: I, J, J1, J2, J3, NODO, II, K, JH, JJ, KK, I1, I2, I5
    integer :: NN, NG_val, JB, JA, KODEB, KODEA
    real(8) :: HT(2,6), GT(2,6)

    NN = 2 * N_nodes
    NG_val = 6 * NE

    ! Initialize H and G
    !$omp parallel do default(shared) private(I, J)
    do I = 1, NN
       do J = 1, NN
          H(I,J) = 0.0d0
       end do
       do J = 1, NG_val
          G(I,J) = 0.0d0
       end do
    end do
    !$omp end parallel do

    CPT = CP * AT * ISTEP

    ! COMPUTE THE COEFFICIENTS OF G AND H MATRICES
    !$omp parallel do default(shared) private(I, J, J1, J2, J3, DIST, GT, HT, NODO, II, K, JH, JJ, KK, C11, C12, C22, I1, I2, I5)
    do I = 1, N_nodes
       do J = 1, NE
          J1 = 2 * J - 1
          J2 = J1 + 1
          J3 = J2 + 1
          if (J3 > N_nodes) J3 = 1
          if ((I - J1) * (I - J2) * (I - J3) == 0 .and. ISTEP == 1) then
             NODO = I - J1 + 1
             if (I == 1 .and. J == NE) NODO = 3
             call LOCIN12(X(I), Y(I), X(J1), Y(J1), X(J2), Y(J2), X(J3), Y(J3), NODO, AT, CP, CS, XNU, GE, RO, GT, HT)
          else
             DIST = DISTP(X(I), Y(I), X(J1), Y(J1), X(J2), Y(J2), X(J3), Y(J3))
             if (DIST <= CPT) then
                call EXTIN12(X(I), Y(I), X(J1), Y(J1), X(J2), Y(J2), X(J3), Y(J3), ISTEP, AT, CP, CS, XNU, GE, RO, GT, HT)
             else
                ! No waves reached this element yet
                cycle
             end if
          end if

          II = 2 * I - 1
          do K = 1, 3
             JH = 4 * J + 2 * K - 5
             if (J == NE .and. K == 3) JH = 1
             JJ = 6 * J + 2 * K - 7
             KK = 2 * K - 1
             ! Using atomic for summation to avoid race conditions if two elements update same H
             ! Actually, the inner loop over elements J updates H(II, JH), JH is independent of I but dependent on J
             ! But H is updated by II which is 2*I-1. Since the parallel loop is over I, II is unique per thread!
             ! So NO race conditions on H and G writes!
             H(II,   JH)   = H(II,   JH)   + HT(1, KK)
             H(II+1, JH)   = H(II+1, JH)   + HT(2, KK)
             H(II,   JH+1) = H(II,   JH+1) + HT(1, KK+1)
             H(II+1, JH+1) = H(II+1, JH+1) + HT(2, KK+1)
             
             G(II,   JJ)   = GT(1, KK)
             G(II+1, JJ)   = GT(2, KK)
             G(II,   JJ+1) = GT(1, KK+1)
             G(II+1, JJ+1) = GT(2, KK+1)
          end do
       end do

       if (ISTEP == 1) then
          if (mod(I, 2) /= 0) then
             ! Odd node
             I1 = I - 2
             if (I == 1) I1 = N_nodes - 1
             I2 = I - 1
             if (I == 1) I2 = N_nodes
             I5 = I + 2
             if (I5 > N_nodes) I5 = 1
             call FTERM12(C11, C12, C22, X(I1), Y(I1), X(I2), Y(I2), X(I), Y(I), X(I+1), Y(I+1), X(I5), Y(I5), CP, CS)
          else
             ! Even node
             C11 = 0.5d0
             C22 = 0.5d0
             C12 = 0.0d0
          end if
          H(2*I-1, 2*I-1) = H(2*I-1, 2*I-1) + C11
          H(2*I-1, 2*I)   = H(2*I-1, 2*I)   + C12
          H(2*I,   2*I-1) = H(2*I,   2*I-1) + C12
          H(2*I,   2*I)   = H(2*I,   2*I)   + C22
       end if
    end do
    !$omp end parallel do

    ! Write the matrices in an auxiliary file
    write(IAU) H(1:NN, 1:NN), G(1:NN, 1:NG_val)

    if (ISTEP > 1) return

    ! Form A1 and B1 for the system of equations
    do J = 1, NN
       JB = J - 2 + 2 * ((J + 1) / 4)
       if (JB < 1) JB = JB + NG_val
       JA = J + 2 * ((J - 1) / 4)
       if (JA > NG_val) JA = JA - 2
       KODEB = KODE(JB)
       KODEA = KODE(JA)

       if (KODEB == 0 .and. KODEA == 0 .and. JB == JA) then
          do I = 1, NN
             A1(I,J) = -GE * G(I,JB)
             B1(I,JB) = -H(I,J)
          end do
       else if (KODEB == 0 .and. KODEA == 0) then
          do I = 1, NN
             A1(I,J) = -GE * (G(I,JB) + G(I,JA))
             B1(I,JB) = -H(I,J)
             B1(I,JA) = 0.0d0
          end do
       else if (KODEB /= 0 .and. KODEA /= 0) then
          do I = 1, NN
             A1(I,J) = H(I,J)
             B1(I,JB) = G(I,JB)
             B1(I,JA) = G(I,JA)
          end do
       else if (KODEB == 0 .and. KODEA == 1) then
          do I = 1, NN
             A1(I,J) = -GE * G(I,JB)
             B1(I,JB) = -H(I,J)
             B1(I,JA) = G(I,JA)
          end do
       else if (KODEB == 1 .and. KODEA == 0) then
          do I = 1, NN
             A1(I,J) = -GE * G(I,JA)
             B1(I,JB) = G(I,JB)
             B1(I,JA) = -H(I,J)
          end do
       end if
    end do
  end subroutine GHMAT12

  subroutine EXTIN12(XP, YP, X1, Y1, X2, Y2, X3, Y3, ISTEP, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val, GT, HT)
    real(8), intent(in) :: XP, YP, X1, Y1, X2, Y2, X3, Y3, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val
    integer, intent(in) :: ISTEP
    real(8), intent(out) :: GT(2,6), HT(2,6)
    real(8) :: GI(10), OME(10)
    integer :: I, J
    real(8) :: A, B, C, D, CPT, F1, F2, F3, XCO, YCO, RA, RD1, RD2, XJA, GH, ETA1, ETA2, RDN
    real(8) :: E1, E2, Z1, Z2, Z3

    GI = [0.9739065285d0, -0.9739065285d0, 0.8650633666d0, -0.8650633666d0, &
          0.6794095682d0, -0.6794095682d0, 0.4333953941d0, -0.4333953941d0, &
          0.1488743389d0, -0.1488743389d0]
    OME = [0.0666713443d0, 0.0666713443d0, 0.1494513491d0, 0.1494513491d0, &
           0.2190863625d0, 0.2190863625d0, 0.2692667193d0, 0.2692667193d0, &
           0.2955242247d0, 0.2955242247d0]

    HT = 0.0d0
    GT = 0.0d0

    A = X3 - 2.0d0*X2 + X1
    B = (X3 - X1) / 2.0d0
    C = Y3 - 2.0d0*Y2 + Y1
    D = (Y3 - Y1) / 2.0d0
    CPT = CP_val * AT_val * ISTEP

    do I = 1, 10
       F1 = GI(I) * (GI(I) - 1.0d0) * 0.5d0
       F2 = 1.0d0 - GI(I)**2
       F3 = GI(I) * (GI(I) + 1.0d0) * 0.5d0
       XCO = (X1 - XP)*F1 + (X2 - XP)*F2 + (X3 - XP)*F3
       YCO = (Y1 - YP)*F1 + (Y2 - YP)*F2 + (Y3 - YP)*F3
       RA = sqrt(XCO**2 + YCO**2)
       
       if (RA > CPT) cycle
       
       RD1 = XCO / RA
       RD2 = YCO / RA
       call FUNDSOLE(E1, E2, Z1, Z2, Z3, RA, ISTEP, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val)
       
       XJA = sqrt((GI(I)*A + B)**2 + (GI(I)*C + D)**2) * OME(I)
       
       GH = XJA * (E1 - E2 + 2.0d0*RD1*RD1*E2)
       GT(1,1) = GT(1,1) + F1 * GH
       GT(1,3) = GT(1,3) + F2 * GH
       GT(1,5) = GT(1,5) + F3 * GH
       
       GH = XJA * (2.0d0*RD1*RD2*E2)
       GT(1,2) = GT(1,2) + F1 * GH
       GT(1,4) = GT(1,4) + F2 * GH
       GT(1,6) = GT(1,6) + F3 * GH
       
       GH = XJA * (E1 - E2 + 2.0d0*RD2*RD2*E2)
       GT(2,2) = GT(2,2) + F1 * GH
       GT(2,4) = GT(2,4) + F2 * GH
       GT(2,6) = GT(2,6) + F3 * GH
       
       ETA1 = (C * GI(I) + D)
       ETA2 = -(A * GI(I) + B)
       RDN = RD1 * ETA1 + RD2 * ETA2
       
       GH = OME(I) * (RDN * (Z1 + 2.0d0*RD1*RD1*Z2) + ETA1*RD1*(Z1 + Z3)) / RA
       HT(1,1) = HT(1,1) + F1 * GH
       HT(1,3) = HT(1,3) + F2 * GH
       HT(1,5) = HT(1,5) + F3 * GH
       
       GH = OME(I) * (RDN * 2.0d0*RD1*RD2*Z2 + ETA1*RD2*Z1 + ETA2*RD1*Z3) / RA
       HT(1,2) = HT(1,2) + F1 * GH
       HT(1,4) = HT(1,4) + F2 * GH
       HT(1,6) = HT(1,6) + F3 * GH
       
       GH = OME(I) * (RDN * 2.0d0*RD1*RD2*Z2 + ETA1*RD2*Z3 + ETA2*RD1*Z1) / RA
       HT(2,1) = HT(2,1) + F1 * GH
       HT(2,3) = HT(2,3) + F2 * GH
       HT(2,5) = HT(2,5) + F3 * GH
       
       GH = OME(I) * (RDN * (Z1 + 2.0d0*RD2*RD2*Z2) + ETA2*RD2*(Z1 + Z3)) / RA
       HT(2,2) = HT(2,2) + F1 * GH
       HT(2,4) = HT(2,4) + F2 * GH
       HT(2,6) = HT(2,6) + F3 * GH
    end do
    
    GT(2,1) = GT(1,2)
    GT(2,3) = GT(1,4)
    GT(2,5) = GT(1,6)
  end subroutine EXTIN12

  subroutine LOCIN12(XP, YP, X1, Y1, X2, Y2, X3, Y3, NODO, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val, GT, HT)
    real(8), intent(in) :: XP, YP, X1, Y1, X2, Y2, X3, Y3, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val
    integer, intent(in) :: NODO
    real(8), intent(out) :: GT(2,6), HT(2,6)
    real(8), parameter :: PI = 3.1415926535d0
    real(8) :: GI(10), OME(10), GIL(10), OMEL(10)
    real(8) :: A, B, C, D, RES1, RES2, RES3, A1, A0, XJAL, CTU, CTT
    real(8) :: GIA, XJA, S1, S2, S3, F1, F2, F3, XCO, YCO, RA, RD1, RD2
    real(8) :: E1, E2, Z1, Z2, Z3, ETA, GH, ETA1, ETA2, RDN, GHA, R1, R3
    integer :: I, J

    GI = [0.9739065285d0, -0.9739065285d0, 0.8650633666d0, -0.8650633666d0, &
          0.6794095682d0, -0.6794095682d0, 0.4333953941d0, -0.4333953941d0, &
          0.1488743389d0, -0.1488743389d0]
    OME = [0.0666713443d0, 0.0666713443d0, 0.1494513491d0, 0.1494513491d0, &
           0.2190863625d0, 0.2190863625d0, 0.2692667193d0, 0.2692667193d0, &
           0.2955242247d0, 0.2955242247d0]
    GIL = [0.0090426309d0, 0.0539712662d0, 0.1353118246d0, 0.2470524162d0, &
           0.3802125396d0, 0.5237923179d0, 0.6657752055d0, 0.7941904160d0, &
           0.8981610912d0, 0.9688479887d0]
    OMEL = [0.1209551319d0, 0.1863635425d0, 0.1956608732d0, 0.1735771421d0, &
            0.1356956729d0, 0.0936467585d0, 0.0557877273d0, 0.0271598109d0, &
            0.0095151826d0, 0.0016381576d0]

    A = X3 - 2.0d0*X2 + X1
    B = (X3 - X1) / 2.0d0
    C = Y3 - 2.0d0*Y2 + Y1
    D = (Y3 - Y1) / 2.0d0
    RES1 = 0.0d0
    RES2 = 0.0d0
    RES3 = 0.0d0

    if (NODO == 1) then
       A1 = 0.5d0; A0 = 0.5d0; XJAL = 2.0d0; RES1 = 1.0d0
    else if (NODO == 2) then
       A1 = 1.0d0; A0 = 0.0d0; XJAL = 1.0d0; RES2 = 1.0d0
    else if (NODO == 3) then
       A1 = -0.5d0; A0 = 0.5d0; XJAL = 2.0d0; RES3 = 1.0d0
    end if

    CTU = (3.0d0 - 4.0d0*XNU_val) / 8.0d0 / PI / GE_val / (1.0d0 - XNU_val)
    CTT = (1.0d0 - 2.0d0*XNU_val) / 4.0d0 / PI / (1.0d0 - XNU_val)

    HT = 0.0d0
    GT = 0.0d0

    do I = 1, 10
       GIA = (GIL(I) - A0) / A1
       XJA = CTU * XJAL * sqrt((GIA*A + B)**2 + (GIA*C + D)**2) * OMEL(I)
       S1 = 0.5d0 * GIA * (GIA - 1.0d0) * XJA
       S2 = (1.0d0 - GIA**2) * XJA
       S3 = 0.5d0 * GIA * (GIA + 1.0d0) * XJA

       if (NODO == 2) then
          GIA = -GIL(I)
          XJA = CTU * sqrt((GIA*A + B)**2 + (GIA*C + D)**2) * OMEL(I)
          S1 = S1 + 0.5d0 * GIA * (GIA - 1.0d0) * XJA
          S2 = S2 + (1.0d0 - GIA**2) * XJA
          S3 = S3 + 0.5d0 * GIA * (GIA + 1.0d0) * XJA
       end if

       F1 = 0.5d0 * GI(I) * (GI(I) - 1.0d0)
       F2 = 1.0d0 - GI(I)**2
       F3 = 0.5d0 * GI(I) * (GI(I) + 1.0d0)
       XCO = (X1 - XP)*F1 + (X2 - XP)*F2 + (X3 - XP)*F3
       YCO = (Y1 - YP)*F1 + (Y2 - YP)*F2 + (Y3 - YP)*F3
       RA = sqrt(XCO*XCO + YCO*YCO)
       RD1 = XCO / RA
       RD2 = YCO / RA
       XJA = sqrt((GI(I)*A + B)**2 + (GI(I)*C + D)**2) * OME(I)
       
       call FUNDSOLE(E1, E2, Z1, Z2, Z3, RA, 1, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val)
       
       ETA = abs(A1 * GI(I) + A0)
       E1 = E1 + CTU * log(ETA)
       
       GH = XJA * (E1 - E2 + 2.0d0*RD1*RD1*E2)
       GT(1,1) = GT(1,1) + F1 * GH + S1
       GT(1,3) = GT(1,3) + F2 * GH + S2
       GT(1,5) = GT(1,5) + F3 * GH + S3
       
       GH = XJA * (2.0d0*RD1*RD2*E2)
       GT(1,2) = GT(1,2) + F1 * GH
       GT(1,4) = GT(1,4) + F2 * GH
       GT(1,6) = GT(1,6) + F3 * GH
       
       GH = XJA * (E1 - E2 + 2.0d0*RD2*RD2*E2)
       GT(2,2) = GT(2,2) + F1 * GH + S1
       GT(2,4) = GT(2,4) + F2 * GH + S2
       GT(2,6) = GT(2,6) + F3 * GH + S3
       
       ETA1 = (C * GI(I) + D)
       ETA2 = -(A * GI(I) + B)
       RDN = RD1 * ETA1 + RD2 * ETA2
       
       GH = OME(I) * (RDN * (Z1 + 2.0d0*RD1*RD1*Z2) + ETA1*RD1*(Z1 + Z3)) / RA
       HT(1,1) = HT(1,1) + F1 * GH
       HT(1,3) = HT(1,3) + F2 * GH
       HT(1,5) = HT(1,5) + F3 * GH
       
       GH = OME(I) * (RDN * (Z1 + 2.0d0*RD2*RD2*Z2) + ETA2*RD2*(Z1 + Z3)) / RA
       HT(2,2) = HT(2,2) + F1 * GH
       HT(2,4) = HT(2,4) + F2 * GH
       HT(2,6) = HT(2,6) + F3 * GH
       
       Z1 = Z1 + CTT
       Z3 = Z3 - CTT
       
       GH = OME(I) * (RDN * 2.0d0*RD1*RD2*Z2 + ETA1*RD2*Z1 + ETA2*RD1*Z3) / RA
       GHA = OME(I) * CTT * (ETA1*RD2 - ETA2*RD1) / RA
       HT(1,2) = HT(1,2) + F1 * GH + (RES1 - F1) * GHA
       HT(1,4) = HT(1,4) + F2 * GH + (RES2 - F2) * GHA
       HT(1,6) = HT(1,6) + F3 * GH + (RES3 - F3) * GHA
       
       GH = OME(I) * (RDN * 2.0d0*RD1*RD2*Z2 + ETA1*RD2*Z3 + ETA2*RD1*Z1) / RA
       HT(2,1) = HT(2,1) + F1 * GH - (RES1 - F1) * GHA
       HT(2,3) = HT(2,3) + F2 * GH - (RES2 - F2) * GHA
       HT(2,5) = HT(2,5) + F3 * GH - (RES3 - F3) * GHA
    end do

    GT(2,1) = GT(1,2)
    GT(2,3) = GT(1,4)
    GT(2,5) = GT(1,6)

    R1 = (X1 - XP)**2 + (Y1 - YP)**2
    R3 = (X3 - XP)**2 + (Y3 - YP)**2

    if (NODO == 1) then
       HT(1,2) = HT(1,2) - CTT * 0.5d0 * log(R3)
       HT(2,1) = HT(2,1) + CTT * 0.5d0 * log(R3)
    else if (NODO == 2) then
       HT(1,4) = HT(1,4) - CTT * 0.5d0 * log(R3 / R1)
       HT(2,3) = HT(2,3) + CTT * 0.5d0 * log(R3 / R1)
    else if (NODO == 3) then
       HT(1,6) = HT(1,6) + CTT * 0.5d0 * log(R1)
       HT(2,5) = HT(2,5) - CTT * 0.5d0 * log(R1)
    end if

  end subroutine LOCIN12
  subroutine FTERM12(C11, C12, C22, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, CP_val, CS_val)
    real(8), intent(out) :: C11, C12, C22
    real(8), intent(in) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, CP_val, CS_val
    real(8), parameter :: PI = 3.1415926535d0
    real(8) :: ZETA1, ZETA2, ZET, CTE

    ZETA1 = atan2(-0.5d0*X1 + 2.0d0*X2 - 1.5d0*X3, 0.5d0*Y1 - 2.0d0*Y2 + 1.5d0*Y3)
    ZETA2 = atan2( 1.5d0*X3 - 2.0d0*X4 + 0.5d0*X5, -1.5d0*Y3 + 2.0d0*Y4 - 0.5d0*Y5)
    ZET = ZETA2 - ZETA1
    
    if (abs(abs(ZET) - PI) < 1.0d-4) stop ' ERROR IN FREETERM12: COINCIDENT ELEMENTS'
    
    if (int(ZET/PI) < 0) then
       ZET = 0.5d0 * (-1.0d0 - ZET/PI)
    else if (int(ZET/PI) == 0) then
       ZET = 0.5d0 * (1.0d0 - ZET/PI)
    else
       ZET = 0.5d0 * (3.0d0 - ZET/PI)
    end if
    
    CTE = (1.0d0 - (CS_val/CP_val)**2) / 4.0d0 / PI
    C11 = ZET - CTE * (sin(2.0d0*ZETA1) - sin(2.0d0*ZETA2))
    C12 = CTE * (cos(2.0d0*ZETA1) - cos(2.0d0*ZETA2))
    C22 = ZET + CTE * (sin(2.0d0*ZETA1) - sin(2.0d0*ZETA2))
  end subroutine FTERM12

  subroutine FUNDSOLE(E1, E2, Z1, Z2, Z3, RA, ISTEP, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val)
    real(8), intent(out) :: E1, E2, Z1, Z2, Z3
    real(8), intent(in) :: RA, AT_val, CP_val, CS_val, XNU_val, GE_val, RO_val
    integer, intent(in) :: ISTEP
    real(8), parameter :: PI = 3.1415926535d0
    real(8) :: P0, P1, P2, S0, S1, S2, CPS, PS0, PS1, PS2
    real(8) :: RAIS0, RAIS1, RAIS2, CTE, CTZ, V1, V2, V3

    E1 = 0.0d0; E2 = 0.0d0; Z1 = 0.0d0; Z2 = 0.0d0; Z3 = 0.0d0
    
    P0 = CP_val * ISTEP * AT_val / RA
    if (P0 <= 1.0d0) return
    
    P1 = max(1.0d0, CP_val * (ISTEP - 1) * AT_val / RA)
    P2 = max(1.0d0, CP_val * (ISTEP - 2) * AT_val / RA)
    S0 = max(1.0d0, CS_val * ISTEP * AT_val / RA)
    S1 = max(1.0d0, CS_val * (ISTEP - 1) * AT_val / RA)
    S2 = max(1.0d0, CS_val * (ISTEP - 2) * AT_val / RA)
    
    CPS = 2.0d0 * (1.0d0 - XNU_val) / (1.0d0 - 2.0d0*XNU_val)
    PS0 = (P0/S0)**2 - 1.0d0
    PS1 = (P1/S1)**2 - 1.0d0
    PS2 = (P2/S2)**2 - 1.0d0
    
    RAIS0 = RAIS(P0) + RAIS(S0)
    RAIS1 = RAIS(P1) + RAIS(S1)
    RAIS2 = RAIS(P2) + RAIS(S2)
    
    CTE = 1.0d0 / 4.0d0 / PI / RO_val / CP_val / CP_val
    E1 = CTE * ((ACOSH(P0) - ACOSH(P1)) + CPS * (ACOSH(S0) - ACOSH(S1)))
    E2 = CTE * PS0 / RAIS0
    if (P1 > 1.0d0) E2 = E2 - CTE * PS1 / RAIS1
    
    V1 = ISTEP / 3.0d0 * (RAIS(P0)*RAIS0 + RAIS(S0)**2) * PS0 / RAIS0
    if (P1 > 1.0d0) V1 = V1 - 2.0d0*(ISTEP - 1)/3.0d0 * (RAIS(P1)*RAIS1 + RAIS(S1)**2) * PS1 / RAIS1
    if (P2 > 1.0d0) V1 = V1 + (ISTEP - 2)/3.0d0 * (RAIS(P2)*RAIS2 + RAIS(S2)**2) * PS2 / RAIS2
    
    V2 = ISTEP * RAIS(P0) / 2.0d0 - (ISTEP - 1) * RAIS(P1) + (ISTEP - 2) * RAIS(P2) / 2.0d0
    V3 = ISTEP * RAIS(S0) / 2.0d0 - (ISTEP - 1) * RAIS(S1) + (ISTEP - 2) * RAIS(S2) / 2.0d0
    
    CTZ = 4.0d0 * GE_val * CTE
    Z1 = CTZ * (V1 - CPS * V3)
    Z2 = CTZ * (-2.0d0 * V1 - (V2 - CPS * V3))
    Z3 = CTZ * (V1 - 2.0d0 * XNU_val / (1.0d0 - 2.0d0*XNU_val) * V2)

  contains
    real(8) function ACOSH(X)
      real(8), intent(in) :: X
      ACOSH = log(X + sqrt(X**2 - 1.0d0))
    end function ACOSH

    real(8) function RAIS(X)
      real(8), intent(in) :: X
      RAIS = sqrt(X*X - 1.0d0) / X
    end function RAIS
  end subroutine FUNDSOLE

  real(8) function DISTP(XP, YP, X1, Y1, X2, Y2, X3, Y3)
    real(8), intent(in) :: XP, YP, X1, Y1, X2, Y2, X3, Y3
    real(8) :: D1, D2, D3
    D1 = (XP - X1)**2 + (YP - Y1)**2
    D2 = (XP - X2)**2 + (YP - Y2)**2
    D3 = (XP - X3)**2 + (YP - Y3)**2
    DISTP = 0.98d0 * sqrt(min(D1, D2, D3))
  end function DISTP

  subroutine PREST12(ISTEP, IAU)
    integer, intent(in) :: ISTEP, IAU
    integer :: J, JB, JA, KODEB, KODEA, NF, I, IMAT, IUT, NN, NG_val
    real(8) :: TIME_VAL, BCVALB, BCVALA, DISP, TRAC

    NN = 2 * N_nodes
    NG_val = 6 * NE

    FI = 0.0d0
    TIME_VAL = AT * ISTEP

    do J = 1, NN
       JB = J - 2 + 2 * ((J + 1) / 4)
       if (JB < 1) JB = JB + NG_val
       JA = J + 2 * ((J - 1) / 4)
       if (JA > NG_val) JA = JA - 2
       
       KODEB = KODE(JB)
       KODEA = KODE(JA)
       
       NF = IFUNC(JB)
       BCVALB = FUNCVAL(NF, TIME_VAL, BC(JB) * FSCALE(NF))
       BCVALA = BCVALB
       
       if (JB /= JA) then
          NF = IFUNC(JA)
          BCVALA = FUNCVAL(NF, TIME_VAL, BC(JA) * FSCALE(NF))
       end if
       
       if (KODEB == 1 .and. KODEA == 1 .and. JB == JA) then
          T(JB, ISTEP) = BCVALB
          BCVALA = 0.0d0
       else if (KODEB == 1 .and. KODEA == 1 .and. JB /= JA) then
          T(JB, ISTEP) = BCVALB
          T(JA, ISTEP) = BCVALA
       else if (KODEB == 0 .and. KODEA == 0) then
          if (abs(BCVALB - BCVALA) > 1.0d-8) stop ' ERROR 1 IN BOUNDARY CONDITIONS'
          U(J, ISTEP) = BCVALB
          BCVALA = 0.0d0
       else if (KODEB == 0 .and. KODEA == 1) then
          U(J, ISTEP) = BCVALB
          T(JA, ISTEP) = BCVALA
       else if (KODEB == 1 .and. KODEA == 0) then
          T(JB, ISTEP) = BCVALB
          U(J, ISTEP) = BCVALA
       else
          stop ' ERROR 2 IN BOUNDARY CONDITIONS'
       end if
       
       if (abs(BCVALB) > 0.0d0 .or. abs(BCVALA) > 0.0d0) then
          do I = 1, NN
             FI(I) = FI(I) + B1(I,JB) * BCVALB + B1(I,JA) * BCVALA
          end do
       end if
    end do

    if (ISTEP == 1) return

    rewind IAU
    read(IAU) ! Skip first step matrices which are not used directly here or handled
    ! Wait, F77:
    ! REWIND IAU
    ! READ(IAU)
    ! DO 140 IMAT=2,ISTEP
    ! READ(IAU) ((H(I,J), I=1,NN), J=1,NN),((G(I,J), I=1,NN), J=1,NG)
    ! IUT=ISTEP-IMAT+1 ...
    ! Yes, skip step 1.
    do IMAT = 2, ISTEP
       read(IAU) H(1:NN, 1:NN), G(1:NN, 1:NG_val)
       IUT = ISTEP - IMAT + 1
       
       do J = 1, NN
          DISP = U(J, IUT)
          if (DISP /= 0.0d0) then
             do I = 1, NN
                FI(I) = FI(I) - H(I,J) * DISP
             end do
          end if
       end do
       
       do J = 1, NG_val
          TRAC = T(J, IUT)
          if (TRAC /= 0.0d0) then
             do I = 1, NN
                FI(I) = FI(I) + G(I,J) * TRAC
             end do
          end if
       end do
    end do
  end subroutine PREST12

  real(8) function FUNCVAL(NF, TT, SCALE)
    integer, intent(in) :: NF
    real(8), intent(in) :: TT, SCALE
    integer :: I, NPOIN

    FUNCVAL = FUNC(NF, 1, 2) * SCALE
    if (SCALE == 0.0d0) return
    
    NPOIN = NPOIF(NF)
    do I = 1, NPOIN - 1
       if (TT < FUNC(NF, I, 1) .or. TT >= FUNC(NF, I+1, 1)) cycle
       FUNCVAL = SCALE * (FUNC(NF, I, 2) + (FUNC(NF, I+1, 2) - FUNC(NF, I, 2)) / &
                 (FUNC(NF, I+1, 1) - FUNC(NF, I, 1)) * (TT - FUNC(NF, I, 1)))
       return
    end do
  end function FUNCVAL

end module assembly_mod
