program main
  use data_mod
  use io_mod
  use assembly_mod
  use solver_mod
  implicit none

  integer :: INP, IPR, IAU, ISTEP, KKK

  INP = 10
  IPR = 11
  IAU = 12

  open(INP, file='str12.dat', status='old')
  open(IPR, file='str12.out', status='replace')
  open(IAU, file='temp.dat', status='replace', access='sequential', form='unformatted')

  call INPUT12(INP, IPR)

  do ISTEP = 1, NSTEP
     print *, 'Step: ', ISTEP

     call GHMAT12(ISTEP, IAU)

     call PREST12(ISTEP, IAU)

     KKK = 0
     if (ISTEP == 1) KKK = 1
     call RSOLVER(A1, FI, 2 * N_nodes, KKK)

     call OUTPUT12(ISTEP, AT * ISTEP, IPR)
  end do

  close(INP)
  close(IPR)
  close(IAU)
end program main
