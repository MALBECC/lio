subroutine write_final_convergence(converged, iterations, energy)
   use fileio_data, only: verbose

   implicit none
   logical         , intent(in) :: converged
   integer         , intent(in) :: iterations
   double precision, intent(in) :: energy

   if (verbose .lt. 1) return;
   if (converged) then
      write(*,7000) iterations, energy
   else
      write(*,7001) iterations, energy
   endif

   return;

7000 FORMAT("  Convergence achieved in ", I6, " iterations. Final energy ", &
            F14.7, " A.U.")
7001 FORMAT("  No convergence achieved in ", I6, " iterations. Final energy ", &
            F14.7, " A.U.")
end subroutine write_final_convergence

subroutine write_energies(E1, E2, En, Ens, Eecp, Exc, ecpmode, E_restrain, &
                          number_restr, nsol)
   use fileio_data, only: style, verbose

   implicit none
   double precision, intent(in) :: E1, E2, En, Ens, Eecp, Exc, E_restrain
   integer         , intent(in) :: number_restr, nsol
   logical         , intent(in) :: ecpmode

   if (verbose .lt. 3) return;
   if (style) then
      write(*,*)
      write(*,7000); write(*,7004); write(*,7001)
      write(*,7005) E1-Eecp
      write(*,7002)
      write(*,7006) E2
      write(*,7002)
      write(*,7007) En
      if (ecpmode) then             ! ECP is active.
         write(*,7002)
         write(*,7008) Eecp
      endif
      if (nsol .gt. 0) then         ! QM/MM calculation is being performed.
         write(*,7002)
         write(*,7013) Ens
      endif
      write(*,7002)
      write(*,7009) Exc
      write(*,7002)
      if (number_restr .gt. 0) then !Restraints are active.
         write(*,7011) E_restrain
         write(*,7002)
      endif
      write(*,7010) E1 + E2 + En + Ens + Exc
      write(*,7003)
      write(*,*)
   else
      write(*,*)
      write(*,'(A)') "Final Energy Contributions in A.U."
      write(*,'(A,F12.6)') "  Total energy = ", E1 + E2 + En + Ens + Exc
      write(*,'(A,F12.6)') "  One electron = ", E1 - Eecp
      write(*,'(A,F12.6)') "  Coulomb      = ", E2
      write(*,'(A,F12.6)') "  Nuclear      = ", En
      write(*,'(A,F12.6)') "  Exch. Corr.  = ", Exc
      if (nsol .gt. 0) write(*,'(A,F12.6)') "  QM-MM energy = ", Ens
      if (ecpmode)     write(*,'(A,F12.6)') "  ECP energy   = ", Eecp
      if (number_restr .gt. 0) &
                       write(*,'(A,F12.6)') "  Restraints   = ", E_restrain
      write(*,*)
   endif

   return

7000 FORMAT(4x,"╔═════════════════════════════════", &
"═══════════╗")
7001 FORMAT(4x,"╠══════════════════╦══════════════", &
"═══════════╣")
7002 FORMAT(4x,"╠══════════════════╬══════════════", &
"═══════════╣")
7003 FORMAT(4x,"╚══════════════════╩══════════════", &
"═══════════╝")
7004 FORMAT(4x,"║     FINAL ENERGY CONTRIBUTIONS IN A.U.     ║")
7005 FORMAT(4x,"║   ONE ELECTRON   ║",4x,F14.7,7x"║")
7006 FORMAT(4x,"║   COULOMB        ║",4x,F14.7,7x"║")
7007 FORMAT(4x,"║   NUCLEAR        ║",4x,F14.7,7x"║")
7008 FORMAT(4x,"║   E. CORE POT    ║",4x,F14.7,7x"║")
7009 FORMAT(4x,"║   EXC. - CORR.   ║",4x,F14.7,7x"║")
7012 FORMAT(4x,"║   E QM-MM        ║",4x,F14.7,7x"║")
7010 FORMAT(4x,"║   TOTAL          ║",4x,F14.7,7x"║")
7011 FORMAT(4x,"║   E. RESTR.      ║",4x,F14.7,7x"║")
7013 FORMAT(4x,"║   NUCLEAR QM-MM  ║",4x,F14.7,7x"║")
end subroutine write_energies

subroutine write_energy_convergence(step, energy, good, told, egood, etold)
   use fileio_data, only: style, verbose

   implicit none
   integer         , intent(in) :: step
   double precision, intent(in) :: energy, good, told, egood, etold

   if (verbose .lt. 2) return;
   if (style) then
      write(*, 8500)
      write(*, 8501) step, energy
      write(*, 8502)
      write(*, 8601)
      write(*, 8602)
      write(*, 8603)
      write(*, 8604) good , told
      write(*, 8605) egood, etold
      write(*, 8606)
   else
      write(*, 8700) step, energy, good, told, egood, etold
   endif

   return;
8500 FORMAT(4x,"╔════════╦═════════════╦══════════", &
 "═╦══════════════════════╗")
8501 FORMAT(4x,"║ iter # ║",2x,I10,1x,"║ QM Energy ║",4x,F14.7,4x,"║")
8502 FORMAT(4x,"╚════════╩═════════════╩══════════", &
 "═╩══════════════════════╝")
8601 FORMAT(4x,"           ╔════════════╦═════════════╗")
8602 FORMAT(4x,"           ║    Value   ║ Conv. Crit. ║")
8603 FORMAT(4x,"╔══════════╬════════════╬═════════", &
"════╣")
8604 FORMAT(4x,"║ Good     ║",1x,E10.3,1x,"║",1x,E10.3,2x,"║")
8605 FORMAT(4x,"║ En. Good ║",1x,E10.3,1x,"║",1x,E10.3,2x,"║")
8606 FORMAT(4x,"╚══════════╩════════════╩═════════", &
"════╝")

8700 FORMAT(2x, "Step = ", I6, 1x, " - QM Energy = ", F12.5, 1x,        &
            "- Rho diff (crit) = ", ES8.2, " (",ES8.2, &
            ") - Energy diff (crit) = ", ES8.2, " (",ES8.2,")")
end subroutine write_energy_convergence
