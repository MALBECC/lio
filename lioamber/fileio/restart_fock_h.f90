interface read_fock_restart
   module procedure read_fock_restart_cd
   module procedure read_fock_restart_cs
   module procedure read_fock_restart_od
   module procedure read_fock_restart_os
end interface read_fock_restart

interface write_fock_restart
   module procedure write_fock_restart_cd
   module procedure write_fock_restart_cs
   module procedure write_fock_restart_od
   module procedure write_fock_restart_os
end interface write_fock_restart
