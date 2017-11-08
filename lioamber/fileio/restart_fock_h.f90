interface read_fock_restart
   module procedure read_fock_restart_c
   module procedure read_fock_restart_o
end interface read_fock_restart

interface write_fock_restart
   module procedure write_fock_restart_c
   module procedure write_fock_restart_o
end interface write_fock_restart
