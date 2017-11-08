interface read_rho_restart
   module procedure read_rho_restart_c
   module procedure read_rho_restart_o
end interface read_rho_restart

interface write_rho_restart
   module procedure write_rho_restart_c
   module procedure write_rho_restart_o
end interface write_rho_restart
