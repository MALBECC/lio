interface read_rho_restart
   module procedure read_rho_restart_ccd
   module procedure read_rho_restart_ccs
   module procedure read_rho_restart_cd
   module procedure read_rho_restart_cs
   module procedure read_rho_restart_ocd
   module procedure read_rho_restart_ocs
   module procedure read_rho_restart_od
   module procedure read_rho_restart_os
end interface read_rho_restart

interface write_rho_restart
   module procedure write_rho_restart_ccd
   module procedure write_rho_restart_ccs
   module procedure write_rho_restart_cd
   module procedure write_rho_restart_cs
   module procedure write_rho_restart_ocd
   module procedure write_rho_restart_ocs
   module procedure write_rho_restart_od
   module procedure write_rho_restart_os
end interface write_rho_restart
