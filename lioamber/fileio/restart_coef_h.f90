interface read_coef_restart
   module procedure read_coef_restart_cd
   module procedure read_coef_restart_cs
   module procedure read_coef_restart_od
   module procedure read_coef_restart_os
end interface read_coef_restart

interface write_coef_restart
   module procedure write_coef_restart_cd
   module procedure write_coef_restart_cs
   module procedure write_coef_restart_od
   module procedure write_coef_restart_os
end interface write_coef_restart
