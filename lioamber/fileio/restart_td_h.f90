interface read_td_restart_verlet
   module procedure read_td_restart_verlet_s
   module procedure read_td_restart_verlet_d
end interface read_td_restart_verlet

interface read_td_restart_magnus
   module procedure read_td_restart_magnus_s
   module procedure read_td_restart_magnus_d
end interface read_td_restart_magnus

interface write_td_restart_verlet
   module procedure write_td_restart_verlet_s
   module procedure write_td_restart_verlet_d
end interface write_td_restart_verlet

interface write_td_restart_magnus
   module procedure write_td_restart_magnus_s
   module procedure write_td_restart_magnus_d
end interface write_td_restart_magnus
