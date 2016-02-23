#####################################################################
#!/bin/bash

clear
make --always-make find_free_unit_ut.x
./find_free_unit_ut.x 1
./find_free_unit_ut.x 2
./find_free_unit_ut.x 3
./find_free_unit_ut.x 4
./find_free_unit_ut.x 0

#####################################################################
