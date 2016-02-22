#####################################################################
#!/bin/bash
ut_prog='safely_open_ut.x'

clear
make --always-make ${ut_prog}
./${ut_prog} 1
./${ut_prog} 2
./${ut_prog} 3
./${ut_prog} 4
echo 'Content of Output-test:'
more Output-test
./${ut_prog} 0

#####################################################################
