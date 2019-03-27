#!/bin/bash

#g09 diazi2.com
#formchk diazi2.chk

echo '&inputgen'                        > inputgen.in
echo '  output_type  = "full_ground",' >> inputgen.in
echo '  output_xyz   = "diazi2.xyz",'  >> inputgen.in
echo '  output_rst   = "diazi2.rst",'  >> inputgen.in
echo '  inpname_fchk = "diazi2.fchk",' >> inputgen.in
echo '&end'                            >> inputgen.in

../inputgen.x inputgen.in > inputgen.log

