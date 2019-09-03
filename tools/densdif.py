#!/usr/bin/python3
#import matplotlib
#import matplotlib.pyplot as plt
#import numpy as npy
#
#
#
################################################################################
def read_dens( file_name ):

    data_file = open( file_name, "r" )

    dens_vals = []
    for data_line in data_file.readlines():
        dens_val1 = float( data_line.split()[0] )
        dens_vals.append( dens_val1 )

    data_file.close()
    return dens_vals
#
#
#
#------------------------------------------------------------------------------#
def print_dens( file_name, dens_vals ):

    data_file = open( file_name, 'w')

    for dens_val1 in dens_vals:
        data_file.write( (" %12.9E  %12.9E \n") % (dens_val1,0.0) )
#        print( '{:.9E}  {:.9E}'.format(dens_val1,0.0), file=data_file )

    data_file.close()
    return 0
#
#
#
################################################################################

dens0 = read_dens("dens_neg.in")
dens1 = read_dens("dens_pos.in")

densf = []
for idx in range( len( dens1 ) ):
    densf.append( dens1[idx] - dens0[idx] )

print_dens("dens_dif.out", densf)

################################################################################
