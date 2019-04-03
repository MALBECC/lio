#!/usr/bin/env python2.7
import re
import os

def obtain_dipole(file_in):
   dip = []
   condition = re.match("\w+td",file_in.name)

   if condition:
      patron = "\s+([0-9.-]+E[-+0-9]\d+)\s+([0-9.-]+E[-+0-9]\d+)\s+([0-9.-]+E[-+0-9]\d+)\s+([0-9.-]+E[-+0-9]\d+)"
   else:
      patron = "\s+([0-9.-E-\d+]+)\s+([0-9.-E-\d+]+)\s+([0-9.-E-\d+]+)\s+([0-9.-E-\d+]+)"

   for line in file_in.readlines():
      m = re.match(patron,line)
      if m:
         dip.append(float(m.group(4)))

   return dip

def error(dip,dip_ok):
   dim1 = len(dip)
   dim2 = len(dip_ok)
   scr = 0
   if dim1 != dim2:
      print "There are differents dipole moment in both outputs"
      return -1

   for num in range(dim1):
      value = abs(dip[num] - dip_ok[num])
      if value > 1e-5:
         src = -1
         print "Error in dipole"
         print "Value of dipole moment",dip[num]
         print "Value of ideal dipole moment",dip_ok[num]

   return scr

def Check(*opt):
   option = "".join(opt)

   # Ouput
   if option == "td":
      file_in = "dipole_moment_td"
   else:
      file_in = "dipole_moment"

   dip = []
   is_file = os.path.isfile(file_in)
   if is_file == False:
      print "The %s file is missing" % file_in
      return -1

   f = open(file_in,"r")
   dip = obtain_dipole(f)
   f.close()
   if not dip:
      print "Error reding in %s" % file_in
      return -1

   # Ideal Ouput
   file_ok = file_in + ".ok"
   dipok = []
   is_file = os.path.isfile(file_ok)
   if is_file == False:
      print "The %s file is missing" % file_ok
      return -1

   f = open(file_ok,"r")
   dipok = obtain_dipole(f)
   f.close()
   if not dipok:
      print "Error reading in %s" % file_ok

   ok_output = error(dip,dipok)
   if ok_output != 0:
      print "Test Dipole:   ERROR"
   else:
      print "Test Dipole:   OK"
