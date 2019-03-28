#!/usr/bin/env python2.7
import re
import os

def obtain_forces(file_in):
   lista = []
   for line in file_in.readlines():
      m = re.match("\s+\d+\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)",line)
      if m:
         lista.append(float(m.group(1)))
         lista.append(float(m.group(2)))
         lista.append(float(m.group(3)))

   return lista

def error(fc,fc_ok):
   dim1 = len(fc)
   dim2 = len(fc_ok)
   scr = 0
   if dim1 != dim2:
      print "La cantidad de forces en ambos outputs es diferente"
      return -1

   for num in range(dim1):
      value = abs(fc[num] - fc_ok[num])
      if value > 1e-3:
         scr = -1
         print "Error detected in Forces"
         print "Value in forces",fc[num]
         print "Value in forces.ok",fc_ok[num]

   return scr

def Check():
   # Output
   fc = []
   is_file = os.path.isfile("forces")
   if is_file == False:
      print "The forces file is missing"
      return -1

   f = open("forces","r")
   fc = obtain_forces(f)
   f.close()
   if not fc:
      print "Error reading in forces"
      return -1

   is_file = os.path.isfile("forces.ok")
   if is_file == False:
      print "The forces.ok file is missing"
      return -1

   fcok = []
   f = open("forces.ok","r")
   fcok = obtain_forces(f)
   f.close()
   if not fcok:
      print "Error reading in forces.ok"
      return -1

   ok_output = error(fc,fcok)
   if ok_output != 0:
      print "Test Forces:   ERROR"
   else:
      print "Test Forces:   OK"
