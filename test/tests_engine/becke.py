#!/usr/bin/env python2.7
import re
import os

def obtain_becke(file_in):
   lista = []
   for line in file_in.readlines():
      m = re.match("\s+\d+\s+\d+\s+([0-9.-]+)",line)
      if m:
         lista.append(float(m.group(1)))

      m = re.match("\s+Total Charge =\s+([0-9.-]+)",line)
      if m:
         lista.append(float(m.group(1)))

   return lista

def error(mull,mull_ok):
   dim1 = len(mull)
   dim2 = len(mull_ok)
   scr = 0
   if dim1 != dim2:
      print "There are diffenrent becke charges in outputs."
      return -1

   for num in range(dim1):
      value = abs(mull[num] - mull_ok[num])
      if value > 1e-2:
         src = -1
         print "Error in becke charges:"
         print "Valor en becke",mull[num]
         print "Valor en becke.ok",mull_ok[num]

   return scr


def Check():
   # Output
   mull = []
   is_file = os.path.isfile("becke")
   if is_file == False:
      print "The becke file is missing."
      return -1

   f = open("becke","r")
   mull = obtain_becke(f)
   f.close()
   if not mull:
      print "Error in reading becke."
      return -1

   # Ideal Output
   is_file = os.path.isfile("becke.ok")
   if is_file == False:
      print "The becke.ok file is missing."
      return -1

   f = open("becke.ok","r")
   mullok = []
   mullok = obtain_becke(f)
   f.close()
   if not mullok:
      print "Error in reading becke.ok."
      return -1

   ok_output = error(mull,mullok)
   if ok_output != 0:
      print "Test Becke:      ERROR"
   else:
      print "Test Becke:      OK"
