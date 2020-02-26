#!/usr/bin/env python2.7
import re
import os

def obtain_mulliken(file_in):
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
      print "There are diffenrent mulliken charges in outputs."
      return -1

   for num in range(dim1):
      value = abs(mull[num] - mull_ok[num])
      if value > 1e-2:
         src = -1
         print "Error in mulliken charges:"
         print "Valor en mulliken",mull[num]
         print "Valor en mulliken.ok",mull_ok[num]

   return scr


def Check():
   # Output
   mull = []
   is_file = os.path.isfile("mulliken")
   if is_file == False:
      print "The mulliken file is missing."
      return -1

   f = open("mulliken","r")
   mull = obtain_mulliken(f)
   f.close()
   if not mull:
      print "Error in reading mulliken."
      return -1

   # Ideal Output
   is_file = os.path.isfile("mulliken.ok")
   if is_file == False:
      print "The mulliken.ok file is missing."
      return -1

   f = open("mulliken.ok","r")
   mullok = []
   mullok = obtain_mulliken(f)
   f.close()
   if not mullok:
      print "Error in reading mulliken.ok."
      return -1

   ok_output = error(mull,mullok)
   if ok_output != 0:
      print "Test Mulliken:   ERROR"
   else:
      print "Test Mulliken:   OK"
