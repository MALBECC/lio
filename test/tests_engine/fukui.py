#!/usr/bin/env python2.7
import re
import os

def obtain_fukui(file_in):
   lista = []
   for line in file_in.readlines():
      m = re.match("\s+\d+\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)",line)
      if m:
         lista.append(float(m.group(1)))
         lista.append(float(m.group(2)))
         lista.append(float(m.group(3)))
         lista.append(float(m.group(4)))

   if len(lista) < 1:
      return -1

   return lista

def error(fuk,fukok):
   dim1 = len(fuk)
   dim2 = len(fukok)
   scr = 0

   if dim1 != dim2:
      print "There are different number of Fukui charges in outputs."
      return -1

   for num in range(dim1):
      value = abs(fuk[num] - fukok[num])
      if value > 1e-2:
          scr = 1
          print "Error in fukui:"
          print "Value of fukui",fuk[num]
          print "Value of fukui.ok",fukok[num]

   return scr

def Check():
   # Output
   fuk = []
   is_file = os.path.isfile("fukui")
   if is_file == False:
      print "The fukui file is missing."
      return -1

   f = open("fukui","r")
   fuk = obtain_fukui(f)
   f.close
   if not fuk:
      print "Error in reading fukui."

   # Ideal Output
   fukok = []
   is_file = os.path.isfile("fukui")
   if is_file == False:
      print "The fukui.ok file is missing."
      return -1

   f = open("fukui.ok","r")
   fukok = obtain_fukui(f)
   f.close
   if not fukok:
      print "Error in reading fukui.ok."

   ok_output = error(fuk,fukok)
  
   if ok_output != 0:
      print "Test Fukui:      ERROR"
   else:
      print "Test Fukui:      OK"

