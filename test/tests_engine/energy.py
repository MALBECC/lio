#!/usr/bin/env python2.7
import re
import os
import subprocess

def obtain_energies(file_in):
   energies = [0.0,0.0,0.0,0.0,0.0,0.0]
   for line in file_in.readlines():
      # Total Energy
      m = re.match("\s+Total energy =\s+([0-9.-]+)",line)
      if m:
         energies[0] = (float(m.group(1)))

      # One Electron Energy
      m = re.match("\s+One electron =\s+([0-9.-]+)",line)
      if m:
         energies[1] = (float(m.group(1)))

      # Coulomb Energy
      m = re.match("\s+Coulomb\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[2] = (float(m.group(1)))

      # Nuclear Energy
      m = re.match("\s+Nuclear\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[3] = (float(m.group(1)))

      # Correlation - Exchange Energy
      m = re.match("\s+Exch. Corr.\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[4] = (float(m.group(1)))

      m = re.match("\s+QM-MM\s+energy =\s+([0-9.-]+)",line)
      if m:
         energies[5] = (float(m.group(1)))

   if len(energies) < 1:
      return -1

   return energies

def error(ene,ene_ok):
   tipo = ["Total energy","One electron","Coulomb","Nuclear","Exch. Corr.","QM-MM energy"]
   dim1 = len(ene)
   dim2 = len(ene_ok)
   scr = 0
   if dim1 != dim2:
      print "There are different energies in both outputs"
      return -1

   for num in range(dim1):
      value = abs(ene[num] - ene_ok[num])
      if value > 1e-4:
         scr = -1
         print "Error in",tipo[num]
         print "Value in salida",ene[num]
         print "Value in salida.ok",ene_ok[num]

   return scr


def Check():
   # Output
   is_file = os.path.isfile("salida")
   if is_file == False:
      print "The salida file is missing"
      return -1

   f = open("salida","r")
   energies = []
   energies = obtain_energies(f)
   f.close()
   if not energies:
      print "Error in reading energies in salida"
      return -1

   # Ideal Output
   is_file = os.path.isfile("salida.ok")
   if is_file == False:
      print "The salida.ok file is missing"
      return -1

   f = open("salida.ok","r")
   energiesok = []
   energiesok = obtain_energies(f)
   f.close()
   if not energiesok:
      print "Error in reading energies in salida.ok"
      return -1

   ok_output = error(energies,energiesok)
   
   if ok_output != 0:
      print "Test Energy:   ERROR"
   else:
      print "Test Energy:   OK"
