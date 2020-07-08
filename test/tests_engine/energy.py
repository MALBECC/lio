#!/usr/bin/env python3
import re
import os
import subprocess

def obtain_energies(file_in):
   energies = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
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

      # Exact exchange for hybrid functionals.
      m = re.match("\s+Exact. Exc.\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[5] = (float(m.group(1)))

      # QM-MM nuclear repulsion.
      m = re.match("\s+QM-MM nuc.\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[6] = (float(m.group(1)))

      # QM-MM electron attraction.
      m = re.match("\s+QM-MM elec.\s+ =\s+([0-9.-]+)",line)
      if m:
         energies[7] = (float(m.group(1)))

      # DFT-D3 energy.
      m = re.match("\s+DFTD3 Energy =\s+([0-9.-]+)",line)
      if m:
         energies[8] = (float(m.group(1)))

   if len(energies) < 1:
      return -1

   return energies

def error(ene,ene_ok):
   tipo = ["Total energy","One electron","Coulomb","Nuclear","Exch. Corr.","Exact Exch.","QM-MM nuclear","QM-MM electronic","DFT-D3"]
   dim1 = len(ene)
   dim2 = len(ene_ok)
   scr = 0
   if dim1 != dim2:
      print("There are different energies in both outputs.")
      return -1

   for num in range(dim1):
      value = abs(ene[num] - ene_ok[num])
      thre = 1.5e-2
      if tipo[num] == "Total energy":
         thre = 1.5e-4
      if value > thre:
         scr = -1
         print("Error in",tipo[num])
         print("Value in output", ene[num])
         print("Value in output.ok", ene_ok[num])

   return scr


def Check():
   # Output
   is_file = os.path.isfile("output")
   if is_file == False:
      print("The output file is missing.")
      return -1

   f = open("output","r")
   energies = []
   energies = obtain_energies(f)
   f.close()
   if not energies:
      print("Error in reading energies in output.")
      return -1

   # Ideal Output
   is_file = os.path.isfile("output.ok")
   if is_file == False:
      print("The output.ok file is missing.")
      return -1

   f = open("output.ok","r")
   energiesok = []
   energiesok = obtain_energies(f)
   f.close()
   if not energiesok:
      print("Error in reading energies in output.ok.")
      return -1

   ok_output = error(energies,energiesok)
   
   if ok_output != 0:
      print("Test Energy:     ERROR")
   else:
      print("Test Energy:     OK")
