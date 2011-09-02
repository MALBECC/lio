#!/usr/bin/ruby1.9.1
require 'gslng'
include GSLng

file1,file2=ARGV[0],ARGV[1]

file1_lines=File.readlines(file1)
file2_lines=File.readlines(file2)

i1=file1_lines.find_index("\n")
i2=file1_lines.find_index("\n")

if (i1 != i2) then raise 'No son archivos comparables' end

file1_lines=file1_lines.first(i1)
file2_lines=file2_lines.first(i2)

file1_lines.map! {|l| l.split.map{|e| e.to_f}}
file2_lines.map! {|l| l.split.map{|e| e.to_f}}

m1 = Matrix.from_array(file1_lines)
m2 = Matrix.from_array(file2_lines)

m = (m1 - m2)
puts (m / m1) * 100
m = m ^ m

v = Vector[0,0,0]
m.each_vec_row {|row| v += row}
v.map! {|e| Math.sqrt(e)}

puts "RMSD: #{v}"
