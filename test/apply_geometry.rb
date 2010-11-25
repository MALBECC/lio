#!/usr/bin/ruby
require 'fileutils'

if (ARGV.include?('-h') || ARGV.size < 2)
  puts "#{__FILE__} <.xyz file> <input file> ... <input file>"
  exit(1)
end

atom_types = {
  'H' => 1, 'O' => 8, 'C' => 7
}

xyz,*files=ARGV

geometry_regex_letter = /^[ \t]*[A-Z]+([ \t]+[0-9-]+.[0-9eE-]*){3,3}[ \t]*$/
geometry_regex_number = /^[ \t]*[0-9]+([ \t]+[0-9-]+.[0-9eE-]*){3,3}[ \t]*$/

# Obtains geometry from .xyz file
geometry = nil
atoms = nil
File.open(xyz) do |f|
  lines = f.readlines
  
  atoms_line = lines.find {|l| l =~ /^ *([0-9]+) *$/}
  if (!atoms_line.nil?) then atoms = $1.to_i
  else raise 'atoms line not found' end
    
  positions = lines.last(atoms)
  if (positions.any? {|p| p !~ geometry_regex_letter}) then
    puts positions.join
    raise 'invalid lines'
  else
    positions.map! do |p|
      atom = p[/[A-Z]/]
      if (!atom_types.include?(atom)) then raise "unknown atom type #{atom}" end
      p[/([A-Z])/,1] = atom_types[atom].to_s
      p
    end
    geometry = positions
  end
end

#puts geometry.join

files.each do |f|
  puts f
  out = nil
  File.open(f) do |io|
    lines = io.readlines
    input_xyz = lines.grep(geometry_regex_number)
    if (input_xyz.size != atoms) then puts input_xyz.join; raise 'couldnt match lines on input' end
    lines[lines.index(input_xyz.first),input_xyz.size] = geometry
    out = lines.join
  end
  
  if (out.nil?) then raise end
  FileUtils.copy(f, "#{f}.bak")
  File.open(f,'w') {|io| io << out }
end
