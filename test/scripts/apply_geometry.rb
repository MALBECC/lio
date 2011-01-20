#!/usr/bin/ruby
require 'fileutils'
require 'geometry'

if (ARGV.include?('-h') || ARGV.size < 2)
  puts "#{__FILE__} <.xyz file> <input file> ... <input file>"
  exit(1)
end

xyz,*files=ARGV

geometry = Geometry.xyz_extract(xyz)

#puts geometry.join

files.each do |f|
  puts f
  out = nil
  File.open(f) do |io|
    lines = io.readlines
    input_xyz = lines.grep(Geometry::REGEX_NUMBER)
    lines[lines.index(input_xyz.first),input_xyz.size] = geometry
    out = lines.join
  end
  
  if (out.nil?) then raise end
  FileUtils.copy(f, "#{f}.bak")
  File.open(f,'w') {|io| io << out }
end
