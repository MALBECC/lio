#!/usr/bin/ruby
require 'geometry'

input=ARGV.first
geometry = Geometry.input_extract(input)
original_size = geometry.size/3
(original_size - 1).times do |i|
  i = original_size - i - 1
  puts "#{i}."
  oxygens = geometry.select {|g| g[0] == 8}.map{|elem| elem.last(3)}
  center = oxygens.inject([0,0,0]) {|c,p| [c[0]+p[0],c[1]+p[1],c[2]+p[2]]}.map {|v| v / oxygens.size.to_f}
  puts "center: #{center}"

  farthest = oxygens.max_by {|o| dif = [center[0] - o[0],center[1] - o[1],center[2] - o[2]]; Math.sqrt(dif[0]**2 + dif[1]**2 + dif[2]**2)}
  idx = geometry.find_index {|g| g.last(3) == farthest}
  puts "farthest: #{farthest} atom: #{idx}" 

  #File.open("agu#{i+1}-old.xyz",'w') {|io| io.puts geometry.size; io.puts; io.puts geometry.map {|g| g.join(' ')}}
  3.times {geometry.delete_at(idx)}
  File.open("agu#{i}.xyz",'w') {|io| io.puts geometry.size; io.puts; io.puts geometry.map {|g| g.join(' ')}}
end
