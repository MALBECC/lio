#!/usr/bin/ruby

module Geometry
  ATOM_TYPES = {
    'H' => 1, 'O' => 8, 'C' => 7
  }

  REGEX_LETTER = /^[ \t]*[A-Z]+([ \t]+[0-9-]+.[0-9eE-]*){3,3}[ \t]*$/
  REGEX_NUMBER = /^[ \t]*[0-9]+([ \t]+[0-9-]+.[0-9eE-]*){3,3}[ \t]*$/

  def self.xyz_extract(xyz)
    # Obtains geometry from .xyz file
    geometry = nil
    File.open(xyz) do |f|
      lines = f.readlines
      puts lines.join("\n")
  
      atoms_line = lines.find {|l| l =~ /^ *([0-9]+) *$/}
      if (!atoms_line.nil?) then atoms = $1.to_i
      else raise 'atoms line not found' end
     
      positions = lines.last(atoms)
      if (positions.any? {|p| p !~ REGEX_LETTER}) then
        puts positions.join
        raise 'invalid lines'
      else
        positions.map! do |p|
          atom = p[/[A-Z]/]
          if (!ATOM_TYPES.include?(atom)) then raise "unknown atom type #{atom}" end
          p[/([A-Z])/,1] = ATOM_TYPES[atom].to_s
          p
        end
        geometry = positions
      end
    end

    return geometry
  end
  
  def self.input_extract(input)
    File.open(input) do |f|
      lines = f.readlines
      lines = lines.grep(REGEX_NUMBER)   
      lines.map! {|l|
        atom,x,y,z = l.split
        [atom.to_i,x.to_f,y.to_f,z.to_f]
      } 
    end
  end
end
