#!/usr/bin/ruby
while s = gets
  if (s =~ /iteration: (?:([0-9]+)s\. )?([0-9]+)us\./) then
    if ($1.nil?) then puts "0.%.6u" % $2
    else puts "%u.%.6u" % [$1,$2] end
  end  
end
