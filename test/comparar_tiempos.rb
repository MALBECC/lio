#!/usr/bin/ruby
def get_tiempos(archivo)
  tiempos = `fgrep TIMER #{archivo}`.split(/\n/)
  tiempos.map! {|t|
    if (t =~ /\d+s\. \d+us\.$/)
      t[/(\d+)s\. \d+us\.$/,1].to_f * 1000.0 + t[/\d+s\. (\d+)us\.$/,1].to_f / 1000.0
    else
      t[/(\d+)us./,1].to_f / 1000.0
    end
  }
  tiempo_total = `fgrep real #{archivo}`[/\d+\.\d+/]
  return tiempo_total, tiempos
end

ARGV.each do |parm|
  puts parm
  [['cpu','c'],['gpu','g']].each do |tipo,archivo|
    tiempo_total, tiempos = get_tiempos("#{parm}/#{parm}-igrid1.#{archivo}")
    puts "#{tipo}: #{tiempo_total}"
    #  primer_tiempo = tiempos.shift
    #  ultimo_tiempo = tiempos.pop
    #  tiempos = tiempos.inject(0) {|n,m| n+m} / tiempos.length
    #  puts "#{tipo}: #{primer_tiempo} #{tiempos} #{ultimo_tiempo} #{total}"

  end
end

