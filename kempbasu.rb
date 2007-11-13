#!/usr/bin/env ruby


require 'tempfile'

def nprocs
  if File.exists? "/proc/cpuinfo"
    output=`cat /proc/cpuinfo | grep "^processor" | wc -l`
    nprocs=output.to_i
  else
    nprocs=2
  end
  return nprocs
end

prgname=File.basename(__FILE__)
dirname=File.dirname(__FILE__)
if ARGV.length < 1
  puts "Missing input matrix filename"
  exit -1
end

inputname=ARGV.shift # pass the other arguments to basu
basename=File.basename(File.basename(inputname,".dat"),".txt");

resultname=basename + "-#{prgname}.txt"
logname=basename + ".log"

puts "Running #{prgname}"
puts "  directory #{dirname}"
puts "  reading #{inputname}"

File.open(inputname,"r") do |file|
  header=file.readline
  samples=header.split("\t").length-1
  line_count=0
  file.each {|line| line_count+=1}

  puts "#{samples} #{line_count}"

  file.rewind

  Tempfile.open("side") do |out|
    out.puts "#{line_count} #{samples} "
    header=file.readline
    file.each do |line| 
      lst=line.split("\t")
      lst.shift
      out.print lst.join("\t")
      #print ">",lst.join("\t")
    end
    out.close

    prog = File.join(dirname,prgname + ".bin")
    tmp_result=out.path + "-#{prgname}"
    
    puts "temp input file is #{out.path}"
    puts "temp result file is #{tmp_result}"

    puts "=" * 20
    $stdout.flush
    cmd = "#{prog} -n #{nprocs} #{ARGV.join(' ')} #{out.path}"
    print "cmd = #{cmd}\n"

    result=system(cmd);

    puts "Program returned #{result}"
    puts "=" * 20
    $stdout.flush

    file.rewind
    2.times { file.readline }
    File.open(resultname,"w+") do |result|
      File.open(tmp_result,"r") do |infile|
        output_header=infile.readline
        result.print header.chomp,"\t",output_header
        infile.each do |output|
          input=file.readline.chomp
          result.print input,"\t",output
        end
      end
    end
    File.unlink(out.path)
    File.unlink(tmp_result)
  end  
end
