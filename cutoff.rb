#!/usr/bin/env ruby

require 'tempfile.rb'

class Array
  def to_i!
    map! {|x| x.to_i}
  end

  def sum
    inject {|sum,x| sum+=x }
  end
end

def mapfile(file,cutoff,&proc) 
  sum=file.readline.chomp.split("\t")
  sum.shift
  sum.to_i!
  $stderr.puts sum.inspect
  min =  sum.min
  $stderr.puts "min = #{min}"

  file.each do |line|
    lst=line.chomp.split("\t")
    tag=lst.shift
    lst.to_i!
    norm=[]
    lst.each_with_index do |x,i| 
      norm[i]=x*1.0/min*sum[i]
    end

    line_sum = norm.sum

    if line_sum >= cutoff 
      yield line_sum,line
    end
    
    #puts norm.inspect
    #puts line_sum
    #puts "======="
  end
end


def usage
  $stderr.puts "Use: cutoff.rb matrix_file cutoff > new_matrix_file"
end

if ARGV.length < 1
  usage()
  exit()
end

#cutoff=ARGV[0].to_i
filename=ARGV[0]
cutoff=ARGV[1].to_i


#File.open(filename,"r") do |file|
#  header=file.readline
#  count=0
#  mapfile(file,cutoff) do |sum,line|
#    count+=1;
#  end
#  $stderr.puts "count = #{count}"  
#end

$stderr.puts "processing matrix"

File.open(filename,"r") do |file|
  header=file.readline
  body=Tempfile.new("side_body")
  new_sum=nil
  mapfile(file,cutoff) do |sum,line|
    lst=line.chomp.split("\t")
    lst.shift

    if new_sum.nil?
      new_sum=Array.new(lst.length).fill(0)
    end

    lst.each_with_index do |x,i| 
      new_sum[i]+=x.to_i
    end
    body.print line
  end
  $stderr.print "new sum "
  $stderr.puts new_sum.inspect

  print header
  print "sum\t",new_sum.join("\t"),"\n"
  body.rewind
  body.each {|line| print line}
  #system("cat #{body.path}");  
end

