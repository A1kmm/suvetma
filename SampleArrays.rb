#!/usr/bin/env ruby

require 'optparse'

numsamples = nil
includeFiles = []
excludeFiles = []

p = OptionParser.new do |opts|
  opts.banner = 'Usage: SampleArrays.rb [options]'
  opts.separator ""
  opts.separator "Specific options:"
  opts.on('-nN', '--num-samples N', OptionParser::DecimalInteger,
          'Specifies the number of arrays to sample') do |n|
    numsamples = n
  end

  opts.on('-iFILE', '--include FILE',
          'Specifies a file containing array names to include in the sample population') do |fn|
    includeFiles << fn
  end

  opts.on('-xFILE', '--exclude FILE',
          'Specifies a file containing array names to exclude from the sample population') do |fn|
    excludeFiles << fn
  end

  opts.separator ""
  opts.separator "General options:"

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

  opts.parse(ARGV)

  if numsamples.nil?
    puts "Need to specify number of samples using -n"
    puts opts
    exit
  end
end
p.parse(ARGV)

population = []

includeFiles.each do |name|
  File.open(name) do |f|
    population.concat f.readlines
  end
  population.uniq!
end

excludeFiles.each do |name|
  File.open(name) do |f|
    population = population - f.readlines
  end
end

# Do sampling without replacement...
(1..numsamples).each do
  index = rand(population.size)
  puts population[index]
  population.delete_at(index)
end
