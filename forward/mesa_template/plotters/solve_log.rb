# use this to create plots of the series of corrections for the newton iterations.

# to get the plot data, set inspectB_flag = .true. in the inspectB routine in mod_hydro_newton.
# if you want to stop after a certain iteration, set inspectB_iter_stop.
# if you need to see other variables, you can add what you want in inspectB.

# the plots show both the corrections and the changes in the corrections.
# if things are failing to converge, this can point out the variable
# and the location in the star that seems to be the source of the trouble.

# to narrow in on a particular region of the star, set first_col and last_col below.
# first_col = 0 starts at the surface; last_col = -1 ends at the center.

# to narrow in on a particular sequence of iterations, set first_row and last_row.
# you may need to skip the initial corrections in order to see what's going on later.

require "scanf"

first_row = 0
last_row = -1

first_col = 0
last_col = -1

#first_col = 1090
#last_col = 1130

load "../../../../utils/image_plot.rb"
load "../../../../utils/image_log.rb"

log_data_dir = '../plot_data/solve_logs'
names_filename = 'names.data'
xlabel = 'grid'
ylabel = 'iteration'
image_Xs = nil
image_Ys = nil

f = File.open('../plot_data/solve_logs/size.data')
sizes = f.scanf('%d %d')
numcols = sizes[0]
numrows = sizes[1]
f.close

puts "numcols #{numcols}"
puts "numrows #{numrows}"

Image_Log_plot.new(                
                log_data_dir,
                names_filename,
                xlabel, ylabel, 
                first_col, last_col, 
                first_row, last_row,
                image_Xs, image_Ys, true, '.log',
                numcols, numrows)
