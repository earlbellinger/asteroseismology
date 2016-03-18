# load.rb for histories
require 'scanf'

load '../mesa_dir.rb'

star_test_path = $MESA_DIR + '/star/test'

load star_test_path + '/star_history/load_star_history.rb'
load '../star_profile/load.rb' # use local star_profile
