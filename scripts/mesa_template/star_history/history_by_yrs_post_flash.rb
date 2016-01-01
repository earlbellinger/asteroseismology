# history_by_age.rb

$star_data_path = '../LOGS'
load 'load.rb'
StarHistory.new(
   'xaxis' => 'by_yrs_post_flash',
        
    # discard log entries with model_number < first_model
    'first_model' => 0, 

    # discard log entries with model_number > last_model
    # if last_model <= 0, then ignore this parameter.
    'last_model' => 0, 

    # shift all ages by this amount (in years)
    # first_age and last_age refer to shifted ages.
    'age_shift' => 0, 

    # discard log entries for models with ages < first_age
    # first_age = -1 means use first log to set first_age
    'first_age' => -0.5,

    # discard log entries for models with ages > last_age
    # last_age = -1 means use final log to set last_age
    'last_age' => 10, 

    # revise last_age by subtracting discard_span
    'discard_span' => 0,

    # discard log entries for model with ages < last_age - age_span
    # age_span = -1 means set age_span = last_age
    'age_span' => -1, 

    # prune set of models if there are more than max_models
    # max_models = -1 mean there is no max.
    'max_models' => 5000, 

    # if need to prune models because exceed max_models,
    # start pruning by discarding models with
    # model number < previous model number + log_cnt.
    # log_cnt = -1 means set log_cnt to ceil(num_models/max_models)
    'log_cnt' => -1, 

    # if still need to discard models after have pruned using log_cnt,
    # then discard the models with the smallest model_numbers until
    # only have max_models remaining.
    
    
    # parameters for positioning legend
    'legend_left_margin' => nil,
    'legend_top_margin' => nil,
    'legend_right' => nil,


   'track_start_profile' => 2, 
   'track_start_param' => 0.9,
   'log_dir' => $star_data_path
   )