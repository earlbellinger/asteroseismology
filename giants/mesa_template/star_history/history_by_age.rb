# history_by_age.rb

$star_data_path = '../LOGS'
#$star_data_path = '../LOGS_Tc_7.2'
#$star_data_path = '../LOGS_mdot_2e-7'

load 'load.rb'
StarHistory.new(
   'xaxis' => 'by_age',
        
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
    'first_age' => -1,

    # discard log entries for models with ages > last_age
    # last_age = -1 means use final log to set last_age
    'last_age' => -1, 

    # revise last_age by subtracting discard_span
    'discard_span' => 0,

    # discard log entries for model with ages < last_age - age_span
    # age_span = -1 means set age_span = last_age
    'age_span' => -1, 

    # if need to prune models because exceed max_models,
    # start pruning by discarding models with
    # model number < previous model number + log_cnt.
    # log_cnt = -1 means set log_cnt to ceil(num_models/max_models)
    'log_cnt' => -1, 

    # if still need to discard models after have pruned using log_cnt,
    # then discard the models with the smallest model_numbers until
    # only have max_models remaining.
    
   
    #'conv_logxm_yaxis' => true,
    #    'min_mass' => -0,
    #    'max_mass' => -14,
    # set range for mass in internal structure plots
    #'min_mass' => 0.99999,
    #'max_mass' => 1.0065,

    
    # parameters for positioning legend

    #'legend_top_margin' => 0.05,# nil,
    #'legend_top_margin' => 0.34,# nil,
    'legend_top_margin' => 0.68,# nil,

    #'legend_left_margin' => 0.6, 'legend_right' => 0.9,
    #'legend_left_margin' => 0.4, 'legend_right' => 0.57,
    'legend_left_margin' => 0.06, 'legend_right' => 0.34,




   'track_start_profile' => 2, 
   'track_start_param' => 0.9,
   'log_dir' => $star_data_path
   )