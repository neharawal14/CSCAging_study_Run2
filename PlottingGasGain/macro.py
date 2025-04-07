import os
current_path="/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/Code_and_checks/verifying_gas_gain_plots/new_reading_all_dependences/cumulative_plots/"
year_period=["2016_full", "2017_full", "2018_full"]    
#year_period=["2016_full"]    
chamber_list = ["ME11a", "ME11b", "ME12HV1", "ME12HV2","ME12HV3", "ME13HV1","ME13HV2","ME13HV3","ME21HV1","ME21HV2", "ME21HV3", "ME22HV1", "ME22HV2", "ME22HV3", "ME22HV4","ME22HV5","ME31HV1","ME31HV2", "ME31HV3", "ME32HV1", "ME32HV2", "ME32HV3", "ME32HV4","ME32HV5","ME41HV1","ME41HV2", "ME41HV3", "ME42HV1", "ME42HV2", "ME42HV3", "ME42HV4","ME42HV5"]
for year_period_num in year_period: 
    os.chdir(f"{current_path}")
    os.system(f"mkdir {year_period_num}")
    os.chdir(f"{year_period_num}")
    os.system("mkdir output_plots")
    os.chdir("output_plots")
    os.system("mkdir gas_gain_intlumi")
    os.chdir("../")
#    os.system("cp ../reading_mean_values.C .")
#    os.system("cp ../reading_mean_graph.C .")
#    os.system("cp ../reading_mean_graph_ratio_same_chamber.C .")
    os.system("mkdir all_channels")
    os.chdir("all_channels")
    for chamber in chamber_list : 
        os.system(f"mkdir {chamber}")
        os.chdir(f"{chamber}")
        os.system(f"mkdir pressure")
        os.system(f"mkdir pressure_second")
        os.system(f"mkdir instlumi")
        os.system(f"mkdir instlumi_second")
        os.system(f"mkdir intlumi_initial")
        os.system(f"mkdir intlumi_final")
        os.system(f"mkdir timesecond_initial")
        os.system(f"mkdir timesecond_final")
        os.chdir("../")
