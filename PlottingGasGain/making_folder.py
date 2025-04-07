import os
path="/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/Code_and_checks/verifying_gas_gain_plots/new_reading_all_dependences/cumulative_plots/AllResults/"
type_channels = ["all_normalised_self","plus_normalised_self", "minus_normalised_self", 
        "odd_chambers_normalised_self", "even_chambers_normalised_self",
        "odd_layers_normalised_self", "even_layers_normalised_self",
        "all_normalised_ME13HV3","plus_normalised_ME13HV3", 
        "odd_layers_normalised_self", "even_layers_normalised_self",
        "upper_plus_normalised_self","upper_minus_normalised_self",
        "lower_plus_normalised_self", "lower_minus_normalised_self", 
      ]
#type_channels = ["plus_minus_endcap"]
#type_channels = ["plus_minus_endcap", "odd_even_chambers" ,  "odd_even_layers"]
#type_channels = [ "lower_plus_minus_endcap","upper_plus_minus_endcap"]
for types in type_channels :
    print(f" type {types}")
    os.chdir(f"{path}")
    os.system(f"mkdir results_gas_gain_{types}")
    os.chdir(f"./results_gas_gain_{types}")
    os.system("mkdir intlumi")
    os.system("mkdir pressure")
    os.system("mkdir instlumi")
    
    os.chdir("pressure")
    os.system("mkdir 2016")
    os.system("mkdir 2017")
    os.system("mkdir 2018")
    os.chdir("2016")
    os.system("mkdir plots_all")
    os.system("mkdir fits")
    os.chdir("../2017")
    os.system("mkdir plots_all")
    os.system("mkdir fits")
    os.chdir("../2018")
    os.system("mkdir plots_all")
    os.system("mkdir fits_all")
    
    os.chdir("../../instlumi")
    os.system("mkdir 2016")
    os.system("mkdir 2017")
    os.system("mkdir 2018")
    os.chdir("2016")
    os.system("mkdir plots_all")
    os.system("mkdir fits")
    os.chdir("../2017")
    os.system("mkdir plots_all")
    os.system("mkdir fits")
    os.chdir("../2018")
    os.system("mkdir plots_all")
    os.system("mkdir fits")
    
    os.chdir("../../intlumi")
    os.system("mkdir 2016")
    os.system("mkdir 2017")
    os.system("mkdir 2018")
    os.system("mkdir run2")
    os.system("mkdir all")
    os.chdir("2016")
    os.system("mkdir fits")
    os.system("mkdir plots_all")
    os.chdir("../2017")
    os.system("mkdir fits")
    os.system("mkdir plots_all")
    os.chdir("../2018")
    os.system("mkdir fits")
    os.system("mkdir plots_all")
    os.chdir("../run2")
    os.system("mkdir fits_all")
    os.system("mkdir fits")
    os.system("mkdir plots_all")
    os.chdir("../all")
    os.system("mkdir fits_all")
    os.system("mkdir fits")
    os.system("mkdir plots_all")

