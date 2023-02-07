import os
import subprocess
import ROOT
# specify the directory path
dataset="D"
input_directory = "/cmsuf/data/store/user/nrawal/rootfiles_2022/SingleMuon_2022/Muon/crab_Muon_Run2022{}-ZMu-PromptReco-v2/".format(dataset)
output_path= "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/output_ntuples_new/2022D_2/"
# specify the subdirectory to ignore
#ignore_subdir = 'log'
subdir_list = []
#string_output = "list_2022_{}.txt".format(dataset)
#output_file = open(string_output,"w")
# iterate over all the directories and files
for subdir, dirs, files in os.walk(input_directory):
		for file in files:
			file_name = os.path.basename(file)
			if(file_name.endswith(".root")):
					#output_file.write(os.path.join(subdir,name_file))
					input_file_path = os.path.join(subdir,file_name)
					final_path = output_path+file_name
#					print("final_file to open ",final_path)
#					print("file to open", input_file_path)

					f = ROOT.TFile.Open(input_file_path)
					if(f.IsZombie()==True) : 
						continue
					process_string = "root -l -b -q ../Src/HistMan_cxx.so ../Src/AnalysisGasGain_cxx.so \'analysisgasgain.C(0,0,\"{}\",\"{}\")\'".format(input_file_path,final_path)
					print("process string ",process_string)
					subprocess.call(process_string,shell=True)
#					print("dir", subdir)
#					print("names ",name_file)
#					print(os.path.join(subdir,name_file))
