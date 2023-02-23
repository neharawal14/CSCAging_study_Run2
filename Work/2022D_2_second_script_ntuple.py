import os
import subprocess
import ROOT
# specify the directory path
dataset="D"
input_directory = "/cmsuf/data/store/user/nrawal/rootfiles_2022/SingleMuon_2022/Muon/crab_Muon_Run2022{}-ZMu-PromptReco-v2/230124_105044/".format(dataset)
output_path= "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/ntuples_output_2022/2022D_2/second/"
# specify the subdirectory to ignore
#ignore_subdir = 'log'
subdir_list = []
#string_output = "list_2022_{}.txt".format(dataset)
#output_file = open(string_output,"w")
# iterate over all the directories and files
log_file=open("log_2022D_2_second.txt","w") 
for subdir, dirs, files in os.walk(input_directory):
		for file in files:
			file_name = os.path.basename(file)
			if(file_name.endswith(".root")):
					input_file_path = os.path.join(subdir,file_name)
					#print("file to open ", input_file_path,"\n")
					final_path = output_path+file_name
					try:
						f = ROOT.TFile.Open(input_file_path)
					except IOError:
						log_file.write("Failed to open file : "+str(input_file_path)+"\n")
					else:
						print(" open file : "+str(input_file_path)+"\n")
						process_string = "root -l -b -q ../Src/HistMan_cxx.so ../Src/AnalysisGasGain_cxx.so \'analysisgasgain.C(0,0,\"{}\",\"{}\")\'".format(input_file_path,final_path)
						log_file.write("process string "+process_string+"\n")
						subprocess.call(process_string,shell=True)
