import os
import subprocess
import ROOT
# specify the directory path
input_directory = "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/ntuples_output_2022"
output_directory = "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/final_files_2022_new"
#output_directory = "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/ntuples_output_new/2022E/final/"
dataset = "C"
output_file_name = "debug_list_2022{}_add_final_ntuple.txt".format(dataset)
output_file = open(output_file_name,"w")
chamber =  ["ME11a", "ME11b", "ME12HV1", "ME12HV2","ME12HV3", "ME13HV1","ME13HV2","ME13HV3","ME21HV1","ME21HV2", "ME21HV3", "ME22HV1", "ME22HV2", "ME22HV3", "ME22HV4","ME22HV5","ME31HV1","ME31HV2", "ME31HV3", "ME32HV1", "ME32HV2", "ME32HV3", "ME32HV4","ME32HV5","ME41HV1","ME41HV2", "ME41HV3", "ME42HV1", "ME42HV2", "ME42HV3", "ME42HV4","ME42HV5"]
# specify the subdirectory to ignore
#ignore_subdir = 'log'
#subdir_list = []
##string_output = "list_2022_{}.txt".format(dataset)
##output_file = open(string_output,"w")
## iterate over all the directories and files
#count=0
#hadd_line = ""
#for subdir, dirs, files in os.walk(input_directory):
#		for file in files:
#			file_name = os.path.basename(file)
#			if(file_name.endswith(chamber[0]+"_tree.root")):
#					#output_file.write(os.path.join(subdir,name_file))
#					input_file_path = os.path.join(subdir,file_name)
#					try:
#						f=ROOT.TFile.Open(input_file_path)
#					except IOError :
#							output_file.write("file not found  "+input_file_path+"\n")
#					else:
#							count=count+1
#							output_file.write("file name "+file_name+"\n")
#							output_file.write(input_file_path+"\n")
#							hadd_line += "{} \n".format(file_name)
#print("count ", count)
#hadd_line = "hadd csc_output_dataset_{}_ME12HV1_tree.root ".format(dataset)+hadd_line
#print("hadd line : ",hadd_line)
##					print("dir", subdir)
##					print("names ",name_file)
##					print(os.path.join(subdir,name_file))
for chamber_name in chamber: 
#	hadd_line = "hadd {}/2022{}/csc_output_dataset_{}_{}_tree.root {}/first/*{}_tree.root {}/second/*{}_tree.root".format(output_directory,dataset,dataset,chamber_name,input_directory,chamber_name,input_directory, chamber_name)
	hadd_line = f"hadd {output_directory}/2022{dataset}/csc_output_dataset_{dataset}_{chamber_name}_tree.root {input_directory}/2022C_1/*{chamber_name}_tree.root {input_directory}/2022C_2/*{chamber_name}_tree.root"
	print(hadd_line)
	subprocess.call(hadd_line,shell=True)

