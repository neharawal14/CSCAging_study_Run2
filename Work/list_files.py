# this program is to list all the files and write in  a directory
import os
dir_path = "/cmsuf/data/store/user/nrawal/SingleMuon/crab_SingleMuon_Run2017B-ZMu-17Nov2017-v1/221027_115247/0000/"
f = open("run_list_2017D_old.txt","w");
for file_1 in os.listdir(dir_path):
   print(file_1)
   f.write("{}{}\n".format(dir_path,file_1))
f.close()
