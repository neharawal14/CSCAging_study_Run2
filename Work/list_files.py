# this program is to list all the files and write in  a directory
import os
dir_path = "/cmsuf/data/store/user/nrawal/rootfiles_2022/SingleMuon_2022/Muon/crab_Muon_Run2022F-ZMu-PromptReco-v1/230124_105133/0002/"
#dir_path='/cmsuf/data/store/user/nrawal/rootfiles_2022/SingleMuon_2022/Muon/crab_Muon_Run2022D-ZMu-PromptReco-v2/230124_105044/0000'
f = open("debug_run_list_2017F.txt","w");
for file_1 in os.listdir(dir_path):
   print(file_1)
   f.write("{}{}\n".format(dir_path,file_1))
f.close()
