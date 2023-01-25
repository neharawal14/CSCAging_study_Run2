# CSCAgeingP5
CSC ageing studies on P5 data

1. Setup 

The following works from the IHEPA machines and requires a CMSSW release to be installed.
It was tested under CMSSW_8_0_27.
 
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_12_4_6
cd CMSSW_12_4_6/src
cmsenv 

git clone git@github.com:neharawal14/CSCAgeing-studies.git

2. Input files 
We will use our produced ntuples as the input file and process them further to get final ntuples

3. Code
The code itself is in the class AnalysisGasGain. 
Before to run the script, compile and link the HistMan and AnalysisGasGain code (if you made changes in it or run it first time) by corresponding macros build_histman.C and build_analysisgasgain.C (being in folder Src):
root -b -q  build_histman.C // to compile and link histogram manager HistMan
root -b -q build_analysisgasgain.C // to compile and link AnalysisGasGain (this takes a few minutes)

You can then test the code by running on a single file: 
root -l -b -q ../Src/HistMan_cxx.so ../Src/AnalysisGasGain_cxx.so analysisgasgain.C\(0,0,\"path/input_rootfile.root\",\"myoutput.root\"\) 

I have used for now Work/single_file_2022.sh script to execute the above command, and this is the way I am processing my ntuples. 
You can modify the procedure accordign to your way (condor or better say slurm scripts) 
 
One file is produced for each station/ring/HV segment separately. They all contain a tree where each event corresponds to a rechit. 
(N.B. one can probably reduce the file size by a significant amount by changing the format and saving multiple hits in the same event...)
Two important variables here are 
_rhsumQRAW: ADC charge as measured during the data taking
_rhsumQ: ADC charge, undoing the HV changes made in the latest data

(By running the above code, or the same code in the script "single_file_2022.sh" we get final 32 output files, one for each HV segment)
.
4. Extracting Gas gain dependency. 
This is done with a C++ class ( ProduceHistosPerChannel) 
Assuming you have the file final_ME21HV1.root in the subfolder Work, do: 
root -l -b 
.L ../Src/ProduceHistosPerChannel.C++
ProduceHistosPerChannel d
d.Loop("final_ME21HV1") //No ".root" extension 

The resulting output file will then contain histograms with all fits, as well as the slope distribution for each variable studied in a single histogram. 
