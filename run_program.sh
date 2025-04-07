#!/bin/bash
#g++ -I $ROOTSYS/include reading_tree_removing_VanderMeer.C `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o executable_VanderMeer
#g++ -I $ROOTSYS/include gas_gain_analysis_code.C `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o executable_gas_gain
#g++ -I $ROOTSYS/include gas_gain_even_odd.C `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o executable_gas_gain_test
g++ -I $ROOTSYS/include AnalysisGasGain.C `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o executable_Gasgain_analysis_time
