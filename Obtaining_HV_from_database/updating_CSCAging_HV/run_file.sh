g++ -I $ROOTSYS/include HV_update.cpp `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o exe_HV
