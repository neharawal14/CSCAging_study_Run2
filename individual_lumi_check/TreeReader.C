#include "TreeReader.h"
TreeReader :: TreeReader(TString input_file_name)
{
TFile *f = TFile::Open(input_file_name);
tree = (TTree*) f->Get("tree_new");

}
void TreeReader :: initialise(){
   tree->SetBranchAddress("_eventNb", &_eventNb, &b_eventNb);
   tree->SetBranchAddress("_runNb", &_runNb, &b_runNb);
   tree->SetBranchAddress("_lumiBlock", &_lumiBlock, &b_lumiBlock);
   tree->SetBranchAddress("_rhid", &_rhid, &b_rhid);
   tree->SetBranchAddress("_stationring", &_stationring, &b_stationring);
   tree->SetBranchAddress("_rhsumQ", &_rhsumQ, &b_rhsumQ);
   tree->SetBranchAddress("_rhsumQ_RAW", &_rhsumQ_RAW, &b_rhsumQ_RAW);
   tree->SetBranchAddress("_HV", &_HV, &b_HV);
   tree->SetBranchAddress("_current", &_current, &b_current);
   tree->SetBranchAddress("_pressure", &_pressure, &b_pressure);
   tree->SetBranchAddress("_temperature", &_temperature, &b_temperature);
   tree->SetBranchAddress("_instlumi", &_instlumi, &b_instlumi);
   tree->SetBranchAddress("_integratelumi", &_integratelumi, &b_integratelumi);
   tree->SetBranchAddress("_timesecond", &_timesecond, &b_timesecond);
//   tree->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
   tree->SetBranchAddress("_bunchcrossing", &_bunchcrossing, &b_bunchcrossing);
}
