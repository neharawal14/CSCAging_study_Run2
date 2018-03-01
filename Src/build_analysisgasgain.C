void  build_analysisgasgain () {
gSystem->Load("HistMan_cxx.so");
gROOT->ProcessLine(".L AnalysisGasGain.cxx++");
}
