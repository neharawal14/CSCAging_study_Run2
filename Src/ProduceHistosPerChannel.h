//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 15 22:47:51 2017 by ROOT version 6.07/03
// from TTree tree/tree
// found on file: TOTME13HV3.root
//////////////////////////////////////////////////////////

#ifndef ProduceHistosPerChannel_h
#define ProduceHistosPerChannel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3.h>
#include <TF1.h>


// Header file for the classes stored in the TTree if any.

class ProduceHistosPerChannel {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       _eventNb;
   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   Int_t           _rhid;
   Int_t           _stationring;
   Double_t        _rhsumQ;
   Double_t        _rhsumQ_RAW;

   Double_t        _pressure;
   Double_t        _temperature;
   Double_t        _instlumi;
   Double_t        _integratelumi;
   UInt_t          _timesecond;
   Int_t           _n_PV;
   Int_t           _bunchcrossing;

   // List of branches
   TBranch        *b__eventNb;   //!
   TBranch        *b__runNb;   //!
   TBranch        *b__lumiBlock;   //!
   TBranch        *b__rhid;   //!
   TBranch        *b__stationring;   //!
   TBranch        *b__rhsumQ;   //!
   TBranch        *b__rhsumQ_RAW;   //!
   TBranch        *b__pressure;   //!
   TBranch        *b__temperature;   //!
   TBranch        *b__instlumi;   //!
   TBranch        *b__integratelumi;   //!
   TBranch        *b__timesecond;   //!
   TBranch        *b__n_PV;   //!
   TBranch        *b__bunchcrossing;   //!


	 TDirectory *dir_var_name;
   ProduceHistosPerChannel(TTree *tree=0);
   virtual ~ProduceHistosPerChannel();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);



   virtual vector <std::pair<double,double > >  GetSlope( TH3D * myh , TString thevar , TString filename, TString title, TFile * outf);

   double ApplyCorrection( double pressure ,TString correctiontype, double p0, double p1 );
   virtual void  Loop(TString detregionstr ="", TString chamber_string = "" );
   
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

	  TString chamber_string_name; 
    TString output_path = "/cmsuf/data/store/user/t2/users/neha.rawal/CSC_2017_data_RAW_RECO/nTuple_output/2017/output_ntuples/";

};

#endif

#ifdef ProduceHistosPerChannel_cxx
ProduceHistosPerChannel::ProduceHistosPerChannel(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
 
}

ProduceHistosPerChannel::~ProduceHistosPerChannel()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ProduceHistosPerChannel::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ProduceHistosPerChannel::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ProduceHistosPerChannel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
   fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
   fChain->SetBranchAddress("_rhid", &_rhid, &b__rhid);
   fChain->SetBranchAddress("_stationring", &_stationring, &b__stationring);
   fChain->SetBranchAddress("_rhsumQ", &_rhsumQ, &b__rhsumQ);
   fChain->SetBranchAddress("_rhsumQ_RAW", &_rhsumQ_RAW, &b__rhsumQ_RAW);
   fChain->SetBranchAddress("_pressure", &_pressure, &b__pressure);
   fChain->SetBranchAddress("_temperature", &_temperature, &b__temperature);
   fChain->SetBranchAddress("_instlumi", &_instlumi, &b__instlumi);
   fChain->SetBranchAddress("_integratelumi", &_integratelumi, &b__integratelumi);
   fChain->SetBranchAddress("_timesecond", &_timesecond, &b__timesecond);
//   fChain->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
   fChain->SetBranchAddress("_bunchcrossing", &_bunchcrossing, &b__bunchcrossing);

   fChain->SetBranchStatus("_temperature",0);
   fChain->SetBranchStatus("_lumiBlock",0);
   fChain->SetBranchStatus("_eventNb",0);
//   fChain->SetBranchStatus("_n_PV",0);
   fChain->SetBranchStatus("_bunchcrossing",0);
   fChain->SetBranchStatus("_timesecond",0);

   Notify();
   
}

Bool_t ProduceHistosPerChannel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ProduceHistosPerChannel::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ProduceHistosPerChannel::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ProduceHistosPerChannel_cxx
