#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
class TreeReader{
	public :
		TTree *tree;
		TreeReader(TString);
		void initialise();
   ULong64_t       _eventNb;
   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   Int_t           _rhid;
   Int_t           _stationring;
   Double_t        _rhsumQ;
   Double_t        _rhsumQ_RAW;
   Double_t        _HV;
   Double_t        _current;

   Double_t        _pressure;
   Double_t        _temperature;
   Double_t        _instlumi;
   Double_t        _integratelumi;
   UInt_t          _timesecond;
   Int_t           _n_PV;
   Int_t           _bunchcrossing;

	 TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_rhid;   //!
   TBranch        *b_stationring;   //!
   TBranch        *b_rhsumQ;   //!
   TBranch        *b_rhsumQ_RAW;   //!
   TBranch        *b_HV;   //!
   TBranch        *b_current;   //!
   TBranch        *b_pressure;   //!
   TBranch        *b_temperature;   //!
   TBranch        *b_instlumi;   //!
   TBranch        *b_integratelumi;   //!
   TBranch        *b_timesecond;   //!
   TBranch        *b_n_PV;   //!
   TBranch        *b_bunchcrossing;   //!



	
};
