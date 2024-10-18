//#define pressure_dependence_removal_instlumi_cxx
#include "pressure_dependence_removal_instlumi.h"
#include "badchannel.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
//#include <stream>
#include <string>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/code_area/files_HVandLumi/nonme11_first.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/code_area/files_HVandLumi/nonme11_second.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/code_area/files_HVandLumi/me11.h"

using namespace std;
enum ring_station_hvsegm {
  me11a,me11b,me12HV1,me12HV2,me12HV3,me13HV1,me13HV2,me13HV3,
  me21HV1,me21HV2,me21HV3,me22HV1,me22HV2,me22HV3,me22HV4,me22HV5,
  me31HV1,me31HV2,me31HV3,me32HV1,me32HV2,me32HV3,me32HV4,me32HV5,
  me41HV1,me41HV2,me41HV3,me42HV1,me42HV2,me42HV3,me42HV4,me42HV5
};


void pressure_dependence_removal_instlumi::defining_bool() { 
 testing_code = false;
// for intlumi corrections, time
 intlumi_initial = true;
 intlumi_final = true;
 time_initial = false;
 time_final = false;

// for anlaysing only 2016
 intlumi_corr_2016 = true;
 intlumi_corr_2017 = true;
 intlumi_corr_2018 = true;

// for pressure corrections
// flags used : pressure_corr, pressure_corr_2016, pressure_corr_2017, pressure_corr_2018, intlumi_instlumi_const
 pressure_corr = true;
 pressure_corr_2016 = true;
 pressure_corr_2017 = true;
 pressure_corr_2018 = true;
 intlumi_instlumi_const = false;
// instlumi corrections
 instlumi_corr = true;
 instlumi_corr_2016 = true;
 instlumi_corr_2017 = true;
 instlumi_corr_2018 = true;
 intlumi_const = false;
 year ="2016";
}

bool instlumicorr =true;
double trimmean =0.85;

bool dropinstlumicorr =false; // Drop inst L correction
bool droppressurecorr =false;// Drop pressure correction
bool seconditer =false;// If you want to perform a second iteration of the fit to pressure and inst L.
bool savehistos =false;
//Offset used in the calculation of the trimmed mean bool istest =false; //For tests, run only on 1/100 of the events
bool debug_statements =false;
bool debug_print = false;
bool iszmumu = false;//If only Zmumu events are used the corrections cannot be derived channel per channel (not enough stat) so we take all channels together

int thestationringhv (TString thename){
  if(thename.Index("ME11a")>0)return me11a ;
  if(thename.Index("ME11b")>0)return me11b ;
  if(thename.Index("ME12HV1")>0)return me12HV1 ;
  if(thename.Index("ME12HV2")>0)return me12HV2 ;
  if(thename.Index("ME12HV3")>0)return me12HV3 ;
  if(thename.Index("ME13HV1")>0)return me13HV1 ;
  if(thename.Index("ME13HV2")>0)return me13HV2 ;
  if(thename.Index("ME13HV3")>0)return me13HV3 ;
  if(thename.Index("ME21HV1")>0)return me21HV1 ;
  if(thename.Index("ME21HV2")>0)return me21HV2 ;
  if(thename.Index("ME21HV3")>0)return me21HV3 ;
  if(thename.Index("ME22HV1")>0)return me22HV1 ;
  if(thename.Index("ME22HV2")>0)return me22HV2 ;
  if(thename.Index("ME22HV3")>0)return me22HV3 ;
  if(thename.Index("ME22HV4")>0)return me22HV4 ;
  if(thename.Index("ME22HV5")>0)return me22HV5 ;
  if(thename.Index("ME31HV1")>0)return me31HV1 ;
  if(thename.Index("ME31HV2")>0)return me31HV2 ;
  if(thename.Index("ME31HV3")>0)return me31HV3 ;
  if(thename.Index("ME32HV1")>0)return me32HV1 ;
  if(thename.Index("ME32HV2")>0)return me32HV2 ;
  if(thename.Index("ME32HV3")>0)return me32HV3 ;
  if(thename.Index("ME32HV4")>0)return me32HV4 ;
  if(thename.Index("ME32HV5")>0)return me32HV5 ;
  if(thename.Index("ME41HV1")>0)return me41HV1 ;
  if(thename.Index("ME41HV2")>0)return me41HV2 ;
  if(thename.Index("ME41HV3")>0)return me41HV3 ;
  if(thename.Index("ME42HV1")>0)return me42HV1 ;
  if(thename.Index("ME42HV2")>0)return me42HV2 ;
  if(thename.Index("ME42HV3")>0)return me42HV3 ;
  if(thename.Index("ME42HV4")>0)return me42HV4 ;
  if(thename.Index("ME42HV5")>0)return me42HV5 ;

  return -1;
}


TString thestationringhv(int i){
  if(i ==me11a)return"ME11a";
  if(i ==me11b)return"ME11b";
  if(i ==me12HV1)return"ME12HV1";
  if(i ==me12HV2)return"ME12HV2";
  if(i ==me12HV3)return"ME12HV3";
  if(i ==me13HV1)return"ME13HV1";
  if(i ==me13HV2)return"ME13HV2";
  if(i ==me13HV3)return"ME13HV3";
  if(i ==me21HV1)return"ME21HV1";
  if(i ==me21HV2)return"ME21HV2";
  if(i ==me21HV3)return"ME21HV3";
  if(i ==me22HV1)return"ME22HV1";
  if(i ==me22HV2)return"ME22HV2";
  if(i ==me22HV3)return"ME22HV3";
  if(i ==me22HV4)return"ME22HV4";
  if(i ==me22HV5)return"ME22HV5";
  if(i ==me31HV1)return"ME31HV1";
  if(i ==me31HV2)return"ME31HV2";
  if(i ==me31HV3)return"ME31HV3";
  if(i ==me32HV1)return"ME32HV1";
  if(i ==me32HV2)return"ME32HV2";
  if(i ==me32HV3)return"ME32HV3";
  if(i ==me32HV4)return"ME32HV4";
  if(i ==me32HV5)return"ME32HV5";
  if(i ==me41HV1)return"ME41HV1";
  if(i ==me41HV2)return"ME41HV2";
  if(i ==me41HV3)return"ME41HV3";
  if(i ==me42HV1)return"ME42HV1";
  if(i ==me42HV2)return"ME42HV2";
  if(i ==me42HV3)return"ME42HV3";
  if(i ==me42HV4)return"ME42HV4";
  if(i ==me42HV5)return"ME42HV5";
  return "";
}


void pressure_dependence_removal_instlumi::Loop(TString input_file_path, TString input_file_name, TString chamber_string, TString output_file_path, TString output_folder_name)
{
  chamber_string_name = chamber_string;
	detregionstr="_dataset";
	output_path = output_file_path;
	output_plots_folder = output_folder_name;
  gStyle->SetOptStat();
  gStyle->SetOptFit(111);
  if(detregionstr.Index("ZMuMu")>=0) iszmumu=true;
  iszmumu = false;
  bool debug = false;

  double nentries;

	time_t initial_time, final_time;
	 
	if(year=="2016"){
	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
	  final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
	}
	if(year=="2017"){
	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
	 final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
	}
	if(year=="2018"){
	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
	 final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
	}


	 gStyle->SetTimeOffset(0.);
   // we are passing name of the input file when calling the macro 
   TString inputfname = input_file_path+ input_file_name;
   TFile * inputf = TFile::Open(inputfname);
   
   if(seconditer)detregionstr =detregionstr+"seconditer"; 
   detregionstr =detregionstr+"_pressure_corrected_"; 

   TTree * tree =(TTree*) inputf->Get("tree");
	 // Init initializes all the branches in this tree
   Init(tree);

   if(testing_code) {
     std::cout<<"running on 1000000 events"<<std::endl;
     nentries = 1000000;
   }
   else {
     nentries = tree->GetEntries();
   }
	 // output root file after processing 
   TFile * outf = new TFile(output_path+"outf"+detregionstr+"_"+chamber_string_name+"_output_run2.root","recreate");  
	 // Making IntegratedLumi vs gas gain slope dependency before starting any pressure and inst lumi corrections
   //outf->cd();
	 // 5 days is 1 bin
  vector< std::pair<double, double > > params_integratelumi_initial;
  vector< std::pair<double, double > > params_timesecond_initial;

  if(intlumi_initial || time_initial){

   TH3D * hchargevsintegratelumi_initial = new TH3D("hchargevsintegratelumi_initial","charge (ADC counts) vs integ lumi (initial)",3000,0,3000, 50, 0,150  ,770,1,771);
   TH3D * hchargevstime_initial = new TH3D("hchargevstime_initial","charge (ADC counts) vs time (initial)",3000,0,3000, 60, initial_time,final_time ,770,1,771);
     for(int i = 0; i < nentries; i++){
  	     LoadTree(i);tree->GetEntry(i);
        if(!_passZmumusel) continue;
        if(i%1000000 ==0)cout << i<<endl;
 			  int rhidreduced = ((int)floor(_rhid/10))%1000;
        if(_rhid> 2000000) rhidreduced +=400;
        int idforcorr = (iszmumu )? 0: rhidreduced ;

        double nominal_HV_value = NominalHV();

        if(debug) std::cout<<" removing the event based on HV"<< " rhid :"<<rhidreduced<<" nominal HV :"<<nominal_HV_value<<" set HV  :"<<_HV<<std::endl;
        if(abs(_HV-nominal_HV_value) >10) continue;

 			  double charge_equalized = _rhsumQ;
        hchargevsintegratelumi_initial->Fill(charge_equalized, _integratelumi, rhidreduced);

        if(time_initial) {
          hchargevstime_initial->Fill(charge_equalized, _timesecond, rhidreduced);
        } // end of only if we analyse time information 
    } // end of entries of tree
  params_integratelumi_initial = GetSlope( hchargevsintegratelumi_initial, "_integratelumi_initial", detregionstr,"",outf, chamber_string_name);
	std::cout<<" done with intlumi information"<<std::endl; 
  if(time_initial) {
       params_timesecond_initial = GetSlope( hchargevstime_initial, "_timesecond_initial", detregionstr,"",outf, chamber_string_name);
  } // end of time_initial


  } // end of intlumi_initial or time_initial


  //Run on all events to extract pressure correction for each channel (rechit) separately. 
// This are the cuts for analysing a particular dataset of time for pressure and inst lumi
// instlumi corrections 
double intlumi_low_cut1_2016 = 12; 
double intlumi_up_cut1_2016 = 32; 
double intlumi_low_cut1_2017 = 55; 
double intlumi_up_cut1_2017 = 75 ; 
double intlumi_low_cut1_2018 = 115; 
double intlumi_up_cut1_2018 = 135 ; 

// pressure corrections cut
double intlumi_low_cut_2016 = 12; 
double intlumi_up_cut_2016 = 32; 
double instlumi_low_cut_2016 = 7000; 
double instlumi_up_cut_2016 = 9000; 
double intlumi_low_cut_2017 = 55; 
double intlumi_up_cut_2017 = 75 ; 
double instlumi_low_cut_2017 = 7000; 
double instlumi_up_cut_2017 = 9000; 
double intlumi_low_cut_2018 = 115; 
double intlumi_up_cut_2018 = 135 ; 
double instlumi_low_cut_2018 = 10000; 
double instlumi_up_cut_2018 = 15000; 

double integratelumi_2016_high = 39.32673126400002;
double integratelumi_2017_high = 83.85340826572357;

vector< std::pair<double, double > > params_pressure_2016;
vector< std::pair<double, double > > params_pressure_2017;
vector< std::pair<double, double > > params_pressure_2018;


std::cout<<"going to pressure information"<<std::endl;
// flags used : pressure_corr, pressure_corr_2016, pressure_corr_2017, pressure_corr_2018, intlumi_instlumi_const
if(pressure_corr) { 
   TH3D * hchargevspressure_2016 = new TH3D("hchargevspressure_2016","charge (ADC counts) vs pressure : 2016",3000,0,3000, 20, 946,986  ,770,1,771);
   TH3D * hchargevspressure_2017 = new TH3D("hchargevspressure_2017","charge (ADC counts) vs pressure : 2017",3000,0,3000, 20, 946,986  ,770,1,771);
   TH3D * hchargevspressure_2018 = new TH3D("hchargevspressure_2018","charge (ADC counts) vs pressure : 2018",3000,0,3000, 20, 946,986  ,770,1,771);
   for(int i = 0; i < nentries; i++){
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl; 
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     double nominal_HV_value = NominalHV();
     // applying HV cuts
     if(abs(_HV-nominal_HV_value) >10) continue;

    if(pressure_corr_2016==true){
		 if(_integratelumi <= integratelumi_2016_high) { 
      // intlumi_instlumi_const is a flag whether to apply the constraints on the instlumi and intlumi
       if(intlumi_instlumi_const == true){
  			 if(_integratelumi<intlumi_low_cut_2016 ) continue;
   	     if(_integratelumi>intlumi_up_cut_2016) continue;
         if(_instlumi<instlumi_low_cut_2016) continue;
         if(_instlumi>instlumi_up_cut_2016 ) continue;
       }
       hchargevspressure_2016->Fill(_rhsumQ, _pressure , rhidreduced);
			 if(debug) std::cout<<" inside the 2016 instlumi loop"<<std::endl;
      }  // end 2016 filling
     } // end of pressure_corr_2016 : fill only when we analyse 2016

      if(pressure_corr_2017==true){
  		 if( integratelumi_2016_high <_integratelumi && _integratelumi<= integratelumi_2017_high) { 
        if(intlumi_instlumi_const == true){
    		 if(_integratelumi<intlumi_low_cut_2017 ) continue;
         if(_integratelumi>intlumi_up_cut_2017) continue;
         if(_instlumi<instlumi_low_cut_2017) continue;
         if(_instlumi>instlumi_up_cut_2017 ) continue;
       }
       hchargevspressure_2017->Fill(_rhsumQ, _pressure , rhidreduced);
			 if(debug) std::cout<<" inside the 2017 instlumi loop"<<std::endl;
      }// end 2017 filling
     } // end of pressure_corr_2017 : fill only when we analyse 2017

     if(pressure_corr_2018==true){
		    if( integratelumi_2017_high <_integratelumi) { 
          if(intlumi_instlumi_const == true){
  			    if(_integratelumi<intlumi_low_cut_2018 ) continue;
   	        if(_integratelumi>intlumi_up_cut_2018) continue;
            if(_instlumi<instlumi_low_cut_2018) continue;
            if(_instlumi>instlumi_up_cut_2018 ) continue;
          }
          // ignore the end part of Run2 when HV was changed and giving us bad fit
          if(_pressure >=976) continue;
          hchargevspressure_2018->Fill(_rhsumQ, _pressure , rhidreduced);
			    if(debug) std::cout<<" inside the 2018 instlumi loop"<<std::endl;
        } // end of 2018 filling  
      } // end of pressure_corr_2018 : fill only when we analyse 2018
	
    } // end of all tree entries for filling pressure histograms

		if(hchargevspressure_2016 != NULL && hchargevspressure_2016->GetEntries() >=100 && pressure_corr_2016){
     std::cout<<" entries  hcharge "<<hchargevspressure_2016->GetEntries()<<std::endl;
    params_pressure_2016 = GetSlope( hchargevspressure_2016, "_pressure", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
    
		if(hchargevspressure_2017 != NULL && hchargevspressure_2017->GetEntries() >=100 && pressure_corr_2017){
     std::cout<<" entries  hcharge "<<hchargevspressure_2017->GetEntries()<<std::endl;
     params_pressure_2017 = GetSlope(hchargevspressure_2017, "_pressure", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
		if(hchargevspressure_2018 != NULL&& hchargevspressure_2018->GetEntries() >=100 && pressure_corr_2018){
     std::cout<<" entries  hcharge "<<hchargevspressure_2018->GetEntries()<<std::endl;
    params_pressure_2018 = GetSlope( hchargevspressure_2018, "_pressure", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}

} // end of pressure_corr

// Now we will correct for the pressure corrections , then only we derive Intlumi corrections
	
   vector< std::pair<double, double > > params_instlumi_2016;
   vector< std::pair<double, double > > params_instlumi_2017;
   vector< std::pair<double, double > > params_instlumi_2018;
    double params_pressure_const_avg_2016;
    double params_pressure_slope_avg_2016;
    double params_pressure_const_avg_2017;
    double params_pressure_slope_avg_2017;
    double params_pressure_const_avg_2018;
    double params_pressure_slope_avg_2018;

   double charge, charge_equalized; 
if(instlumi_corr) {

  if(pressure_corr_2016){
    params_pressure_const_avg_2016 =  (params_pressure_2016[0].first + params_pressure_2016[771].first)/2. ;
    params_pressure_slope_avg_2016 =  (params_pressure_2016[0].second + params_pressure_2016[771].second)/2. ;
  }
  if(pressure_corr_2017){
    params_pressure_const_avg_2017 =  (params_pressure_2017[0].first + params_pressure_2017[771].first)/2. ;
    params_pressure_slope_avg_2017 =  (params_pressure_2017[0].second + params_pressure_2017[771].second)/2. ;
  }
  if(pressure_corr_2018){
    params_pressure_const_avg_2018 =  (params_pressure_2018[0].first + params_pressure_2018[771].first)/2. ;
    params_pressure_slope_avg_2018 =  (params_pressure_2018[0].second + params_pressure_2018[771].second)/2. ;
  }

   // Instlumi bins 
   if(debug) std::cout<<" going to instlumi information"<<std::endl;
   TH3D * hchargevsinstlumi_2016 = new TH3D("hchargevsinstlumi_2016","charge (ADC counts) vs instlumi : 2016",3000,0,3000, 42, 0, 21000, 770, 1, 771);
   TH3D * hchargevsinstlumi_2017 = new TH3D("hchargevsinstlumi_2017","charge (ADC counts) vs instlumi : 2017",3000,0,3000, 42, 0, 21000, 770, 1, 771);
   TH3D * hchargevsinstlumi_2018 = new TH3D("hchargevsinstlumi_2018","charge (ADC counts) vs instlumi : 2018",3000,0,3000, 42, 0, 21000, 770, 1, 771);
   
    for(int i = 0; i < nentries; i++){
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl; 
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;//First (second) endcap have rechit ID < (>) 2000000
    
     double nominal_HV_value = NominalHV();
     if(abs(_HV-nominal_HV_value) >10) continue;

     if(instlumi_corr_2016==true){
		  if(_integratelumi <= integratelumi_2016_high) { 
       if(intlumi_const == true){
  			 if(_integratelumi<intlumi_low_cut1_2016 ) continue;
   	     if(_integratelumi>intlumi_up_cut1_2016) continue;
       }
  		  double charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 );
        double charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 ) ; 
        hchargevsinstlumi_2016->Fill(charge_equalized, _instlumi , rhidreduced);
			 if(debug) std::cout<<" inside the 2016 instlumi loop"<<std::endl;
      }// end 2016 filling
     } // end of instlumi_corr_2016 : fill only when we analyse 2016
     
     if(instlumi_corr_2017==true){
      if( integratelumi_2016_high <_integratelumi && _integratelumi<= integratelumi_2017_high) { 
       if(intlumi_const == true){
    		 if(_integratelumi<intlumi_low_cut1_2017 ) continue;
         if(_integratelumi>intlumi_up_cut1_2017) continue;
       }
	     double charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 );
        double charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 ) ; 

       hchargevsinstlumi_2017->Fill(charge_equalized, _instlumi , rhidreduced);
			 if(debug) std::cout<<" inside the 2017 instlumi loop"<<std::endl;
      } // end 2017 filling
     } // end of instlumi_corr_2017 : fill only when we analyse 2017

      if(instlumi_corr_2018==true){ 
       if( integratelumi_2017_high <_integratelumi) { 
        if(intlumi_const == true){
  			 if(_integratelumi<intlumi_low_cut_2018 ) continue;
   	     if(_integratelumi>intlumi_up_cut_2018) continue;
        }
       double charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 );
        double charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 ) ; 

        hchargevsinstlumi_2018->Fill(charge_equalized, _instlumi , rhidreduced);
			  if(debug) std::cout<<" inside the 2018 instlumi loop"<<std::endl;
       } // end 2018 filling
      } //end of instlumi_corr_2018 : fill only when we analyse 2018  
	  } // end of all tree entries for filling instlumi histograms
		
    if(hchargevsinstlumi_2016 != NULL && hchargevsinstlumi_2016->GetEntries() >=100 && instlumi_corr_2016){
    params_instlumi_2016 = GetSlope( hchargevsinstlumi_2016, "_instlumi", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		} 
    
    if(hchargevsinstlumi_2017 != NULL && hchargevsinstlumi_2017->GetEntries() >=100 && instlumi_corr_2017){
			std::cout<<" entries  Hcharge "<<hchargevsinstlumi_2017->GetEntries()<<std::endl;
      params_instlumi_2017 = GetSlope(hchargevsinstlumi_2017, "_instlumi", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
    if(hchargevsinstlumi_2018 != NULL && hchargevsinstlumi_2018->GetEntries() >=100 && instlumi_corr_2018){
     params_instlumi_2018 = GetSlope( hchargevsinstlumi_2018, "_instlumi", detregionstr,"",outf, chamber_string_name);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
   
  } // end of loop of instlumi corrections derivations

  //Now final charge after pressure, instlumi correction
     vector< std::pair<double, double > > params_integratelumi;
     vector< std::pair<double, double > > params_timesecond;
    double params_instlumi_const_avg_2016;
    double params_instlumi_slope_avg_2016;
    double params_instlumi_const_avg_2017;
    double params_instlumi_slope_avg_2017;
    double params_instlumi_const_avg_2018;
    double params_instlumi_slope_avg_2018;



 if(intlumi_final || time_final) { 
  TH3D * hchargevsintegratelumi = new TH3D("hchargevsintegratelumi","charge (ADC counts) vs integ lumi",3000,0,3000, 50, 0,150  ,770,1,771);
  TH3D * hchargevstime = new TH3D("hchargevstime","charge (ADC counts) vs time",3000,0,3000, 60, initial_time,final_time  ,770,1,771);
 
  if(instlumicorr){
      if(instlumi_corr_2016){
          params_instlumi_const_avg_2016 =  (params_instlumi_2016[0].first + params_instlumi_2016[771].first)/2. ;
          params_instlumi_slope_avg_2016 =  (params_instlumi_2016[0].second + params_instlumi_2016[771].second)/2. ;
          params_pressure_const_avg_2016 =  (params_pressure_2016[0].first + params_pressure_2016[771].first)/2. ;
          params_pressure_slope_avg_2016 =  (params_pressure_2016[0].second + params_pressure_2016[771].second)/2. ;
       }
      if(instlumi_corr_2017){
          params_instlumi_const_avg_2017 =  (params_instlumi_2017[0].first + params_instlumi_2017[771].first)/2. ;
          params_instlumi_slope_avg_2017 =  (params_instlumi_2017[0].second + params_instlumi_2017[771].second)/2. ;
          params_pressure_const_avg_2017 =  (params_pressure_2017[0].first + params_pressure_2017[771].first)/2. ;
          params_pressure_slope_avg_2017 =  (params_pressure_2017[0].second + params_pressure_2017[771].second)/2. ;

      }
      if(instlumi_corr_2018){
        params_instlumi_const_avg_2018 =  (params_instlumi_2018[0].first + params_instlumi_2018[771].first)/2. ;
        params_instlumi_slope_avg_2018 =  (params_instlumi_2018[0].second + params_instlumi_2018[771].second)/2. ;
        params_pressure_const_avg_2018 =  (params_pressure_2018[0].first + params_pressure_2018[771].first)/2. ;
        params_pressure_slope_avg_2018 =  (params_pressure_2018[0].second + params_pressure_2018[771].second)/2. ;
      }
   }
  else{
      if(pressure_corr_2016){
          params_pressure_const_avg_2016 =  (params_pressure_2016[0].first + params_pressure_2016[771].first)/2. ;
          params_pressure_slope_avg_2016 =  (params_pressure_2016[0].second + params_pressure_2016[771].second)/2. ;
        }
      if(pressure_corr_2017){
        params_pressure_const_avg_2017 =  (params_pressure_2017[0].first + params_pressure_2017[771].first)/2. ;
        params_pressure_slope_avg_2017 =  (params_pressure_2017[0].second + params_pressure_2017[771].second)/2. ;
      }
      if(pressure_corr_2018){
        params_pressure_const_avg_2018 =  (params_pressure_2018[0].first + params_pressure_2018[771].first)/2. ;
        params_pressure_slope_avg_2018 =  (params_pressure_2018[0].second + params_pressure_2018[771].second)/2. ;
      }
  }

    // When we apply correction we should take average dependence for pressure : avg of plus and minus endcap
    for(int i=0 ; i<tree->GetEntries(); i++){
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;
     double charge, charge_equalized; 
     double nominal_HV_value = NominalHV();
     if(debug) std::cout<<" removing the event based on HV"<< " rhid :"<<rhidreduced<<" nominal HV :"<<nominal_HV_value<<" set HV  :"<<_HV<<std::endl;
     if(abs(_HV-nominal_HV_value) >10) continue;
		 //applying pressure correction by taking average of slope dependence in plus and minus endcap 
      if(intlumi_corr_2016 ==true) {
        if(_integratelumi <= integratelumi_2016_high){

          // If instlumi correction which means pressure corrections are already applied, if not then we need to just apply pressure corrections
          if(instlumicorr){
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 ) ; 

      		  charge  = charge * ApplyCorrection( _instlumi ,"instlumi",   params_instlumi_const_avg_2016 ,params_instlumi_slope_avg_2016 );
            charge_equalized  = charge_equalized * ApplyCorrection( _instlumi ,"instlumi", params_instlumi_const_avg_2016 ,params_instlumi_slope_avg_2016 ) ; 

          }
          else{
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2016 ,params_pressure_slope_avg_2016 ) ; 
          }
		    }// to apply 2016 corrections  only on events of 2016
      } // end of 2016 corrections

      if(intlumi_corr_2017==true) {
         if( integratelumi_2016_high < _integratelumi && _integratelumi<= integratelumi_2017_high) { 
          if(instlumicorr){
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 ) ; 

      		  charge  = charge * ApplyCorrection( _instlumi ,"instlumi",   params_instlumi_const_avg_2017 ,params_instlumi_slope_avg_2017 );
            charge_equalized  = charge_equalized * ApplyCorrection( _instlumi ,"instlumi", params_instlumi_const_avg_2017 ,params_instlumi_slope_avg_2017 ) ; 

          }
          else{
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2017 ,params_pressure_slope_avg_2017 ) ; 
          }
         } // end of 2017 corrections for intlumi
      } // end of 2017 corrections

     if(intlumi_corr_2018==true) {
        if(_integratelumi > integratelumi_2017_high){
          if(instlumicorr){
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 ) ; 

      		  charge  = charge * ApplyCorrection( _instlumi ,"instlumi",   params_instlumi_const_avg_2018 ,params_instlumi_slope_avg_2018 );
            charge_equalized  = charge_equalized * ApplyCorrection( _instlumi ,"instlumi", params_instlumi_const_avg_2018 ,params_instlumi_slope_avg_2018 ) ; 

          }
          else{
      		  charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 );
            charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg_2018 ,params_pressure_slope_avg_2018 ) ; 
          }

        }    // end of 2018 corrections for intlumi
      } // end of 2018 corrections 
    hchargevsintegratelumi->Fill(charge_equalized, _integratelumi, rhidreduced);
   } // end of going through each tree entry

    if(intlumi_final) {
     params_integratelumi = GetSlope( hchargevsintegratelumi, "_integratelumi", detregionstr,"",outf, chamber_string_name); 
	  std::cout<<"working for time second information after correction "<<std::endl; 
    } // end of intlumi_final

    if(time_final) {
      params_timesecond = GetSlope( hchargevstime, "_timesecond", detregionstr,"",outf, chamber_string_name);
    } // end of time_final

 } // end of both intlumi and time final corrections
	 
std::cout<<" Now we will fill the final root file with all this information and close it "<<std::endl;
outf->Close();
} // end of Loop function

vector < std::pair<double, double > >  pressure_dependence_removal_instlumi::GetSlope( TH3D * myh , TString thevar , TString filename, TString title, TFile * outf, TString chamber_name_string){

  TString name_histogram = myh->GetName();
  if(debug_print) std::cout<<"inside the slope function for the var "<<thevar<<std::endl;
  //Will store all the fit results in a TTree (one entry per channel)
  TString treename_goodchannels = thevar; 
  std::cout<<" name of histgoram "<<name_histogram<<std::endl;
	if(name_histogram=="hchargevspressure_2016") treename_goodchannels = treename_goodchannels+"_2016";
	else if(name_histogram=="hchargevspressure_2017") treename_goodchannels = treename_goodchannels+"_2017";
	else if(name_histogram=="hchargevspressure_2018") treename_goodchannels = treename_goodchannels+"_2018";
  else if(name_histogram=="hchargevsinstlumi_2016") treename_goodchannels = treename_goodchannels+"_2016";
  else if(name_histogram=="hchargevsinstlumi_2017") treename_goodchannels = treename_goodchannels+"_2017";
  else if(name_histogram=="hchargevsinstlumi_2018") treename_goodchannels = treename_goodchannels+"_2018";
  else treename_goodchannels = "tree_all_goodchannels"+thevar ; 

	std::cout<<" name of the tree good channels ************************ pressure **********"<<treename_goodchannels<<std::endl;
//  TTree * theouttree_goodchannels = new TTree(treename_goodchannels,"");
  TString treename ; 	
	if(name_histogram=="hchargevspressure_2016") treename = "tree_"+thevar+"_2016";
	else if(name_histogram=="hchargevspressure_2017") treename ="tree_"+thevar+ "_2017";
	else if(name_histogram=="hchargevspressure_2018") treename = "tree_"+thevar+"_2018";
	else if(name_histogram=="hchargevsinstlumi_2016") treename = "tree_"+thevar+"_2016";
	else if(name_histogram=="hchargevsinstlumi_2017") treename = "tree_"+thevar+"_2017";
	else if(name_histogram=="hchargevsinstlumi_2018") treename = "tree_"+thevar+"_2018";
  else treename= "tree_"+thevar ; 

//  TTree * theouttree = new TTree(treename,"");
//  Float_t _slope(-1),_slope_error(-1),_chi2(-1);
//	Float_t _slope_goodchannels_plus(-1), _slope_error_goodchannels_plus(-1);
//	Float_t _slope_goodchannels_minus(-1), _slope_error_goodchannels_minus(-1);
//  Int_t _ndof(-1),_layer(-1),_chamber(-1),_endcap(-1),_stationringHVseg(-1); 
//  Bool_t _isbadchannel(false);
//  theouttree_goodchannels->Branch("_slope_goodchannels_plus",&_slope_goodchannels_plus,"_slope_goodchannels_plus/F");
//  theouttree_goodchannels->Branch("_slope_goodchannels_minus",&_slope_goodchannels_minus,"_slope_goodchannels_minus/F");
//  theouttree_goodchannels->Branch("_slope_error_goodchannels_plus",&_slope_error_goodchannels_plus,"_slope_error_goodchannels_plus/F");
//  theouttree_goodchannels->Branch("_slope_error_goodchannels_minus",&_slope_error_goodchannels_minus,"_slope_error_goodchannels_minus/F");
//  theouttree->Branch("_slope",&_slope,"_slope/F");
//  theouttree->Branch("_slope_error",&_slope_error,"_slope_error/F");
//  theouttree->Branch("_chi2",&_chi2 ,"_chi2/F");
//  theouttree->Branch("_ndof",&_ndof,"_ndof/I");
//  theouttree->Branch("_layer",&_layer,"_layer/I");
//  theouttree->Branch("_chamber",&_chamber,"_chamber/I");
//  theouttree->Branch("_endcap",&_endcap,"_endcap/I");
//  theouttree->Branch("_stationringHVseg",&_stationringHVseg,"_stationringHVseg/I");
//  theouttree->SetAutoSave(1000000);
  TString xtitle ;
  if(thevar.Index("pressure")>=0) xtitle = "Pressure (hPa)";
  if(thevar.Index("instlumi")>=0) xtitle = "Inst Lumi (#mub s)^{-1}";
  if(thevar.Index("integrate")>=0) xtitle = "Integrated luminosity (fb^{-1})";
  if(thevar.Index("time")>=0) xtitle = "time";
 
 
	outf->cd();

	TString dir_name_var ;
	if(name_histogram=="hchargevspressure_2016") dir_name_var = thevar+"_2016";
	else if(name_histogram=="hchargevspressure_2017") dir_name_var = thevar+"_2017";
	else if(name_histogram=="hchargevspressure_2018") dir_name_var = thevar+"_2018";
	else if(name_histogram=="hchargevsinstlumi_2016") dir_name_var = thevar+"_2016";
	else if(name_histogram=="hchargevsinstlumi_2017") dir_name_var = thevar+"_2017";
	else if(name_histogram=="hchargevsinstlumi_2018") dir_name_var = thevar+"_2018";
	else dir_name_var = thevar;

	TDirectoryFile *dir_var_name =  (TDirectoryFile*) outf->mkdir(dir_name_var);
	dir_var_name->cd();
	// result function is used to store the value of the slope and constant after fitting 
  vector < std::pair<double, double >  > result ; //Assume that fitted function has two parameters
  
  for(int i = 0; i<772;i++){
    std::pair<double, double > initpair(0,0); 
    result.push_back(initpair);
  }

  if(thevar.Index("pressure")>=0 &&droppressurecorr ) return result; 
  if(thevar.Index("instlumi")>=0 &&dropinstlumicorr ) return result; 

  // this lowedge and highedge are for integrated luminosity  slope 
  double lowedge =  -0.01; 
  double highedge  = 0.01;
  
  if( thevar.Index("instlumi")>=0 )  lowedge = -0.00005;
  if( thevar.Index("instlumi")>=0 )  highedge = 0.00005;
  if( thevar.Index("instlumi_2")>=0 )  lowedge = -0.000005;
  if( thevar.Index("instlumi_2")>=0 )  highedge = 0.000005;
 
	if( thevar.Index("pressure")>=0 )  lowedge = -0.05;
  if( thevar.Index("pressure")>=0 )  highedge = 0.03;
 	if( thevar.Index("time")>=0 )  lowedge = 1462060800;
  if( thevar.Index("time")>=0)  highedge = 1477958400;
 	//if( thevar.Index("time")>=0 )  lowedge = 1656633600;
  //if( thevar.Index("time")>=0)  highedge = 1669852800;
 	// this histogram store the slope values for each channel and it will be a guassian distribution  
//	TH1D * h_slope = new TH1D ("h_slope"+thevar+"_"+filename,"",200, lowedge, highedge );
//  TH1D * h_slopeuncty = new TH1D ("h_slopeuncty"+thevar+"_"+filename,"",200, lowedge/100., highedge/100. );
//  TH1D * h_chi2 = new TH1D ("h_chi2"+thevar+"_"+filename,"",200,0,1000);
	// Loop over reduced rechit ID
  // Reduced rechit ID has the following format: (A+B)*10+C, where A=1,..36 (chamber nb), B =0 (endcap 1) or 40 (endcap 2)  and C =1,..., 6 (layer). 
	// For calculating normalized charge distribution for all layers together
	// h=0 corresponds to the plus endcap, h=771 corresponds to minus endcap
	//ofstream mean_values;
	//mean_values.open("mean_values.txt");
  TString rhidshort;

  // For debugging purposes I am going to save total number of entries in each bin of the histograms
  TH2D *num_entries_2D_hist_12th_bin = new TH2D("num_entries_2D_hist_12th_bin", " # rechit entries : 12th bin  ", 36,1, 37, 6,1,7);
  TH2D *num_entries_2D_hist_15th_bin = new TH2D("num_entries_2D_hist_15th_bin", " # rechit entries : 15th bin ", 36,1, 37, 6,1,7);

 TH2D *num_entries_2D_hist_all_bins_plus = new TH2D("num_entries_2D_hist_all_bins_plus", " # rechit entries : cumulative channels : plus-endcap ", 36,1, 37, 6,1,7);
 TH2D *num_entries_2D_hist_all_bins_minus = new TH2D("num_entries_2D_hist_all_bins_minus", " # rechit entries : cumulative channels : minus-endcap ", 36,1 , 37, 6,1,7);

 TH1D *num_entries_1D_hist_all_bins_plus; 
 TH1D *num_entries_1D_hist_all_bins_minus; 
 int nbins; 
 int up_limit; 
 if(thevar=="_pressure" || thevar=="_instlumi"){
   nbins = 300;
   up_limit = 60000;
 }
 else if(thevar=="_integratelumi_initial" || thevar=="_integratelumi"){
   nbins = 1000;
   up_limit = 100000;
 }
 num_entries_1D_hist_all_bins_plus = new TH1D("num_entries_1D_hist_all_bins_plus", " # rechit entries : cumulative channels : plus-endcap ", nbins, 0, up_limit);
 num_entries_1D_hist_all_bins_minus = new TH1D("num_entries_1D_hist_all_bins_minus", " # rechit entries : cumulative channels : minus-endcap ", nbins, 0, up_limit);

  bool flag_12 = false;
  bool flag_15 = false;

  // finding the list of channels for which the channels have less entries
  // and excluding them from the cumulative plots
   int numChannels = myh->GetNbinsZ();
   std::vector<int> num_entries(numChannels, 0);
   std::vector<int> low_entry_channels;
	 int nb_channels= 0;
   for(int h = 0; h < 772 ; h++){
		// This is to check thing if(h!=0 &&h!=771) continue;
    //Skipping empty entries
    if( (1<= h && h<= 6) || (401<= h && h<=406))continue;
    //if(h!=0 && h%10>=7)continue;
    if(h!=0 && h%10>=7)continue;
	 	// h >7 implies that the one which are extra layers after reduced rechit Id = 366 for chamber36, endcap1, layer6 :
		// and reducedrechitId = 411 for chamber 1 , endcap 2, layer 1  , we dont want to count them 
    if(h!=0 && h%10==0)continue;
    if(h>370&& h<=400)continue;
    //N.B.: h=0 takes all rechits together 
		TString endcap = (h<=400)? "_Endcap1":"_Endcap2";
    int chambernb = (h<=400)?  (int)floor(h/10)  : (int)floor( (h-400) /10) ;
		// this condition will remove if there are extra chambers between 400 and 410 

    if(h!=0 && h!=771) rhidshort ="chamber"+ (TString) Form("%d", chambernb)  +"_layer"+ (TString)Form("%d",h%10) + endcap;
    int layernb = h%10; 
    if(h==0) rhidshort = "allgoodchannels_plus";
    if(h==771) rhidshort = "allgoodchannels_minus";

    for (int j = 1; j <= myh->GetNbinsY(); j++) {
    TH1D* proj1 = (h == 0 || h == 771) ? nullptr : (TH1D*)(myh->ProjectionX("_px", j, j, h, h))->Clone();
    if (proj1) {  // Only add if the projection is not null
        num_entries[h] += proj1->GetEntries();
        //num_entries_all_bins += proj1->GetEntries();
        delete proj1;  // Clean up the memory to avoid memory leaks
      }
    }
    if(num_entries[h]!=0) nb_channels++;
 
   if(h!=0 && h!=771){ 
     //std::cout<<"*******************entries in the channel "<<h<<" : "<<num_entries[h]<<std::endl;
     if(h<=400){
     num_entries_2D_hist_all_bins_plus->SetBinContent(chambernb, layernb, num_entries[h]);
     num_entries_1D_hist_all_bins_plus->Fill(num_entries[h]);
     }
     if(h>400){
     num_entries_2D_hist_all_bins_minus->SetBinContent(chambernb, layernb, num_entries[h]);
     num_entries_1D_hist_all_bins_minus->Fill(num_entries[h]);
     }
    }
   }

		std::cout<<" nb channels "<<nb_channels<<std::endl;
		double sum = 0;
    double sum_sq = 0;
    for (int n : num_entries) {
        sum += n;
        sum_sq += n * n;
    }

    // do this only for pressure and instlumi but not integrate lumi since number of entries are a lot, and the program just looses its mind
    // Defining nrechits = 10000 as threshold for integrated luminosity
  double threshold;
  if(thevar=="_integratelumi_initial" || thevar=="_integratelumi"){
   threshold = 20000;
//   if(chamber_string=="ME11b") threshold = 15000;
//   if(chamber_string=="ME12HV1") threshold = 15000;
//   if(chamber_string=="ME12HV2") threshold = 15000;
//   if(chamber_string=="ME12HV3") threshold = 15000;
  }
  else{
//	double mean = sum /nb_channels;
//  double std_dev = std::sqrt(sum_sq / nb_channels - mean * mean);
//  // Define threshold and identify low entry channels
//   double k = 1.5;
//   std::cout<<" sum values "<<sum<<" sum_sq "<<sum_sq<<std::endl;
//   std::cout<<" mean values "<<mean<<" std dev "<<std_dev<<std::endl;
//   threshold = mean - k * std_dev;
   threshold = 20000;
  }
  for (int h = 0; h < 772; ++h) {
       if(num_entries[h]==0) continue;
        if (num_entries[h] < threshold) {
            low_entry_channels.push_back(h);
        }
    }
    std::cout << "Channels with low entries (threshold = " << threshold << "):" << std::endl;
    for (int h : low_entry_channels) {
       std::cout << "Channel " << h << " with " << num_entries[h] << " entries" << std::endl;
      }

      std::cout<<" number of channels with low entries "<<low_entry_channels.size()<<std::endl;

      
    // Starting processing each channel
    for(int h = 0; h < 772 ; h++){
      // This is to check thing if(h!=0 &&h!=771) continue;
      //Skipping empty entries
      if( (1<= h && h<= 6) || (401<= h && h<=406))continue;
      //if(h!=0 && h%10>=7)continue;
      if(h!=0 && h%10>=7)continue;
      // h >7 implies that the one which are extra layers after reduced rechit Id = 366 for chamber36, endcap1, layer6 :
      // and reducedrechitId = 411 for chamber 1 , endcap 2, layer 1  , we dont want to count them 
      if(h!=0 && h%10==0)continue;
      if(h>370&& h<=400)continue;
      //N.B.: h=0 takes all rechits together 
      TString endcap = (h<=400)? "_Endcap1":"_Endcap2";
      int chambernb = (h<=400)?  (int)floor(h/10)  : (int)floor( (h-400) /10) ;
      // this condition will remove if there are extra chambers between 400 and 410 

      if(h!=0 && h!=771) rhidshort ="chamber"+ (TString) Form("%d", chambernb)  +"_layer"+ (TString)Form("%d",h%10) + endcap;
      int layernb = h%10; 
      if(h==0) rhidshort = "allgoodchannels_plus";
      if(h==771) rhidshort = "allgoodchannels_minus";

      //std::cout<<" testing only  good channels plus and minus, h value"<<h<<" rhid values "<<rhidshort<<" endcap "<<endcap<<std::endl;
      
      //Declare a new histo to store trimmed mean for different values of the variable of interest (pressure, inst L,...)
      TString htrimmeanvsXname = "htrimmean"+filename+title+thevar+"_"+rhidshort;
      TH1D * htrimmeanvsX = new TH1D(htrimmeanvsXname,"", myh->GetNbinsY() , myh->GetYaxis()->GetBinLowEdge(1) , myh->GetYaxis()->GetBinLowEdge( myh->GetNbinsY()+1) );
      double renormalfactor = 1;    
      double error_renormalfactor =0 ;    
      double error_value=0 ; 
      bool flag= false;		
      double gas_gain;
      double gas_gain_error;



      //Loop over the bins of the variable of interest
      //Get the rechit ADC charge distribution for a given bin of the variable of interest

    // Two 2D histograms for filling : 
      
      //TCanvas *c_individual;
      for(int j = 1; j <= myh->GetNbinsY(); j++){
        TH1D* proj = (h==0 || h==771)? (TH1D*) (myh->ProjectionX("_px",j,j,1,1))->Clone() : (TH1D*)(myh->ProjectionX("_px",j,j,h,h))->Clone();  

        // saving all the distributions 
        
  /*			if(rhidshort.Contains("chamber1") || rhidshort.Contains("chamber4") ||rhidshort.Contains("chamber5") || rhidshort.Contains("chamber6") || rhidshort.Contains("chamber7") || rhidshort.Contains("chamber8") || rhidshort.Contains("chamber9") || rhidshort.Contains("chamber20") ||rhidshort.Contains("chamber21") || rhidshort.Contains("chamber24") || rhidshort.Contains("chamber25") || rhidshort.Contains("chamber36") ) {  
        c_individual = new TCanvas("c_individual", " for a  bin , charge distribution ");
        c_individual->cd();
        gStyle->SetOptStat(111111211);
        proj->Draw();
        TString bin_nb = TString::Format("%d",j);
        proj->SetTitle("charge : "+rhidshort+" bin : "+bin_nb);
        c_individual->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"_new_2016/int_lum_bin_channel_"+rhidshort+"_bin_"+bin_nb+"_"+treename+".pdf");  
        } */
  //		  std::cout<<"Before adding the channels the integral for "<<rhidshort<<" bin "<<bin_nb<<" var "<<thevar<<" integral "<<proj->Integral()<<" entries "<<proj->GetEntries()<<std::endl;
        if(h==0 || h==771) {
        proj->Reset("ICESM");
        proj->ResetStats();
      }	
        double normal = 0;
        double chargemeantrimm = 0; double  integral = 0;
        int added_events=0;
        int chan_initial=0;
        int chan_final=0;
        if(h == 0){chan_initial=1; chan_final =400;}
        if(h == 771){chan_initial=401; chan_final =771;}
        //Special cases: all good channels in a single histo
        // you need to add charges from all the bins of rhid
        if(h == 0 || h==771){

            TCanvas *c;
            for(int ichan = chan_initial; ichan < chan_final ; ichan++){
        
              if( (1<= ichan && ichan<= 6) || (401<= ichan && ichan<=406))continue;
              if(ichan%10>=7 ||ichan%10==0) continue;	  
              if(ichan>370 &&ichan<400) continue;

              if (std::find(low_entry_channels.begin(), low_entry_channels.end(), h) != low_entry_channels.end()) {
                std::cout<<" we are going to skip thses entries from our cumulative calculations : channel : "<<h<<std::endl;
              continue;  // Skip the rest of the loop if h is in low_entry_channels
              }

//              std::cout<<"counted the entries to form our cumulative calculations : channel :"<<h<<std::endl;
              
              int theendcap = (ichan<=400)? 1:2;
              int thechamber = (ichan<=400)?  (int)floor(ichan/10)  : (int)floor( (ichan-400) /10) ;
              TString rhidshort ="chamber"+ (TString) Form("%d", thechamber)  +"_layer"+ (TString)Form("%d",ichan%10) + "_Endcap"+ theendcap;

              TH1D * h_prov = (TH1D*) ( myh->ProjectionX("_px",j,j, ichan,ichan) )->Clone();
  /*

        if(rhidshort.Contains("chamber1") || rhidshort.Contains("chamber4") ||rhidshort.Contains("chamber5") || rhidshort.Contains("chamber6") || rhidshort.Contains("chamber7") || rhidshort.Contains("chamber8") || rhidshort.Contains("chamber9") || rhidshort.Contains("chamber20") ||rhidshort.Contains("chamber21") || rhidshort.Contains("chamber24") || rhidshort.Contains("chamber25") || rhidshort.Contains("chamber36")  ) {  
        c_individual = new TCanvas("c_individual", " for 1st bin , charge distribution ");
        c_individual->cd();
        gStyle->SetOptStat(111111111);
        //gStyle->SetOptStat("ksiourmen");
        h_prov->Draw();
        TString s = TString::Format("%d",ichan);
        h_prov->SetTitle("charge : "+rhidshort+ " channel : "+s + " correspoinding :  "+h);
        c_individual->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/int_lum_bin_"+rhidshort+"_channel_"+s+"_bin_"+bin_nb+".pdf"); 
        }  */

              TH1D * h_prov_new = (TH1D*) ( myh->ProjectionX("_px",j,j, ichan,ichan) )->Clone();
              h_prov_new->Reset();
              h_prov_new->ResetStats();

              h_prov->Sumw2();
              h_prov_new->Sumw2();

              // atleast we will add only those channels for which we have more than 50 entries
              if(h_prov->Integral()<50) continue;

              normal = h_prov->Integral();
              // first trim the histogram then normalize with integral
              int last_bin = 0;
              double integral =0;
              for(int it = 1; it<=  h_prov->GetNbinsX() ;it++){
                  if(integral < trimmean*normal){
                  integral+=h_prov->GetBinContent(it);	
                  last_bin =it;
                  }
              }
              double new_integral=0; 
              int entries_last_bin;
              for(int it=1; it<last_bin; it++){
                  new_integral += h_prov->GetBinContent(it); 
              }
              entries_last_bin = (trimmean*normal - new_integral); 
                
              //if(debug_statements) std::cout<<"entries in last bin : rhid"<<rhidshort<<" var"<<thevar<<" bin number"<<j<<" before : "<<h_prov->GetBinContent(last_bin)<<" after :"<<entries_last_bin<<std::endl;
                
                  for(int it=1; it<last_bin ; it++) {
                    h_prov_new->SetBinContent(it,h_prov->GetBinContent(it));
                    h_prov_new->SetBinError(it,h_prov->GetBinError(it));
                  }
              
                  h_prov_new->SetBinContent(last_bin, entries_last_bin);
                  h_prov_new->SetBinError(last_bin,h_prov->GetBinError(last_bin) * (entries_last_bin / h_prov->GetBinContent(last_bin)));

                // h_prov_new->SetBinError(last_bin, h_prov->GetBinError(last_bin));
      
                  for(int it=last_bin+1; it<=h_prov->GetNbinsX() ; it++) {
                    h_prov_new->SetBinContent(it,0);
                    h_prov_new->SetBinError(it,0);
                  } 
                  
                  proj->Add(h_prov_new);// to each channel we normalize with respect to total integral acroos the channel , not wrt only trimmed mean integral
                  //normal = h_prov_new->Integral();

                //	proj->Add(h_prov_new,1./normal);// to each channel we normalize with respect to total integral acroos the channel , not wrt only trimmed mean integral
                  // new addition for evaluating mean based on individual channel
          

              // adding one more if condition to make plots for plus endcap and minus endcap separately 
                  delete h_prov;
                  delete h_prov_new;
                  // To make sure the memory is free after deleting, so that it do not give any issue
                  h_prov = nullptr;
                  h_prov_new = nullptr;
          } // End of loop of channels
        }// End of special case
        

        normal = proj->Integral();
        chargemeantrimm = 0; integral = 0;


        TH1D * h_trim = (TH1D * )proj->Clone();     
        gStyle->SetOptStat(111111211);



        TH1D * h_trim_new = (TH1D*) proj->Clone();

        // before truncating lets plot charge distribution 
        //Previously we were truncating only individual channel not all good channels 
        // truncation for all the channels, along with allgoodchannels
        // I am reseting stats, so new mean and integral is not interferred with old one
        // Storing the number of events in 1st bin of all the chambers 

      if(h!=0 && h!=771 && (j==12 || j==15)){ 
          int n_entries_bin = proj->Integral();
          if(j==12){ 
          num_entries_2D_hist_12th_bin->SetBinContent(chambernb, layernb, n_entries_bin);
          }
          if(j==15){ 
          num_entries_2D_hist_15th_bin->SetBinContent(chambernb, layernb, n_entries_bin);
          }
        } 
        if(h!=0 && h!=771){//Do the truncation, if h ==0 and h==771 the truncation is already done
            // The next feew stepas are to find till which bin we need to chopp off tail
            h_trim_new->Reset();
            h_trim_new->ResetStats();
            float normal = h_trim->Integral();
            int last_bin = 0;
            for(int it = 1; it<=   proj->GetNbinsX() ;it++) {
              if(integral < trimmean * normal){
                integral+=proj->GetBinContent(it);	
                last_bin =it;
              }
            }
					double new_integral =0; 
					int entries_last_bin;
					for(int it=1; it<last_bin; it++){
						new_integral += proj->GetBinContent(it); 
					}
					entries_last_bin = (int) (trimmean*normal - new_integral); 
					if(debug_statements) std::cout<<"entries in last bin : rhid"<<rhidshort<<" var"<<thevar<<" bin number"<<j<<" before : "<<h_trim->GetBinContent(last_bin)<<" after :"<<entries_last_bin<<std::endl;
					for(int it=1; it<last_bin ; it++) {
						h_trim_new->SetBinContent(it,h_trim->GetBinContent(it));
						h_trim_new->SetBinError(it,h_trim->GetBinError(it));
					}
					h_trim_new->SetBinContent(last_bin, entries_last_bin);
				  h_trim_new->SetBinError(last_bin, h_trim->GetBinError(last_bin) * (entries_last_bin / h_trim->GetBinContent(last_bin)));
					for(int it=last_bin+1; it<=proj->GetNbinsX() ; it++) {
						h_trim_new->SetBinContent(it,0);
						h_trim_new->SetBinError(it,0);
					}
					float final_integral = new_integral+entries_last_bin;
				} 

			 gStyle->SetOptStat("kKsSiourRmMen");
			 TString title_for_proj = rhidshort+" : "+thevar+" : bin : "+j;
			 proj->SetTitle(title_for_proj);
       // Saving before trim and after trim histograms for few chambers	
/*       if(proj->Integral()>0 && savehistos && (rhidshort.Contains("chamber4") || rhidshort.Contains("allgoodchannels")) ){
          std::cout<<"  entered in plots *******************************************************************"<<std::endl;
         
				 TCanvas *canvas_charge = new TCanvas();
 			 	 canvas_charge->Divide(2,1);
 			 	 canvas_charge->cd(1);

	    	 proj->Draw("E");
			   gStyle->SetOptStat(111111211);

 			 	 canvas_charge->cd(2);
	    	 //h_trim_new->Draw("E");
	    	 h_trim_new->Draw("E");
 			 	 canvas_charge->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"_test/charge_distribution"+title+"_"+rhidshort+"vs"+"_era"+dir_name_var+"_"+thevar+"_bin"+j+"_channel"+h+"_before_trim.pdf"); 

          //myfile<<j<<"\t"<<rhidshort<<"\t"<<h_trim_new->GetMean()<<"\t"<<h_trim_new->GetMeanError()<<std::endl;
			 }  */
       if(proj->Integral() < 50) continue;
			 // question is if already the histogram have more than 50 entries, then how the trimmed histogram can have just one entry. Not possible , but don't know why I applied additional cut. Are these two things not same ?  Integral is basically total bin height * mean value, however , 
  	   if(h_trim_new->Integral() ==1) continue;
       htrimmeanvsX->SetBinContent(j, h_trim_new->GetMean() ); 
     	 htrimmeanvsX->SetBinError(j, h_trim_new->GetMeanError() ) ; 


        delete proj;
        delete h_trim_new;
				delete h_trim;

        proj = nullptr;
        h_trim_new = nullptr;
        h_trim = nullptr;
       } //end of loop over  bins of variable of intereset

    if(htrimmeanvsX->Integral()==0) continue;
    TCanvas * c3 = new TCanvas;
    c3->cd();
		//std::cout<<"entering to save htrimmeanvsx "<<std::endl;
    //Some cosmetic stuff now...
    htrimmeanvsX->SetTitle(chamber_string_name+" : "+filename+title+"_"+rhidshort);
    //    htrimmeanvsX->GetYaxis()->SetRangeUser(0,600);
    if(thevar.Index("pressure")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("integratelumi")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("time")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 && rhidshort =="allgoodchannels_plus" )  htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 && rhidshort =="allgoodchannels_minus" )  htrimmeanvsX->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 )  htrimmeanvsX->GetXaxis()->SetRangeUser(0,20000);
    if(thevar.Index("time")>=0 ){
		 htrimmeanvsX->GetXaxis()->SetTimeDisplay(1);
		 htrimmeanvsX->GetXaxis()->SetLabelSize(0.02);
		 htrimmeanvsX->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    }
    htrimmeanvsX->GetXaxis()->SetTitle(xtitle);
		htrimmeanvsX->GetYaxis()->SetTitle("Trimmed mean charge");
		//if(thevar.Index("integratelumi")>=0 ) htrimmeanvsX->GetYaxis()->SetTitle("Normalized Trimmed mean charge");
		if(thevar.Index("integratelumi")>=0 ) htrimmeanvsX->GetYaxis()->SetTitle("Trimmed mean charge");
  
		htrimmeanvsX->SetMarkerStyle(20); htrimmeanvsX->SetMarkerSize(0.7);
    htrimmeanvsX->SetName(filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar);
    gStyle->SetOptStat("001111111");
		//Defines the range for the fit,    
    double fitlowedge (0), fithighedge(44);
    if(dir_name_var.Index("pressure_2016") >=0  ) fitlowedge = 951; 
    if(dir_name_var.Index("pressure_2017") >=0  ) fitlowedge = 949; 
    if(dir_name_var.Index("pressure_2018") >=0  ) fitlowedge = 953; 
    if(dir_name_var.Index("pressure_2016") >=0  ) fithighedge = 981; 
    if(dir_name_var.Index("pressure_2017") >=0  ) fithighedge = 979; 
    if(dir_name_var.Index("pressure_2018") >=0  ) fithighedge = 985; 

    if( thevar.Index("instlumi")>=0 ) fitlowedge = 1000; 
    if( thevar.Index("instlumi")>=0 ) fithighedge =20000 ;


			TF1 *fa1 =nullptr; 
      if(thevar.Index("pressure") >=0) {
		  	fa1= new TF1("fa1","exp([0]) * exp([1]*(x-967))", fitlowedge,fithighedge);
        fa1->SetParameters(5, -0.005); 
        if(htrimmeanvsX==NULL) continue;
		    htrimmeanvsX->Fit(fa1, "R");  
      }
     if(thevar.Index("instlumi")>=0){
			fa1= new TF1("fa1","exp([0]) * exp([1]*(x))", fitlowedge,fithighedge );
			fa1= new TF1("fa1","expo", fitlowedge,fithighedge );

     if(htrimmeanvsX==NULL) continue;
  	  htrimmeanvsX->Fit(fa1, "R"); }
      gStyle->SetOptStat("001111111");
      c3->SetName("c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar);



    
			if(htrimmeanvsX->Integral()>0) {

      bool flag_low_entry_channel=false;
      // to provide the flag for events with low rechits 
      if (std::find(low_entry_channels.begin(), low_entry_channels.end(), h) != low_entry_channels.end()) {
              std::cout<<" we are going to skip thses entries from our cumulative calculations : channel : "<<h<<std::endl;
              flag_low_entry_channel=true;
      }
        htrimmeanvsX->SetBinContent(htrimmeanvsX->GetNbinsX() + 1, flag_low_entry_channel ? 1 : 0);
        //std::cout<<" writing this "<<std::endl;
				htrimmeanvsX->Write();
			}
      // Calling the python function for the fit by passing the histogram 
     

			std::pair <double, double> parampair ;
			if(thevar.Index("integratelumi")>=0 ){ 
			 parampair.first = 0;
			 parampair.second =  0;		
			}
			if(thevar.Index("instlumi")>=0 ||thevar.Index("pressure")>=0   ){ 
			parampair.first =  fa1->GetParameter(0) ;
			parampair.second =  fa1->GetParameter(1) ;
			}
			//result[h] = parampair;
			// modifying result so that it contains slope from all good channels together
			result[h] = parampair;
			
			if(thevar.Index("integratelumi")>=0 ){ 
				double theslope =0; 
			}
  	  delete htrimmeanvsX;
      delete c3;

      htrimmeanvsX = nullptr;
      c3 = nullptr;

      delete fa1; 
      fa1= nullptr;

	} //end of loop for all the channelss 
 			if(num_entries_2D_hist_12th_bin!=0){
      num_entries_2D_hist_12th_bin->GetXaxis()->SetTitle("Chmaber nb");
      num_entries_2D_hist_12th_bin->GetYaxis()->SetTitle("layer nb");
      num_entries_2D_hist_12th_bin->Write();
      }
      if(num_entries_2D_hist_15th_bin!=0){
      num_entries_2D_hist_15th_bin->GetXaxis()->SetTitle("Chmaber nb");
      num_entries_2D_hist_15th_bin->GetYaxis()->SetTitle("layer nb");
      num_entries_2D_hist_15th_bin->Write();
      }
      num_entries_2D_hist_all_bins_plus->GetXaxis()->SetTitle("Chmaber nb");
      num_entries_2D_hist_all_bins_plus->GetYaxis()->SetTitle("layer nb");
      num_entries_2D_hist_all_bins_plus->Write();
      num_entries_2D_hist_all_bins_minus->GetXaxis()->SetTitle("Chmaber nb");
      num_entries_2D_hist_all_bins_minus->GetYaxis()->SetTitle("layer nb");
      num_entries_2D_hist_all_bins_minus->Write();   

      num_entries_1D_hist_all_bins_plus->GetXaxis()->SetRangeUser(0, num_entries_1D_hist_all_bins_plus->GetBinLowEdge(num_entries_1D_hist_all_bins_plus->GetNbinsX()+1));   
      num_entries_1D_hist_all_bins_minus->GetXaxis()->SetRangeUser(0, num_entries_1D_hist_all_bins_minus->GetBinLowEdge(num_entries_1D_hist_all_bins_minus->GetNbinsX()+1));   
      num_entries_1D_hist_all_bins_plus->Write();   
      num_entries_1D_hist_all_bins_minus->Write();   
  //myfile.close();
// 	delete h_slope;
//	delete h_chi2;
//	delete h_slopeuncty;
//  h_slope = nullptr;
//  h_chi2 = nullptr;
//  h_slopeuncty = nullptr;

	std::cout<<" done with one type of variable "<<thevar<<std::endl;
  return result;    
 
} //end of GetSlope function

double pressure_dependence_removal_instlumi::ApplyCorrection( double X ,TString correctiontype, double p0, double p1 ){
  double refvalue = 0;
  if(correctiontype =="pressure"&&!droppressurecorr){
    refvalue =967;
    double thecorr =exp(p1*(refvalue-X)); 

    return thecorr; 
  }

  
  if(correctiontype =="instlumi"&& !dropinstlumicorr){
    refvalue =15000; 
    double thecorr = exp(p1*(refvalue-X)); 

    return thecorr; 
    
  }
  
  return 1;
}

double pressure_dependence_removal_instlumi::NominalHV(){

    std::pair<double,double> chargeandHV(0,0);
    double dHV_(0),HV_(0);
  if(_stationring ==11 || _stationring ==14){

    if( _stationring ==14) _rhid -=30000 ;
    _rhid -= _rhid%10; _rhid/=10;

    dHV_ = (dHV_HV_ME11(_rhid) ).first;
    HV_ = (dHV_HV_ME11(_rhid) ).second;
    
  }
  else if(_runNb>= 281613){
    dHV_ = (dHV_HV_NonME11_v2(_rhid) ).first;
    HV_ = (dHV_HV_NonME11_v2(_rhid) ).second;
  } 
  else{
    dHV_ = (dHV_HV_NonME11_v1(_rhid) ).first;
    HV_ = (dHV_HV_NonME11_v1(_rhid) ).second;
  }  

  // all ME11 set to 2900 V and nonME11 to 3600 V
  if( _runNb <277792 ){HV_ -=  dHV_; dHV_=0;}
  
  // all ME11 set to 2900 V , non ME11 have new values
  if( _runNb < 281613&& (_stationring==11|| _stationring ==14)  ){HV_ -=  dHV_; dHV_=0;}

  // all ME11 above 281613 have new values and all non ME11 above 277792 have new values
  else if(_runNb >= 324077 && (_stationring==11|| _stationring ==14))  { HV_= HV_ -32 ; dHV_= dHV_-32; }
  // all ME11 in 2018 now have correct values 
   // Non ME11 outer rings voltage was lowered by 35V
  // ME12, ME13, ME22, ME32, ME42 
  else if(_runNb >= 324077 && (_stationring==12|| _stationring ==13|| _stationring ==22|| _stationring ==32|| _stationring ==42))  { HV_= HV_ -    35 ; dHV_= dHV_-35; }

  return HV_;
}

