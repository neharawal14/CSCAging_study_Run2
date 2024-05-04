#define pressure_dependence_removal_cxx
#include "pressure_dependence_removal.h"
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
using namespace std;
enum ring_station_hvsegm {
  me11a,me11b,me12HV1,me12HV2,me12HV3,me13HV1,me13HV2,me13HV3,
  me21HV1,me21HV2,me21HV3,me22HV1,me22HV2,me22HV3,me22HV4,me22HV5,
  me31HV1,me31HV2,me31HV3,me32HV1,me32HV2,me32HV3,me32HV4,me32HV5,
  me41HV1,me41HV2,me41HV3,me42HV1,me42HV2,me42HV3,me42HV4,me42HV5
};
bool dropinstlumicorr =false; // Drop inst L correction
bool droppressurecorr =false;// Drop pressure correction
bool dropcurrentcorr = true;
bool dropHVequalisation = false;
bool seconditer =false;// If you want to perform a second iteration of the fit to pressure and inst L.
bool savehistos =false;

double trimmean =0.7;
//Offset used in the calculation of the trimmed mean bool istest =false; //For tests, run only on 1/100 of the events
bool debug_plots = true;
bool debug_statements =false;
bool seconditer_plots =false;// If you want to perform a second iteration of the fit to pressure and inst L.

bool corrfrom25to36fbonly = true; // Extract pressure dependency from the range 15 to 20 /fb, 7e33 to 9e33, extract inst L. from the range 15 to 20 /fb 

bool debug_print = false;
bool timesecond_plot = false;
bool iszmumu = false;//If only Zmumu events are used the corrections cannot be derived channel per channel (not enough stat) so we take all channels together

TH3D * hchargevspressure_2; 
TH3D* hchargevsinstlumi_2;
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


bool badrun(int runnb){// A few runs that turned out to be bad or for which some info (inst L, pressure) is missing
  if(runnb==273150)return true;
  if(runnb==273426)return true;
  if(runnb==274157)return true;
  if(runnb==274443)return true;
  if(runnb==275064)return true;
  if(runnb==275285)return true;
  if(runnb==275286)return true;
  if(runnb==275289)return true;
  if(runnb==275757)return true;
  if(runnb==275758)return true;
  if(runnb==275759)return true;
  if(runnb==275781)return true;
  if(runnb==275841)return true;
  if(runnb==275846)return true;
  if(runnb==275887)return true;
  if(runnb==275922)return true;
  if(runnb==276064)return true;
  if(runnb==276095)return true;
  if(runnb==276832)return true;
  if(runnb==278309)return true;
  if(runnb==278821)return true;
  if(runnb==279995)return true;
  if(runnb==280002)return true;
  if(runnb==280006)return true;
  if(runnb==280007)return true;
  if(runnb==281691)return true;
  return false;

}

void pressure_dependence_removal::Loop(TString input_file_path, TString input_file_name, TString chamber_string, TString output_file_path, TString output_folder_name)
{
  chamber_string_name = chamber_string;
	detregionstr="_dataset";
	output_path = output_file_path;
	output_plots_folder = output_folder_name;
  gStyle->SetOptStat();
  gStyle->SetOptFit(111);
  if(detregionstr.Index("ZMuMu")>=0) iszmumu=true;
  iszmumu = false;


	time_t initial_time, final_time;
	 
//	if(year=="2016"){
//	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
//	  final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
//	}
//	if(year=="2017"){
//	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
//	 final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
//	}
//	if(year=="2018"){
//	 	initial_time	= 1462060800 ; // 1 May 2016 : 00 : 00 : 00
//	 final_time = 1477958400; // 31 Dec 2016 : 00 : 00 : 00
//	}
//

	gStyle->SetTimeOffset(0.);
	//  TH1F * chargepresscorrhighinstlumi = new TH1F("chargepresscorrhighinstlumi","",10000,0,10000);
	//  TH1F * chargepresscorrlowinstlumi = new TH1F("chargepresscorrlowinstlumi","",10000,0,10000);
   // we are passing name of the input file when calling the macro 
   TString inputfname = input_file_path+ input_file_name;
   TFile * inputf = TFile::Open(inputfname);
	 //   all these conditions are checked to know what type of cuts you are applying
   //   if(istest)detregionstr =detregionstr+"test"; 
   if(seconditer)detregionstr =detregionstr+"seconditer"; 
   if(corrfrom25to36fbonly)detregionstr =detregionstr+"_pressure_corrected_"; 
   //if(check7to9e33only) detregionstr+="check7to9e33only";
   // if(dropinstlumicorr)detregionstr+="noinstlumicorr";
 	 // if(droppressurecorr)detregionstr+="nopressurecorr";

   TTree * tree =(TTree*) inputf->Get("tree");
	 // Init initializes all the branches in this tree
   Init(tree);
	 // output root file after processing 
   TFile * outf = new TFile(output_path+"outf"+detregionstr+"_"+chamber_string_name+"_output_everything_new.root","recreate");  
	 // Making IntegratedLumi vs gas gain slope dependency before starting any pressure and inst lumi corrections
   TH3D * hchargevsintegratelumi_initial = new TH3D("hchargevsintegratelumi_initial","charge (ADC counts) vs integ lumi (initial)",3000,0,3000, 160, 0,160  ,770,1,771);
	 // 5 days is 1 bin
//   TH3D * hchargevstime_initial = new TH3D("hchargevstime_initial","charge (ADC counts) vs time (initial)",3000,0,3000, 60, initial_time,final_time ,770,1,771);
	 // Loop over tree and reach charge for each rechit
   for(int i = 0; i < tree->GetEntries(); i++){
  // for(int i = 0; i < 10000000; i++){
         //if(_integratelumi==0) continue;
         //if(i%100 !=0 &&istest) continue;
	 	     LoadTree(i);tree->GetEntry(i);
         if(i%1000000 ==0)cout << i<<endl;
         // if(badrun(_runNb) ) continue;
         //   if(check7to9e33only  &&_instlumi>9000 ) continue;
         //  if(check7to9e33only  &&_instlumi<7000 ) continue;
         // what is significance of this ?
		     // if(check7to9e33only) _integratelumi -= 17000;
				 int rhidreduced = ((int)floor(_rhid/10))%1000;
         if(_rhid> 2000000) rhidreduced +=400;
         int idforcorr = (iszmumu )? 0: rhidreduced ;
				 double charge = _rhsumQ_RAW;
				 std::pair<float, float> 	 HV_to_equalise = UncorrGasGain_HVInitial2016(sumq,fRun,_stationring,_rhid);
				 float HV_to_equalise = HV_to_equalise.first;
				 double charge_equalised = _rhsumQ_RAW * HV_to_equalise;
 				 hchargevsintegratelumi_initial->Fill(charge, _integratelumi, rhidreduced);
 //        hchargevstime_initial->Fill(charge, _timesecond, rhidreduced);
   }
	 std::cout<<"working for information with integrate luminoisty"<<hchargevsintegratelumi_initial->GetEntries()<<std::endl;
  //vector< std::pair<double, double > > params_integratelumi_initial;
//  params_integratelumi_initial = GetSlope( hchargevsintegratelumi_initial, "_integratelumi_initial", detregionstr,"",outf);

	std::cout<<" done with intlumi information"<<std::endl;
    //outf->cd();
		//hchargevsintegratelumi_initial->Write();
		//outf->Close();
//	 std::cout<<"working for time second information"<<std::endl;
//	 vector< std::pair<double, double > > params_timesecond_initial;
//   params_timesecond_initial = GetSlope( hchargevstime_initial, "_timesecond_initial", detregionstr,"",outf);

   //Run on all events to extract pressure correction for each channel (rechit) separately. 
   //This is done with the help of a 3d histogram, storing charge, pressure and local rechit ID (inside a given station/ring/HV chamber)
 
double intlumi_low_cut_2016 = 12; 
double intlumi_up_cut_2016 = 16; 
double instlumi_low_cut_2016 = 7000; 
double instlumi_up_cut_2016 = 9000; 
double intlumi_low_cut_2017 = 55; 
double intlumi_up_cut_2017 = 60 ; 
double instlumi_low_cut_2017 = 7000; 
double instlumi_up_cut_2017 = 9000; 
double intlumi_low_cut_2018 = 115; 
double intlumi_up_cut_2018 = 133 ; 
double instlumi_low_cut_2018 = 10000; 
double instlumi_up_cut_2018 = 15000; 

double integratelumi_2016_high = 39.32673126400002;
double integratelumi_2017_high = 83.85340826572357;

bool debug = false;
int n_entries_2017 = 0;
	 std::cout<<"going to pressure information"<<std::endl;
   TH3D * hchargevspressure_2016 = new TH3D("hchargevspressure_2016","charge (ADC counts) vs pressure : 2016",3000,0,3000, 40, 943,984  ,770,1,771);
   TH3D * hchargevspressure_2017 = new TH3D("hchargevspressure_2017","charge (ADC counts) vs pressure : 2017",3000,0,3000, 40, 943,984  ,770,1,771);
   TH3D * hchargevspressure_2018 = new TH3D("hchargevspressure_2018","charge (ADC counts) vs pressure : 2018",3000,0,3000, 40, 943,984  ,770,1,771);
   for(int i = 0; i < tree->GetEntries(); i++){
   //for(int i = 0; i < 10000000; i++){
      //     if(droppressurecorr)break;
      //     if(i%100 !=0&&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl; if(badrun(_runNb) ) continue;
		 if(_integratelumi <= integratelumi_2016_high) { 
			 if(_integratelumi<intlumi_low_cut_2016 ) continue;
 	     if(_integratelumi>intlumi_up_cut_2016) continue;
       if(_instlumi<instlumi_low_cut_2016) continue;
       if(_instlumi>instlumi_up_cut_2016 ) continue;
     	 int rhidreduced = ((int)floor(_rhid/10))%1000;
       if(_rhid> 2000000) rhidreduced +=400;//First (second) endcap have rechit ID < (>) 2000000
     //     Reduced rechit ID has the following format: (A+B)*10+C, where A=1,..36 (chamber nb), B =0 (endcap 1) or 40 (endcap 2)  and C =1,..., 6 (layer).
		 
//		 if(isbadchannel("ME11",rhidreduced,_runNb, _integratelumi)  ) continue;
       hchargevspressure_2016->Fill(_rhsumQ_RAW, _pressure , rhidreduced);
			 if(debug) std::cout<<" inside the 2016 instlumi loop"<<std::endl;
     }
		 if( integratelumi_2016_high <_integratelumi && _integratelumi<= integratelumi_2017_high) { 
			 if(_integratelumi<intlumi_low_cut_2017 ) continue;
 	     if(_integratelumi>intlumi_up_cut_2017) continue;
       if(_instlumi<instlumi_low_cut_2017) continue;
       if(_instlumi>instlumi_up_cut_2017 ) continue;
     	 int rhidreduced = ((int)floor(_rhid/10))%1000;
       if(_rhid> 2000000) rhidreduced +=400;//First (second) endcap have rechit ID < (>) 2000000

			 n_entries_2017 = n_entries_2017+1;
       hchargevspressure_2017->Fill(_rhsumQ_RAW, _pressure , rhidreduced);
			 ///if(debug) std::cout<<" inside the 2017 instlumi loop"<<std::endl;
     }
		 if( integratelumi_2017_high <_integratelumi) { 
			 if(_integratelumi<intlumi_low_cut_2018 ) continue;
 	     if(_integratelumi>intlumi_up_cut_2018) continue;
       if(_instlumi<instlumi_low_cut_2018) continue;
       if(_instlumi>instlumi_up_cut_2018 ) continue;
     	 int rhidreduced = ((int)floor(_rhid/10))%1000;
       if(_rhid> 2000000) rhidreduced +=400;//First (second) endcap have rechit ID < (>) 2000000
       hchargevspressure_2018->Fill(_rhsumQ_RAW, _pressure , rhidreduced);

			 if(debug) std::cout<<" inside the 2018 instlumi loop"<<std::endl;
     }
	 }

    vector< std::pair<double, double > > params_pressure_2016;
    vector< std::pair<double, double > > params_pressure_2017;
    vector< std::pair<double, double > > params_pressure_2018;
		if(hchargevspressure_2016 != NULL && hchargevspressure_2016->GetEntries() >=100){
    params_pressure_2016 = GetSlope( hchargevspressure_2016, "_pressure", detregionstr,"",outf);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
		if(hchargevspressure_2017 != NULL && hchargevspressure_2017->GetEntries() >=100){
			std::cout<<" entries "<<n_entries_2017<<std::endl;
			std::cout<<" entries  Hcharge "<<hchargevspressure_2017->GetEntries()<<std::endl;
     params_pressure_2017 = GetSlope(hchargevspressure_2017, "_pressure", detregionstr,"",outf);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
		if(hchargevspressure_2018 != NULL&& hchargevspressure_2018->GetEntries() >=100){
    params_pressure_2018 = GetSlope( hchargevspressure_2018, "_pressure", detregionstr,"",outf);  // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   
		}
//	outf->Close();
	std::cout<<" after closing the file "<<std::endl;


////	 // For now, I am just applying pressuer corrections and studying gas gain wrt integrated luminosity
//// //Now, run on all events, apply pressure correction and extract inst L correction
////  /*   TH3D * hchargevsinstlumi = new TH3D("hchargevsinstlumi","charge (ADC counts) vs inst lumi",3000,0,3000, 42, 0,21000  ,770,1,771);
////  double charge;   
////    for(int i = 0; i < tree->GetEntries(); i++){
////      //     if(dropinstlumicorr)break;
////      //    if(i%100 !=0 &&istest) continue;
////     LoadTree(i);tree->GetEntry(i);
////     //   if(i%1000000 ==0)cout << i<<endl;
////     if(badrun(_runNb) ) continue;
////     if(corrfrom25to36fbonly &&_integratelumi>45) continue;
////     if(corrfrom25to36fbonly &&_integratelumi<25) continue;
//// 
////     int rhidreduced = ((int)floor(_rhid/10))%1000;
////     if(_rhid> 2000000) rhidreduced +=400;
////     
////     int idforcorr = (iszmumu )? 0: rhidreduced ;
////		   
////		 if(rhidreduced <400 ) { charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second ); }
////		 if(rhidreduced >=400 ) { charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second );}
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  ;
////		 // modifying the corections, so that the same correction applies to all channels, and that correction is derived from all good channels slope and that's why idforcorr = 0 for this case
////     hchargevsinstlumi->Fill(charge, _instlumi, rhidreduced);
////
////     if(debug_print)		std::cout<<"charge after pressure correction "<<charge<<std::endl;
////
//////		 bool badchannel = isbadchannel(detregionstr , rhidreduced);
//////     if(!badchannel&&_instlumi>8000&&_instlumi<10000&&_integratelumi>12000&&_integratelumi<16000) 
////			 //chargepresscorrhighinstlumi->Fill(charge);
//////		 if(!badchannel&&_instlumi>3000&&_instlumi<5000&&_integratelumi>12000 &&_integratelumi<16000 )
////			 //chargepresscorrlowinstlumi->Fill(charge);
////     
////   }
////
////   vector< std::pair<double, double > > params_instlumi;
////   params_instlumi = GetSlope( hchargevsinstlumi, "_instlumi", detregionstr,"",outf); 
////	*/ 
////  //
/////*
////
//// std::cout<<"issue after 1st set of iteration "<<std::endl;
////   //Second iteration (optional): same game
////
////   vector< std::pair<double, double > > params_pressure_2;
////   vector< std::pair<double, double > > params_instlumi_2;
////  
////   if(seconditer_plots){
////      hchargevspressure_2 = new TH3D("hchargevspressure_2","charge (ADC counts) vs pressure (2nd iteration)",3000,0,3000, 40, 944,984  ,770,1,771);
////
////   //Pressure correction
////   for(int i = 0; i < tree->GetEntries(); i++){
////
//////     if(i%100 !=0&&istest) continue;
////     LoadTree(i);tree->GetEntry(i);
////     if(i%1000000 ==0)cout << i<<endl;
////     if(badrun(_runNb) ) continue;
////     if(corrfrom25to36fbonly &&_integratelumi>36 ) continue;
////     if(corrfrom25to36fbonly &&_integratelumi<25) continue;
////     if(corrfrom25to36fbonly &&_instlumi<14000 ) continue;
////     if(corrfrom25to36fbonly &&_instlumi>19000) continue;
////
////     int rhidreduced = ((int)floor(_rhid/10))%1000;
////     if(_rhid> 2000000) rhidreduced +=400;
////     int idforcorr = (iszmumu )? 0: rhidreduced ;
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second );
////		 // modifying the corections, so that the same correction applies to all channels, and that correction is derived from all good channels slope and that's why idforcorr = 0 for this case
////     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
////
////	 //	 * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
////     hchargevspressure_2->Fill(charge, _pressure , rhidreduced);
////     }
////    params_pressure_2 = GetSlope( hchargevspressure_2, "_pressure_2", detregionstr,"",outf);
////   } 
////   
////   //Now charge after pressure correction
////   
////    hchargevsinstlumi_2 = new TH3D("hchargevsinstlumi_2","charge (ADC counts) vs inst lumi",150,0,3000,42, 0,21000  ,770,1,771);
////   
////   for(int i = 0; i < tree->GetEntries(); i++){
////     
////   //  if(i%100 !=0 &&istest) continue;
////     LoadTree(i);tree->GetEntry(i);
////     
////     if(i%1000000 ==0)cout << i<<endl;
////    // if(badrun(_runNb) ) continue;
////     if(corrfrom25to36fbonly &&_integratelumi<25) continue;
////     if(corrfrom25to36fbonly &&_integratelumi>36) continue;
////     
////     int rhidreduced = ((int)floor(_rhid/10))%1000;
////     if(_rhid> 2000000) rhidreduced +=400;
////     int idforcorr = (iszmumu )? 0: rhidreduced ;
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second )  ;
////		 // modifying the corections, so that the same correction applies to all channels, and that correction is derived from all good channels slope and that's why idforcorr = 0 for this case
////     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[0]).first , (params_instlumi[0]).second ) * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[0]).first , (params_pressure_2[0]).second )  ;
////     hchargevsinstlumi_2->Fill( charge, _instlumi, rhidreduced);
////     
////   }
////   params_instlumi_2 = GetSlope( hchargevsinstlumi_2, "_instlumi_2", detregionstr,"",outf);
//// */   
////  // 
////	 std::cout<<" pressure corrections derived , going for application "<<std::endl;
//////	 TFile *file_output = new TFile(output_path+"output_after_pressure_depenedence_removal_"+chamber_string_name +"_everything_separate.root","RECREATE");
//////   tree_new = new TTree("tree_new", "tree");
//////	 tree_new->SetDirectory(0); 
//////   Setup_new_tree(); 
////   //Now final charge after pressure, instlumi correction
   TH3D * hchargevsintegratelumi = new TH3D("hchargevsintegratelumi","charge (ADC counts) vs integ lumi",3000,0,3000, 160, 0,160  ,770,1,771);
 //  TH3D * hchargevstime = new TH3D("hchargevstime","charge (ADC counts) vs time",3000,0,3000, 60, initial_time,final_time  ,770,1,771);

	// std::cout<<" pressure corrections derived , filling the root file "<<std::endl;
   //for(int i = 0; i < tree->GetEntries(); i++){
//	 int j;
//	 int save_entries =  100000;

   	//tree_new->SetAutoFlush(-30000000); //("", TObject::kOverwrite);
 //  int num_cycles = tree->GetEntries()/save_entries;
//	 if(num_cycles!=0) {
//   for(int j = 0; j <= num_cycles -1; j++){
//   for(int i = j*1000000; i < (j+1)* 1000000; i++){
     // std::cout<<"entered into integratedlumi loop correction "<<i<<std::endl; 
		 //if(_integratelumi==0) continue;
     // std::cout<<"case for which integratedlumi correction exists "<<i<<std::endl; 
     //if(i%100 !=0 &&istest) continue;
		 //tree_new->SetMaxTreeSize(Long64_t size)
		 //int half_entries = (tree->GetEntries()/2);
		 // When we apply correction we should take average dependence for pressure : avg of plus and minus endcap
		 for(int i=0 ; i<tree->GetEntries(); i++){
		 //for(int i=0; i<10000000; i++){

     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     // what is significance of this ?
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;
	//	 if(isbadchannel("ME11",rhidreduced,_runNb,_integratelumi)  ) continue;
     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
     double charge, charge_equalized; 
      
		//applying pressure correction by taking average of slope dependence in plus and minus endcap 
		 //double params_pressure_const_avg =  (params_pressure[0].first + params_pressure[771].first)/2. ;
		 //double params_pressure_slope_avg =  (params_pressure[0].second + params_pressure[771].second)/2. ;
		 //charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   params_pressure_const_avg ,params_pressure_slope_avg );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure", params_pressure_const_avg ,params_pressure_slope_avg ) ; 
    if(_integratelumi <= integratelumi_2016_high){
		 if(rhidreduced <400 ) {
			 if(debug) std::cout<<" inside the pressure correction for 2016"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2016[0]).first , (params_pressure_2016[0]).second );
//    charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
		 }
     if(rhidreduced >=400 ) { 
			 if(debug) std::cout<<" inside the pressure correction for 2016"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2016[771]).first , (params_pressure_2016[771]).second );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second ); 
		 }
		}
	if( integratelumi_2016_high < _integratelumi && _integratelumi<= integratelumi_2017_high) { 
		 if(rhidreduced <400 ) {
			 if(debug) std::cout<<" inside the pressure correction for 2017"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2017[0]).first , (params_pressure_2017[0]).second );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
		 }
     if(rhidreduced >=400 ) { 
			 if(debug) std::cout<<" inside the pressure correction for 2017"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2017[771]).first , (params_pressure_2017[771]).second );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second ); 
		 }
		}
    if(_integratelumi > integratelumi_2017_high){
		 if(rhidreduced <400 ) {
			 if(debug) std::cout<<" inside the pressure correction for 2018"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2018[0]).first , (params_pressure_2018[0]).second );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
		 }
     if(rhidreduced >=400 ) { 
			 if(debug) std::cout<<" inside the pressure correction for 2018"<<" event "<<i<<std::endl;
			 charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2018[771]).first , (params_pressure_2018[771]).second );
     //charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second ); 
		 }
		}

	 // ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
   //		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;
  //		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second );
//		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[0]).first , (params_pressure_2[0]).second );
 //	 * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;

 hchargevsintegratelumi->Fill(charge, _integratelumi, rhidreduced);
// hchargevsintegratelumi->Fill(_rhsumQ_RAW, _integratelumi, rhidreduced);
 }
//////   hchargevstime->Fill(charge, _timesecond, rhidreduced);
//////   new_eventNb       = _eventNb       ; 
//////   new_runNb         = _runNb         ;  
//////   new_lumiBlock     = _lumiBlock     ;
//////   new_rhid          = _rhid          ;  
//////   new_stationring   = _stationring   ;
//////   new_rhsumQ        = charge_equalized        ; 
//////	 new_HV           =  _HV;
//////	 new_current           =  _current;
//////   new_rhsumQ_RAW    = charge    ;
//////   new_pressure      = _pressure      ; 
//////   new_temperature   = _temperature   ;
//////   new_instlumi      = _instlumi      ; 
//////   new_integratelumi = _integratelumi ;
//////   new_timesecond    = _timesecond    ;
//////   new_n_PV          = _n_PV          ;
//////   new_bunchcrossing = _bunchcrossing ;
//////	 tree_new->Fill();
////
//////	 if(i%save_entries ==0) {
//////   	tree_new->Write("", TObject::kOverwrite);
//////	 }
////	// }
////  // tree_new->Write("", TObject::kOverwrite);
////	// }
////    /*if(num_cycles==0) { j=0;} 
////    if(j==num_cycles || num_cycles==0){
////	 		for(int i = j*1000000; i < tree->GetEntries(); i++){
////     // std::cout<<"entered into integratedlumi loop correction "<<i<<std::endl; 
////		 //if(_integratelumi==0) continue;
////     // std::cout<<"case for which integratedlumi correction exists "<<i<<std::endl; 
////     //if(i%100 !=0 &&istest) continue;
////     LoadTree(i);tree->GetEntry(i);
////     if(i%1000000 ==0)cout << i<<endl;
////     if(badrun(_runNb) ) continue;
////     if(check7to9e33only  &&_instlumi<8000 ) continue;
////     if(check7to9e33only  &&_instlumi>14000 ) continue;
////     // what is significance of this ?
////  	 if(check7to9e33only) _integratelumi -= 17000;
////     int rhidreduced = ((int)floor(_rhid/10))%1000;
////     if(_rhid> 2000000) rhidreduced +=400;
////     int idforcorr = (iszmumu )? 0: rhidreduced ;
////  //	 if(isbadchannel("ME11",rhidreduced,_runNb,_integratelumi)  ) continue;
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
////     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
////     double charge, charge_equalized; 
////
////  	 if(rhidreduced <400 ) { charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second );
////     charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[0]).first , (params_pressure[0]).second ); }
////     if(rhidreduced >=400 ) { charge  = _rhsumQ_RAW * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second );
////     charge_equalized  = _rhsumQ * ApplyCorrection( _pressure ,"pressure",   (params_pressure[771]).first , (params_pressure[771]).second ); }
////
////  	 // ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
////     //		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;
////    //		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second );
//////		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[0]).first , (params_pressure_2[0]).second );
////   //	 * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;
////
//////   hchargevsintegratelumi->Fill(charge, _integratelumi, rhidreduced);
//////   hchargevstime->Fill(charge, _timesecond, rhidreduced);
//// //  new_eventNb       = _eventNb       ; 
//// //  new_runNb         = _runNb         ;  
//// //  new_lumiBlock     = _lumiBlock     ;
////   new_rhid          = _rhid          ;  
//// //  new_stationring   = _stationring   ;
////   new_rhsumQ        = charge_equalized        ; 
////   new_HV           =  _HV;
////   new_rhsumQ_RAW    = charge    ;
////   new_pressure      = _pressure      ; 
//// //  new_temperature   = _temperature   ;
//// //  new_instlumi      = _instlumi      ; 
//// //  new_integratelumi = _integratelumi ;
////   new_timesecond    = _timesecond    ;
//// //  new_n_PV          = _n_PV          ;
//// //  new_bunchcrossing = _bunchcrossing ;
////   tree_new->Fill();
////  }
////   tree_new->Write("", TObject::kOverwrite);
////  }*/ 
////	 //
	 std::cout<<"after filling th eluminoisty distribution "<<std::endl; 
   vector< std::pair<double, double > > params_integratelumi;
   params_integratelumi = GetSlope( hchargevsintegratelumi, "_integratelumi", detregionstr,"",outf); 
	 std::cout<<"working for time second information after correction "<<std::endl; 
//	 vector< std::pair<double, double > > params_timesecond;
 //  params_timesecond = GetSlope( hchargevstime, "_timesecond", detregionstr,"",outf);
//
	 //file_output->cd();
	 std::cout<<" error in root file  filling the root file "<<std::endl;
	 }
////	 //tree_new->Write();
////	 //file_output->Close(); 
////	
/////*
////	 int nbins_time = hchargevstime->GetNbinsX();
////	 for(int i=1; i<= nbins_time ; i++){
////		 double bin_content = hchargevstime->GetBinContent(i);
////		 double binLowEdge = hchargevstime->GetXaxis()->GetBinLowEdge(i);
////		 double binUpEdge = hchargevstime->GetXaxis()->GetBinLowEdge(i+1);
//// 
////		 if (bin_content == 0 || TMath::IsNaN(bin_content)) continue;
////	   outfile<<bin_content<<"\t"<<binLowEdge<<"\t"<<binUpEdge<<std::endl;
////	 } */
//// // outf->Close(); 
////   std::cout<<"issue realized after coming to end"<<std::endl;
  
//}

vector < std::pair<double, double > >  pressure_dependence_removal::GetSlope( TH3D * myh , TString thevar , TString filename, TString title, TFile * outf){

	TString name_histogram = myh->GetName();

   if(debug_print) std::cout<<"inside the slope function for the var "<<thevar<<std::endl;
  //Will store all the fit results in a TTree (one entry per channel)
  TString treename_goodchannels = "tree_all_goodchannels"+thevar ; 
  std::cout<<" name of histgoram "<<name_histogram<<std::endl;
	if(name_histogram=="hchargevspressure_2016") treename_goodchannels = treename_goodchannels+"_2016";
	else if(name_histogram=="hchargevspressure_2017") treename_goodchannels = treename_goodchannels+"_2017";
	else if(name_histogram=="hchargevspressure_2018") treename_goodchannels = treename_goodchannels+"_2018";

	std::cout<<" name of the tree good channels ************************ pressure **********"<<treename_goodchannels<<std::endl;
  TTree * theouttree_goodchannels = new TTree(treename_goodchannels,"");
  TString treename ; 	
	if(name_histogram=="hchargevspressure_2016") treename = "tree_"+thevar+"_2016";
	else if(name_histogram=="hchargevspressure_2017") treename ="tree_"+thevar+ "_2017";
	else if(name_histogram=="hchargevspressure_2018") treename = "tree_"+thevar+"_2018";
  else
  treename= "tree_"+thevar ; 


  TTree * theouttree = new TTree(treename,"");
  Float_t _slope(-1),_slope_error(-1),_chi2(-1);
	Float_t _slope_goodchannels_plus(-1), _slope_error_goodchannels_plus(-1);
	Float_t _slope_goodchannels_minus(-1), _slope_error_goodchannels_minus(-1);
  Int_t _ndof(-1),_layer(-1),_chamber(-1),_endcap(-1),_stationringHVseg(-1); 
  Bool_t _isbadchannel(false);
  theouttree_goodchannels->Branch("_slope_goodchannels_plus",&_slope_goodchannels_plus,"_slope_goodchannels_plus/F");
  theouttree_goodchannels->Branch("_slope_goodchannels_minus",&_slope_goodchannels_minus,"_slope_goodchannels_minus/F");
  theouttree_goodchannels->Branch("_slope_error_goodchannels_plus",&_slope_error_goodchannels_plus,"_slope_error_goodchannels_plus/F");
  theouttree_goodchannels->Branch("_slope_error_goodchannels_minus",&_slope_error_goodchannels_minus,"_slope_error_goodchannels_minus/F");
  theouttree->Branch("_slope",&_slope,"_slope/F");
  theouttree->Branch("_slope_error",&_slope_error,"_slope_error/F");
  theouttree->Branch("_chi2",&_chi2 ,"_chi2/F");
  theouttree->Branch("_ndof",&_ndof,"_ndof/I");
  theouttree->Branch("_layer",&_layer,"_layer/I");
  theouttree->Branch("_chamber",&_chamber,"_chamber/I");
  theouttree->Branch("_endcap",&_endcap,"_endcap/I");
  theouttree->Branch("_stationringHVseg",&_stationringHVseg,"_stationringHVseg/I");
  theouttree->SetAutoSave(1000000);
  TString xtitle ;
  if(thevar.Index("pressure")>=0) xtitle = "Pressure (hPa)";
  if(thevar.Index("instlumi")>=0) xtitle = "Inst Lumi (#mub s)^{-1}";
  if(thevar.Index("integrate")>=0) xtitle = "Integrated luminosity (fb^{-1})";
  if(thevar.Index("time")>=0) xtitle = "time";

	// Taking projection of charge distribution on the variable on interest axis
  //	TH1D* projglobal = (TH1D*) myh-> ProjectionX("_px",0, myh->GetNbinsY());
  // projglobal->Sumw2();
  /* TCanvas * c = new TCanvas;
  c->cd();
  projglobal->SetTitle(filename+title);
  projglobal->SetMarkerStyle(20);
  projglobal->SetMarkerSize(0.7);
  projglobal->GetXaxis()->SetTitle("Charge");
  projglobal->GetYaxis()->SetTitle("Normalized nb of entries");
  projglobal->GetYaxis()->SetTitleSize(0.05);
  projglobal->DrawNormalized("lep");
  projglobal->SetName("charge"+filename+title); */
  //  if(thevar.Index("integrate")>=0 && projglobal->Integral()>0 &&savehistos)
  //    c->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/charge_"+thevar+"_"+filename+".pdf");
  
	outf->cd();

	TString dir_name_var ;
	if(name_histogram=="hchargevspressure_2016") dir_name_var = thevar+"_2016";
	else if(name_histogram=="hchargevspressure_2017") dir_name_var = thevar+"_2017";
	else if(name_histogram=="hchargevspressure_2018") dir_name_var = thevar+"_2018";
	else 
		dir_name_var = thevar;

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
	TH1D * h_slope = new TH1D ("h_slope"+thevar+"_"+filename,"",200, lowedge, highedge );
  TH1D * h_slopeuncty = new TH1D ("h_slopeuncty"+thevar+"_"+filename,"",200, lowedge/100., highedge/100. );
  TH1D * h_chi2 = new TH1D ("h_chi2"+thevar+"_"+filename,"",200,0,1000);
	// Loop over reduced rechit ID
  // Reduced rechit ID has the following format: (A+B)*10+C, where A=1,..36 (chamber nb), B =0 (endcap 1) or 40 (endcap 2)  and C =1,..., 6 (layer). 
	// For calculating normalized charge distribution for all layers together
	// h=0 corresponds to the plus endcap, h=771 corresponds to minus endcap
	ofstream mean_values;
	mean_values.open("mean_values.txt");
  TString rhidshort;
	//TH1D *new_mean_positive;
	//TH1D *new_mean_negative;
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
    if(h==0) rhidshort = "allgoodchannels_plus";
    if(h==771) rhidshort = "allgoodchannels_minus";

		if(debug_statements) std::cout<<" testing only  good channels plus and minus, h value"<<h<<" rhid values "<<rhidshort<<" endcap "<<endcap<<std::endl;
    //Declare a new histo to store trimmed mean for different values of the variable of interest (pressure, inst L,...)
    TString htrimmeanvsXname = "htrimmean"+filename+title+thevar+"_"+rhidshort;
    TString htrimmeanvsXname_special = "htrimmean"+filename+title+thevar+"_"+rhidshort+"_special";
    TH1D * htrimmeanvsX = new TH1D(htrimmeanvsXname,"", myh->GetNbinsY() , myh->GetYaxis()->GetBinLowEdge(1) , myh->GetYaxis()->GetBinLowEdge( myh->GetNbinsY()+1) );

		TH1D * htrimmeanvsX_special = nullptr; 
	  if(h==0 || h==771){	
		htrimmeanvsX_special = new TH1D(htrimmeanvsXname_special,"", myh->GetNbinsY() , myh->GetYaxis()->GetBinLowEdge(1) , myh->GetYaxis()->GetBinLowEdge( myh->GetNbinsY()+1)) ; }
    /*
		  //Get the rechit ADC charge distribution integrated over all values of the variable of interest
			//    TH1D* projall = myh->ProjectionX("_px",0, myh->GetNbinsY(), h, h) ;  
			//    projall->SetName("charge" +filename+title+"_"+rhidshort );

		// previously we were writing this normalized distribution to the root file also but now just saving them in the folder
		// outf->cd();
		// only good for allgoodchannels and should be done later
    // if(thevar.Index("pressure")>=0 && projall->Integral()>0) projall->Write();
		 TCanvas *canvas_charge_per_chamber = new TCanvas;
       if(thevar.Index("integratelumi_initial")>=0 && projall->Integral()>0 &&savehistos &&rhidshort=="allgoodchannels" ){
  		 canvas_charge_per_chamber->cd();
  	   projall->GetYaxis()->SetTitle("Normalized no. of rechit entries");
  	   projall->GetXaxis()->SetTitle("Raw Charge (ADC counts)");
  	   projall->SetTitle("");
       gStyle->SetOptStat();
  	   projall->DrawNormalized();	 
  		 canvas_charge_per_chamber->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/normalized_charge_"+rhidshort+"_integratelumi_initial.pdf");
			 canvas_charge_per_chamber->Close();
	     } */

    double renormalfactor = 1;    
    double error_renormalfactor =0 ;    
    double error_value=0 ; 
    bool flag= false;		
    //new_mean_positive=nullptr;
		//new_mean_negative=nullptr;
		double gas_gain;
		double gas_gain_error;
		//Loop over the bins of the variable of interest
    //Get the rechit ADC charge distribution for a given bin of the variable of interest
    //for(int j = 1; j <= myh->GetNbinsY(); j++){
		// Lets test for 1 bin
		
/*	  TCanvas *c_individual;
		  TH1D *my_hist = (TH1D*) (myh->ProjectionX("_px",1,1,16,16))->Clone() ;
		  c_individual = new TCanvas("c_individual", " for 1st bin , charge distribution ");
			c_individual->cd();
			gStyle->SetOptStat(111111211);
			my_hist->Draw();
			my_hist->SetTitle("charge : channel "+rhidshort);
			c_individual->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/attempt_first_int_lum_bin_channel_"+rhidshort+".pdf"); */

    for(int j = 1; j <= myh->GetNbinsY(); j++){
      TH1D* proj = (h==0 || h==771)? (TH1D*) (myh->ProjectionX("_px",j,j,1,1))->Clone() : (TH1D*)(myh->ProjectionX("_px",j,j,h,h))->Clone();  
			TString bin_nb = TString::Format("%d",j);

			
			// saving all the distributions 
/*		c_individual = new TCanvas("c_individual", " for a  bin , charge distribution ");
			c_individual->cd();
			gStyle->SetOptStat(111111211);
			proj->Draw();
			TString bin_nb = TString::Format("%d",j);
			proj->SetTitle("charge : "+rhidshort+" bin : "+bin_nb);
			c_individual->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/first_int_lum_bin_channel_"+rhidshort+"_bin_"+bin_nb+".pdf");*/
//		  std::cout<<"Before adding the channels the integral for "<<rhidshort<<" bin "<<bin_nb<<" var "<<thevar<<" integral "<<proj->Integral()<<" entries "<<proj->GetEntries()<<std::endl;
			if(h==0 || h==771) {
			 proj->Reset("ICESM");
			 proj->ResetStats();
		 }	
//			if(debug_statements) std::cout<<"Before adding the channels the integral for allgoodchannels after clearing  for bin "<<j<<" variable "<<thevar<<" integral "<<proj->Integral()<<"mean "<<proj->GetMean()<<std::endl;
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
      //   if(h==0)			 new_mean_positive =new TH1D("new_mean_positive", " Mean of positive endcap",50,100,600);	
		  //		 if(h==771)		 new_mean_negative =new TH1D("new_mean_negative", " Mean of negative endcap",50,100,600);	

          TCanvas *c;
      	  for(int ichan = chan_initial; ichan < chan_final ; ichan++){
       
            if( (1<= ichan && ichan<= 6) || (401<= ichan && ichan<=406))continue;
       	    if(ichan%10>=7 ||ichan%10==0) continue;	  
       	    if(ichan>370 &&ichan<400) continue;
       	    int theendcap = (ichan<=400)? 1:2;
       	    int thechamber = (ichan<=400)?  (int)floor(ichan/10)  : (int)floor( (ichan-400) /10) ;
    				TString rhidshort ="chamber"+ (TString) Form("%d", thechamber)  +"_layer"+ (TString)Form("%d",ichan%10) + "_Endcap"+ theendcap;

            //std::cout<<"just channel "<<ichan<<std::endl; 
       	    //if(isbadchannel("ME11",ichan,_ruNb)  ) continue;
       	    TH1D * h_prov = (TH1D*) ( myh->ProjectionX("_px",j,j, ichan,ichan) )->Clone();
			// saving all the distributions 
/*		  c_individual = new TCanvas("c_individual", " for 1st bin , charge distribution ");
			c_individual->cd();
			gStyle->SetOptStat(111111111);
			//gStyle->SetOptStat("ksiourmen");
			h_prov->Draw();
			TString s = TString::Format("%d",ichan);
			h_prov->SetTitle("charge : "+rhidshort+ " channel : "+s + " correspoinding :  "+h);
			c_individual->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/int_lum_bin_"+rhidshort+"_channel_"+s+"_bin_"+bin_nb+".pdf"); */

       	    TH1D * h_prov_new = (TH1D*) ( myh->ProjectionX("_px",j,j, ichan,ichan) )->Clone();
						h_prov_new->Reset();
						h_prov_new->ResetStats();

						h_prov->Sumw2();
						h_prov_new->Sumw2();
						//h_prov_new->Sumw2();
						// atleast we will add only those channels for which we have more than 50 entries
       	    if(h_prov->Integral()<50) continue;

					// this next command was just for some test , and should not be run in actual program	
						//						if(rhidshort.Contains("chamber1") || rhidshort.Contains("chamber4") ||rhidshort.Contains("chamber5") || rhidshort.Contains("chamber6") || rhidshort.Contains("chamber7") || rhidshort.Contains("chamber8") || rhidshort.Contains("chamber9") || rhidshort.Contains("chamber20") ||rhidshort.Contains("chamber21") || rhidshort.Contains("chamber24") || rhidshort.Contains("chamber25") || rhidshort.Contains("chamber36")  ) { std::cout<<" yes it did come to other points bin "<<j<<std::endl; continue; }
            // testing if we need to store results for just one endcap : endcap =1 -> +endcap
						// if(theendcap==1) continue;					
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
							  h_prov_new->SetBinError(last_bin, h_prov->GetBinError(last_bin));
		
								/*	for(int it=last_bin+1; it<=h_prov->GetNbinsX() ; it++) {
							  	h_prov_new->SetBinContent(it,0);
							  	h_prov_new->SetBinError(it,0);
						   	}*/
								normal = h_prov_new->Integral();

								// new addition for evaluating mean based on individual channel
								if(h==0) {  mean_values<<" mean of the trimmed histogram before : rhid "<<rhidshort<<" mean "<<h_prov->GetMean()<<" error "<<h_prov->GetMeanError()<<"after trimming "<<h_prov_new->GetMean()<<" error "<<h_prov_new->GetMeanError()<<std::endl; }
														//new_mean_positive->Fill(h_prov_new->GetMean(),1./h_prov_new->GetMeanError());}
								if(h==771) {mean_values<<" mean of the trimmed histogram before : rhid "<<rhidshort<<"mean "<<h_prov->GetMean()<<" error "<<h_prov->GetMeanError()<<"after trimming "<<h_prov_new->GetMean()<<" error "<<h_prov_new->GetMeanError()<<std::endl; }
									//new_mean_negative->Fill(h_prov_new->GetMean(),1./h_prov_new->GetMeanError());}
//						   std::cout<<" integral individual channel before trimming  "<<rhidshort<<"for channel "<<ichan<<" : "<<h_prov->Integral()<<" after trimming "<<h_prov_new->Integral()<<std::endl;
							//	if(debug_statements) std::cout<<" integral of projection adding to allgoodchannels "<<rhidshort<<"for channel "<<ichan<<" : "<<normal<<" final events "<<added_events<<std::endl;
					
			   				proj->Add(h_prov_new,1./normal);// to each channel we normalize with respect to total integral acroos the channel , not wrt only trimmed mean integral

            // adding one more if condition to make plots for plus endcap and minus endcap separately 
  		          delete h_prov;
			          delete h_prov_new;
         } // End of loop of channels
		 	}// End of special case
      

		// for h=0 and h==771 , lets find the average mean values in new way
//		TF1 *f1 = new TF1("f1","gaus",100,600);
//		TF1 *f2 = new TF1("f2","gaus",100,600);
		
//		if(h==0) {new_mean_positive->Fit(f1,"R+"); gas_gain=f1->GetParameter(1);	 gas_gain_error=f1->GetParError(1); mean_values<<"average trimmed value in goodchannels plus "<<gas_gain<<"erorr "<<gas_gain_error<<std::endl;
	  // plotting the average plot 
/*		TCanvas *c_1 = new TCanvas("c_1","average ");
	  c_1->cd();
		gStyle->SetOptFit();
	  new_mean_positive->Draw();
		new_mean_positive->SetTitle(" Average charge : "+rhidshort+" : "+thevar+" : bin : "+bin_nb);
*/	
//		c_1->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/average_gas_gain_positive_bin_"+bin_nb+"_"+thevar+".pdf");
		
//		}
//		if(h==771) {new_mean_negative->Fit(f2,"R+");gas_gain=f2->GetParameter(1);	 gas_gain_error=f2->GetParError(1);mean_values<<"average trimmed value in goodchannels minus "<<gas_gain<<"erorr "<<gas_gain_error<<std::endl;
			  // plotting the average plot 
/*		TCanvas *c_1 = new TCanvas("c_1","average ");
	  c_1->cd();
		gStyle->SetOptFit();
		new_mean_negative->SetTitle(" Average charge : "+rhidshort+" : "+thevar+" : bin : "+bin_nb);
	  new_mean_negative->Draw();*/
//		c_1->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/average_gas_gain_negative_bin_"+bin_nb+"_"+thevar+".pdf");
//		}
//	    delete new_mean_positive; 
 //     delete new_mean_negative; 
//
      normal = proj->Integral();
      chargemeantrimm = 0; integral = 0;
      TH1D * h_trim = (TH1D * )proj->Clone();     
			gStyle->SetOptStat(111111211);
      TH1D * h_trim_new = (TH1D*) proj->Clone();

			// before truncating lets plot charge distribution 
      //Previously we were truncating only individual channel not all good channels 
      // truncation for all the channels, along with allgoodchannels
			// I am reseting stats, so new mean and integral is not interferred with old one
	    if(h!=0 && h!=771){//Do the truncation, if h ==0 and h==771 the truncation is already done

			// The next feew stepas are to find till which bin we need to chopp off tail
				  h_trim_new->Reset();
				  h_trim_new->ResetStats();

       	  //if(isbadchannel("ME11",ichan,_runNb)  ) continue;
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
				  h_trim_new->SetBinError(last_bin, h_trim->GetBinError(last_bin));
					for(int it=last_bin+1; it<=proj->GetNbinsX() ; it++) {
						h_trim_new->SetBinContent(it,0);
						h_trim_new->SetBinError(it,0);
					}
					float final_integral = new_integral+entries_last_bin;
				} 


			 gStyle->SetOptStat("kKsSiourRmMen");
			 TString title_for_proj = rhidshort+" : "+thevar+" : bin : "+j;
			 proj->SetTitle(title_for_proj);
//		  if(debug_statements) std::cout<<" rhid channel "<<rhidshort<<" bin number "<<j<<" Integral "<<h_trim->Integral()<<" value of mean after trimming "<<h_trim->GetMean()<<" Error on mean "<<h_trim->GetMeanError()<<std::endl;
/*    	 if(proj->Integral()>0 && (thevar.Index("pressure")>=0 || thevar.Index("integratelumi")>=0) && savehistos && (rhidshort.Contains("chamber1") || rhidshort.Contains("allgoodchannels")) && (j==11 || j==12 || j==8) ){
//    	 if(h_trim->Integral()>0 && (thevar.Index("timesecond")>=0) && savehistos && (rhidshort="chamber1"||rhidshort="chamber1" ||rhidshort="chamber558" ||rhidshort="chamber556") || rhidshort.Contains("allgoodchannels")) && (j==11 || j==12) ){
          std::cout<<"  entred in plots *******************************************************************"<<std::endl;
         
				 TCanvas *canvas_charge = new TCanvas();
 			 	 canvas_charge->cd();
			   gStyle->SetOptStat(111111211);
	    	 proj->Draw("E");
 			 	 canvas_charge->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/charge_distribution"+title+"_"+rhidshort+"vs"+thevar+"_bin"+j+"_channel"+h+"_before_trim.pdf"); 
			 }

    	 if(proj->Integral()>0 && (thevar.Index("pressure")>=0 || thevar.Index("integratelumi") >=0) && savehistos && (rhidshort.Contains("chamber1") || rhidshort.Contains("allgoodchannels")) && (j==11 || j==12 || j==8) ){
    	 //if(h_trim_new->Integral()>0 && (thevar.Index("integratelumi")>=0 || thevar.Index("pressure")>=0) && savehistos ){
			 	 TCanvas *canvas_charge = new TCanvas();
 			 	 canvas_charge->cd();
			   gStyle->SetOptStat(111111211);
	    	 h_trim_new->Draw("E");
 			 	 canvas_charge->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/charge_distribution"+title+"_"+rhidshort+"vs"+thevar+"_bin"+j+"_channel"+h+"_after_trim.pdf"); 
       } 
*/
     // we are fillling entries only if atleast a distribution has more than 50 entries and its not just 1 entry after trimming
		 // Previously we were filling integrated luminosity variable by normalizing with respect to 1st bin, now we are filling the bin with raw value of charge
           /* 	 if(thevar.Index("integratelumi")>=0 ){
       	     //Special case for charge vs lumi: one sets Q(first bin) = 1;   
       	     if( j== 1 ){
       	       htrimmeanvsX->SetBinContent(j, 1 );
       	       renormalfactor = h_trim->GetMean();
       	       htrimmeanvsX->SetBinError(j, h_trim->GetRMS()/sqrt(h_trim->Integral())/renormalfactor);
       	     }
       	     else {
       	       htrimmeanvsX->SetBinContent(j, h_trim->GetMean()/renormalfactor );
       	       htrimmeanvsX->SetBinError(j, h_trim->GetRMS()/sqrt(h_trim->Integral() )/renormalfactor );
       	     }
       	     }   */
       	     //Special case for charge vs lumi: one sets Q(first bin) = 1;   
					   // now we are treating the uncertainty in renomal factor also
					   //
/*       	 if(thevar.Index("integratelumi")>=0 ){
//					  if(debug_statements)  std::cout<<"coming to  all  channels : rhid :var  "<<rhidshort<<" : "<<thevar<<std::endl;
//				  	if(debug_statements)  std::cout<<" particular rhid bin number in integratelumi "<<j<<" entries "<<h_trim->Integral()<<"mean "<<h_trim->GetMean()<<std::endl;

       	    if(j==1){
						// Integral greater than 50 ensures that we have enough entries
	       	    renormalfactor =h_trim_new->GetMean();
//							if(debug_statements) std::cout<<"1st bin int lumi : rhidshort :"<<rhidshort<<" avlue "<<renormalfactor<<std::endl;
							error_renormalfactor = h_trim_new->GetMeanError();
							//error_renormalfactor = h_trim_new->GetMeanError();
 						  error_value = sqrt(pow(h_trim_new->GetMeanError()* 1/renormalfactor ,2) +pow((h_trim_new->GetMean()/(renormalfactor* renormalfactor))*error_renormalfactor,2) );
						 // storing only for all good channels	
							// propogation of errors
//							error_value = sqrt( pow( (h_trim_new->GetRMS()/sqrt(h_trim_new->Integral())) * 1/renormalfactor,2) +pow((h_trim->GetMean()/(renormalfactor* renormalfactor))*error_renormalfactor,2) );
//							flag = true;
//
	       	    htrimmeanvsX->SetBinContent(j, 1 );
	       	    htrimmeanvsX->SetBinError(j,error_value);
	       	    //htrimmeanvsX->SetBinError(j,error_renormalfactor);
						}
       	  else {
//					if(debug_statements) std::cout<<"int luim  bin : "<<j<<" rhidshort :"<<rhidshort<<" mean "<<h_trim->GetMean()<<std::endl;
						// propogation of errors
						error_value = sqrt(pow(h_trim_new->GetMeanError()* 1/renormalfactor ,2) +pow((h_trim_new->GetMean()/(renormalfactor* renormalfactor))*error_renormalfactor,2) );
       	    htrimmeanvsX->SetBinContent(j, h_trim_new->GetMean()/renormalfactor );
       	    htrimmeanvsX->SetBinError(j,error_value);
       	  }
 				} */

       if(proj->Integral() < 50) continue;
			 // question is if already the histogram have more than 50 entries, then how the trimmed histogram can have just one entry. Not possible , but don't know why I applied additional cut. Are these two things not same ?  Integral is basically total bin height * mean value, however , 
  	   if(h_trim_new->Integral() ==1) continue;
       htrimmeanvsX->SetBinContent(j, h_trim_new->GetMean() ); 
     	 htrimmeanvsX->SetBinError(j, h_trim_new->GetMeanError() ) ; 

			 if(h==0 || h==771){
			 htrimmeanvsX_special->SetBinContent(j, gas_gain ); 
     	 htrimmeanvsX_special->SetBinError(j, gas_gain_error ) ; 
			 }

        delete proj;
        delete h_trim_new;
				delete h_trim;
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
    if(thevar.Index("instlumi")>=0 )  htrimmeanvsX->GetXaxis()->SetRangeUser(8000,20000);
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
    if( thevar.Index("pressure")>=0 ) fitlowedge = 950; 
    if( thevar.Index("pressure")>=0 ) fithighedge = 972 ;
    if( thevar.Index("instlumi")>=0 ) fitlowedge = 2000; 
    if( thevar.Index("instlumi")>=0 ) fithighedge =12000 ;

    //Defines the type of function to fit
/*    if(thevar.Index("integratelumi")>=0 || thevar.Index("time")>=0){ 
		  gStyle->SetOptStat();
      htrimmeanvsX->Draw("LEP");
		}
*/
			TF1 *fa1 =nullptr; 
    if(thevar.Index("pressure") >=0) {
			fa1= new TF1("fa1","expo", fitlowedge,fithighedge );
      if(htrimmeanvsX==NULL) continue;
		  htrimmeanvsX->Fit(fa1, "R");  }
    if(thevar.Index("instlumi")>=0){
			fa1= new TF1("fa1","expo", fitlowedge,fithighedge );

      if(htrimmeanvsX==NULL) continue;
		  htrimmeanvsX->Fit(fa1, "R"); }
    //else   fa1= new TF1("f1","pol1", fitlowedge,fithighedge );
    //Now, let's fit
//		htrimmeanvsX->Fit(fa1, "R"); 

    gStyle->SetOptStat("001111111");
    c3->SetName("c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar);
			if(htrimmeanvsX->Integral()>0) {
				htrimmeanvsX->Write();
			}
      htrimmeanvsX->Draw();

			// Jus tthe same thing for a special plot with mean taken from gaussian distribution 
/*		if(h==0 || h==771){
    c3 = new TCanvas;
    htrimmeanvsX_special->SetTitle(chamber_string_name+" : "+filename+title+"_"+rhidshort+"_special");
    if(thevar.Index("pressure")>=0 )   htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("integratelumi")>=0 )   htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("time")>=0 )   htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 )   htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 && rhidshort =="allgoodchannels_plus" )  htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 && rhidshort =="allgoodchannels_minus" )  htrimmeanvsX_special->GetYaxis()->SetRangeUser(150,500);
    if(thevar.Index("instlumi")>=0 )  htrimmeanvsX_special->GetXaxis()->SetRangeUser(8000,20000);
    if(thevar.Index("time")>=0 ){
		 htrimmeanvsX_special->GetXaxis()->SetTimeDisplay(1);
		 htrimmeanvsX_special->GetXaxis()->SetLabelSize(0.02);
		 htrimmeanvsX_special->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    }
    htrimmeanvsX_special->GetXaxis()->SetTitle(xtitle);
		htrimmeanvsX_special->GetYaxis()->SetTitle("Trimmed mean charge");
		if(thevar.Index("integratelumi")>=0 ) htrimmeanvsX_special->GetYaxis()->SetTitle("Trimmed mean charge");
  
		htrimmeanvsX_special->SetMarkerStyle(20); htrimmeanvsX_special->SetMarkerSize(0.7);
    htrimmeanvsX_special->SetName(filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+"_special");
    gStyle->SetOptStat("001111111");
		//Defines the range for the fit,    
     fitlowedge = 0, fithighedge = 44;
    if( thevar.Index("pressure")>=0 ) fitlowedge = 948; 
    if( thevar.Index("pressure")>=0 ) fithighedge = 981 ;
    if( thevar.Index("instlumi")>=0 ) fitlowedge = 10000; 
    if( thevar.Index("instlumi")>=0 ) fithighedge =19000 ;

    //Defines the type of function to fit
    if(thevar.Index("integratelumi")>=0 || thevar.Index("time")>=0){ 
		  gStyle->SetOptStat();
      htrimmeanvsX_special->Draw("LEP");
		}

    if(thevar.Index("pressure") >=0) {
			fa1= new TF1("fa1","expo", fitlowedge,fithighedge );
		  htrimmeanvsX_special->Fit(fa1, "R");  }
    if(thevar.Index("instlumi")>=0){
			fa1= new TF1("fa1","expo", fitlowedge,fithighedge );
		  htrimmeanvsX_special->Fit(fa1, "R"); }

    gStyle->SetOptStat("001111111");
    c3->SetName("c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+"_special");
		if(htrimmeanvsX_special->Integral()>0) htrimmeanvsX_special->Write();
			c3->cd();
      htrimmeanvsX_special->Draw();
    } */


			// this is to save the plots for one chamber
			//if(htrimmeanvsX->Integral()>0 && (thevar.Index("integratelumi")>=0 || thevar.Index("pressure") >=0 || thevar.Index("instlumi")>=0) && savehistos && (rhidshort=="allgoodchannels" || rhidshort=="chamber15_layer3_Endcap2")){c3->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
//			if(htrimmeanvsX->Integral()>0 && (thevar.Index("integratelumi")>=0 || thevar.Index("pressure") >=0 || thevar.Index("instlumi")>=0) && (rhidshort =="allgoodchannels" || rhidshort=="chamber36_layer6_Endcap2" || rhidshort =="chamber9_layer1_Endcap1") ){c3->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
			//if(htrimmeanvsX->Integral()>0 && (thevar.Index("integratelumi")>=0 || thevar.Index("pressure") >=0 || thevar.Index("instlumi")>=0) && savehistos){c3->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }

			 /*  
				 if(thevar.Index("integratelumi") >=0){
				 std::pair <double, double> parampair_1 ;
				 parampair_1.first =  0 ;
				 parampair_1.second = 0;
				 result[h] = parampair_1;
				}
			 */

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
				//std::cout<<" issue while savig the output text file :  "<<thevar<<std::endl;
        if(h==0){	
      	 //outfile<<"charge\tttime_initial\ttime_final for all good channels positive \n"<<std::endl;
      	 double charge_value = 0 ;
      	 for(int i=1 ; i<=htrimmeanvsX->GetNbinsX(); i++) { 
           charge_value = htrimmeanvsX->GetBinContent(i);
					 if(charge_value==0) continue;
            std::cout << std::fixed << std::setprecision(0);
      	   //outfile<<charge_value<<"\t"<<std::setprecision(10)<<htrimmeanvsX->GetXaxis()->GetBinLowEdge(i)<<"\t"<<htrimmeanvsX->GetXaxis()->GetBinLowEdge(i+1)<<std::endl;
      	 }

        }
       if(h==771){	
      	 //outfile<<"charge\tttime_initial\ttime_final for all good channels negative \n"<<std::endl;
      	 double charge_value = 0 ;
      	 for(int i=1 ; i<=htrimmeanvsX->GetNbinsX(); i++) { 
           charge_value = htrimmeanvsX->GetBinContent(i);
					 if(charge_value==0) continue;
            std::cout << std::fixed << std::setprecision(0);
      	   //outfile<<charge_value<<"\t"<<std::setprecision(10)<<htrimmeanvsX->GetXaxis()->GetBinLowEdge(i)<<"\t"<<htrimmeanvsX->GetXaxis()->GetBinLowEdge(i+1)<<std::endl;
	       }
        }
			 //outfile.close();

		//	std::cout<<" issue before ending the loop :  "<<thevar<<std::endl;
  	delete htrimmeanvsX;
  	delete htrimmeanvsX_special;
      delete c3;

		
			if(thevar.Index("instlumi")>=0 ||thevar.Index("pressure")>=0   ){ 
			double theslope = fa1->GetParameter(1);
			double unc_slope = fa1->GetParError(1);
			double chi2overN = (fa1->GetNDF () >0)? (fa1->GetChisquare ()/fa1->GetNDF () ) : -10 ; 
			
/*			if( thevar.Index("instlumi") >=0&& theslope>0.00001 )    {
     if( thevar.Index("instlumi") >=0 ) {std::cout<<" **************************************tifying if outlier exists "<<thevar<<" slopes "<<theslope<<" rhid "<<rhidshort<<std::endl;}
		c3->SaveAs(output_plots_folder+"plotfolder_"+chamber_string_name+"/outliers_trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); } */

		if(h==0){
   		_slope_goodchannels_plus=theslope;
   		_slope_error_goodchannels_plus=unc_slope;
   		theouttree_goodchannels->Fill();
		}

		if(h==771){
   		_slope_goodchannels_minus=theslope;
   		_slope_error_goodchannels_minus=unc_slope;
   		theouttree_goodchannels->Fill();
		}
    
   if( theslope  !=0 && h!=0 && h!=771){
      h_slope->Fill(theslope);
      h_slopeuncty->Fill(unc_slope);
      h_chi2->Fill(chi2overN) ;
   }

    if( theslope  !=0 && h!=0 && h!=771 ){
      int chambandlay = (h>400)? (h-400) : h; 
      _slope =theslope;
			//_slope_error =fa1->GetParError(1);
      _slope_error =unc_slope;
      _chi2 =fa1->GetChisquare () ;
      _ndof = fa1->GetNDF ();
      _chamber = (int) floor( chambandlay/10 ) ; 
      _layer = h%10 ;
      _endcap = (int)floor( h/400 ) +1 ;
      _stationringHVseg =  thestationringhv(filename);
      _isbadchannel = false; // isbadchannel("ME11",h,_runNb);
      theouttree->Fill();
     } 
		}//end of if condition for  writing tree for instlumi and intlumi parameters  
  

	} //end of loop for all the channelss 
  
  if( thevar.Index("instlumi")>=0  ) h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs instant Luminosity (/#mu b^{-1} s^{-1})");
  if( thevar.Index("time")>=0  ) h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs time");
  if( thevar.Index("integrate")>=0  ) h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs #int L (/fb^{-1})");
  if( thevar.Index("pressure")>=0  )  h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs P (hPa^{-1})");

//	if(thevar.Index("integratelumi_initial")<0)  {
//	outf = TFile::Open();
//	outf = TFile::Open();
	std::cout<<" error in loop filling root file "<<std::endl; 
/*	outf->cd();
  h_slope->Write();
  h_slopeuncty->Write();
  h_chi2->Write();
  theouttree->Write();
  theouttree_goodchannels->Write(); */
//  outf->Close();

 	delete h_slope;
	delete h_chi2;
	delete h_slopeuncty;

	std::cout<<" done with one type of variable "<<thevar<<std::endl;
  return result;    

} //end of GetSlope function

double pressure_dependence_removal::ApplyCorrection( double X ,TString correctiontype, double p0, double p1 ){
  double refvalue = 0;
  if(correctiontype =="pressure"&&!droppressurecorr){
    refvalue =960;
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
