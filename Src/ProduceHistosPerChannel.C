#define ProduceHistosPerChannel_cxx
#include <iostream>
#include <stdio.h>

#include <fstream>
#include <iomanip>
#include <string>

#include <sstream>


#include "ProduceHistosPerChannel.h"
#include "badchannel.h"
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

double trimmean =0.7;//Offset used in the calculation of the trimmed mean
bool istest =false; //For tests, run only on 1/100 of the events
bool debug_trimmed_plots = false;
bool savehistos =true;
bool seconditer_plots =true;// If you want to perform a second iteration of the fit to pressure and inst L.
bool seconditer =false;// If you want to perform a second iteration of the fit to pressure and inst L.

bool corrfrom15to23fbonly = true; // Extract pressure dependency from the range 15 to 20 /fb, 7e33 to 9e33, extract inst L. from the range 15 to 20 /fb 

bool check7to9e33only = false; // Only study data in the 7e33 to 9e33 lumi range
bool dropinstlumicorr =false; // Drop inst L correction
bool droppressurecorr =false;// Drop pressure correction

bool debug_print = false;
bool iszmumu = false;//If only Zmumu events are used the corrections cannot be derived channel per channel (not enough stat) so we take all channels together
TH3D * hchargevspressure_2, * hchargevsinstlumi_2;
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
void ProduceHistosPerChannel::Loop(TString detregionstr, TString chamber_string)
{
  chamber_string_name = chamber_string;
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(111);
  if(detregionstr.Index("ZMuMu")>=0) iszmumu=true;
  iszmumu = false;
  TH1F * chargepresscorrhighinstlumi = new TH1F("chargepresscorrhighinstlumi","",10000,0,10000);
  TH1F * chargepresscorrlowinstlumi = new TH1F("chargepresscorrlowinstlumi","",10000,0,10000);


 
   //TString inputfname = "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_newfiles/2017_again/final_files/2017/input_root_files/"+ detregionstr+".root" ;
   TString inputfname = "/cmsuf/data/store/user/t2/users/neha.rawal/CSC_2017_data_RAW_RECO/nTuple_output/2017/remove_HV_trips/2017_delievered_files_final/output_"+ chamber_string +".root" ;
   TFile * inputf = new TFile(inputfname,"open");
   if(istest)detregionstr =detregionstr+"test"; 
   if(seconditer)detregionstr =detregionstr+"seconditer"; 
   if(corrfrom15to23fbonly)detregionstr =detregionstr+"15to23fbonly_7e33to10e33"; 
   if(dropinstlumicorr)detregionstr+="noinstlumicorr";
   if(droppressurecorr)detregionstr+="nopressurecorr";
   if(check7to9e33only) detregionstr+="check7to9e33only";
   //TFile * outf  = new TFile("/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_newfiles/2017_again/final_files/2017/output_root_files/outf"+detregionstr+"_third.root","recreate");  
   TFile * outf  = new TFile("/cmsuf/data/store/user/t2/users/neha.rawal/CSC_2017_data_RAW_RECO/nTuple_output/2017/output_ntuples/output_files_removed_trips_again/outf"+detregionstr+"_output.root","recreate");  
   TTree * tree =(TTree*) inputf->Get("tree");
   Init(tree);
   cout << tree->GetEntries()<<endl;
 
	 // Making IntegratedLumi vs gas gain slope dependency before starting any pressure and inst lumi corrections
   TH3D * hchargevsintegratelumi_initial = new TH3D("hchargevsintegratelumi_initial","charge (ADC counts) vs integ lumi (initial)",3000,0,3000, 46, 0,46  ,770,1,771);

   for(int i = 0; i < tree->GetEntries(); i++){
     // std::cout<<"entered into integratedlumi loop "<<i<<std::endl; 
		 //if(_integratelumi==0) continue;
     // std::cout<<"case for which integratedlumi exists "<<i<<std::endl; 
     if(i%100 !=0 &&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     if(check7to9e33only  &&_instlumi>9000 ) continue;
     if(check7to9e33only  &&_instlumi<7000 ) continue;
     // what is significance of this ?
		 if(check7to9e33only) _integratelumi -= 17000;
     
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;

     double charge  = _rhsumQ_RAW;

     hchargevsintegratelumi_initial->Fill(charge, _integratelumi/1000, rhidreduced);
     
   }
   
   
   vector< std::pair<double, double > > params_integratelumi_initial;
   params_integratelumi_initial = GetSlope( hchargevsintegratelumi_initial, "_integratelumi_initial", detregionstr,"",outf);
 
   //Run on all events to extract pressure correction for each channel (rechit) separately. 
   //This is done with the help of a 3d histogram, storing charge, pressure and local rechit ID (inside a given station/ring/HV chamber)
   TH3D * hchargevspressure = new TH3D("hchargevspressure","charge (ADC counts) vs pressure",300,0,3000, 40, 940,980  ,770,1,771);

   for(int i = 0; i < tree->GetEntries(); i++){
     
     if(droppressurecorr)break;
     if(i%100 !=0&&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl; if(badrun(_runNb) ) continue;
//		 cout<<"entering to integrate lumi loop"<<endl;
//		 cout<<"integrate lumi : "<<_integratelumi<<"\t"<<"instant lumi : "<<_instlumi<<endl;

		 if(corrfrom15to23fbonly &&_integratelumi>23000 ) continue;
     if(corrfrom15to23fbonly &&_integratelumi<15000) continue;
		 
//		 cout<<"entering to inst lumi loop"<<endl;
     if(corrfrom15to23fbonly &&_instlumi<7000 ) continue;
     if(corrfrom15to23fbonly &&_instlumi>10000) continue;
     if(check7to9e33only  &&_instlumi<7000 ) continue;
     if(check7to9e33only  &&_instlumi>10000 ) continue;
     
		 //cout<<"exiting to integrate lumi loop"<<endl;
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;//First (second) endcap have rechit ID < (>) 2000000
     //     Reduced rechit ID has the following format: (A+B)*10+C, where A=1,..36 (chamber nb), B =0 (endcap 1) or 40 (endcap 2)  and C =1,..., 6 (layer).
 //    std::cout<<" pressure information gettig filled"<<std::endl;
     hchargevspressure->Fill(_rhsumQ_RAW, _pressure , rhidreduced);
   }

   vector< std::pair<double, double > > params_pressure;
   params_pressure = GetSlope( hchargevspressure, "_pressure", detregionstr,"",outf);

/*   	 if(debug_print){
       std::cout<<" slope vector for pressure correction"<<std::endl; 
			 for(auto it = std::begin(params_pressure) ; it != std::end(params_pressure); ++it){
				 std::cout<<*it<<"  ";
	
			}
			 std::cout<<std::endl;
  	 }*/
   // The function on the line above fits the trim mean charge vs pressure for each rechit and returns the fitted parameters.   



   //Now, run on all events, apply pressure correction and extract inst L correction

   TH3D * hchargevsinstlumi = new TH3D("hchargevsinstlumi","charge (ADC counts) vs inst lumi",300,0,3000, 50, 0,15000  ,770,1,771);
   int prevrunnb(0), prevlumiblock(0);
   
   for(int i = 0; i < tree->GetEntries(); i++){
     if(dropinstlumicorr)break;
     if(i%100 !=0 &&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     if(corrfrom15to23fbonly &&_integratelumi>23000) continue;
     if(corrfrom15to23fbonly &&_integratelumi<15000) continue;
 
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
//     std::cout<<" pressure correction getting applied "<<std::endl;
     
     int idforcorr = (iszmumu )? 0: rhidreduced ;
     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  ;
     hchargevsinstlumi->Fill( charge, _instlumi, rhidreduced);
 

     if(debug_print)		std::cout<<"charge after pressure correction "<<charge<<std::endl;

		 bool badchannel = isbadchannel(detregionstr , rhidreduced);
     if(!badchannel&&_instlumi>8000&&_instlumi<10000&&_integratelumi>12000&&_integratelumi<16000) chargepresscorrhighinstlumi->Fill(charge); if(!badchannel&&_instlumi>3000&&_instlumi<5000&&_integratelumi>12000 &&_integratelumi<16000 ) chargepresscorrlowinstlumi->Fill(charge);
     
   }

   vector< std::pair<double, double > > params_instlumi;
   params_instlumi = GetSlope( hchargevsinstlumi, "_instlumi", detregionstr,"",outf);
   


 std::cout<<"issue after 1st set of iteration "<<std::endl;
   //Second iteration (optional): same game

   vector< std::pair<double, double > > params_pressure_2;
   vector< std::pair<double, double > > params_instlumi_2;
  
   if(seconditer_plots){
      hchargevspressure_2 = new TH3D("hchargevspressure_2","charge (ADC counts) vs pressure (2nd iteration)",3000,0,3000, 40, 940,980  ,770,1,771);

   //Pressure correction
   for(int i = 0; i < tree->GetEntries(); i++){

     if(i%100 !=0&&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     if(corrfrom15to23fbonly &&_integratelumi>23000 ) continue;
     if(corrfrom15to23fbonly &&_integratelumi<15000) continue;
     if(corrfrom15to23fbonly &&_instlumi<7000 ) continue;
     if(corrfrom15to23fbonly &&_instlumi>10000) continue;

     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;
     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second );
	 //	 * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
     hchargevspressure_2->Fill(charge, _pressure , rhidreduced);
     
   }
   params_pressure_2 = GetSlope( hchargevspressure_2, "_pressure_2", detregionstr,"",outf);
   
   
 std::cout<<"issue before inst luminoisty 2 "<<std::endl;
   //Now charge after pressure correction
   
    hchargevsinstlumi_2 = new TH3D("hchargevsinstlumi_2","charge (ADC counts) vs inst lumi",300,0,3000, 50, 0,15000  ,770,1,771);
   
   int prevrunnb(0), prevlumiblock(0);
   for(int i = 0; i < tree->GetEntries(); i++){
     
     if(i%100 !=0 &&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     if(corrfrom15to23fbonly &&_integratelumi<15000) continue;
     if(corrfrom15to23fbonly &&_integratelumi>23000) continue;
     
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;
     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) * ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second )  ;
     hchargevsinstlumi_2->Fill( charge, _instlumi, rhidreduced);
     
   }
   params_instlumi_2 = GetSlope( hchargevsinstlumi_2, "_instlumi_2", detregionstr,"",outf);
   
   }
   
   
 std::cout<<"issue before integrated luminoisty "<<std::endl;
   
   //Now final charge after pressure, instlumi correction
   prevrunnb=0;
   prevlumiblock=0;
   TH3D * hchargevsintegratelumi = new TH3D("hchargevsintegratelumi","charge (ADC counts) vs integ lumi",3000,0,3000, 46, 0,46  ,770,1,771);

   for(int i = 0; i < tree->GetEntries(); i++){
     // std::cout<<"entered into integratedlumi loop correction "<<i<<std::endl; 
		 //if(_integratelumi==0) continue;
     // std::cout<<"case for which integratedlumi correction exists "<<i<<std::endl; 


     if(i%100 !=0 &&istest) continue;
     LoadTree(i);tree->GetEntry(i);
     if(i%1000000 ==0)cout << i<<endl;
     if(badrun(_runNb) ) continue;
     if(check7to9e33only  &&_instlumi<6000 ) continue;
     if(check7to9e33only  &&_instlumi>10000 ) continue;
     // what is significance of this ?
		 if(check7to9e33only) _integratelumi -= 17000;
     
     int rhidreduced = ((int)floor(_rhid/10))%1000;
     if(_rhid> 2000000) rhidreduced +=400;
     int idforcorr = (iszmumu )? 0: rhidreduced ;

     //double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
     double charge  = _rhsumQ_RAW* ApplyCorrection( _pressure ,"pressure",   (params_pressure[idforcorr]).first , (params_pressure[idforcorr]).second );
			 //* ApplyCorrection(_instlumi, "instlumi", (params_instlumi[idforcorr]).first , (params_instlumi[idforcorr]).second ) ;
    // we need to compute the integratelumi again, since we do not need the iterations till second histograms, 
		// We can make SaveHistos false, but we can just save the integratelumi plots and write all the trees 
//		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second )  * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;
		 if(seconditer) charge = charge *  ApplyCorrection( _pressure ,"pressure",   (params_pressure_2[idforcorr]).first , (params_pressure_2[idforcorr]).second );
	 //	 * ApplyCorrection(_instlumi, "instlumi", (params_instlumi_2[idforcorr]).first , (params_instlumi_2[idforcorr]).second ) ;

     hchargevsintegratelumi->Fill(charge, _integratelumi/1000, rhidreduced);
     
   }
   
   
   vector< std::pair<double, double > > params_integratelumi;
   params_integratelumi = GetSlope( hchargevsintegratelumi, "_integratelumi", detregionstr,"",outf);
   

  hchargevsintegratelumi_initial->Write();
  hchargevspressure_2->Write();
  hchargevspressure->Write();
  hchargevsinstlumi->Write();
  hchargevsintegratelumi->Write();
  hchargevsinstlumi_2->Write();

  //// chargepresscorrhighinstlumi->Write();
  // chargepresscorrlowinstlumi->Write();

   outf->Close();
 std::cout<<"issue realized after coming to end"<<std::endl;
   
}

  

vector < std::pair<double, double > >  ProduceHistosPerChannel::GetSlope( TH3D * myh , TString thevar , TString filename, TString title, TFile * outf){

  //Will store all the fit results in a TTree (one entry per channel)
  TString treename = "tree_"+thevar ; 
  TTree * theouttree = new TTree(treename,"");
  Float_t _slope(-1),_slope_error(-1),_chi2(-1);
  Int_t _ndof(-1),_layer(-1),_chamber(-1),_endcap(-1),_stationringHVseg(-1); 
  Bool_t _isbadchannel(false);
  theouttree->Branch("_slope",&_slope,"_slope/F");
  theouttree->Branch("_slope_error",&_slope_error,"_slope_error/F");
  theouttree->Branch("_chi2",&_chi2 ,"_chi2/F");
  theouttree->Branch("_ndof",&_ndof,"_ndof/I");
  theouttree->Branch("_layer",&_layer,"_layer/I");
  theouttree->Branch("_chamber",&_chamber,"_chamber/I");
  theouttree->Branch("_endcap",&_endcap,"_endcap/I");
  theouttree->Branch("_stationringHVseg",&_stationringHVseg,"_stationringHVseg/I");

  //Will also print the slope values in a text files
//  ofstream textfile;
//  textfile.open ("slopes_"+filename+"_"+thevar+".txt");

 
  TString xtitle ;
  if(thevar.Index("pressure")>=0) xtitle = "Pressure (hPa)";
  if(thevar.Index("instlumi")>=0) xtitle = "Inst Lumi (#mub s)^{-1}";
  if(thevar.Index("integrate")>=0) xtitle = "Integrated luminosity (fb^{-1})";


  TH1D* projglobal = (TH1D*) myh-> ProjectionX("_px",0, myh->GetNbinsY());

  projglobal->Sumw2();
  TCanvas * c = new TCanvas;
  c->cd();
  projglobal->SetTitle(filename+title);
  projglobal->SetMarkerStyle(20);
  projglobal->SetMarkerSize(0.7);
  projglobal->GetXaxis()->SetTitle("Charge");
  projglobal->GetYaxis()->SetTitle("Normalized nb of entries");
  projglobal->GetYaxis()->SetTitleSize(0.05);
  projglobal->DrawNormalized("lep");
  projglobal->SetName("charge"+filename+title);
/*  if(thevar.Index("pressure")>=0 && projglobal->Integral()>0 &&savehistos){ 
    //projglobal->SaveAs("plotfolder/charge_"+filename+".pdf");
   c->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/charge_"+thevar+"_"+filename+".pdf");
  }
  if(thevar.Index("integrate")>=0 && projglobal->Integral()>0 &&savehistos){ 
    //projglobal->SaveAs("plotfolder/charge_"+filename+".pdf");
    c->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/charge_"+thevar+"_"+filename+".pdf");
  }
*/
    outf->cd();
		dir_var_name = outf->mkdir(thevar);
		dir_var_name->cd();
		//std::cout<<dir_var_name<<std::endl;



  vector < std::pair<double, double >  > result ; //Assume that fitted function has two parameters
  
  for(int i = 0; i<771;i++){
    std::pair<double, double > initpair(0,0); 
    result.push_back(initpair);
  }

  if(thevar.Index("pressure")>=0 &&droppressurecorr ) return result; 
  if(thevar.Index("instlumi")>=0 &&dropinstlumicorr ) return result; 
  

  // this lowedge and highedge are for integrated luminosity 
  double lowedge =  -0.005; 
  double highedge  = 0.005;
  
  if( thevar.Index("instlumi")>=0 )  lowedge = -0.00003;
  if( thevar.Index("instlumi")>=0 )  highedge = 0.00003;
  if( thevar.Index("instlumi_2")>=0 )  lowedge = -0.00001;
  if( thevar.Index("instlumi_2")>=0 )  highedge = 0.00001;
 
	if( thevar.Index("pressure")>=0 )  lowedge = -0.02;
  if( thevar.Index("pressure")>=0 )  highedge = 0.01;
  
 // I am  trying to make plots for all layers together for a chamber  
//  TH1D * h_trim_layers = (TH1D * )proj->Clone();     
  TH1D * h_slope = new TH1D ("h_slope"+thevar+"_"+filename,"",200, lowedge, highedge );
  TH1D * h_slopeuncty = new TH1D ("h_slopeuncty"+thevar+"_"+filename,"",200, lowedge/100., highedge/100. );
  TH1D * h_chi2 = new TH1D ("h_chi2"+thevar+"_"+filename,"",200,0,1000);


  for(int h = 0; h <771 ; h++){//Loop over reduced rechit ID
    //Reduced rechit ID has the following format: (A+B)*10+C, where A=1,..36 (chamber nb), B =0 (endcap 1) or 40 (endcap 2)  and C =1,..., 6 (layer). 
    //Skipping empty entries
    if(h!=0 && h%10>=7)continue;
	 	// h >7 implies that the one which are extra layers after reduced rechit Id = 366 for chamber36, endcap1, layer6 :
		// and reducedrechitId = 411 for chamber 1 , endcap 2, layer 1  , we dont want to count them 
    if(h!=0 && h%10==0)continue;
    if(h>370&& h<=400)continue;
    //N.B.: h=0 takes all rechits together 
    
		TString endcap = (h<=400)? "_Endcap1":"_Endcap2";
    int chambernb = (h<=400)?  (int)floor(h/10)  : (int)floor( (h-400) /10) ;  // this condition will remove if there are extra chambers between 400 and 410 

		// layer all together 
		/*if(i==12) i=0 
		 * while(i!=12) {
		 * add all the charge 
      TH1D * htrimmeanvsX = new TH1D(htrimmeanvsXname,"", myh->GetNbinsY() , myh->GetYaxis()->GetBinLowEdge(1) , myh->GetYaxis()->GetBinLowEdge( myh->GetNbinsY()+1)  );
		 * i=i+1
		 * }
		 *
		 * for(int layer_nb = 0 ; layer_nb <6; layer_nb++){

		} */
    TString rhidshort ="chamber"+ (TString) Form("%d", chambernb)  +"_layer"+ (TString)Form("%d",h%10) + endcap;
    if(h==0) rhidshort = "allgoodchannels";

    //Declare histo to store trimmed mean for different values of the variable of interest (pressure, inst L,...)
    TString htrimmeanvsXname = "htrimmean"+filename+title+thevar+"_"+rhidshort;
    TH1D * htrimmeanvsX = new TH1D(htrimmeanvsXname,"", myh->GetNbinsY() , myh->GetYaxis()->GetBinLowEdge(1) , myh->GetYaxis()->GetBinLowEdge( myh->GetNbinsY()+1)  );

		TCanvas *canvas_charge_per_chamber = new TCanvas;
    //Get the rechit ADC charge distribution integrated over all values of the variable of interest
    TH1D* projall = myh->ProjectionX("_px",0, myh->GetNbinsY() ,h,h) ;  
    projall->SetName("charge" +filename+title+"_"+rhidshort );

		// previously we were writing this normalized distribution to the root file also but now just saving them in the folder
		// outf->cd();
    // if(thevar.Index("pressure")>=0 && projall->Integral()>0) projall->Write();
    // myh is a 3D histogram 
    if(thevar.Index("integratelumi_initial")>=0 && projall->Integral()>0 &&savehistos &&rhidshort=="allgoodchannels"){
  		  canvas_charge_per_chamber->cd();
  	  	projall->GetYaxis()->SetTitle("Normalized no. of rechit entries");
  	  	projall->GetXaxis()->SetTitle("Raw Charge (ADC counts)");
  	  	projall->SetTitle("");
  	    projall->DrawNormalized();	 
  			canvas_charge_per_chamber->SaveAs(output_path+"plotfolders_removed_trips_again/plotfolder_"+chamber_string_name+"/normalized_charge_"+rhidshort+"_integratelumi_initial.pdf");
				//Only write the histo if the variable of interest is pressure i.e. before any correction is applied the normalized charge distribution for each bin
	  	}
 		// saving the normalized distribution for each chamber and layer for all the pressure values
/* 		if(thevar.Index("pressure")>=0 && projall->Integral()>0 &&savehistos  && (rhidshort=="chamber10_layer2_Endcap1" || rhidshort=="chamber35_layer4_Endcap1") ){
		  canvas_charge_per_chamber->cd();
	  	projall->GetYaxis()->SetTitle("Normalized no. of entries");
	  	projall->GetXaxis()->SetTitle("Raw Charge (ADC counts)");
	    projall->DrawNormalized();	 
			canvas_charge_per_chamber->SaveAs(output_path+"plotfolder/plotfolders_"+chamber_string_name+"/normalized_charge_"+rhidshort+".pdf");//Only write the histo if the variable of interest is pressure i.e. before any correction is applied the normalized charge distribution for each bin
	  	}
  		if(thevar.Index("integrate")>=0 && projall->Integral()>0 &&savehistos  && (rhidshort=="chamber10_layer2_Endcap1" || rhidshort=="chamber35_layer4_Endcap1")){
		  canvas_charge_per_chamber->cd();
	  	projall->GetYaxis()->SetTitle("Normalized no. of entries");
	  	projall->GetXaxis()->SetTitle("Raw Charge (ADC counts)");
	    projall->DrawNormalized();	 
			canvas_charge_per_chamber->SaveAs(output_path+"plotfolder/plotfolders_"+chamber_string_name+"/normalized_charge_"+rhidshort+"_integratelumi.pdf");//Only write the histo if the variable of interest is pressure i.e. before any correction is applied the normalized charge distribution for each bin
	  	}
 */   
    double renormalfactor = 1;    

    for(int j = 1; j <= myh->GetNbinsY(); j++){//Loop over the bins of the variable of interest
    
      //Get the rechit ADC charge distribution for a given bin of the variable of interest
      TH1D* proj = (h==0)? (TH1D*) (myh->ProjectionX("_px",j,j,0,0))->Clone()  : (TH1D*)(myh->ProjectionX("_px",j,j,h,h))->Clone();  
      double normal = 0;

      double chargemeantrimm = 0; double  integral = 0;
      if(h == 0){//Special cases: all good channels in a single histo

          TCanvas *c;
       	 for(int ichan = 1; ichan < 771 ; ichan++){
       
       	    if(ichan%10>=7 ||ichan%10==0) continue;	  
       	    if(ichan>370 &&ichan<400) continue;
       	    int theendcap = (ichan<400)? 1:2;
       	    int thechamber = (theendcap==1)?  (int)floor(ichan/10)  : (int)floor( (ichan-400) /10) ;
       	    if(isbadchannel(filename,ichan)  ) continue;
       
       	    TH1F * h_prov = (TH1F*)( myh->ProjectionX("_px",j,j, ichan,ichan)  )->Clone();
       
       	    if(h_prov->GetEntries() <=100){
       	      //cout <<"Skipping, small nb of entries: " <<  h_prov->GetEntries() <<endl; 
       	      //cout <<"ichan: "<< ichan<<", "<<endl;
       	      continue;
       	    }
	
		        // saving the trimmed distributions to the folder for some of the chambers, it will be trimmed plots for every pressure bin
            // plots for allgoodchannels (h=0)	
		        //	if(thevar.Index("pressure") > 0  || thevar.Index("instlumi") >0) std::cout<<"rhidshort : "<<rhidshort<<std::endl;
  
      	    integral = 0;
      	    normal = h_prov->Integral();
      	    for(int it = 1; it<=   h_prov->GetNbinsX() ;it++) {
      	      if(integral < trimmean * normal){ integral+= h_prov->GetBinContent(it);}
      	      else h_prov->SetBinContent(it,0); // Builds a truncated ("trimmed") histo
	          }

						std::cout<<"integral for weight "<<integral<<std::endl;
						std::cout<<"mean of hist of individual channel "<<h_prov->GetMean()<<std::endl;
					c = new TCanvas("c", "canvas charge individual all good channels", 600,800);
      		c->cd();
      	  h_prov->GetXaxis()->SetTitle("Charge");
      		h_prov->GetYaxis()->SetTitle("Entries");
      		h_prov->Draw();
      		c->SaveAs(output_path+"plotfolders_removed_trips_again/plotfolder_"+chamber_string_name+"_new/trimmed_mean_integratelumi_allgoodchannels_aftertrim.pdf");
          c->Clear();
	          proj->Add( h_prov,1./integral );//Each channel is added to the total with an equal weight. 
	          delete h_prov;
         }

				 //plotting charge distribution y for all good channesl
/*    		 if(thevar.Index("integratelumi_initial") > 0 ) {
      		TCanvas * c = new TCanvas("c", "canvas chargemean for instlumi ", 600,800);
      		c->cd();
      	  proj->GetXaxis()->SetTitle("Charge");
      		proj->GetYaxis()->SetTitle("Entries");
      		proj->Draw();
      		c->SaveAs(output_path+"plotfolders_removed_trips/plotfolder_"+chamber_string_name+"/trimmed_mean_integratelumi_allgoodchannels_aftertrim.pdf");
	       } */
      }// End of special case

      normal = proj->Integral();
      chargemeantrimm = 0; integral = 0;
      TH1D * h_trim = (TH1D * )proj->Clone();     
	
		  // trimmed mean plot for few chambers for all the pressure bins	
    	/*	if(thevar.Index("pressure") > 0 && (rhidshort=="chamber10_layer2_Endcap1" || rhidshort =="chamber35_layer4_Endcap")) {
      		TCanvas * c = new TCanvas("c", "canvas for trimmed mean for pressure ", 600,800);
	      	c->cd();
	        h_trim->GetXaxis()->SetTitle("Charge");
	      	h_trim->GetYaxis()->SetTitle("Entries");
	      	h_trim->Draw();
	      	TString stri;
	      	stri.Form("%d",j);
	      	c->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/trimmed_mean_pressure_"+stri+"_"+rhidshort+"_beforetrim.pdf");
	       }
      	 if(thevar.Index("instlumi") > 0 && (rhidshort=="chamber10_layer2_Endcap1"|| rhidshort =="chamber35_layer4_Endcap1")) {
      		TCanvas * c = new TCanvas("c", "canvas for trimmed mean for instlumi ", 600,800);
      		c->cd();
      	  h_trim->GetXaxis()->SetTitle("Charge");
      		h_trim->GetYaxis()->SetTitle("Entries");
      		h_trim->Draw();
      		TString stri;
      		stri.Form("%d",j);
      		c->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/trimmed_mean_instlumi_"+stri+"_"+rhidshort+"_beforetrim.pdf");
      	 } */
     	 /*if(thevar.Index("integratelumi_initial") > 0 && (rhidshort=="chamber10_layer2_Endcap1"|| rhidshort =="chamber35_layer4_Endcap1")) {
     		TCanvas * c = new TCanvas("c", "canvas for trimmed mean for instlumi ", 600,800);
     		c->cd();
     	  h_trim->GetXaxis()->SetTitle("Charge");
     		h_trim->GetYaxis()->SetTitle("Entries");
     		h_trim->Draw();
     		TString stri;
     		stri.Form("%d",j);
     		c->SaveAs(output_path+"plotfolders/plotfolder_ME12HV1_other_plots/trimmed_mean_integratelumi_"+stri+"_"+rhidshort+"_beforetrim.pdf");
     	 }*/

       if(h!=0){ //Do the truncation, if h ==0 the truncation is already done
        	for(int it = 1; it<=   proj->GetNbinsX() ;it++) {
	         if(integral < trimmean * normal){ integral+= proj->GetBinContent(it);}
	         else h_trim->SetBinContent(it,0);
        	}
        }
        // truncation done 
      
       if(h_trim->Integral() >0){
       	 if(thevar.Index("integratelumi")>=0 ){
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
       	 }
     	   else{//if charge vs anything else than integrated lumi
     	     htrimmeanvsX->SetBinContent(j, h_trim->GetMean() ); 
     	     htrimmeanvsX->SetBinError(j, h_trim->GetRMS()/sqrt(h_trim->Integral() ) ); 
     	   }
       }
           
        delete proj;
        delete h_trim;
       } //end of loop over  bins of variable of intereset
     
    if(htrimmeanvsX->Integral()==0) continue;


    TCanvas * c3 = new TCanvas;

    //Some cosmetic stuff now...
    htrimmeanvsX->SetTitle(filename+title+"_"+rhidshort);
    htrimmeanvsX->GetYaxis()->SetRangeUser(0,700);
    if(thevar.Index("pressure")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(0,1000);
    if(thevar.Index("integratelumi")>=0 )   htrimmeanvsX->GetYaxis()->SetRangeUser(0.5,1.5);
    htrimmeanvsX->GetXaxis()->SetTitle(xtitle);
		htrimmeanvsX->GetYaxis()->SetTitle("Trimmed mean charge");
		if(thevar.Index("integratelumi")>=0 ) htrimmeanvsX->GetYaxis()->SetTitle("Normalized Trimmed mean charge");
		if(thevar.Index("integratelumi_initial")>=0 ) htrimmeanvsX->GetYaxis()->SetTitle("Normalized Trimmed mean charge");
    
		htrimmeanvsX->SetMarkerStyle(20); htrimmeanvsX->SetMarkerSize(0.7);
    htrimmeanvsX->SetName(filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar);
    c3->SetName("c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar);
	 	if(htrimmeanvsX->Integral()>0) htrimmeanvsX->Write();
    c3->cd(); 
		//Defines the range for the fit,    
    double fitlowedge (0), fithighedge(44);
    if( thevar.Index("pressure")>=0 ) fitlowedge = 953; 
    if( thevar.Index("pressure")>=0 ) fithighedge = 973 ;
    if( thevar.Index("instlumi")>=0 ) fitlowedge = 5000; 
    if( thevar.Index("instlumi")>=0 ) fithighedge =12000 ;

    //Defines the type of function to fit
    if(thevar.Index("integratelumi")>=0){ 

      htrimmeanvsX->Draw("LEP");
		}
			TF1 *fa1; 
    if(thevar.Index("pressure") >=0) {
			fa1= new TF1("f1","expo", fitlowedge,fithighedge );

		  htrimmeanvsX->Fit(fa1, "R");  }
    if(thevar.Index("instlumi")>=0){
			fa1= new TF1("f1","expo", fitlowedge,fithighedge );

		  htrimmeanvsX->Fit(fa1, "R"); }
    //else   fa1= new TF1("f1","pol1", fitlowedge,fithighedge );
    //Now, let's fit
//		htrimmeanvsX->Fit(fa1, "R"); 


    // this is to save the plots for 1 chamber 1 layer and 1 endcap every time
		// taking so much time so stopped it for now
    // generally the below line is commented out and used was commented out and used
		//if(htrimmeanvsX->Integral()>0 && (thevar.Index("integratelumi")>=0 || thevar.Index("pressure") >=0 || thevar.Index("instlumi")>=0) && savehistos){c3->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
/*	 	 if(htrimmeanvsX->Integral()>0 &&( thevar.Index("integratelumi")>=0  || thevar.Index("pressure")>=0  || thevar.Index("instlumi")>=0 ) )
	  	{
				c3->SaveAs(output_path+"plotfolders_removed_trips_again/plotfolder_"+chamber_string_name+"/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); 
			}
		if(htrimmeanvsX->Integral()>0 && (thevar.Index("integratelumi")>=0 ))
		{
			c3->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"_other_plots/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); 
		}*/
   // std::cout<<"the var "<<thevar<<" rhid "<<rhidshort<<std::endl;
   ///		if(htrimmeanvsX->Integral()>0&& thevar=="_integratelumi_initial" &&rhidshort=="chamber10_layer2_Endcap1"){c3->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
   //		if(htrimmeanvsX->Integral()>0&& thevar=="_integratelumi" &&rhidshort == "chamber10_layer2_Endcap1"){c3->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
   //		if(htrimmeanvsX->Integral()>0&& thevar=="_pressure_2" &&rhidshort == "chamber10_layer2_Endcap1"){c3->SaveAs(output_path+"plotfolders/plotfolder_"+chamber_string_name+"/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); }
   // for now I just need to save 1 plot not all to minimize
   //		if(htrimmeanvsX->Integral()>0&& savehistos && (rhidshort == "chamber10_layer2_Endcap1" || rhidshort == "chamber32_layer2_Endcap1" ||rhidshort == "chamber20_layer4_Endcap1") ){c3->SaveAs("plotfolder/c_"+filename+"trimmean"+title+"_"+rhidshort+"vs"+thevar+".pdf"); 
    //}
    
      delete htrimmeanvsX;
      delete c3;

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
    result[h] = parampair;
    
    if(thevar.Index("integratelumi")>=0 ){ 
			double theslope =0; 
		}
		
    if(thevar.Index("instlumi")>=0 ||thevar.Index("pressure")>=0   ){ 
    double theslope = fa1->GetParameter(1);
    double chi2overN = (fa1->GetNDF () >0)? (fa1->GetChisquare ()/fa1->GetNDF () ) : -10 ; 
		
    if(  theslope  !=0 && h!=0){
      h_slope->Fill(theslope);
      h_slopeuncty->Fill(fa1->GetParError(1));
      h_chi2->Fill(chi2overN) ;
    }

    if( theslope  !=0 ){
      int chambandlay = (h>400)? (h-400) : h; 
      _slope =theslope;
      _slope_error =fa1->GetParError(1);
      _chi2 =fa1->GetChisquare () ;
      _ndof = fa1->GetNDF ();
      _chamber = (int) floor( chambandlay/10 ) ; 
      _layer = h%10 ;
      _endcap = (int)floor( h/400 ) +1 ;
      _stationringHVseg =  thestationringhv(filename);
      _isbadchannel = isbadchannel(filename,h);
      theouttree->Fill();
     }
		}//end of if condition for  writing tree for instlumi and intlumi parameters  
	} //end of loop for all the channelss 
/*	      if( thevar.Index("pressure")>=0 )textfile <<"pressure information"<<std::endl;	  	
      if( thevar.Index("pressure")>=0 )textfile << filename<<"\t"   <<(int)floor( h/400 ) +1 << "\t" <<   (int) floor( chambandlay/10 )   <<"\t"  <<h%10  << "\t"<< fixed<<setprecision(11)<< theslope <<"\t" <<  fa1->GetParError(1) << "\t" << chi2overN  <<endl;

			if( thevar.Index("instlumi")>=0 )textfile <<"instlumi information"<<std::endl;	  	
      if( thevar.Index("instlumi")>=0 )textfile << filename<<"\t"   <<(int)floor( h/400 ) +1 << "\t" <<   (int) floor( chambandlay/10 )   <<"\t"  <<h%10  << "\t"<< fixed<<setprecision(11)<< theslope <<"\t" <<  fa1->GetParError(1) << "\t" << chi2overN  <<endl;


	      if( thevar.Index("integrate")>=0 )textfile <<"integrate information"<<std::endl;	  	
      if( thevar.Index("integrate")>=0 )textfile << filename<<"\t"   <<(int)floor( h/400 ) +1 << "\t" <<   (int) floor( chambandlay/10 )   <<"\t"  <<h%10  << "\t"<< fixed<<setprecision(11)<< theslope <<"\t" <<  fa1->GetParError(1) << "\t" << chi2overN  <<endl; */


  if( thevar.Index("instlumi")>=0  ) h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs instant Luminosity (/#mu b^{-1} s^{-1})");
  if( thevar.Index("integrate")>=0  ) h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs #int L (/fb^{-1})");
  if( thevar.Index("pressure")>=0  )  h_slope->GetXaxis()->SetTitle("Slope of #Delta Q vs P (hPa^{-1})");
  outf->cd();
  h_slope->Write();
  h_slopeuncty->Write();
  h_chi2->Write();
  theouttree->Write();

  return result;    
	

} //end of GetSlope function



double ProduceHistosPerChannel::ApplyCorrection( double X ,TString correctiontype, double p0, double p1 ){
  double refvalue = 0;
  if(correctiontype =="pressure"&&!droppressurecorr){
    refvalue =960;
    double thecorr =exp(p1*(refvalue-X)); 

    return thecorr; 
  }

  
  if(correctiontype =="instlumi"&& !dropinstlumicorr){
    refvalue =10000; 
    double thecorr = exp(p1*(refvalue-X)); 

    return thecorr; 
    
  }
  
  return 1;

}
