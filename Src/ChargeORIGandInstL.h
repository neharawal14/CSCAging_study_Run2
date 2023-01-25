//#include <string> 
#include <sstream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>



#include "../files_HVandLumi/nonme11_first.h"
#include "../files_HVandLumi/nonme11_second.h"
#include "../files_HVandLumi/me11.h"

double instlumi(int runnb= 284029, int lumis= 100){
  TString runnb_tstr = (TString) Form("%d",runnb);
  TString ls_tstr = (TString)Form("%d",lumis);
 
  TString tstrtest = "../files_HVandLumi/InstLumiPerRun/2022/"+runnb_tstr+".csv";
  ifstream lumifile;
  lumifile.open(tstrtest.Data());
  string value;

  runnb_tstr = (TString) Form("%d",runnb)+":";
  ls_tstr = (TString)Form("%d",lumis)+":";

  while ( lumifile.good() )
    {
      getline ( lumifile, value ); 
 
      size_t posrunnb = value.find(runnb_tstr.Data() );
      if(posrunnb!=string::npos){
	size_t posLS = value.find( ls_tstr.Data() ,posrunnb+10);
	if(posLS!=string::npos){

	size_t pos2 =  value.find( "STABLE BEAMS");
	pos2 =  value.find( ",",pos2+1);pos2 =  value.find( ",",pos2+1);pos2 =  value.find( ",",pos2+1);
	string thelumistring = string( value, pos2+1, 8);
	double thelumi = stod(thelumistring) ; 
	
	return thelumi/23.31;
	}
      }
    }
  return 0;
} 




std::pair<double,double> UncorrGasGain_HVInitial2016(double charge = 400, int runnb = 277792, int stationring=12, int rhid =  1120151){
    std::pair<double,double> chargeandHV(0,0);
  double dHV_(0),HV_(0);
  if(stationring ==11 || stationring ==14){

    if( stationring ==14) rhid -=30000 ;
    rhid -= rhid%10; rhid/=10;

    dHV_ = (dHV_HV_ME11(rhid) ).first;
    HV_ = (dHV_HV_ME11(rhid) ).second;
    
  }
  else if(runnb>= 281613){
    dHV_ = (dHV_HV_NonME11_v2(rhid) ).first;
    HV_ = (dHV_HV_NonME11_v2(rhid) ).second;
  } 
  else{
    dHV_ = (dHV_HV_NonME11_v1(rhid) ).first;
    HV_ = (dHV_HV_NonME11_v1(rhid) ).second;
  }  

  if( runnb <277792 ){HV_ -=  dHV_; dHV_=0;}
  if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV_ -=  dHV_; dHV_=0;}

  
  /*Recipe:
	
    dHV=-ln(G/G0)/B where
    G - gas gain in HV segment of the CSC layer (as Trimmed Mean of
    distribution SumQ=3X3 ADC
    with parameter 0.7)
    G0 - reference global gas gain
    B - slope (from HV scan data in 2015)
    = for inner CSC (ME2/1,ME3/1, ME4/1)     G0=263, B=5.193*10**(-3) 1/V
    for outer CSC (ME1/2, ME1/3, ME234/2)  G0=373, B=5.463*10**(-3) 1/V
    for ME11: G0= 330, B= 6.262*10**(-3) 1/V
    S7 of https://indico.cern.ch/event/555989/contributions/2242293/attachments/1309116/1957949/terentiev_2015_CSC_HV_corrections_talk_v4.pdf
  */

  double const_B_ = 0;
  if( (stationring==11|| stationring ==14) ) const_B_ = 6.26e-3;
  else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B_ = 5.193e-3;
  else  const_B_ = 5.463e-3;

  double chargeold_ = charge * exp( -const_B_* dHV_ ) ;

  chargeandHV.first = chargeold_;
  chargeandHV.second = HV_;

  return chargeandHV ;

}


