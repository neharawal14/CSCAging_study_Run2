//#include <string> 
#include <sstream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>



#include "../files_HVandLumi/nonme11_first.h"
#include "../files_HVandLumi/nonme11_second.h"
#include "../files_HVandLumi/me11.h"
ifstream ifile_me11( "/cms/data/store/user/t2/users/laurentthomas/CSCHVStudies/CMSSW_8_0_27/src/custogasgain/files_HVandLumi/CAEN_ME11_0.7_2015.txt" );
ifstream ifile_nonme11_v1( "/cms/data/store/user/t2/users/laurentthomas/CSCHVStudies/CMSSW_8_0_27/src/custogasgain/files_HVandLumi/Ch_dhvG0_nonME11_inner_outer_0.7_2015.txt" );
ifstream ifile_nonme11_v2( "/cms/data/store/user/t2/users/laurentthomas/CSCHVStudies/CMSSW_8_0_27/src/custogasgain/files_HVandLumi/Ch_dhvG0_nonME11_inner_outer_0.7_2015_v2.txt" );
//ifstream lumifile("/cms/raid/raid9/laurentthomas/CSCHVStudies/CMSSW_8_0_23/src/custogasgain/files_HVandLumi/instlumi_bu.csv") ;




double instlumi(int runnb= 284029, int lumis= 100){
  TString runnb_tstr = (TString) Form("%d",runnb);
  TString ls_tstr = (TString)Form("%d",lumis);
 
  TString tstrtest = "/home/laurentthomas/CSCHVStudies/CMSSW_8_0_27/src/custogasgain_new/files_HVandLumi/InstLumiPerRun/"+runnb_tstr+".csv";
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

std::pair<double,double> UncorrGasGain_HVInitial2016_nonME11_v2(double charge = 400, int runnb = 277792, int stationring=12, int rhid =  1120151){
  std::pair<double,double> chargeandHV(0,0);
  
  double dHV_ = (dHV_HV_NonME11_v2(rhid) ).first;
  double HV_ = (dHV_HV_NonME11_v2(rhid) ).second;
  if( runnb <277792 ){HV_ -=  dHV_; dHV_=0;}
  if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV_ -=  dHV_; dHV_=0;}


  double const_B_ = 0; 
  if( (stationring==11|| stationring ==14) ) const_B_ = 6.26e-3; 
  else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B_ = 5.193e-3; 
  else  const_B_ = 5.463e-3;
  
  double chargeold_ = charge * exp( -const_B_* dHV_ ) ;

  chargeandHV.first = chargeold_;
  chargeandHV.second = HV_;
  
  //    cout << charge<<", " <<runnb<< ", "<< stationring << ", "<< rhid<<endl;
    //  cout << "hv, dhv, charge, chargeold " << HV_<<", " <<dHV_<<", "<< charge<<", "<< chargeold_ <<endl;
  return chargeandHV ;
  
  //  ifstream file = ifile_nonme11_v2;
//This function unsets the HV change that happened during the middle of the 2016 data taking and retrieves the corresponding uncorrected charge.
  // The HV corrections were applied in two steps:

  //  1) first version of non-ME1/1 corrections: 15:00 26th of July 2016 . First stable beam run after that: 277792
  //  2) second version of non-ME1/1 corrections (with cap at 150V) and also the ME1/1 corrections were applied at 16th of September 2016 i.e.  First stable beam run after that: 281613

  //  cout << charge<<", " <<runnb<< ", "<< stationring << ", "<< rhid<<endl;

  int posforread = 37;

  if( stationring ==14) rhid -=30000 ;
  if( stationring==11|| stationring ==14) {rhid -= rhid%10; rhid/=10; posforread-=5;}
  
  

  /*

  double dHV(0), HV(0);
  if(stationring==11|| stationring ==14){
    dHV =  dHV_HV_NonME11_first(rhid);
    HV =   dHV_HV_NonME11_second(rhid);
  }
  else if(runnb < 281613){
    dHV_HV_NonME11_second(rhid);
  else  dHV_HV_ME11(rhid);



  	string dHVstring = string( value, posforread,5);
	double dHV = stod(dHVstring) ;
	string HVstring = string( value, posforread+5 ,6);
	double HV = stod(HVstring) ;
	if( runnb <277792 ){HV -=  dHV; dHV=0;}
	if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV -=  dHV; dHV=0;}
	
  */	

	//Now calculating the uncorrected charge
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
  /*
	double const_B = 0; 
	if( (stationring==11|| stationring ==14) ) const_B = 6.26e-3; 
	else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B = 5.193e-3; 
	else  const_B = 5.463e-3;
		
	double chargeold = charge * exp( -const_B* dHV ) ;

	chargeandHV.first = chargeold;
	chargeandHV.second = HV;
}
  

  */




  string value;
  
  //  return chargeandHV ;
  //  ifstream ifile_nonme11_v2;
  ifile_nonme11_v2.clear();
  ifile_nonme11_v2.seekg (0, ios::beg);


  /*  if( stationring==11|| stationring ==14) ifile_nonme11_v2.open( "/cms/raid/raid9/laurentthomas/CSCHVStudies/CMSSW_8_0_23/src/custogasgain/ifile_nonme11_v2s_HVandLumi/CAEN_ME11_0.7_2015.txt" );
  else if(runnb<281613) ifile_nonme11_v2.open( "/cms/raid/raid9/laurentthomas/CSCHVStudies/CMSSW_8_0_23/src/custogasgain/ifile_nonme11_v2s_HVandLumi/Ch_dhvG0_nonME11_inner_outer_0.7_2015.txt" );
  else ifile_nonme11_v2.open( "/cms/raid/raid9/laurentthomas/CSCHVStudies/CMSSW_8_0_23/src/custogasgain/ifile_nonme11_v2s_HVandLumi/Ch_dhvG0_nonME11_inner_outer_0.7_2015_v2.txt" );
  string value;*/
  while ( ifile_nonme11_v2.good() )
    {
      getline ( ifile_nonme11_v2, value ); 
      TString rhid_tstr = (TString) Form("%d",rhid);
      size_t posrhid = value.find(rhid_tstr.Data() );
      if(posrhid!=string::npos){
	//	cout << string( value, 0, value.length()-1 )<<endl;
	
	//cout << string( value, posforread, 5 )<<endl;
	string dHVstring = string( value, posforread,5);
	double dHV = stod(dHVstring) ;
	string HVstring = string( value, posforread+5 ,6);
	double HV = stod(HVstring) ;
	if( runnb <277792 ){HV -=  dHV; dHV=0;}
	if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV -=  dHV; dHV=0;}
	
	

	//Now calculating the uncorrected charge
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

	double const_B = 0; 
	if( (stationring==11|| stationring ==14) ) const_B = 6.26e-3; 
	else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B = 5.193e-3; 
	else  const_B = 5.463e-3;
		
	double chargeold = charge * exp( -const_B* dHV ) ;

	chargeandHV.first = chargeold;
	chargeandHV.second = HV;

	//	cout << "hv, dhv, charge, chargeold " << HV<<", " <<dHV<<", "<< charge<<", "<< chargeold <<endl;

	return chargeandHV;	  
	}
    }
  
  double dummy ; 
  //cout << "Should not reach that point nonme11v2" << endl;
  cin >>dummy;
  return chargeandHV; 

}
































std::pair<double,double> UncorrGasGain_HVInitial2016_nonME11_v1(double charge = 400, int runnb = 277792, int stationring=12, int rhid =  1120151){
  
  return( dHV_HV_NonME11_v1(rhid) );
  int posforread = 37;
  std::pair<double,double> chargeandHV(0,0);
  if( stationring ==14) rhid -=30000 ;
  if( stationring==11|| stationring ==14) {rhid -= rhid%10; rhid/=10; posforread-=5;}
  

  string value;
  

  ifile_nonme11_v1.clear();
  ifile_nonme11_v1.seekg (0, ios::beg);

  while ( ifile_nonme11_v1.good() )
    {
      getline ( ifile_nonme11_v1, value ); 
      TString rhid_tstr = (TString) Form("%d",rhid);
      size_t posrhid = value.find(rhid_tstr.Data() );
      if(posrhid!=string::npos){
	//	cout << string( value, 0, value.length()-1 )<<endl;
	
	//cout << string( value, posforread, 5 )<<endl;
	string dHVstring = string( value, posforread,5);
	double dHV = stod(dHVstring) ;
	string HVstring = string( value, posforread+5 ,6);
	double HV = stod(HVstring) ;
	if( runnb <277792 ){HV -=  dHV; dHV=0;}
	if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV -=  dHV; dHV=0;}
	
	

	//Now calculating the uncorrected charge
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

	double const_B = 0; 
	if( (stationring==11|| stationring ==14) ) const_B = 6.26e-3; 
	else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B = 5.193e-3; 
	else  const_B = 5.463e-3;
		
	double chargeold = charge * exp( -const_B* dHV ) ;

	chargeandHV.first = chargeold;
	chargeandHV.second = HV;

	//	cout << "hv, dhv, charge, chargeold " << HV<<", " <<dHV<<", "<< charge<<", "<< chargeold <<endl;

	return chargeandHV;	  
	}
    }
  
  double dummy ; 
  cout << "Should not reach that point nonme11v1" << endl;
  cin >>dummy;
  return chargeandHV; 

}












std::pair<double,double> UncorrGasGain_HVInitial2016_ME11(double charge = 400, int runnb = 277792, int stationring=12, int rhid =  1120151){
  
  return( dHV_HV_ME11(rhid) );
  int posforread = 37;
  std::pair<double,double> chargeandHV(0,0);
  if( stationring ==14) rhid -=30000 ;
  if( stationring==11|| stationring ==14) {rhid -= rhid%10; rhid/=10; posforread-=5;}
  

  string value;
  

  ifile_me11.clear();
  ifile_me11.seekg (0, ios::beg);

  while ( ifile_me11.good() )
    {
      getline ( ifile_me11, value ); 
      TString rhid_tstr = (TString) Form("%d",rhid);
      size_t posrhid = value.find(rhid_tstr.Data() );
      if(posrhid!=string::npos){
	//	cout << string( value, 0, value.length()-1 )<<endl;
	
	//cout << string( value, posforread, 5 )<<endl;
	string dHVstring = string( value, posforread,5);
	double dHV = stod(dHVstring) ;
	string HVstring = string( value, posforread+5 ,6);
	double HV = stod(HVstring) ;
	if( runnb <277792 ){HV -=  dHV; dHV=0;}
	if( runnb < 281613&& (stationring==11|| stationring ==14)  ){HV -=  dHV; dHV=0;}
	
	

	//Now calculating the uncorrected charge
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

	double const_B = 0; 
	if( (stationring==11|| stationring ==14) ) const_B = 6.26e-3; 
	else if( (stationring==21|| stationring ==31 || stationring ==41) ) const_B = 5.193e-3; 
	else  const_B = 5.463e-3;
		
	double chargeold = charge * exp( -const_B* dHV ) ;

	chargeandHV.first = chargeold;
	chargeandHV.second = HV;

	//	cout << "hv, dhv, charge, chargeold " << HV<<", " <<dHV<<", "<< charge<<", "<< chargeold <<endl;

	return chargeandHV;	  
	}
    }
  
  double dummy ; 
  cout << "Should not reach that point me11" << endl;
  cin >>dummy;
  return chargeandHV; 

}

