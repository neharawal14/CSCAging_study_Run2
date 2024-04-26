#include "TreeReader.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include <fstream>
using namespace std;

//Open the root file
// fill : 8102, 8118, 8233,  8236 , 8247
	int fillNb = 8102;
	TString few_chambers = "no";
	TString pressure_equalised = "yes";
void draw_hist(TH1D *h, TString, TString , bool);
std::pair<double, double> trimmed_mean(TH1D *h, TString, TString,bool);
	//TFile *f = TFile::Open();
void draw_hist(TH1D *h, TString saving_name, TString chamber, bool instlumi){
	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(001121211);
	c->cd();
	h->Draw();
	TString saving_name_new;
	TString pressure_file, nb_chambers;
	if(pressure_equalised=="yes") pressure_file = "pressure_equalised"; 
	else pressure_file = "not_equalised";

	if(few_chambers=="yes"){
			nb_chambers = "small_division";
	}
		else {
		nb_chambers = "more_division";}

	if(instlumi) saving_name_new = "./results_ME21HV1/fill_"+TString::Format("%d",fillNb)+"/"+pressure_file+"/"+nb_chambers+"/histogram_charge_"+saving_name+"_chamber_"+chamber+"instlumi.pdf";
	else saving_name_new = "./results_ME21HV1/fill_"+TString::Format("%d",fillNb)+"/"+pressure_file+"/"+nb_chambers+"/histogram_charge_"+saving_name+"_chamber_"+chamber+"time.pdf";
	c->SaveAs(saving_name_new);
	c->Close();
}


int main(){
	TString chamber ="ME21HV1";
	bool ifinstlumi = false;
	TString input_file;
	double y_low, y_up;
	//if(fillNb == 8233 || fillNb==8236 || fillNb ==8247){
	if(chamber=="ME21HV1"){
   y_low = 240; 
		y_up = 290;
	}
/*	else if( fillNb==8247){
    y_low = 265; 
		y_up = 320;
	} */
	else{
    y_low = 260; 
		y_up = 320; 
	}
if(pressure_equalised=="yes"){ input_file=	"/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/2024_HV_equalisation/making_pressure_equalised_ntuples/pressure_equalised_files/output_after_pressure_depenedence_removal_"+chamber+"_everything.root";}
else { input_file = "/cmsuf/data/store/user/t2/users/neha.rawal/CSCAgeing_2022/final_files_2022_new/2022_all/csc_output_2022_"+chamber+"_tree.root";
}

// Fill nb : 8102a
// Start Time
//2022-08-06 17:19:14 --- UTC: 15:19:14
//Stable Beams Start
//2022-08-06 19:25:35 --- UTC: 17:25:35
//Stable Beams End
//2022-08-07 10:07:11 --- UTC: 08:07:11
//End Time
//2022-08-07 10:08:50 --- UTC: 08:08:50
//Delivered Lumi Total
//518.97 pb−1
// Run 356939 - 356941
// Run 356943 - 356956

	// Start time -  2022-08-06 19:25:25 UTC  //1659813925 000
//	// time 2 = 2022-08-06 20:30: 16 UTC 1659817816
	// time 3 = 2022-08-06 21:29: 45 UTC 1659821385
	// Middle time -  2022-08-07 00:03:00 UTC // 1659830580
	// time 4 = 2022-08-07 03:12:25 UTC 1659841945 
	// end time - 2022-08-07 07:09:54 UTC //1659856180 000
	std::vector<int> runNb ; 
	std::vector<ULong64_t> time_my; 
	std::vector<double> instlumi_my ;
if(fillNb==8102){	
	runNb= {356939,356940,356941,356942,356943,356944,356945,356946,
 356947,356948,356949,356950,356951,356952,356953,	356953,356954,356955,356956	};
if(few_chambers=="yes"){
	time_my = {1659813925, 1659830580, 1659856180};
	instlumi_my = {13748.88,  9871.87,  6850.11};
}
else{
	time_my = {1659813925, 1659821385, 1659830580,1659841945, 1659856180};
	instlumi_my = {13748.88, 11744.03,  9871.87, 8406.88,  6850.11};
	}
}

if(fillNb==8118){
	runNb = {357325, 357326,357327,357328,357329,
		357330,357331,357332, 357333};

  if(few_chambers!="yes"){
		// time =  {2022-08-11 21:46:41, 23:06:34, 2022-08-12 01:28:52, 03:36:50, 05:20:32};
	 time_my = {1660254401 ,1660259194,1660267732 ,1660275410 ,  1660281632};
	instlumi_my = {15809.7 , 14028.2,  11789.63, 10268.4, 9362.2};
 }
	else{
	time_my = {1660254401 ,1660267732,   1660281632};
	instlumi_my = {15809.7,  11789.63,  9362.2};
	}
 }
if(fillNb==8233){
   runNb= {359899,359901,359902,359903,359904,359905,359906,359907,359908};

  if(few_chambers!="yes"){
  // time = { 2022-10-06 10:30:26, 11:09:46, 12:05:08, 2022-10-06 13:16:34};
	 time_my = { 1665052226, 1665054586,1665057908, 1665062194}; 
	 instlumi_my = {1392.82, 1273.41, 1157.50, 1068.95};
	}
	else{
	 time_my = { 1665052226, 1665054586, 1665062194}; 
	 instlumi_my = {1392.82, 1273.41,1068.95};
	}
}
if(fillNb==8236){
   runNb= {359985,359986,359987,359988,359989,359990,359991,359992,359993,359994,359995,359996,359997,359998,359999,36000, 360001,360002, 360003,360004,360005};

  if(few_chambers!="yes"){
  // time = { 2022-10-07 17:49:58, 19:56:43, 22:37:13, 2022-10-08 1:36:25, 04:11:21};
	 time_my = {1665164998, 1665172603,  1665182233,  1665192985, 1665202281}; 
	 instlumi_my = {19052.02, 15776.09, 13071.51, 10963.99, 9526.38};
	}
	else{
  // time = { 2022-10-07 17:49:58, 22:37:13, 04:11:21};
	 time_my = {1665164998,   1665182233,   1665202281}; 
	 instlumi_my = {19052.02,  13071.51,  9526.38};

	}
}
if(fillNb==8247){
   runNb= {360116,360120,360121,360122,360123,360124,360125,360126,360127,360128,360129,360130,360131,360132};

  if(few_chambers!="yes"){
  // time = { 2022-10-10 14:42:14, 16:34:54, 19:58:41,23:01:35 , 2022-10-11 03:22:39 };
	 time_my = {1665412934, 1665419694, 1665431921, 1665442895, 1665458559}; 
	 instlumi_my = {18323.64, 18375.29, 143234.36,  11920.71, 9413};
	}
	else{
  // time = { 2022-10-07 17:49:58, 22:37:13, 04:11:21};
	 time_my = {1665412934,  1665431921, 1665458559}; 
	 instlumi_my = {18323.64, 143234.36,  9413};

	}
}
// Fill nb : 8118
	TreeReader *tree = new TreeReader(input_file);
	tree->initialise();

	std::vector<TH1D*> hcharge; 


	 for(int i=0; i<time_my.size()-1; i++){
    TH1D * h = new TH1D("h","charge (ADC counts)  : time "+chamber + " : "+i+" fill : instlumi",3000,0,3000);
     if(ifinstlumi) h->SetTitle("charge (ADC counts)  : time "+chamber + " : "+i+" fill : instlumi ");
		 else h->SetTitle("charge (ADC counts)  : time "+chamber + " : "+i+" fill : time");
     hcharge.push_back(h);
	 }
  bool flag;
  std::vector<int> nb_entries = 
    {0, 0,0,0,0,
		0,0,0,0,0,0,0,
		0,0, 0,0,0,0};
//ofstream tex_file;
std::vector<TString> num= {"first", "second", "third", "fourth", "fifth"};
//tex_file.open("star"+chamber+"_time.txt");
for(int i = 0; i < tree->tree->GetEntries(); i++){
				 flag = false;
	 	     tree->tree->LoadTree(i);
				 tree->tree->GetEntry(i);
         if(i%1000000 ==0)cout << i<<endl;

				 for(int j=0; j<runNb.size();j++){
				 if(tree->_runNb == runNb[j]) {
					 nb_entries[j] = nb_entries[j]+1;
					 flag=true;
				  }
				 }

			if(flag==false) continue;
			else{
				for(int j=0; j<hcharge.size(); j++){
				//for(int j=0; j<2; j++){
				//
			if(ifinstlumi){
				 if(tree->_instlumi > instlumi_my[j+1]  && tree->_instlumi <= instlumi_my[j]) {
				  double charge = tree->_rhsumQ_RAW;
				  //double pressure_charge = tree->_rhsumQ_RAW_new;
          hcharge[j]->Fill(charge);
         }
			}
			else{
				 if(tree->_timesecond <time_my[j+1]  && tree->_timesecond >= time_my[j]) {
				   //if(tree->_instlumi <13748.887  && tree->_instlumi >= 9871.873) {
					 //tex_file<<" instlumi values first half : run:  "<<tree->_runNb<<" : "<<tree->_instlumi<<tree->_instlumi<<std::endl;
				//std::cout<<" : half : run:  "<<j<<" runNb "<<tree->_runNb<<" : "<<tree->_instlumi<<" time "<<tree->_timesecond<<std::endl;
				  double charge = tree->_rhsumQ_RAW;
				  //double pressure_charge = tree->_rhsumQ_RAW_new;
          hcharge[j]->Fill(charge);
         }
				} // another else end here
			}
			// else end here
     }
}
// end of tree loop

				 for(int j=0; j<runNb.size();j++){
					 std::cout<<" number of entries in run "<<runNb[j]<<" num "<<nb_entries[j]<<std::endl;
				 }
std::vector<std::pair<double, double>> trimmed_mean_vec;

				for(int j=0; j<hcharge.size(); j++){
					std::cout<<"title "<<hcharge.at(j)->GetTitle()<<" : "<<hcharge.at(j)->GetEntries()<<std::endl;
				}
				for(int j=0; j<hcharge.size(); j++){
		 	    std::pair<double, double> trimmed_mean_value;
					std::cout<<"title "<<hcharge.at(j)->GetTitle()<<" : "<<hcharge.at(j)->GetEntries()<<std::endl;
			    trimmed_mean_value = trimmed_mean(hcharge.at(j), num.at(j), chamber, ifinstlumi);
					trimmed_mean_vec.push_back(trimmed_mean_value);
				}
 
				for(int j=0; j<hcharge.size(); j++){
				std::cout<<" trimmed mean first "<<trimmed_mean_vec[j].first<<std::endl;
				}

TH1D *h_new = new TH1D("h_new", "hist new",time_my.size()-1, time_my[0], time_my[time_my.size()-1]); 
//time_my = {1659813925, 1659821385, 1659830580,1659841945, 1659856180};
for(int i=1; i<time_my.size(); i++){
 h_new->SetBinContent(i, trimmed_mean_vec[i-1].first);
 h_new->SetBinError(i, trimmed_mean_vec[i-1].second);
}
TString title_m;
if(pressure_equalised=="yes") { 
 if(ifinstlumi) title_m = " Gas gain with instlumi : fill : "+TString::Format("%d", fillNb) + " pressure_equalised ";
 else title_m = " Gas gain with time : fill : "+TString::Format("%d", fillNb) + " pressure_equalised";
}
else { 
 if(ifinstlumi) title_m = " Gas gain with instlumi : fill : "+TString::Format("%d", fillNb) + " not pressure_equalised ";
 else title_m = " Gas gain with time : fill : "+TString::Format("%d", fillNb) + " not pressure_equalised";
}

h_new->SetTitle(title_m);
h_new->GetYaxis()->SetTitle(" Gas gain");// time = { 2022-10-07 17:49:58, 22:37:13, 04:11:21};
     h_new->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
		 h_new->GetXaxis()->SetTimeOffset(0,"gmt");
		 h_new->GetXaxis()->SetTimeDisplay(1);
		 h_new->GetXaxis()->SetLabelSize(0.02);
		 h_new->GetYaxis()->SetRangeUser(y_low, y_up);
TCanvas * c  = new TCanvas();
gStyle->SetOptStat(0);
c->cd();
h_new->Draw();
TString save_name_1, nb_chambers;
 if(few_chambers=="yes"){
         nb_chambers = "small_division";
     }
       else {
       nb_chambers = "more_division";}
if(pressure_equalised=="yes") 
save_name_1= "./results_ME21HV1/fill_"+TString::Format("%d", fillNb)+"/pressure_equalised/"+nb_chambers+"/gas_gain_fill_"+TString::Format("%d", fillNb)+".pdf";
else
save_name_1= "./results_ME21HV1/fill_"+TString::Format("%d", fillNb)+"/not_equalised/"+nb_chambers+"/gas_gain_fill_"+TString::Format("%d", fillNb)+".pdf";
c->SaveAs(save_name_1);

// Draw the slope
 
return 0 ; 
}

// Find the trimmed mean of the histogram 
std::pair<double, double> trimmed_mean(TH1D * hist, TString name, TString chamber, bool instlumi_bool){

	double trimmean = 0.7;
	double normal; 

	TH1D *h_prov = (TH1D *) hist->Clone();
	TH1D *h_prov_new = (TH1D *) hist->Clone();

	h_prov_new->Reset();
	h_prov_new->ResetStats();
	h_prov_new->SetTitle((TString) h_prov->GetTitle()+" after trimming");
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
							
				for(int it=1; it<last_bin ; it++) {
							  	h_prov_new->SetBinContent(it,h_prov->GetBinContent(it));
							  	h_prov_new->SetBinError(it,h_prov->GetBinError(it));
						   	}
						
								h_prov_new->SetBinContent(last_bin, entries_last_bin);
							  h_prov_new->SetBinError(last_bin, h_prov->GetBinError(last_bin));

							TString save_name_before = name+"_before_trim"	;
							TString save_name_after = name+"_after_trim"	;

				draw_hist(h_prov, save_name_before, chamber, instlumi_bool);
				draw_hist(h_prov_new, save_name_after, chamber, instlumi_bool);
				double	mean_trimmed = h_prov_new->GetMean();
				double	mean_trimmed_error = h_prov_new->GetMeanError();
				std::pair<double, double> 	trim_mean = std::make_pair(mean_trimmed,mean_trimmed_error);
				return trim_mean;
}


