#include "SystematicRemoval.h"

void SystematicRemoval::initialise(TString input_path_string, TString saving_path_string, TString normalising_chamber_string){

  gStyle->SetOptFit(2211);
   flag= true;
  input_path  = input_path_string;
 // "outf_dataset_pressure_corrected__"+chamber[i]+"_output_run2_const.root";
  TString input_file_name_ref =  input_path+"outf_dataset_pressure_corrected__"+normalising_chamber_string+"_output_run2_const.root";
//"outf_dataset_pressure_corrected__"+normalising_chamber_string+"_output_run2.root";
  input_file_ref = TFile::Open(input_file_name_ref,"READ");
  saving_path = saving_path_string;
 // chamber_name = chamber_name_string;
  normalising_chamber = normalising_chamber_string;

  TDirectoryFile *dir_ref = (TDirectoryFile*) input_file_ref->Get("_integratelumi");
  if(!dir_ref){
  return;}
  if (dir_ref->GetListOfKeys()->GetSize() <= 0){
  return;
  } 
  std::cout<<" access directory"<<std::endl;
//  TString  histogram_plus_name = "_dataset_pressure_corrected_trimmean_allgoodchannels_plusvs_integratelumi";
  TString  histogram_name = "_dataset_pressure_corrected_trimmean_allgoodchannelsvs_integratelumi";
  //TString  histogram_minus_name = "_dataset_pressure_corrected_trimmean_allgoodchannels_minusvs_integratelumi";

  h_goodchannels_avg_ref = (TH1D*) dir_ref->Get(histogram_name);
  h_goodchannels_avg_ref->GetYaxis()->SetRangeUser(200,500);
}

void SystematicRemoval::adding_plus_minus(TString chamber_name){
  
  TString input_file_name_org = input_path+"outf_dataset_pressure_corrected__"+chamber_name+"_output_run2_const.root";
  std::cout<<" input file trying to open "<<input_file_name_org<<std::endl;
  input_file_org = TFile::Open(input_file_name_org,"READ");
 
  TDirectoryFile *dir_org = (TDirectoryFile*) input_file_org->Get("_integratelumi");
  if(!dir_org){
     flag=false;
    return;}
  if (dir_org->GetListOfKeys()->GetSize() <= 0){
  flag = false; 
  return;
  } 

  std::cout<<" access directory"<<std::endl;
  TString  histogram_name = "_dataset_pressure_corrected_trimmean_allgoodchannelsvs_integratelumi";

  std::cout<<" histogram name  "<<histogram_name<<std::endl;
  h_goodchannels_avg = (TH1D*) dir_org->Get(histogram_name);
  h_goodchannels_avg->GetYaxis()->SetRangeUser(200,500);

  std::cout<<" will draw histogram  "<<histogram_name<<std::endl;
 TCanvas *c = new TCanvas();
 c->cd();
 h_goodchannels_avg->Draw();
 h_goodchannels_avg->SetTitle("cumulative : "+chamber_name);
 //h_goodchannels_avg->GetYaxis()->SetTitle("integrated luminosity (fb^{-1})");
 c->SaveAs(saving_path+"integratelumi_"+chamber_name+".pdf"); 


}

void SystematicRemoval::ratio(TString chamber_name){

  if(flag==false) return;

  std::cout<<" will start ratio module "<<std::endl;
  h_normalised_avg = (TH1D*) h_goodchannels_avg->Clone();
  h_normalised_avg->Reset();
  h_normalised_avg->ResetStats();
  h_normalised_avg_first = (TH1D*) h_goodchannels_avg->Clone();
  h_normalised_avg_first->Reset();
  h_normalised_avg_first->ResetStats();

  TCanvas *c5 = new TCanvas();
  c5->cd();
  h_goodchannels_avg->Draw();
  h_goodchannels_avg->SetMarkerColor(kBlue);
  h_goodchannels_avg_ref->SetMarkerColor(kRed);
  h_goodchannels_avg_ref->Draw("SAME");
  c5->SaveAs(saving_path+"ratio_normalised_"+chamber_name+"_avg_ME13HV3.pdf");

  double first_non_zero_value =0;
  double first_non_zero_err_value =0;
  double first_non_zero_ref_value =0;
  double first_non_zero_ref_err_value =0;

  //std::cout<<" inside chamber place "<<chamber_name<<std::endl; 

for (int i = 1; i <= h_goodchannels_avg->GetNbinsX(); i++) {
    double A_i = h_goodchannels_avg->GetBinContent(i);
    double B_i = h_goodchannels_avg_ref->GetBinContent(i);
    double sigma_A_i = h_goodchannels_avg->GetBinError(i);
    double sigma_B_i = h_goodchannels_avg_ref->GetBinError(i);

  //std::cout<<" inside chamber place "<<chamber_name<<" A "<<A_i<<" B "<<B_i<<std::endl; 


    if (B_i != 0) {
        double ratio_value_avg = A_i / B_i;
        double sigma_R_i = ratio_value_avg* sqrt((sigma_A_i / A_i) * (sigma_A_i / A_i) + (sigma_B_i / B_i) * (sigma_B_i / B_i));

        h_normalised_avg->SetBinContent(i, ratio_value_avg);
        h_normalised_avg->SetBinError(i, sigma_R_i);

    }
    else {
        h_normalised_avg->SetBinContent(i, 0);
        h_normalised_avg->SetBinError(i, 0);
        h_normalised_avg_first->SetBinContent(i, 0);
        h_normalised_avg_first->SetBinError(i, 0);
    }

}
    h_normalised_avg->GetYaxis()->SetRangeUser(0.6,1.4);
// Taking ratio of entries with first bin entry
bool flag_first_bin=false;
double value_first_bin, value_first_bin_err;
for (int i = 1; i <= h_normalised_avg->GetNbinsX(); i++) {
    double A_i = h_normalised_avg->GetBinContent(i);

    if(A_i==0) {
        h_normalised_avg_first->SetBinContent(i, 0);
        h_normalised_avg_first->SetBinError(i, 0);
    }
    if(A_i!=0 && flag_first_bin==false){
        flag_first_bin = true;
        value_first_bin = h_normalised_avg->GetBinContent(i);
        value_first_bin_err = h_normalised_avg->GetBinError(i);
    }
    if(A_i!=0 && flag_first_bin==true){
        double value = h_normalised_avg->GetBinContent(i); 
        double value_err = h_normalised_avg->GetBinError(i); 

        double ratio = value/value_first_bin;
        double ratio_err = ratio * sqrt( (value_first_bin_err/value_first_bin) *  (value_first_bin_err/value_first_bin)+
                 (value_err/value) * (value_err/value) ) ;
        h_normalised_avg_first->SetBinContent(i, ratio);
        h_normalised_avg_first->SetBinError(i, ratio_err);
    }
} // end of histogram bins ratios
    h_normalised_avg_first->GetYaxis()->SetRangeUser(0.95,1.05);
  TCanvas* c = new TCanvas();
  c->cd();
  h_normalised_avg->Draw("P");
  c->SaveAs(saving_path+"integratelumi_normalised_"+chamber_name+"_avg.pdf");

  TCanvas* c1 = new TCanvas();
  c1->cd();
  h_normalised_avg_first->Draw("P");
  c1->SaveAs(saving_path+"integratelumi_normalised_"+chamber_name+"_avg_normalised_first_bin.pdf");

  std::cout<<" done with ratio"<<std::endl;
} // end of ratio fucntion module
//Normalising the fit
 std::pair<float, float> SystematicRemoval::normalised_fit(TString chamber_name){
 // std::cout<<" fitting the  ratio"<<std::endl;
  std::pair<float, float> mean_slope;
 if(flag==false){
  mean_slope.first= 0; 
  mean_slope.second= 0; 
  return mean_slope;
 }
 
double bin_low_edge, bin_up_edge;
bool flag_first=false;
for (int i = 1; i <= h_normalised_avg_first->GetNbinsX(); i++) {
    double A_i = h_normalised_avg_first->GetBinContent(i);

    if (A_i != 0 && !flag_first) {
        // Set the lower edge to the first non-zero bin
        bin_low_edge = h_normalised_avg_first->GetBinLowEdge(i);
        flag_first = true;
    }

    if (A_i != 0) {
        // Continuously update the upper edge to the upper edge of the last non-zero bin
        bin_up_edge = h_normalised_avg_first->GetBinLowEdge(i + 1); // Use next bin's lower edge
    }
}

// Ensure the range is valid
if (bin_low_edge >= bin_up_edge) {
    std::cerr << "Error: Invalid fit range. No non-zero entries found.\n";
}

  //TF1 *fa1 = new TF1("fa1","[0]*x+[1]", 0,145);
  TF1 *fa1 = new TF1("fa1","[0]*x+[1]",bin_low_edge, bin_up_edge );
  fa1->SetParName(0, "slope");
  fa1->SetParName(1, "constant");
  fa1->SetParameters(0.00001,1);

  h_normalised_avg_first->Fit(fa1, "R");
  //h_normalised_plus->Fit(fa1);

      // Manually add fit results using TPaveText
    TPaveText *pt = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText("Fit Results:");
    pt->AddText(Form("Slope: %.4f #pm %.4f", fa1->GetParameter(0), fa1->GetParError(0)));
    pt->AddText(Form("Intercept: %.4f #pm %.4f", fa1->GetParameter(1), fa1->GetParError(1)));

  TCanvas *c = new TCanvas();
  c->cd();
  gStyle->SetOptFit(1111);
  h_normalised_avg_first->Draw("P");
  //h_normalised_plus->Draw("P");
  fa1->Draw("same");
  //pt->Draw("same");
  c->Update();
  c->SaveAs(saving_path+"integratelumi_normalised_"+chamber_name+"_first_bin_avg_fit.pdf");

//  std::cout<<" done fitting the  ratio"<<std::endl;
  //std::pair<float, float> mean_plus;
	float slope, slope_error;
	slope = fa1->GetParameter(0);	
	slope_error =  fa1->GetParError(0);	

	//mean_plus.first= slope_plus;
	//mean_plus.second= slope_error_plus;

	//mean_minus.first= slope_minus;
	//mean_minus.second= slope_error_minus;

  mean_slope.first = slope;
  mean_slope.second = slope_error ;
  std::cout<<" mean values "<<mean_slope.first<<" errror "<<mean_slope.second<<std::endl;
	return mean_slope;
}

