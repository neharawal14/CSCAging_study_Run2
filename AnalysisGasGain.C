#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TF1.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TPaveStats.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// Declaring all the functions 

class GasDependence {

		public : 
    // variables
    // Map year strings to condition functions
    TString area_name;
    TString year = "2016";
    TString chamber_name = "ME12HV1";
    TString type="all_channels";
    std::vector<std::pair<int, int>> run_ranges;
    std::vector<std::pair<int, int>> time_ranges;
    bool debug = false;

    double _rhsumQ_equalised_HV_data, _pressure, _rhsumQ_RAW;
    double _rhsumQ;
    Int_t _rhid;
    double _integratelumi, _instlumi;
    ULong64_t _lumiBlock;
    double _HV, _HV_nominal;
    ULong64_t _runNb; 
    ULong64_t _eventNb; 
    UInt_t _timesecond;
    TTree *tree;

    void initialise(TString chamber_name_string , TString year_value ,TString area_name_string, TTree* ); 
    double analysing_dependence(const std::map<std::pair<int, int>, std::tuple<TH1D*, double, int> > histograms, TString var, TFile *);
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> reading_tree_making_maps(TString variable, bool, double, double) ;
    std::string get_channel_name(int rhid);
    TH1D *  trimmed_mean(TH1D * h);
    double ApplyCorrection(double X ,TString correctiontype, double slope );
    double Draw_cumulative_summary_histogram(const std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms ,std::map<int, TH1D*> cumulativeHistogram, TString type, TString var, double binWidth);

    void createCumulativeHistogram(std::map<int, TH1D*>& histMap, 
                               int bin, 
                               const TString var, 
                               const TString condition, 
                               double bin_value, 
                               double binWidth, 
                               const TString chamber_name, 
                               const TString year) ;


};

void GasDependence :: initialise(TString chamber_name_string , TString year_value ,TString area_name_string, TTree *tree_here){

  area_name = area_name_string;
  chamber_name = chamber_name_string;
  year = year_value;
  tree  = tree_here;
  // Set tree branches
  tree->SetBranchAddress("_rhsumQ_equalised_HV_data", &_rhsumQ_equalised_HV_data);
  tree->SetBranchAddress("_rhsumQ", &_rhsumQ);
  tree->SetBranchAddress("_rhsumQ_RAW", &_rhsumQ_RAW);
  tree->SetBranchAddress("_pressure", &_pressure);
  tree->SetBranchAddress("_integratelumi", &_integratelumi);
  tree->SetBranchAddress("_timesecond", &_timesecond);
  tree->SetBranchAddress("_instlumi", &_instlumi);
  tree->SetBranchAddress("_rhid", &_rhid);
  tree->SetBranchAddress("_HV", &_HV);
  tree->SetBranchAddress("_HV_nominal", &_HV_nominal);
  tree->SetBranchAddress("_runNb", &_runNb);
  tree->SetBranchAddress("_eventNb", &_eventNb);
  tree->SetBranchAddress("_lumiBlock", &_lumiBlock);
}


std::string GasDependence :: get_channel_name(int rhid){
     int rhid_reduced = static_cast<int>(std::floor(rhid / 10)) % 1000;
     if (rhid > 2000000) {
         rhid_reduced += 400;
     }
     std::string endcap = (rhid_reduced <= 400) ? "_Endcap1" : "_Endcap2";
     int chamber_nb;
     if (rhid_reduced <= 400) {
         chamber_nb = static_cast<int>(std::floor(rhid_reduced / 10));
     } else {
         chamber_nb = static_cast<int>(std::floor((rhid_reduced - 400) / 10));
     }
     std::string channel_name;
     if (rhid_reduced != 0 && rhid_reduced != 771) {
         channel_name = "chamber" + std::to_string(chamber_nb) +
                        "_layer" + std::to_string(rhid_reduced % 10) +
                        endcap;
     } else {
         channel_name = "undefined";
     }
 
     return channel_name;
 }


std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> GasDependence :: reading_tree_making_maps(TString variable, bool correction, double slope_pressure, double slope_instlumi) {
     // Map to store histograms and mean values
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> result;

    // Map to track cumulative sums and counts for calculating mean values
    std::map<std::pair<int, int>, std::pair<double, int>> binStats;

    double min_var, max_var;
    int nBins;

    if (variable == "pressure") {
        min_var = 940;
        max_var = 990;
        nBins = 50;
    } 
    else if (variable == "instlumi") {
        min_var = 0;
        max_var = 25000;
        nBins = 50;

    }
    else if (variable == "intlumi") {
        min_var = 0;
        max_var = 150;
        nBins = 150;
    }
    else if (variable == "time") {
      int bin_width = 86400;
      if(year=="2016") { 
        min_var  = 1462838400 ; // 10 May 2016 : 00 : 00 : 00
        max_var= 1477871999; // 30 Oct 2016 : 23 : 59 : 59
       }
      else if(year=="2017"){
        min_var  = 1497484800 ; // 15 June 2017 : 00 : 00 : 00 
        max_var= 1510790399; // 15 Nov 2017 : 23 : 59 : 59
      }
      else if(year=="2018"){
        min_var  = 1527206400 ; // 25 April 2018 : 00 : 00 : 00
        max_var =  1540511999; // 25 Oct 2018 : 23 : 59 : 59
      }
      nBins = ((max_var - min_var + 1) / bin_width);
    }


    double binWidth = (max_var - min_var) / nBins;
    int nEntries = tree->GetEntries();

    int skipped_entries = 0;
    int skipped_entries_time = 0;

    int total_inst_entries = 0;
    int accepted_inst_entries = 0;
    int rejected_inst_entries = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // removing the recuperated and wrong gas composition period 
        // Skip invalid entries based on conditions
        if (abs(_HV - _HV_nominal) > 10) continue;
        // This was for testing purposes :
       int Bin;
        if(variable=="pressure")  Bin = static_cast<int>( ((_pressure - min_var) / binWidth)+1);
        else if(variable=="instlumi")  Bin = static_cast<int>( ((_instlumi - min_var) / binWidth) +1);
        else if(variable=="intlumi")  Bin = static_cast<int>( ((_integratelumi - min_var) / binWidth) +1);
        else if(variable=="time")  Bin = static_cast<int>( ((_timesecond - min_var) / binWidth) +1);
        if (Bin < 0 || Bin >= nBins) continue; // Skip out-of-range pressures

        // Create a reduced key for this bin
        std::pair<int, int> reducedKey = std::make_pair(_rhid, Bin);

        // Ensure the histogram exists for the reduced key
        if (result.find(reducedKey) == result.end()) {
            // Create a new histogram for this bin ; name of the histogram depends on whether before or after correction
            // This is to potentially not have same name histogram and craeting memory leak
            TH1D *hist = nullptr;
            if(correction==false){
            hist = new TH1D(
                Form("h_rhid_%d_%s_%d", _rhid, variable.Data(), Bin),
                Form("Histogram for RHID %d, %s Bin %d", _rhid, variable.Data(), Bin),
                3000, 0, 3000
            );
            }
            else if(correction==true){
            hist = new TH1D(
                Form("h_rhid_%d_%s_second_%d", _rhid, variable.Data(), Bin),
                Form("Histogram for RHID %d, %s Bin %d : after correction", _rhid, variable.Data(), Bin),
                3000, 0, 3000
            );
            }
            result[reducedKey] = std::make_tuple(hist, 0.0, 0); // Initialize the histogram and mean
        }

        double equalised_charge =0;
        equalised_charge = _rhsumQ_equalised_HV_data * ApplyCorrection(_pressure ,"pressure", slope_pressure) * ApplyCorrection(_instlumi, "instlumi", slope_instlumi);

        // Fill the histogram with charge data
        std::get<0>(result[reducedKey])->Fill(equalised_charge);

        // Update cumulative sum and count for calculating the mean pressure
        if(variable=="pressure"){        
        binStats[reducedKey].first += _pressure; // Cumulative sum of pressure values
        binStats[reducedKey].second += 1;        // Count of pressure value
        }

        else if(variable=="instlumi"){        
        binStats[reducedKey].first += _instlumi; // Cumulative sum of pressure values
        binStats[reducedKey].second += 1;        // Count of pressure values
        }
        else if(variable=="intlumi"){        
        binStats[reducedKey].first += _integratelumi; // Cumulative sum of pressure values
        binStats[reducedKey].second += 1;        // Count of pressure values
        }
        else if(variable=="time"){        
        binStats[reducedKey].first += _timesecond; // Cumulative sum of pressure values
        binStats[reducedKey].second += 1;        // Count of pressure values
        }

    } // end of going through each tree entry
//    std::cout<<" total entries "<<nEntries<<" skipped entries "<<skipped_entries<<std::endl;
    std::cout<<" total entries inst "<<total_inst_entries<<" skipped entries "<<rejected_inst_entries<<" accepted "<<accepted_inst_entries<<std::endl;

    // Calculate the mean pressure for each bin and store it in the result map
    for (auto& [key, stats] : binStats) {
        double sumPressure = stats.first;
        int count = stats.second;
        double meanPressure = (count > 0) ? (sumPressure / count) : 0.0;

         // Update the mean value in the result map
         std::get<1>(result[key]) = meanPressure;
         std::get<2>(result[key]) = count;
         TH1D* h1 = std::get<0>(result[key]);
         int nentries = h1->Integral();
         double mean = h1->GetMean();
         if(debug) std::cout<<" key "<<key.first<<" bin "<<key.second<< " mean "<<meanPressure<<" histogram value "<<nentries<<" count "<<count<<" mean value "<<mean<<std::endl;
     } // end of reading mean value for each bin

    return result;
}
// Trimming histogram and providing trimmed histogram
TH1D * GasDependence ::  trimmed_mean(TH1D * h){
 
   TH1D *h_trim = (TH1D*) h->Clone();
   TH1D * h_trim_new = (TH1D*) h->Clone();
 
   float integral = 0, trimmean= 0.85;
   // trim the histogram now
              h_trim_new->Reset();
              h_trim_new->ResetStats();
              // Counting the overflow entry also, for trimming purposes
              float normal = h_trim->Integral() + h_trim->GetBinContent(h_trim->GetNbinsX()+1);
              int last_bin = 0;
              for(int it = 1; it<=  h_trim->GetNbinsX() ;it++) {
                if(integral < trimmean * normal){
                  integral+=h_trim->GetBinContent(it);
                  //std::cout<<" old bin entry "<<it<<" entry "<<h_trim->GetBinContent(it)<<std::endl;
                  last_bin =it;
                }
              }
            double new_integral =0;
            int entries_last_bin;
            for(int it=1; it<last_bin; it++){
              new_integral += h_trim->GetBinContent(it);
            }
 
            entries_last_bin = (int) (trimmean*normal - new_integral);
            for(int it=1; it<last_bin ; it++) {
              h_trim_new->SetBinContent(it,h_trim->GetBinContent(it));
              //std::cout<<" new bin entry "<<it<<" entry "<<h_trim_new->GetBinContent(it)<<std::endl;
              h_trim_new->SetBinError(it,h_trim->GetBinError(it));
            }
            h_trim_new->SetBinContent(last_bin, entries_last_bin);
            if(entries_last_bin!=0) {
            h_trim_new->SetBinError(last_bin, h_trim->GetBinError(last_bin) * (entries_last_bin / h_trim->GetBinContent(last_bin)));
            }
            for(int it=last_bin+1; it<=h_trim->GetNbinsX() ; it++) {
              h_trim_new->SetBinContent(it,0);
              h_trim_new->SetBinError(it,0);
            }
            float final_integral = new_integral+entries_last_bin;
            float check_integral = normal * trimmean;
           // std::cout<<" final integral "<<final_integral<<" normal "<<check_integral<<std::endl;
           std::pair<float, float> trimmed_mean_value;
           trimmed_mean_value.first = h_trim_new->GetMean();
           trimmed_mean_value.second = h_trim_new->GetMeanError();
           return h_trim_new;
 }

 
int main(int argc, char *argv[]) {

    TString chamber_name = TString::Format("%s", argv[1]) ;
    TString year = TString::Format("%s", argv[2]) ;
    TString area_name = TString::Format("%s", argv[3]) ;

    bool debug = false;
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/home-n/nrawal/CSCAgeing/Run2_combine_new_selections/"+year+"_updated_after_removal/csc_output_"+year+"_"+chamber_name+"_tree_updated.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return 0 ;
    }
    // Get the tree
    TTree *tree = (TTree *)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Cannot find tree named 'tree' in the file!" << std::endl;
        return  0;
    }

    TString output_file = "./cumulative_plots/"+area_name + "/"+"dataset_output_"+chamber_name+"_"+year+".root";
    TFile*  outputFile = new TFile(output_file,"RECREATE");
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_intlumi_initial;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_time_initial;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_intlumi_final;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_time_final;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_pressure;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_pressure_second;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_instlumi;
    std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms_instlumi_second;
    double slope_pressure = 0, slope_instlumi = 0;
    double slope_pressure_second = 0, slope_instlumi_second = 0;
    double slope_intlumi_initial = 0 , slope_intlumi_final = 0;
    double slope_time_initial = 0 , slope_time_final = 0;

    // Map to store cumulative histograms for each pressure bin
   GasDependence *obj_intlumi_initial = new GasDependence();
   obj_intlumi_initial->initialise(chamber_name, year, area_name, tree);
   histograms_intlumi_initial = obj_intlumi_initial->reading_tree_making_maps("intlumi", false, 0, 0);
   slope_intlumi_initial = obj_intlumi_initial->analysing_dependence(histograms_intlumi_initial, "intlumi_initial", outputFile);
   histograms_intlumi_initial.clear();
   delete obj_intlumi_initial;

    GasDependence *obj_time_initial = new GasDependence();
    obj_time_initial->initialise(chamber_name, year, area_name, tree);
    histograms_time_initial = obj_time_initial->reading_tree_making_maps("time", false, 0, 0);
    slope_time_initial = obj_time_initial->analysing_dependence(histograms_time_initial, "timesecond_initial", outputFile);
    delete obj_time_initial;

   // Pressure corrections now
    GasDependence *obj_pressure = new GasDependence();
    obj_pressure->initialise(chamber_name, year, area_name, tree);
    histograms_pressure = obj_pressure->reading_tree_making_maps("pressure", false, slope_pressure, slope_instlumi);
    if(debug){
    for(const auto &[key, tuple]: histograms_pressure){
     std::pair<int, int> key_here = key;
     double mean = std::get<1>(tuple);
     TH1D * h  = std::get<0>(tuple);
     std::cout<<" key "<<key.first<<" second "<<key.second<<std::endl;
     std::cout<<" mean "<<mean<<" entries "<<h->Integral()<<std::endl;
    }
    }
    slope_pressure = obj_pressure->analysing_dependence(histograms_pressure, "pressure", outputFile);
    histograms_pressure.clear();
    delete obj_pressure;

    GasDependence *obj_pressure_second = new GasDependence();
    obj_pressure_second->initialise(chamber_name, year, area_name, tree);
    histograms_pressure_second = obj_pressure_second->reading_tree_making_maps("pressure", true, slope_pressure, slope_instlumi);
    slope_pressure_second = obj_pressure_second->analysing_dependence(histograms_pressure_second, "pressure_second", outputFile);
    histograms_pressure_second.clear();
    delete obj_pressure_second;

    // instlumi dependence 
    GasDependence *obj_instlumi = new GasDependence();
    obj_instlumi->initialise(chamber_name, year, area_name, tree);
    histograms_instlumi = obj_instlumi->reading_tree_making_maps("instlumi", false, slope_pressure, slope_instlumi);
    slope_instlumi = obj_instlumi->analysing_dependence(histograms_instlumi, "instlumi", outputFile);
    histograms_instlumi.clear();
    delete obj_instlumi;

    GasDependence *obj_instlumi_second = new GasDependence();
    obj_instlumi_second->initialise(chamber_name, year, area_name, tree);
    histograms_instlumi_second = obj_instlumi_second->reading_tree_making_maps("instlumi", true, slope_pressure, slope_instlumi);
    slope_instlumi_second = obj_instlumi_second->analysing_dependence(histograms_instlumi_second, "instlumi_second", outputFile);
    histograms_instlumi_second.clear();
    delete obj_instlumi_second;

///    // final intlumi dependence
    GasDependence *obj_intlumi_final = new GasDependence();
    obj_intlumi_final->initialise(chamber_name, year, area_name, tree);
    histograms_intlumi_final = obj_intlumi_final->reading_tree_making_maps("intlumi", true, slope_pressure, slope_instlumi);
    slope_intlumi_final = obj_intlumi_final->analysing_dependence(histograms_intlumi_final, "intlumi_final", outputFile);
    histograms_intlumi_final.clear();
    delete obj_intlumi_final;

   GasDependence *obj_time_final = new GasDependence();
   obj_time_final->initialise(chamber_name, year, area_name, tree);
   histograms_time_final = obj_time_final->reading_tree_making_maps("time", true, slope_pressure, slope_instlumi);
   slope_time_final = obj_time_final->analysing_dependence(histograms_time_final, "timesecond_final", outputFile);
   delete obj_time_final;

 
    file->Close();
    outputFile->Close();

    return 0;
}
double GasDependence :: analysing_dependence(const std::map<std::pair<int, int>, std::tuple<TH1D*, double, int> > histograms, TString var, TFile * output_file){
///    std::map<int, TH1D*> cumulativeHistograms;
///    std::map<int, TH1D*> cumulativeHistograms_plus;
///    std::map<int, TH1D*> cumulativeHistograms_minus;
    std::map<std::string, std::map<int, TH1D*>> cumulativeHistograms;
    bool single_channel_check=false;

     TString dir_name = var;
     TDirectoryFile *dir_var =  (TDirectoryFile*) output_file->mkdir(dir_name);
     dir_var->cd();
     double min_var, max_var;
     int nBins;
     if(var=="pressure" || var=="pressure_second") {
       min_var = 940; 
       max_var = 990;
       nBins = 50; 
    
     }
    else if(var=="instlumi" || var=="instlumi_second") {
       min_var = 0; 
       max_var = 25000;
       nBins = 50; 
     }
    else if(var=="intlumi_initial" || var=="intlumi_final") {
       min_var = 0; 
       max_var = 150;
       nBins = 150; 
     }
    else if (var == "timesecond_initial" || var=="timesecond_final") {
      int bin_width = 86400;
      if(year=="2016") { 
        min_var  = 1462838400 ; // 10 May 2016 : 00 : 00 : 00
        max_var= 1477871999; // 30 Oct 2016 : 23 : 59 : 59
       }
      else if(year=="2017"){
        min_var  = 1497484800 ; // 15 June 2017 : 00 : 00 : 00 
        max_var= 1510790399; // 15 Nov 2017 : 23 : 59 : 59
      }
      else if(year=="2018"){
        min_var  = 1527206400 ; // 25 April 2018 : 00 : 00 : 00
        max_var =  1540511999; // 25 Oct 2018 : 23 : 59 : 59
      }
      nBins = ((max_var - min_var + 1) / bin_width);
    }


    double binWidth = (max_var - min_var) / nBins;
    //    int rhid_value = 2121921;
    //    TString rhid_string = "Chamber 19 layer 2 Endcap2";
    
    // Individual trimmed mean TGraphs for individual channel 
     std::map<int, TGraphAsymmErrors*> h_summary ; 
     for (const auto& [key, value] : histograms){
          // one new summary Histogram for each rhid
          int rhid = key.first;
          std::string channel_name_string =  get_channel_name(rhid);
          TString channel_string(channel_name_string);
          if(h_summary.find(rhid)==h_summary.end()){
            TGraphAsymmErrors *hist = new TGraphAsymmErrors(); 
            TString graph_name = TString::Format("dataset_trimmed_%s_%s", channel_string.Data(), var.Data());
            TString graph_title = chamber_name+" : "+year+" : "+channel_string+" : "+var;

            hist->SetName(graph_name);
            hist->SetTitle(graph_title);
            hist->GetYaxis()->SetTitle("Trimmed mean charge");
            if(var=="pressure" || var=="pressure_second")  hist->GetXaxis()->SetTitle("Pressure (hPa)");
            if(var=="instlumi" || var=="instlumi_second")  hist->GetXaxis()->SetTitle("Instlumi (*10^{30} cm^{2} s^{-1})");
            if(var=="intlumi_initial" || var=="intlumi_final")  hist->GetXaxis()->SetTitle("Integrated lumi  (fb^{-1})");
            if(var=="timesecond_initial" || var=="timesecond_final"){
                hist->GetXaxis()->SetTimeDisplay(1);
                hist->GetXaxis()->SetLabelSize(0.02);
                hist->GetXaxis()->SetTimeFormat("%Y/%m/%d"); 
                hist->GetXaxis()->SetTitle("time");
            }

            gStyle->SetOptStat(111112211);
            gStyle->SetOptFit(1111);

            h_summary[rhid] = hist;
            
          }
     }
     std::map<int, int> n_points;
     double x_var_value;
     double y_charge_value;
     double y_charge_err;

   
     std::map<TString, int> map_chamber = 
     { {"ME11a", 36},  {"ME11b", 36}, {"ME12HV1", 36}, {"ME12HV2", 36}, {"ME12HV3", 36},
       {"ME13HV1", 36}, {"ME13HV2", 36}, {"ME13HV3", 36},
       {"ME21HV1", 18}, {"ME21HV2", 18}, {"ME21HV3", 18},
       {"ME31HV1", 18}, {"ME31HV2", 18}, {"ME31HV3", 18},
       {"ME41HV1", 18}, {"ME41HV2", 18}, {"ME41HV3", 18},
       {"ME22HV1", 36}, {"ME22HV2", 36}, {"ME22HV3", 36},{"ME22HV4", 36},{"ME22HV5", 36},
       {"ME32HV1", 36}, {"ME32HV2", 36}, {"ME32HV3", 36},{"ME32HV4", 36},{"ME32HV5", 36},
       {"ME42HV1", 36}, {"ME42HV2", 36}, {"ME42HV3", 36},{"ME42HV4", 36},{"ME42HV5", 36},
     };
    int upper_nb = map_chamber[chamber_name] ;  
    int lower_nb = upper_nb/2;
    // Trim histograms and make cumulative distribution for each bin of the variable
     for (const auto& [key, value] : histograms){
         // key.first will be chamber name; key.second will be BinNb
         TH1D* h = std::get<0>(value);
         double var_value = std::get<1>(value);
         int entries = std::get<2>(value);
         int rhid = key.first;
         int Bin = key.second;
         int bin_value = (key.second-1) *binWidth + min_var;
         TString Bin_string = TString::Format("%d",Bin);

         // finding the endcap 
         int rhid_reduced = static_cast<int>(std::floor(rhid / 10)) % 1000;
          if (rhid > 2000000) {
           rhid_reduced += 400;
          }
          TString endcap = (rhid_reduced <= 400) ? "positive" : "negative";

 
         TH1D *h_trimmed = trimmed_mean(h);
         if(debug) std::cout << "Trimming histogram for RHID: " << key.first
                   << ", "<<var<<" Bin: " << key.second << "bin value "<<bin_value<< " hist entries "<<h->Integral()<< " after trim "<<h_trimmed->Integral()<<" value  "<<var_value<<" entries "<<entries<<std::endl;

         // _dataset_pressure_corrected_trimmean_chamber35_layer4_Endcap2vs_pressure
         // Converting rhid into string //2143411
         std::string channel_name_string =  get_channel_name(rhid);
         int chamber_nb;
         int layer_nb;
         if (rhid_reduced <= 400) {
          chamber_nb = static_cast<int>(std::floor(rhid_reduced / 10));
         } else {
          chamber_nb = static_cast<int>(std::floor((rhid_reduced - 400) / 10));
         }
         layer_nb = rhid_reduced % 10;
         TString channel_string(channel_name_string);
  
         h_trimmed->SetName("dataset_trimmed_"+channel_string+"_bin_"+Bin_string+"_vs_"+var);
         if(h_trimmed!=NULL || h_trimmed->GetEntries()>0) {
           if(n_points.find(rhid)==n_points.end()){
             n_points[rhid]= 0;
           }
          // TO avoid filling zero values to graph

           if(h_trimmed->Integral()<=0) { 
           continue; 
           }

          // To test if I use center of the bin as my gravity 
          // x_var_value = (Bin-1) * binWidth+ min_var+ 0.5 * binWidth;
          // double low_x_err = 0.5 * binWidth;
          // double up_x_err = 0.5 * binWidth;
          //  h_summary[rhid]->SetPoint(n_points[rhid], x_var_value, y_charge_value);
           
           x_var_value = var_value ;
           double low_x_err = x_var_value - bin_value ;
           double up_x_err = bin_value+binWidth - x_var_value;

           y_charge_value = h_trimmed->GetMean();
           y_charge_err = h_trimmed->GetMeanError();

           h_summary[rhid]->SetPoint(n_points[rhid], x_var_value, y_charge_value);
//           h_summary[rhid]->SetPointError(n_points[rhid], low_x_err, up_x_err, y_charge_err, y_charge_err);
           h_summary[rhid]->SetPointError(n_points[rhid], 0,0 , y_charge_err, y_charge_err);
           n_points[rhid] = n_points[rhid]+1;
         }
        // Define condition categories dynamically
        std::vector<std::string> conditions = {"all"}; // Always include "all"
        if (endcap == "positive") conditions.push_back("plus");
        if (endcap == "negative") conditions.push_back("minus");
        if (chamber_nb % 2 == 0) conditions.push_back("even_chambers"); // Example: Add more conditions
        if (chamber_nb % 2 != 0) conditions.push_back("odd_chambers"); // Example: Add more conditions
        if (layer_nb % 2 != 0) conditions.push_back("odd_layers"); // Example: Add more conditions
        if (layer_nb % 2 == 0) conditions.push_back("even_layers"); // Example: Add more conditions
        if (_lumiBlock % 2 == 0) conditions.push_back("even_lumiBlock"); // Example: Add more conditions
        if (_lumiBlock % 2 != 0) conditions.push_back("odd_lumiBlock"); // Example: Add more conditions
   
        if(chamber_nb >=lower_nb+1 && chamber_nb <=upper_nb && endcap=="positive") conditions.push_back("upper_plus");  
        if(chamber_nb >=lower_nb+1 && chamber_nb <=upper_nb && endcap=="negative") conditions.push_back("upper_minus");  
        if(chamber_nb >=1 && chamber_nb <=lower_nb && endcap =="positive") conditions.push_back("lower_plus");  
        if(chamber_nb >=1 && chamber_nb <=lower_nb && endcap =="negative") conditions.push_back("lower_minus");  
        for (const auto& cond : conditions) {
            createCumulativeHistogram(cumulativeHistograms[cond], Bin, var, cond, bin_value, binWidth, chamber_name, year);
            cumulativeHistograms[cond][Bin]->Add(h_trimmed);
        }

         delete h;
         delete h_trimmed;
    }
    // end of making cumulative Histograms

    // Drawing individual histogram and write them to root file
    for(const auto &[rhid, hist] : h_summary){
       if(hist){
         TF1 *expFit2 = nullptr;

         if(var=="pressure" || var=="pressure_second" || var=="instlumi" || var=="instlumi_second"){

          // Founding first and last point of the edges         
           double fitlowedge = 0, fithighedge = 0; // Initialize edges
           int nPoints = hist->GetN(); // Get the number of points in the graph
           bool foundFirst = false;
           // Loop through all points to find the first and last valid points
           for (int i = 0; i < nPoints; i++) {
               double x, y;
               hist->GetPoint(i, x, y);
               if (y > 0) { // Replace with your condition (e.g., y > threshold)
                   if (!foundFirst) {
                       fitlowedge = x-binWidth; // First valid x-coordinate
                       foundFirst = true;
                   }
                   fithighedge = x+binWidth; // Continuously update with the last valid x-coordinate
               }
           }
           
           std::cout << "Fit range: [" << fitlowedge << ", " << fithighedge << "]" << std::endl;
               if(var=="pressure" || var=="pressure_second"){
               expFit2 = new TF1("expFit2", "exp([0]) * exp([1]*(x-967))", fitlowedge, fithighedge);
               //expFit2 = new TF1("expFit2", "exp([0]) * exp([1]*(x-967))", min_var, max_var);
               expFit2->SetParameters(5, -0.005); // Initial guesses
               } 
               if(var=="pressure_second"){
               expFit2 = new TF1("expFit2", "exp([0]) * exp([1]*(x-967))", fitlowedge, fithighedge);
               } 

               else if(var=="instlumi" || var=="instlumi_second"){
               //expFit2 = new TF1("expFit2", "exp([0]) * exp([1]*(x-10000))", min_var, max_var);
               expFit2 = new TF1("expFit2", "exp([0]) * exp([1]*(x-10000))", fitlowedge, fithighedge);
               }
               hist->Fit(expFit2, "R");
          }
         delete expFit2;
         hist->Write();
       }
     }
    std::map<std::string, double> slope_values;

    for (const auto& [cond, histMap] : cumulativeHistograms) {
    slope_values[cond] = this->Draw_cumulative_summary_histogram(histograms, histMap, cond, var, binWidth);
    }
// std::vector<double> slope_value_vector;
// slope_value_vector.push_back(slope_value, slope_value_plus, slope_value_minus);
// return slope_value_vector;
 double slope_value = slope_values["all"];
 return slope_value;
}
void GasDependence ::createCumulativeHistogram(std::map<int, TH1D*>& histMap, 
                               int bin, 
                               const TString var, 
                               const TString condition, 
                               double bin_value, 
                               double binWidth, 
                               const TString chamber_name, 
                               const TString year) {
    if (histMap.find(bin) == histMap.end()) {
        TString histName = Form("cumulative_%s_%d_%s", var.Data(), bin, condition.Data());
        TString histTitle = Form("cumulative: %s: %s : %.0f <= %s < %.0f : %s", 
                                 chamber_name.Data(), year.Data(), bin_value, var.Data(), bin_value + binWidth, condition.Data());

        histMap[bin] = new TH1D(histName, histTitle, 3000, 0, 3000); // Example: adjust binning as needed
        histMap[bin]->Sumw2();
        histMap[bin]->GetXaxis()->SetTitle("charge (ADC)");
        histMap[bin]->GetYaxis()->SetTitle("Nb. of entries");
        histMap[bin]->SetTitle(histTitle); 
    }
}
// To draw summary histogram from cumulativeHistogrmas
double GasDependence :: Draw_cumulative_summary_histogram(const std::map<std::pair<int, int>, std::tuple<TH1D*, double, int>> histograms , std::map<int, TH1D*> cumulativeHistograms, TString type_cumulative, TString var, double binWidth) {
  double slope_value = 0;
    TGraphAsymmErrors* summaryHistogram = new TGraphAsymmErrors();
    TString summary_title;
    if(type_cumulative=="all") { 
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var;
    }
    else if(type_cumulative=="plus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : +endcap";
    }
    else if(type_cumulative=="minus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : -endcap";
    }
    else if(type_cumulative=="odd_chambers"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : odd chambers";
    }
    else if(type_cumulative=="even_chambers"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : even chambers";
    }
    else if(type_cumulative=="lower_plus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : lower +endcap";
    }
    else if(type_cumulative=="lower_minus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : lower -endcap";
    }
    else if(type_cumulative=="upper_plus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : upper +endcap";
    }
    else if(type_cumulative=="upper_minus"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : upper -endcap";
    }
    else if(type_cumulative=="odd_layers"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : odd layers";
    }
    else if(type_cumulative=="even_layers"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : even layers";
    }
    else if(type_cumulative=="odd_lumiBlock"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : odd lumiBlock";
    }
    else if(type_cumulative=="even_lumiBlock"){
     summary_title = "Cumulative : "+chamber_name+" : "+year+ " : "+var+" : even lumiBlock";
    }

   summaryHistogram->SetTitle(summary_title);
   // Finding the mean x value for summary Histogram
   std::map<int, double> sum_var_bin;
   std::map<int, double> mean_var_bin;
   std::map<int, double> entries_var_bin;

   std::map<int, int > entries_total;
   for (const auto& [key, value] : histograms){
         int Bin = key.second;
         double var_value = std::get<1>(value);
         int rhid = key.first;
         int entries_total = std::get<2>(value);
      
         if(mean_var_bin.find(Bin)==mean_var_bin.end()){
           sum_var_bin[Bin] = 0; 
           mean_var_bin[Bin] = 0; 
           entries_var_bin[Bin] = 0; 
          } // end of initialising mean, sum and entries
         sum_var_bin[Bin] = sum_var_bin[Bin]+(var_value * entries_total); 
         entries_var_bin[Bin] = entries_var_bin[Bin]+entries_total; 
     }
     for(auto &[key, value] : sum_var_bin){
         mean_var_bin[key] = value/entries_var_bin[key];
     }
    int n_points_summary = 0;
   
    // Reading cumulative Histogram into a summary plot 
    for (auto& pair : cumulativeHistograms) {
        int pressureBin = pair.first;
        TH1D* cumulativeHistogram = pair.second;
         if(cumulativeHistogram->Integral()<=0) { 
           continue; 
          }
        // Calculate the mean of the cumulative histogram
        double mean = cumulativeHistogram->GetMean();
        double mean_error = cumulativeHistogram->GetMeanError();

        summaryHistogram->SetPoint(n_points_summary, mean_var_bin[pressureBin], mean);
        summaryHistogram->SetPointError(n_points_summary, 0, 0, mean_error, mean_error);
        if(debug) std::cout<<" point "<<n_points_summary<<" mean "<<mean_var_bin[pressureBin]<<" bin actual "<<pressureBin<<" mean value "<<mean<<std::endl;
        n_points_summary++;

        // Write the cumulative histogram to the output file
        TCanvas *c1 = new TCanvas();
        c1->cd();
        gStyle->SetOptStat(111112211);
        cumulativeHistogram->Draw();

        TString saving_name = "cumulative_plots/"+area_name+"/all_channels/"+chamber_name+"/"+var+TString::Format("/cumulative_charge_%d_%s_",pressureBin, var.Data())+type_cumulative+".pdf";
        //c1->SaveAs(saving_name);
        //cumulativeHistogram->Write();
        delete c1;
        delete cumulativeHistogram; // Clean up cumulative histograms
    }
     summaryHistogram->GetYaxis()->SetRangeUser(250,600);
     summaryHistogram->GetYaxis()->SetTitle("Trimmed mean charge");
      if(var=="pressure" || var=="pressure_second")  {
        summaryHistogram->GetXaxis()->SetTitle("Pressure (hPa)");
        summaryHistogram->GetXaxis()->SetRangeUser(940, 990);
      }
      if(var=="instlumi" || var=="instlumi_second"){
        summaryHistogram->GetXaxis()->SetTitle("Instlumi (*10^{30} cm^{2} s^{-1})");
        summaryHistogram->GetXaxis()->SetRangeUser(0,23000);
      }
      if(var=="intlumi_initial" || var=="intlumi_final") {
        summaryHistogram->GetXaxis()->SetRangeUser(0,145);
        summaryHistogram->GetXaxis()->SetTitle("Integrated lumi  (fb^{-1})");
      }
      if(var=="timesecond_initial" || var=="timesecond_final"){
          summaryHistogram->GetXaxis()->SetTimeDisplay(1);
          summaryHistogram->GetXaxis()->SetLabelSize(0.02);
          summaryHistogram->GetXaxis()->SetTimeFormat("%Y/%m/%d"); 
          summaryHistogram->GetXaxis()->SetTitle("time");
       }

     if(debug) std::cout<<" number of points "<<summaryHistogram->GetN();
     summaryHistogram->SetMarkerStyle(20);
     summaryHistogram->SetMarkerSize(0.5);

     TCanvas *c = new TCanvas();
     gStyle->SetOptStat(111112211);
     gStyle->SetOptFit(1111);
    
     if(debug) std::cout<<" before plotting and fitting summary"<<std::endl; 
    
     if(var=="pressure" || var=="pressure_second" || var=="instlumi" || var=="instlumi_second"){
      double fitlowedge = 0, fithighedge = 0; // Initialize edges
      int nPoints = summaryHistogram->GetN(); // Get the number of points in the graph
      bool foundFirst = false;
      
      // Loop through all points to find the first and last valid points
      for (int i = 0; i < nPoints; i++) {
          double x, y;
          summaryHistogram->GetPoint(i, x, y);
      
          if (y > 0) { // Replace with your condition (e.g., y > threshold)
              if (!foundFirst) {
                  fitlowedge = x-binWidth; // First valid x-coordinate
                  foundFirst = true;
              }
              fithighedge = x+binWidth; // Continuously update with the last valid x-coordinate
          }
      }
     if(debug) std::cout << "Fit range: [" << fitlowedge << ", " << fithighedge << "]" << std::endl;
     TF1 *expFit = nullptr;
     if(var=="pressure"){
     expFit = new TF1("expFit", "exp([0]) * exp([1]*(x-967))", fitlowedge, fithighedge);
     expFit->SetParameters(5, -0.005); // Initial guesses
     summaryHistogram->Fit(expFit, "R");
     slope_value = expFit->GetParameter(1);
     } 
     else if(var=="pressure_second"){
     expFit = new TF1("expFit", "exp([0]) * exp([1]*(x-967))", fitlowedge, fithighedge);
     summaryHistogram->Fit(expFit, "R");
     slope_value = expFit->GetParameter(1);
     } 
     else if(var=="instlumi" || var=="instlumi_second"){
     expFit = new TF1("expFit", "exp([0]) * exp([1]*(x-1000))", fitlowedge, fithighedge);
     //expFit = new TF1("expFit", "exp([0]) * exp([1]*(x-10000))", min_var, max_var);
     summaryHistogram->Fit(expFit, "R");
     slope_value = expFit->GetParameter(1);
     }
    }
     c->cd();
     summaryHistogram->Draw("AP");
     gPad->Update();
     c->Update();
     c->Update();
 
     c->SaveAs("cumulative_plots/"+area_name+"/all_channels/"+chamber_name+"/"+ var+"/"+chamber_name+"_cumulative_trimmed_mean_"+var+"_"+type_cumulative+".pdf");
     summaryHistogram->SetName("dataset_trimmed_"+chamber_name+"_allgoodchannelsvs_"+var+"_"+type_cumulative);

     summaryHistogram->Write();
     delete summaryHistogram;
     return slope_value;
}

double GasDependence :: ApplyCorrection(double X ,TString correctiontype,  double slope ){
  double refvalue = 0;
  if(correctiontype=="pressure"){
    refvalue =967 ;
  }
  else if (correctiontype=="instlumi"){
    refvalue = 10000 ;
  }
 double thecorr =exp(slope*(refvalue-X));
 //double thecorr = p1*(refvalue-X);
 return thecorr;
  }
