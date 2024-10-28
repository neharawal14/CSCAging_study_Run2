#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include "tree_class.h"
// Create a unique key for rhid and runNb pair
// Create a unique key for rhid and runNb pair
 std::string make_key(int rhid, int runNb) {
     return std::to_string(rhid) + "_" + std::to_string(runNb);
     }
double equalised_charge(double raw_charge, double HV_read,  int stationring);

int main(int argc, char *argv[]) {

    std::string chamber_name = argv[1];
    std::string year = argv[2];
    std::string chamber_value = argv[3];
    // Use a map with only averageHV as the value
    // Use a map with only averageHV as the value
    std::unordered_map<std::string, float> hv_map;
    std::string input_path = "/eos/home-n/nrawal/CSCAgeing/test_files/debug_HV_values/";
if(chamber_name!="ME11a" && chamber_name!="ME11b")
 {
    std::string input_file_plus = "average_hv_by_rhid_run_"+chamber_name+"_plus_"+year+".txt";
    std::string input_file_minus = "average_hv_by_rhid_run_"+chamber_name+"_minus_"+year+".txt";
    std::string input_file_string_plus = input_path+input_file_plus;
    std::string input_file_string_minus = input_path+input_file_minus;
    std::cout<<" input file "<<input_file_string_plus<<std::endl;
    std::ifstream infile_plus(input_file_string_plus);
    std::ifstream infile_minus(input_file_string_minus);
    int rhid, runNb;
    float avg_hv;
    int rhid_minus, runNb_minus;
    float avg_hv_minus;


    if (!infile_plus) {
        std::cerr << "Error opening averageHV plus.txt" << std::endl;
        return 1;
    }
    if (!infile_minus) {
        std::cerr << "Error opening averageHV minus.txt" << std::endl;
        return 1;
    }

    // Read the text file and populate the map
    std::string line;
    // Skip the header line
    std::getline(infile_plus, line);
    // Read each subsequent line
    while (std::getline(infile_plus, line)) {
       std::istringstream iss(line);
       std::string rhid_str, runNb_str, avg_hv_str;
       // Read values separated by commas
       if (!std::getline(iss, rhid_str, ',') ||
            !std::getline(iss, runNb_str, ',') ||
            !std::getline(iss, avg_hv_str, ',')) {
            std::cerr << "Error: Malformed line in file." << std::endl;
            continue;
        }

        // Convert string values to integers and float
        rhid = std::stoi(rhid_str);
        runNb = std::stoi(runNb_str);
        avg_hv = std::stof(avg_hv_str);

//        std::cout<<" hid "<<rhid<<" run "<<runNb<<std::endl;
        // Store in the map
        hv_map[make_key(rhid, runNb)] = avg_hv;
    }


    // Read the text file and populate the map
    std::string line_minus;
    // Skip the header line
    std::getline(infile_minus, line_minus);
    // Read each subsequent line
    while (std::getline(infile_minus, line_minus)) {
       std::istringstream iss(line_minus);
       std::string rhid_str_minus, runNb_str_minus, avg_hv_str_minus;
       // Read values separated by commas
       if (!std::getline(iss, rhid_str_minus, ',') ||
            !std::getline(iss, runNb_str_minus, ',') ||
            !std::getline(iss, avg_hv_str_minus, ',')) {
            std::cerr << "Error: Malformed line in file." << std::endl;
            continue;
        }

        // Convert string values to integers and float
        rhid_minus = std::stoi(rhid_str_minus);
        runNb_minus = std::stoi(runNb_str_minus);
        avg_hv_minus = std::stof(avg_hv_str_minus);

//        std::cout<<" hid "<<rhid<<" run "<<runNb<<std::endl;
        // Store in the map
        hv_map[make_key(rhid_minus, runNb_minus)] = avg_hv_minus;
    }

}

else{
    std::string input_file = "average_hv_by_rhid_run_"+chamber_name+"_"+year+".txt";
    std::string input_file_string = input_path+input_file;
    std::ifstream infile(input_file_string);
    int rhid, runNb;
    float avg_hv;
    int rhid_minus, runNb_minus;
    float avg_hv_minus;


    if (!infile) {
        std::cerr << "Error opening averageHV plus.txt" << std::endl;
        return 1;
    }

    // Read the text file and populate the map
    std::string line;
    // Skip the header line
    std::getline(infile, line);
    // Read each subsequent line
    while (std::getline(infile, line)) {
       std::istringstream iss(line);
       std::string rhid_str, runNb_str, avg_hv_str;
       // Read values separated by commas
       if (!std::getline(iss, rhid_str, ',') ||
            !std::getline(iss, runNb_str, ',') ||
            !std::getline(iss, avg_hv_str, ',')) {
            std::cerr << "Error: Malformed line in file." << std::endl;
            continue;
        }

        // Convert string values to integers and float
        rhid = std::stoi(rhid_str);
        runNb = std::stoi(runNb_str);
        avg_hv = std::stof(avg_hv_str);

//        std::cout<<" hid "<<rhid<<" run "<<runNb<<std::endl;
        // Store in the map
        hv_map[make_key(rhid, runNb)] = avg_hv;
    }
}
//      for (auto i : hv_map){
//      std::cout<<i.first<<" second runNb "<<i.second<<std::endl;
//    }

    // Open the ROOT file and access the tree
    TString path = "/eos/home-n/nrawal/CSCAgeing/Run2_combine_new_selections/"+year+"/";
    TString new_path = "/eos/home-n/nrawal/CSCAgeing/Run2_combine_new_selections/"+year+"_updated/";
    TString input_file_name  = "csc_output_"+year+"_"+chamber_value+"_tree.root";
    TString output_file_name = "csc_output_"+year+"_"+chamber_value+"_tree_HVupdated.root";

    TString input_file = path+input_file_name;
    TString output_file = new_path+output_file_name;
//    std::cout<<" file will be opening "<<path<<std::endl;

    TreeStructure *tree = new TreeStructure();
    tree->initialise(input_file, output_file);
    std::cout<<" done iwth initialised my tree"<<std::endl;

//    Double_t _HV_new;
//    Int_t _rhid_new;
//    ULong64_t  _runNb_new;
//    // Access branches
//    Double_t _HV;
//    Int_t _rhid;
//    ULong64_t  _runNb;
//    tree->SetBranchAddress("_HV", &_HV);
//    tree->SetBranchAddress("_rhid", &_rhid);
//    tree->SetBranchAddress("_runNb", &_runNb);
//
//    tree_new.Branch("_HV",&_HV_new,"_HV_new/D");
//    tree_new.Branch("_rhid",&_rhid_new,"_rhid_new/I");
//    tree_new.Branch("_runNb",&_runNb_new,"_runNb_new/l");

    // Loop over all entries in the tree
    //Long64_t nEntries = tree->GetEntries();
    //Long64_t nEntries = 1000;
    Double_t HV_value;
    Long64_t entries = tree->tree->GetEntries();
    std::cout<<" total entries "<<entries<<std::endl;
    tree->Setup_new_tree();
    for (Long64_t i = 0; i < entries; ++i) {
     // std::cout<<" getting first entry "<<std::endl;

      //std::cout<<" entries "<<tree->tree->GetEntries()<<std::endl;
      tree->tree->GetEntry(i);
    //std::cout<<" read tree entry  pass muons   "<<std::endl;
    tree->new_passZmumusel = tree->old_passZmumusel;
    //std::cout<<" after pass muons   "<<std::endl;
    tree->new_passisomuondzdxy = tree->old_passisomuondzdxy;
    tree->new_eventNb = tree->old_eventNb;
    tree->new_runNb  = tree->old_runNb;
    tree->new_lumiBlock = tree->old_lumiBlock;
    tree->new_rhid    = tree->old_rhid;

    //std::cout<<" after rhid   "<<std::endl;
    tree->new_stationring   = tree->old_stationring;
    tree->new_current = 0;
    tree->new_pressure = tree->old_pressure;
    tree->new_temperature = tree->old_temperature;
    tree->new_timesecond = tree->old_timesecond;
    tree->new_n_PV = tree->old_n_PV;
    tree->new_bunchcrossing = tree->old_bunchcrossing;
    tree->new_etamuon = tree->old_etamuon;
    tree->new_phimuon = tree->old_phimuon;
    tree->new_ptmuon = tree->old_ptmuon;

    //std::cout<<" after pt muons   "<<std::endl;
    tree->new_z_pt = tree->old_z_pt;
    tree->new_z_eta = tree->old_z_eta;
    tree->new_z_phi = tree->old_z_phi;
    tree->new_z_mass = tree->old_z_mass;
    tree->new_isolation1 = tree->old_isolation1;
    tree->new_isolation2 = tree->old_isolation2;
    tree->new_iso_PF_first = tree->old_iso_PF_first;
    tree->new_iso_PF_second = tree->old_iso_PF_second;
    tree->new_instlumi = tree->old_instlumi;
    tree->new_integratelumi = tree->old_integratelumi;
    tree->new_rhsumQ = tree->old_rhsumQ;
    tree->new_rhsumQ_RAW = tree->old_rhsumQ_RAW;


   //   std::cout<<" reading key  "<<std::endl;
    std::string key = make_key(tree->old_rhid,tree->old_runNb);
   //   std::cout<<" after reading key  "<<std::endl;
    //std::cout<<" entry "<<i<<" runNb "<<_runNb<<" rhid "<<_rhid<<std::endl;
    // Check if the key exists in the map and update HV
    if (hv_map.find(key) != hv_map.end()) {
            HV_value = hv_map[key];  // Update HV with averageHV from the map
           // std::cout<<" entry "<<i<<" updated HV "<<HV_value<<" runNb "<<tree->old_runNb<<" rhid "<<tree->old_rhid<<std::endl;
            tree->new_HV_nominal = tree->old_HV;
            tree->new_HV = HV_value;
            // equalising charge according to HV in database
            tree->new_rhsumQ_equalised_HV_data = equalised_charge(tree->old_rhsumQ_RAW, HV_value, tree->old_stationring);
            tree->tree_new->Fill();  // Update the tree with the new HV value
        }

        else{
            std::cout<<" entry not exist "<<i<<" updated HV "<<tree->old_HV<<" runNb "<<tree->old_runNb<<" rhid "<<tree->old_rhid<<std::endl;
//             std::cerr << "Error: not finding the run " << std::endl;
//             return 1;
           
        }
    }

    // Write the updated tree to the ROOT file
    tree->WriteNew();

    return 0;
}

double equalised_charge(double raw_charge, double HV_read,  int stationring){
    double dHV_;
    double const_B_;
    if( (stationring==11|| stationring ==14) ) {const_B_ = 6.26e-3; dHV_ = HV_read - 2900;}
    else if( (stationring==21|| stationring ==31 || stationring ==41) ) {const_B_ = 5.193e-3;dHV_ = HV_read - 3600;}
    else  {const_B_ = 5.463e-3; dHV_ = HV_read - 3600;}
   
    // equalised charge from voltage in data 
    double charge_equalised = raw_charge * exp(-1* const_B_* dHV_ ) ;
    return charge_equalised ;

}


