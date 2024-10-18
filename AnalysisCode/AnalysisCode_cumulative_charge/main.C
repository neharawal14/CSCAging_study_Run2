#include "CumulativePlots.h"
#include "individual_channels.h"
#include "MeanPlots.h"
#include "SystematicRemoval.h"
int main(int argc, char *argv[]){

  std::vector<TString> chamber ;
  for(int i=0; i<argc-1 ; i++){
    TString chamber_name = TString::Format("%s", argv[i+1]);
    std::cout<<" arg "<<chamber_name<<std::endl;	
    chamber.push_back(chamber_name);
  }
  CumulativePlot plot_obj;
  individual_channels plot_obj_indie;
	std::cout<<" started program"<<std::endl;
  std::pair<float, float> mean_values_pressure_pair_2016;
  std::pair<float, float> mean_values_pressure_pair_2017;
  std::pair<float, float> mean_values_pressure_pair_2018;
  std::pair<float, float> mean_values_instlumi_pair_2016;
  std::pair<float, float> mean_values_instlumi_pair_2017;
  std::pair<float, float> mean_values_instlumi_pair_2018;

  std::vector<std::pair<float, float>> mean_values_pressure_2016;
  std::vector<std::pair<float, float>> mean_values_pressure_2017;
  std::vector<std::pair<float, float>> mean_values_pressure_2018;
  std::vector<std::pair<float, float>> mean_values_instlumi_2016;
  std::vector<std::pair<float, float>> mean_values_instlumi_2017;
  std::vector<std::pair<float, float>> mean_values_instlumi_2018;

 // std::vector<std::pair<float, float>> mean_values_intlumi;

    	TString input_path = "/eos/home-n/nrawal/CSCAgeing/2024_new_plots/new_ntuples_result_run2/";
//      outf_dataset_pressure_corrected__ME11a_output_run2.root
	for(int i=0; i<chamber.size();i++){
    	std::cout<<" about to initialise the chamber progra:"<<std::endl;
    	TString input_file_name = input_path+"outf_dataset_pressure_corrected__"+chamber[i]+"_output_run2.root";
    	TString saving_path = "/afs/cern.ch/work/n/nrawal/CSCAgeing_code_study/results_AnalysisCode/New_selections_run2/cumulative_plots/";

    	std::cout<<" about to initialise the chamber progra:"<<std::endl;
  	  plot_obj.initialise(chamber[i], input_file_name, saving_path);
  	  plot_obj_indie.initialise(chamber[i], input_file_name, saving_path);
  	  //plot_obj.initialise(chamber[i], input_file_name, saving_path);
    	std::cout<<" started progra:"<<std::endl;
    	//plot_obj_indie.plot_individual_channels("_pressure_2016");
    	//plot_obj_indie.plot_individual_channels("_pressure_2017");
    	//plot_obj_indie.plot_individual_channels("_pressure_2018");
    	//plot_obj_indie.plot_individual_channels("_instlumi_2016");
    	//plot_obj_indie.plot_individual_channels("_instlumi_2017");
    	//plot_obj_indie.plot_individual_channels("_instlumi_2018");
      //plot_obj_indie.plot_individual_channels("_integratelumi_initial");
    	//plot_obj_indie.plot_individual_channels("_integratelumi");

    	std::cout<<" done with individual channels "<<std::endl;
    	plot_obj.plot_goodchannels("_pressure_2016");
	    plot_obj.plot_goodchannels("_pressure_2017");
	    plot_obj.plot_goodchannels("_pressure_2018");
    	std::cout<<" started progra:"<<chamber[i]<<std::endl;

    	plot_obj.plot_goodchannels("_instlumi_2016");
    	plot_obj.plot_goodchannels("_instlumi_2017");
    	plot_obj.plot_goodchannels("_instlumi_2018");

////      std::cout<<" done with plotting channels"<<chamber[i]<<std::endl;
////    	//mean_values_pressure_pair_2016 = plot_obj.fit_goodchannels("_pressure_2016");
////      mean_values_pressure_pair_2017 = plot_obj.fit_goodchannels("_pressure_2017");
////      //mean_values_pressure_pair_2018 = plot_obj.fit_goodchannels("_pressure_2018");
////    	//mean_values_instlumi_pair_2016 = plot_obj.fit_goodchannels("_instlumi_2016");
////      mean_values_instlumi_pair_2017 = plot_obj.fit_goodchannels("_instlumi_2017");
////      //mean_values_instlumi_pair_2018 = plot_obj.fit_goodchannels("_instlumi_2018");
////
////      //mean_values_pressure_2016.push_back(mean_values_pressure_pair_2016);
////      mean_values_pressure_2017.push_back(mean_values_pressure_pair_2017);
////      //mean_values_pressure_2018.push_back(mean_values_pressure_pair_2018);
////      //mean_values_instlumi_2016.push_back(mean_values_instlumi_pair_2016);
////      mean_values_instlumi_2017.push_back(mean_values_instlumi_pair_2017);
////      //mean_values_instlumi_2018.push_back(mean_values_instlumi_pair_2018);
////
////      std::cout<<" started progra:"<<std::endl;
////  	  plot_obj.plot_goodchannels("_integratelumi_initial");
////   	  plot_obj.plot_goodchannels("_integratelumi");
//// }
////
//// // after obtaining the mean of slope for all the respective pressure or int lumi files 
//// // Lets make the plot for mean of slope vs chamber 
//// // First for pressure
////  std::cout<<" will plot mean values now on  canvas"<<std::endl;
////  for(int i=0;i<chamber.size(); i++){
////      std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2016[i].first<<" error "<<mean_values_pressure_2016[i].second<<std::endl;
////      //std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2017[i].first<<std::endl;
////  }
/////*
////  for(int i=0;i<chamber.size(); i++){
////      std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2017[i].first<<" error "<<mean_values_pressure_2017[i].second<<std::endl;
////      //std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2017[i].first<<std::endl;
////  }
////  for(int i=0;i<chamber.size(); i++){
////      std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2018[i].first<<" error "<<mean_values_pressure_2018[i].second<<std::endl;
////      //std::cout<<" mean :  chamber "<<chamber[i]<<" : "<<mean_values_pressure_2017[i].first<<std::endl;
////  }
////*/
//////	MeanPlot mean_plot_pressure_2016;
////	MeanPlot mean_plot_pressure_2017;
//////	MeanPlot mean_plot_pressure_2018;
//////	MeanPlot mean_plot_instlumi_2016;
////	MeanPlot mean_plot_instlumi_2017;
//////	MeanPlot mean_plot_instlumi_2018;
////
////  std::cout<<" initialising"<<std::endl;
//////	mean_plot_pressure_2016.initialise(mean_values_pressure_2016, chamber, "2016", "_pressure"); 
////	mean_plot_pressure_2017.initialise(mean_values_pressure_2017, chamber, "2017","_pressure"); 
//////	mean_plot_pressure_2018.initialise(mean_values_pressure_2018, chamber, "2018","_pressure"); 
//////  std::cout<<" mean values "<<std::endl;
//////	mean_plot_instlumi_2016.initialise(mean_values_instlumi_2016, chamber, "2016","_instlumi"); 
////	mean_plot_instlumi_2017.initialise(mean_values_instlumi_2017, chamber, "2017","_instlumi"); 
//////	mean_plot_instlumi_2018.initialise(mean_values_instlumi_2018, chamber, "2018","_instlumi"); 
////
//////	mean_plot_pressure_2016.plot_mean();
////  mean_plot_pressure_2017.plot_mean();
//////	mean_plot_pressure_2018.plot_mean();
//////	mean_plot_instlumi_2016.plot_mean();
////  mean_plot_instlumi_2017.plot_mean();
//////	mean_plot_instlumi_2018.plot_mean();
////
////  //We will reduce systematics by taking ratio with ME31HV2
////
///////  std::cout<<" here the normalisation "<<std::endl;
///////  std::pair<float, float> mean_values_integratelumi_value;
///////  std::vector<std::pair<float, float> > mean_values_integratelumi;
///////  TString saving_path = "/afs/cern.ch/work/n/nrawal/CSCAgeing_code_study/results_AnalysisCode/New_selections_2017/systematic_normalised_plots/";
///////  SystematicRemoval sys_obj;
///////  sys_obj.initialise(input_path, saving_path, "ME13HV3");
///////  for(int i=0; i<chamber.size(); i++){
///////    sys_obj.adding_plus_minus(chamber[i]);
///////    std::cout<<" added now going to ratio "<<std::endl;
///////    sys_obj.ratio(chamber[i]);
///////    std::cout<<" done ratio ones  values "<<std::endl;
///////    mean_values_integratelumi_value = sys_obj.normalised_fit(chamber[i]);
///////    std::cout<<" mean values "<<mean_values_integratelumi_value.first<<std::endl;
///////    mean_values_integratelumi.push_back(mean_values_integratelumi_value);
  }
///////	MeanPlot mean_plot_intlumi;
///////  std::cout<<" initiaiseing mean plot intlumi"<<std::endl;
///////	mean_plot_intlumi.initialise(mean_values_integratelumi, chamber, "run2", "_integratelumi"); 
///////	mean_plot_intlumi.plot_mean();
///////
////	std::cout<<" started program"<<std::endl;

return 0;
}
