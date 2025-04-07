std::vector<TString> years = {"2016", "2017", "2018"};
//TString ratio_type = "plus_minus_endcap"; 
//TString ratio_type = "odd_even_chambers"; 
//TString ratio_type = "odd_even_layers"; 
//TString ratio_type = "lower_plus_minus_endcap"; 
//TString ratio_type = "upper_plus_minus_endcap"; 

std::map<TString, std::pair<TString, TString> > map_ratio_type = { 
  {"plus_minus_endcap", {"minus", "plus"}} , 
  {"odd_even_chambers", {"odd_chambers", "even_chambers"}} , 
  {"odd_even_layers", {"odd_layers", "even_layers"}},  
  {"lower_plus_minus_endcap", {"lower_minus", "lower_plus"}},  
  {"upper_plus_minus_endcap", {"upper_minus", "upper_plus"}}  
  };

class GasGain{
  public : 
    TString ratio_type ; 
    TFile * output_file ;
   
    void initialise_variable(TString); 
    void reading_gas_gain_cumulative();
    TGraphAsymmErrors* combining_different_years(std::vector<TGraphAsymmErrors*> yearly_normalised_graphs, TString);
    void draw_mean_plot(std::vector<float> mean_values_vector, std::vector<float> mean_error_values_vector, TString var, std::vector<TString> chamber_name, TString year);
    std::pair<float, float> normalised_fit(TGraphAsymmErrors* graph, TString chamber_name,  TString year ) ;
    std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> normalising(TGraphAsymmErrors* graph_original, TString chamber_name, TString year,  TGraphAsymmErrors* graph_reference) ;
};
TGraphAsymmErrors* GasGain :: combining_different_years(std::vector<TGraphAsymmErrors*> yearly_normalised_graphs, TString chamber){
        // Combine yearly normalized histograms
        int n_bins = yearly_normalised_graphs[0]->GetN();
        int total_n_bins = 0;
        for(int i=0;  i<yearly_normalised_graphs.size(); i++){
				  total_n_bins = total_n_bins +  yearly_normalised_graphs[i]->GetN() ;
        }
				std::vector<double> x_final(total_n_bins), y_final(total_n_bins), exl_final(total_n_bins), exh_final(total_n_bins), eyl_final(total_n_bins), eyh_final(total_n_bins);
        double first_bin_value_2016 = 1.0; // Default normalization
        double first_bin_err_value_2016 = 0.0; // Default normalization
        bool first_bin_found = false;
        for (int i = 0; i < n_bins; i++) {
            double x_ref, y_ref;
            yearly_normalised_graphs[0]->GetPoint(i, x_ref, y_ref);
            if (!first_bin_found && y_ref > 0) {
               std::cout<<" chamber "<<chamber<<" normalising  value "<<first_bin_value_2016<<" error "<<first_bin_err_value_2016<<std::endl;
                first_bin_value_2016 = y_ref; // Set first bin as normalization reference
                first_bin_err_value_2016 = yearly_normalised_graphs[0]->GetErrorYhigh(i); // Set first bin as normalization reference
                first_bin_found = true;
            }
          }
           double total_error_ratio = 0;
           int j =0;
            for (auto& graph : yearly_normalised_graphs) {
                double x, y;
                for(int i=0 ; i<graph->GetN(); i++){
                 graph->GetPoint(i, x, y);
                 if (y > 0) {
                    double normalized_value = y / first_bin_value_2016; // Normalize to first bin of 2016
                    total_error_ratio = normalized_value * sqrt(pow(graph->GetErrorYhigh(i)/y,2) +  pow(first_bin_err_value_2016/first_bin_value_2016, 2));
										y_final[j] = normalized_value; 
										x_final[j] = x; 
										exl_final[j] = graph->GetErrorXlow(i); 
										exh_final[j] = graph->GetErrorXhigh(i); 
										eyl_final[j] = total_error_ratio;
										eyh_final[j] = total_error_ratio; 
                  }
                 j = j+1;
               }
            }
        TGraphAsymmErrors * final_normalised_histograms = new TGraphAsymmErrors(total_n_bins, x_final.data(), y_final.data(), exl_final.data(), exh_final.data(), eyl_final.data(), eyh_final.data());

    std::cout << "Final normalised histograms saved!" << std::endl;
		return final_normalised_histograms;
}

std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> GasGain ::  normalising(TGraphAsymmErrors* graph_original, TString chamber_name, TString year,  TGraphAsymmErrors* graph_reference) {

  TGraphAsymmErrors * graph_original_clone = (TGraphAsymmErrors*) graph_original->Clone();
    int n = graph_original->GetN();
    std::vector<double> x(n), y(n), exl(n), exh(n), eyl(n), eyh(n);
    
    TString string_first = map_ratio_type[ratio_type].first;
    TString string_second = map_ratio_type[ratio_type].second;

		// Draw both the original and reference on same canvas for a comparison before taking ratio
    graph_original->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    graph_original->SetTitle(chamber_name + " : Gas gain after pressure/instlumi corrections ("+year+")");
    graph_original_clone->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    graph_original_clone->SetTitle(chamber_name + " : Gas gain after pressure/instlumi corrections ("+year+")");

    graph_reference->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    graph_reference->SetTitle(chamber_name+" : Gas gain after pressure/instlumi corrections ("+year+")");


		TCanvas *c1 = new TCanvas();
		c1->cd();
		//graph_original->Draw("AP");	
		graph_original_clone->Draw("AP");	
		graph_original->SetMarkerColor(kBlue);
		graph_original->SetMarkerStyle(20);
		graph_original->SetMarkerSize(0.5);
		graph_original_clone->SetMarkerColor(kBlue);
		graph_original_clone->SetMarkerStyle(20);
		graph_original_clone->SetMarkerSize(0.5);
		graph_reference->SetMarkerColor(kRed);
		graph_reference->SetMarkerStyle(20);
		graph_reference->SetMarkerSize(0.5);
		graph_reference->Draw("P same");
      
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();
    
    // Loop over all points in the graph
    for (int i = 0; i < graph_original_clone->GetN(); ++i) {
        double x, y;
        graph_original_clone->GetPoint(i, x, y);
        if (y > 0) { // Ignore points with zero or negative Y
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
    }
    
    if (minY == std::numeric_limits<double>::max()) {
        minY = 1e-6; // Default small value if all points are zero or skipped
    }
    
    // Optionally, add a margin to Y-axis range
    double marginFactor = 0.02; // 10% margin
    double rangeMin = minY - marginFactor * fabs(minY);
    double rangeMax = maxY + marginFactor * fabs(maxY);
 
    double minY1 = std::numeric_limits<double>::max();
    double maxY1 = -std::numeric_limits<double>::max();
    
    // Loop over all points in the graph
    for (int i = 0; i < graph_reference->GetN(); ++i) {
        double x, y;
        graph_reference->GetPoint(i, x, y);
        if (y > 0) { // Ignore points with zero or negative Y
            if (y < minY1) minY1 = y;
            if (y > maxY1) maxY1 = y;
        }
    }
    
    if (minY1 == std::numeric_limits<double>::max()) {
        minY1 = 1e-6; // Default small value if all points are zero or skipped
    }
    
    // Optionally, add a margin to Y-axis range
    double marginFactor1 = 0.05; // 10% margin
    double rangeMin1 = minY1 - marginFactor * fabs(minY1);
    double rangeMax1 = maxY1 + marginFactor * fabs(maxY1);
   
    double rangeMinall , rangeMaxall;
    if(rangeMin <rangeMin1) rangeMinall  = rangeMin;
    else rangeMinall = rangeMin1; 
    if(rangeMax >rangeMax1) rangeMaxall =  rangeMax;
    else rangeMaxall = rangeMax1; 

    // Set the Y-axis range
    //graph_original_clone->GetYaxis()->SetRangeUser(rangeMinall, rangeMaxall);
    graph_original_clone->GetYaxis()->SetRangeUser(300, 480);
    TLegend *leg = new TLegend(0.1, 0.8, 0.30, 0.9); 
    leg->AddEntry(graph_original_clone, chamber_name+" "+string_first, "lep");
    leg->AddEntry(graph_reference, chamber_name + " "+string_second, "lep");

    leg->Draw("same");
    c1->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+ratio_type+"_before_ratio_"+year+".pdf");

  
    for (int i = 0; i < n; i++) {
        double xo, yo, xr, yr;
        graph_original->GetPoint(i, xo, yo);
        graph_reference->GetPoint(i, xr, yr);
        
        //std::cout<<" error in y "<<graph_original->GetErrorYhigh(i)<<" value "<<yr<<" x "<<xr<<std::endl;
        if (yr != 0) {
            y[i] = yo / yr;
            eyh[i] = y[i] * sqrt(pow(graph_original->GetErrorYhigh(i) / yo, 2) + pow(graph_reference->GetErrorYhigh(i) / yr, 2));
            eyl[i] = y[i] * sqrt(pow(graph_original->GetErrorYlow(i) / yo, 2) + pow(graph_reference->GetErrorYlow(i) / yr, 2));
        } else {
            y[i] = 0;
            eyh[i] = 0;
            eyl[i] = 0;
        }
        x[i] = xo;
        exl[i] = graph_original->GetErrorXlow(i);
        exh[i] = graph_original->GetErrorXhigh(i);
        //std::cout<<" error in y "<<eyh[i]<<std::endl;
    }
 
    TGraphAsymmErrors* h_normalised = new TGraphAsymmErrors(n, x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
    h_normalised->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    h_normalised->SetTitle(chamber_name + " "+string_first+" endcap :  Gas gain after pressure/instlumi corrections (normalised wrt "+chamber_name+" "+string_second+" endcap) ("+year+")");
		h_normalised->SetMarkerStyle(20);
		h_normalised->SetMarkerSize(0.5);

		// Draw the histogram after taking ratio, and save the normalised hist 
		TCanvas * c= new TCanvas();
		c->cd();
		h_normalised->Draw("AP");
    c->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+ratio_type+"_after_ratio_"+year+".pdf");
	
  	// Normalised the norm histogram with respect to the first bin 
    bool flag_first_bin = false;
    double value_first_bin = 0, value_first_bin_err = 0;
    
    int n1 = h_normalised->GetN();
    std::vector<double> x1(n1), y1(n1), exl1(n1), exh1(n1), eyl1(n1), eyh1(n1);
    
    for (int i = 0; i < n1; i++) {
        double xi, yi;
        h_normalised->GetPoint(i, xi, yi);
        
        if (yi == 0) {
            y1[i] = 0;
            eyl1[i] = 0;
            eyh1[i] = 0;
        }
        if (yi != 0 && !flag_first_bin) {
            flag_first_bin = true;
            value_first_bin = yi;
            value_first_bin_err = (h_normalised->GetErrorYlow(i) + h_normalised->GetErrorYhigh(i)) / 2;
        }
        if (yi != 0 && flag_first_bin) {
            double value_err = (h_normalised->GetErrorYlow(i) + h_normalised->GetErrorYhigh(i)) / 2;
            
            double ratio = yi / value_first_bin;
            double ratio_err = ratio * sqrt(pow(value_first_bin_err / value_first_bin, 2) + pow(value_err / yi, 2));
            
            y1[i] = ratio;
            eyl1[i] = ratio_err;
            eyh1[i] = ratio_err;
        }
        
        x1[i] = xi;
        exl1[i] = h_normalised->GetErrorXlow(i);
        exh1[i] = h_normalised->GetErrorXhigh(i);
    }
    
  TGraphAsymmErrors* h_normalised_first = new TGraphAsymmErrors(n1, x1.data(), y1.data(), exl1.data(), exh1.data(), eyl1.data(), eyh1.data());
    h_normalised_first->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    h_normalised_first->SetTitle(chamber_name + " "+string_first+" endcap : Gas gain after pressure/instlumi corrections (normalised wrt "+string_second+" endcap) ("+year+")");
		h_normalised_first->SetMarkerStyle(20);
		h_normalised_first->SetMarkerSize(0.5);

  h_normalised_first->GetYaxis()->SetRangeUser(0.95,1.05);

  TCanvas* c2 = new TCanvas();
  c2->cd();
  h_normalised_first->Draw("AP");
  c2->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+ratio_type+"_after_ratio_normalised_"+year+".pdf");
    
  std::pair<TGraphAsymmErrors*,TGraphAsymmErrors* > hist_normalised;
  hist_normalised.first = h_normalised;
  hist_normalised.second = h_normalised_first;

  return hist_normalised;
  
}

std::pair<float, float> GasGain :: normalised_fit(TGraphAsymmErrors* graph, TString chamber_name,  TString year) {
   
    std::cout<<" number of points to fit "<<graph->GetN()<<std::endl; 
    double bin_low_edge, bin_up_edge;
    bool flag_first = false;
    int n = graph->GetN();
    
    for (int i = 0; i < n; i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (y != 0 && !flag_first) {
            bin_low_edge = x;
            flag_first = true;
        }
        if (y != 0) bin_up_edge = x;
    }
    
    if (bin_low_edge >= bin_up_edge) {
        std::cerr << "Error: Invalid fit range. No non-zero entries found.\n";
        return {0, 0};
    }

    std::cout<<" bin low "<<bin_low_edge<<" bin up "<<bin_up_edge<<std::endl; 
    TF1* fa1 = new TF1("fa1", "[0]*x+[1]", bin_low_edge, bin_up_edge);
    fa1->SetParameters(0.00001, 1);
    graph->Fit(fa1, "R");
    
    float slope = fa1->GetParameter(0);
    float slope_error = fa1->GetParError(0);
 
  TCanvas *c = new TCanvas();
  c->cd();
  gStyle->SetOptFit(1111);
  graph->Draw("AP");
  fa1->Draw("same");
  c->Update();

  c->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/"+year+"/fits/"+chamber_name+"_"+ratio_type+"_after_ratio_normalised_fit_"+year+".pdf");
  return {slope, slope_error};
}


void GasGain :: draw_mean_plot(std::vector<float> mean_values_vector, std::vector<float> mean_error_values_vector, TString var, std::vector<TString> chamber_name, TString year){
    
    std::cout<<" started in plot mean plot "<<std::endl;
    std::vector<float> value ={1,2,3,4,5,6,7,8,9,10,11,12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

    int size_chamber = chamber_name.size();
    value.resize(size_chamber);
    //std::vector<float> value = {1,2,3, 4,5, 6, 7, 8};
    TGraphErrors *graph_mean_values = new TGraphErrors(size_chamber, value.data(), mean_values_vector.data(), 0, mean_error_values_vector.data());
    std::cout<<" after declaring graph "<<std::endl;
    graph_mean_values->GetXaxis()->SetTickLength(1);
    TString graph_title;
    if(var=="pressure" || var=="instlumi"){
      graph_title = " Slope of gas gain dependence on "+var+" : "+year;
    }
    else if(var=="pressure_second"){
      graph_title = " Slope of gas gain dependence on "+var+" after pressure corr : "+year;
    }
    else if(var=="instlumi_second"){
      graph_title = " Slope of gas gain dependence on "+var+" after pressure/instlumi corr : "+year;
    }
    else if(var=="intlumi"){
      graph_title = " Slope of gas gain dependence on integrated luminosity (after pressure/instlumi cor.) : "+year;

    }

    if(var=="pressure"){
      graph_mean_values->GetYaxis()->SetRangeUser(-0.010,0.002); 
    }
    else if(var=="pressure_second") graph_mean_values->GetYaxis()->SetRangeUser(-0.0005,0.0005);
    else if(var=="instlumi") graph_mean_values->GetYaxis()->SetRangeUser(-0.00001,0.00001);
    else if(var=="instlumi_second") graph_mean_values->GetYaxis()->SetRangeUser(-0.00001,0.00001);
    else if(var=="intlumi") graph_mean_values->GetYaxis()->SetRangeUser(-0.0003,0.0003);
    std::cout<<" before declaring axes "<<std::endl;
    TAxis *axis = graph_mean_values->GetXaxis();
    axis->Draw();
   
    for(int i=0; i<chamber_name.size(); i++){
       graph_mean_values->GetXaxis()->SetBinLabel(graph_mean_values->GetXaxis()->FindBin(i + 1.), chamber_name[i]); // Find out     which bin on the x-axis the point corresponds to and set the bin label
    }
    graph_mean_values->GetXaxis()->SetTitleOffset(0.1); 
    TCanvas *canv1 = new TCanvas();
    canv1->cd();
    canv1->SetLeftMargin(0.12);
    canv1->SetGrid();
    gPad->SetGrid();
    graph_mean_values->Draw("AP");
    graph_mean_values->SetTitle(graph_title);
    Double_t *gr_xarray = graph_mean_values->GetX();
    Double_t *gr_yarray = graph_mean_values->GetY();

    std::cout<<" before starting marker "<<std::endl;
    std::map <TString, int> marker_colour = { {"ME11a", 8}, {"ME11b", 8}, 
      {"ME12HV1", 2}, {"ME12HV2", 2}, {"ME12HV3", 2},
      {"ME13HV1", 2}, {"ME13HV2", 2}, {"ME13HV3", 2},
      {"ME21HV1", 4}, {"ME21HV2", 4}, {"ME21HV3", 4},
      {"ME22HV1", 2},{"ME22HV2", 2},{"ME22HV3", 2},{"ME22HV4", 2},{"ME22HV5", 2},
      {"ME31HV1", 4}, {"ME31HV2", 4}, {"ME31HV3", 4},
      {"ME32HV1", 2},{"ME32HV2", 2},{"ME32HV3", 2},{"ME32HV4", 2},{"ME32HV5", 2},
      {"ME41HV1", 4}, {"ME41HV2", 4}, {"ME41HV3", 4},
      {"ME42HV1", 2},{"ME42HV2", 2},{"ME42HV3", 2},{"ME42HV4", 2},{"ME42HV5", 2}
    };
    //int marker_colour[32] = {8, 8, 2, 2, 2, 2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    //int marker_colour[32] = {8, 8, 2, 2, 2, 2,2,2,  4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    //int marker_colour[32] = {8, 2, 2, 2, 2, 2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    TLegend *legend_1 = new TLegend(0.7,0.7,0.9,0.9);
    for (Int_t j=0; j<chamber_name.size(); j++) {
    std::cout<<" declaring gr x arrayr "<<std::endl;
        TMarker *m = new TMarker(gr_xarray[j], gr_yarray[j], 20);
        m->SetMarkerColor(marker_colour[chamber_name[j]]);
        m->Draw();

      std::cout<<"after declaring gr x arrayr "<<std::endl;
			if(j==0) {
				legend_1->AddEntry(m," ME11a, ME11b (10^{0})","p");
			}
			if(j==3) {
			legend_1->AddEntry(m," Outer Chambers (10^{0})","p"); }
			if(j==8)
			legend_1->AddEntry(m," Inner Chambers (20^{0})","p");
   }
   std::cout<<" before drawing mulitgraph "<<std::endl;
	 legend_1->Draw("SAME");

   TLatex* cmslabel_1;
   TLatex* text1,*text2;
   cmslabel_1 = new TLatex(0.18,0.82, "CMS #bf{#it{Preliminary}}");
   cmslabel_1->SetNDC(kTRUE);
   cmslabel_1->SetTextSize(0.06);
   cmslabel_1->SetTextFont(42);
   cmslabel_1->Draw("same");

	 canv1->SaveAs("results_gas_gain_"+ratio_type+"/mean_slope_values_"+var+"_"+year+"_fit.pdf"); 
  }
 void GasGain:: initialise_variable(TString ratio_type_string){

   ratio_type = ratio_type_string;
 }

void  reading_mean_graph_allyears_ratio_self(TString ratio_type){

  GasGain *obj = new GasGain();
  obj->initialise_variable(ratio_type);
  obj->reading_gas_gain_cumulative();
  delete obj;
  
}

void GasGain:: reading_gas_gain_cumulative(){
    std::vector<float> mean_intlumi_run2_value_list;
    std::vector<float> mean_intlumi_run2_error_value_list;

    std::vector<float> mean_intlumi_2016_value_list;
    std::vector<float> mean_intlumi_2016_error_value_list;
    std::vector<float> mean_intlumi_2017_value_list;
    std::vector<float> mean_intlumi_2017_error_value_list;
    std::vector<float> mean_intlumi_2018_value_list;
    std::vector<float> mean_intlumi_2018_error_value_list;
  
    std::vector<TString> chambers = {"ME11a", "ME11b", 
    "ME12HV1" ,
    "ME12HV2" , "ME12HV3",
    "ME13HV1" ,"ME13HV2" , "ME13HV3",
    "ME21HV1", "ME21HV2", "ME21HV3", 
    "ME22HV1", 
    "ME22HV2" ,"ME22HV3","ME22HV4", 
    "ME22HV5",
    "ME31HV1", "ME31HV2" ,"ME31HV3",
    "ME32HV1", "ME32HV2" ,"ME32HV3","ME32HV4",
    "ME32HV5",
    "ME41HV1", "ME41HV2" ,"ME41HV3",
    "ME42HV1", "ME42HV2" ,"ME42HV3","ME42HV4",
    "ME42HV5",
  };

 //end of deciding reference
 
    // reading ME13HV3 histogram as reference histogram
    std::cout<<" will start looking into different chambers"<<std::endl;

    output_file = new TFile("output_file_"+ratio_type+".root","RECREATE");
    std::map<TString, TGraphAsymmErrors*> final_normalised_histograms;

    TString string_first = map_ratio_type[ratio_type].first;
    TString string_second = map_ratio_type[ratio_type].second;
    std::cout<<" string first "<<string_first<<" string second "<<string_second<<std::endl;

   for (auto& chamber : chambers) {
      std::vector<TGraphAsymmErrors*> yearly_normalised_graphs;
      TGraphAsymmErrors* final_normalised_histogram;
      for (auto& year : years) {

          TDirectoryFile  *dir_intlumi ; 
          TGraphAsymmErrors*  h_reference, *h_intlumi_final; 
          
          TString input_file = year+"_full/dataset_output_"+chamber+"_"+year+".root";
    		  std::cout<<" input file "<<input_file<<std::endl;
				  TFile * file = TFile::Open(input_file, "READ");
            dir_intlumi = (TDirectoryFile*) file->Get("intlumi_final");
            h_reference = (TGraphAsymmErrors*) dir_intlumi->Get("dataset_trimmed_"+chamber+"_allgoodchannelsvs_intlumi_final_"+string_first);
            h_intlumi_final = (TGraphAsymmErrors*) dir_intlumi->Get("dataset_trimmed_"+chamber+"_allgoodchannelsvs_intlumi_final_"+string_second);
		   
           std::pair< TGraphAsymmErrors *  , TGraphAsymmErrors * > normalised_histograms; 
          // Normalising gas gain wrt intlumi of final chambers with reference chamber
           normalised_histograms = this->normalising(h_intlumi_final, chamber, year,  h_reference);  
           TGraphAsymmErrors * normalised_hist = normalised_histograms.first ;
           TGraphAsymmErrors * normalised_ratio_hist = normalised_histograms.second;

           output_file->cd();
           normalised_hist->Write();
          yearly_normalised_graphs.push_back(normalised_hist);
          std::pair<float, float> mean_intlumi_pair = this->normalised_fit(normalised_ratio_hist, chamber, year);
		    	 if(year=="2016") {
		       mean_intlumi_2016_value_list.push_back(mean_intlumi_pair.first);
		       mean_intlumi_2016_error_value_list.push_back(mean_intlumi_pair.second);
          } 
		    	else if(year=="2017"){
		       mean_intlumi_2017_value_list.push_back(mean_intlumi_pair.first);
		       mean_intlumi_2017_error_value_list.push_back(mean_intlumi_pair.second);
          }
		      else if(year=="2018"){
		       mean_intlumi_2018_value_list.push_back(mean_intlumi_pair.first);
		       mean_intlumi_2018_error_value_list.push_back(mean_intlumi_pair.second);
          }
       }
 	  final_normalised_histogram = 	this->combining_different_years(yearly_normalised_graphs, chamber);

    final_normalised_histogram->SetMarkerStyle(20);
    final_normalised_histogram->SetMarkerSize(0.5);
 
//    output_file->cd();
    final_normalised_histogram->Write(); 

    TCanvas *c = new TCanvas();
    c->cd();
    final_normalised_histogram->Draw("AP");
    final_normalised_histogram->SetTitle("Gas gain after pressure/instlumi corrections : "+chamber + " "+string_first +" (normalised wrt "+string_second + ")");
    final_normalised_histogram->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    final_normalised_histogram->GetYaxis()->SetTitle("Gas gain");
    final_normalised_histogram->GetYaxis()->SetRangeUser(0.95,1.05);

    c->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/all/plots_all/final_normalised_"+chamber+"_"+ratio_type+".pdf");
    double bin_low_edge, bin_up_edge;
    bool flag_first = false;
    int n = final_normalised_histogram->GetN();
    for (int i = 0; i < n; i++) {
         double x, y;
         final_normalised_histogram->GetPoint(i, x, y);
         if (y != 0 && !flag_first) {
             bin_low_edge = x;
             flag_first = true;
         }
         if (y != 0) bin_up_edge = x;
     }
 
     if (bin_low_edge >= bin_up_edge) {
         std::cerr << "Error: Invalid fit range. No non-zero entries found.\n";
     }
    TF1 *fa1 = new TF1("fa1","[0] *x+[1]", bin_low_edge, bin_up_edge);
    final_normalised_histogram->Fit(fa1,"R");
    TCanvas *c1 = new TCanvas();
    c1->cd();
    gStyle->SetOptFit(1111);
    final_normalised_histogram->Draw("AP");
    c1->SaveAs("results_gas_gain_"+ratio_type+"/intlumi/all/fits/final_normalised_"+chamber+"_"+ratio_type+"_fit.pdf");

    mean_intlumi_run2_value_list.push_back(fa1->GetParameter(0));
    mean_intlumi_run2_error_value_list.push_back(fa1->GetParError(0));

  } // end of all chambers

   if (std::find(years.begin(), years.end(), "2016") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2016_value_list, mean_intlumi_2016_error_value_list, "intlumi", chambers, "2016");
   }
   if (std::find(years.begin(), years.end(), "2017") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2017_value_list, mean_intlumi_2017_error_value_list, "intlumi", chambers, "2017");
   }
   if (std::find(years.begin(), years.end(), "2018") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2018_value_list, mean_intlumi_2018_error_value_list, "intlumi", chambers, "2018");
   }
     this->draw_mean_plot(mean_intlumi_run2_value_list, mean_intlumi_run2_error_value_list, "intlumi", chambers, "run2");

     output_file->Close();
}
