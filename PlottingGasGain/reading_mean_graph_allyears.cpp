std::vector<TString> years = {"2016", "2017", "2018"};

class GasGain{
  public : 
    TString reference;
    TString normalising_term; 
    TString type_of_channels; 

    void initialise_variables(TString reference_name, TString normalising_name, TString type_name);
    void reading_gas_gain_cumulative(TString type_channels);
    void reading_pressure_instlumi_cumulative(TString , TString type_channels);
    void draw_mean_plot(std::vector<float> mean_values_vector, std::vector<float> mean_error_values_vector, TString var, std::vector<TString> chamber_name, TString year, TString type_channels);
    std::pair<float, float> normalised_fit(TGraphAsymmErrors* graph, TString chamber_name,TString,  TString year , TString type_channels) ;
    TGraphAsymmErrors* combining_different_years(std::vector<TGraphAsymmErrors*> yearly_normalised_graphs, TString, TString );
    std::pair<TGraphAsymmErrors*, TGraphAsymmErrors*>  normalising(TGraphAsymmErrors* graph_original, TString chamber_name, TString year,  TGraphAsymmErrors* graph_reference, TString, TString type_channels) ;
};

TGraphAsymmErrors* GasGain :: combining_different_years(std::vector<TGraphAsymmErrors*> yearly_normalised_graphs, TString chamber, TString type_channels){
        // Combine yearly normalized histograms
        int n_bins = yearly_normalised_graphs[0]->GetN();
        int total_n_bins = 0;
        for(int i=0;  i<yearly_normalised_graphs.size(); i++){
				  total_n_bins = total_n_bins +  yearly_normalised_graphs[i]->GetN() ;
        }
        std::cout<<" total bins "<<total_n_bins<<std::endl;
				std::vector<double> x_final(total_n_bins), y_final(total_n_bins), exl_final(total_n_bins), exh_final(total_n_bins), eyl_final(total_n_bins), eyh_final(total_n_bins);
        int j=0;
            for (auto& graph : yearly_normalised_graphs) {
                for(int i=0 ; i<graph->GetN(); i++){
                double x, y;
                 graph->GetPoint(i, x, y);
                 if (y > 0) {
										y_final[j] = y; 
										x_final[j] = x; 
										exl_final[j] = graph->GetErrorXlow(i); 
										exh_final[j] = graph->GetErrorXhigh(i); 
										eyl_final[j] =  graph->GetErrorYlow(i);
										eyh_final[j] = graph->GetErrorYhigh(i);
                 j=j+1;
                  }
            }
          }// endof filling graph
        TGraphAsymmErrors * final_normalised_histograms = new TGraphAsymmErrors(total_n_bins, x_final.data(), y_final.data(), exl_final.data(), exh_final.data(), eyl_final.data(), eyh_final.data());

            TCanvas *c1 = new TCanvas();
            c1->cd();
            final_normalised_histograms->Draw("AP");
            final_normalised_histograms->GetYaxis()->SetRangeUser(250,600);
            c1->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/run2/plots_all/"+chamber+"_individual_before_ratio_run2.pdf");
            
            
///
///				std::vector<double> x_final(total_n_bins), y_final(total_n_bins), exl_final(total_n_bins), exh_final(total_n_bins), eyl_final(total_n_bins), eyh_final(total_n_bins);
///        double first_bin_value_2016 = 1.0; // Default normalization
///        double first_bin_err_value_2016 = 0.0; // Default normalization
///        bool first_bin_found = false;
//        for (int i = 0; i < n_bins; i++) {
//            double x_ref, y_ref;
//            yearly_normalised_graphs[0]->GetPoint(i, x_ref, y_ref);
//            if (!first_bin_found && y_ref > 0) {
//               std::cout<<" chamber "<<chamber<<" normalising  value "<<first_bin_value_2016<<" error "<<first_bin_err_value_2016<<std::endl;
//                first_bin_value_2016 = y_ref; // Set first bin as normalization reference
//                first_bin_err_value_2016 = yearly_normalised_graphs[0]->GetErrorYhigh(i); // Set first bin as normalization reference
//                first_bin_found = true;
//            }
//          }
//           double total_error_ratio = 0;
//           int j =0;
//            for (auto& graph : yearly_normalised_graphs) {
//                double x, y;
//                for(int i=0 ; i<graph->GetN(); i++){
//                 graph->GetPoint(i, x, y);
//                 if (y > 0) {
//                    double normalized_value = y / first_bin_value_2016; // Normalize to first bin of 2016
//                    total_error_ratio = normalized_value * sqrt(pow(graph->GetErrorYhigh(i)/y,2) +  pow(first_bin_err_value_2016/first_bin_value_2016, 2));
//										y_final[j] = normalized_value; 
//										x_final[j] = x; 
//										exl_final[j] = graph->GetErrorXlow(i); 
//										exh_final[j] = graph->GetErrorXhigh(i); 
//										eyl_final[j] = total_error_ratio;
//										eyh_final[j] = total_error_ratio; 
//                  }
//                 j = j+1;
//               }
//            }
    std::cout << "Final normalised histograms saved!" << std::endl;
		return final_normalised_histograms;
}

std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> GasGain ::  normalising(TGraphAsymmErrors* graph_original, TString chamber_name, TString year,  TGraphAsymmErrors* graph_reference, TString reference_chamber, TString type_channels) {
    int n = graph_original->GetN();
    std::vector<double> x(n), y(n), exl(n), exh(n), eyl(n), eyh(n);

   TString type_string; 
   if(type_channels=="all_channels") type_string = "all channels";
   else if(type_channels=="plus_channels") type_string = "+ endcap";
   else if(type_channels=="minus_channels") type_string = "- endcap";

		// Draw both the original and reference on same canvas for a comparison before taking ratio
    graph_original->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    graph_original->SetTitle(chamber_name + " : "+year+" : Gas gain after pressure/instlumi corrections ("+year+") : " +type_string);
    graph_reference->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    graph_reference->SetTitle(reference_chamber+" : Gas gain after pressure/instlumi corrections ("+year+") : " +type_string);

    graph_original->GetYaxis()->SetRangeUser(280, 600);
    graph_reference->GetYaxis()->SetRangeUser(280, 600);

		TCanvas *c1 = new TCanvas();
		c1->cd();
		graph_original->Draw("AP");	
		graph_original->SetMarkerColor(kBlue);
		graph_original->SetMarkerStyle(20);
		graph_original->SetMarkerSize(0.5);
		graph_reference->SetMarkerColor(kRed);
		graph_reference->SetMarkerStyle(20);
		graph_reference->SetMarkerSize(0.5);

		graph_reference->Draw("P same");

    TLegend *leg = new TLegend(0.1, 0.8, 0.30, 0.9); 
    leg->AddEntry(graph_original, chamber_name, "lep");
    leg->AddEntry(graph_reference, reference_chamber, "lep");

    leg->Draw("same");
    c1->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+reference_chamber+"_before_ratio_"+year+".pdf");

  
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
    h_normalised->SetTitle(chamber_name + " :  Gas gain after pressure/instlumi corrections (normalised wrt "+reference_chamber+") ("+year+") : "+type_string);
		h_normalised->SetMarkerStyle(20);
		h_normalised->SetMarkerSize(0.5);

		// Draw the histogram after taking ratio, and save the normalised hist 
		TCanvas * c= new TCanvas();
		c->cd();
		h_normalised->Draw("AP");
    c->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+reference_chamber+"_after_ratio_"+year+".pdf");
	
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
    h_normalised_first->SetTitle(chamber_name + " : "+year+" : Gas gain after pressure/instlumi corrections (normalised wrt "+reference_chamber+") ("+year+") : "+type_string);
		h_normalised_first->SetMarkerStyle(20);
		h_normalised_first->SetMarkerSize(0.5);

  h_normalised_first->GetYaxis()->SetRangeUser(0.95,1.05);

  TCanvas* c2 = new TCanvas();
  c2->cd();
  h_normalised_first->Draw("AP");
  c2->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/"+year+"/plots_all/"+chamber_name+"_"+reference_chamber+"_after_ratio_normalised_"+year+".pdf");
  
  std::pair<TGraphAsymmErrors*,TGraphAsymmErrors* > hist_normalised;
  hist_normalised.first = h_normalised;
  hist_normalised.second = h_normalised_first;

  return hist_normalised;
}

std::pair<float, float> GasGain :: normalised_fit(TGraphAsymmErrors* graph, TString chamber_name, TString reference_chamber, TString year, TString type_channels) {
   
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

  c->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/"+year+"/fits/"+chamber_name+"_"+reference_chamber+"_after_ratio_normalised_fit_"+year+".pdf");
  return {slope, slope_error};
}


void GasGain :: draw_mean_plot(std::vector<float> mean_values_vector, std::vector<float> mean_error_values_vector, TString var, std::vector<TString> chamber_name, TString year, TString type_channels){
    
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
      graph_title = " Slope of gas gain dependence on "+var+" : "+year+" : "+type_channels+" channels";
    }
    else if(var=="pressure_second"){
      graph_title = " Slope of gas gain dependence on "+var+" after pressure corr : "+year+" : "+type_channels+" channels";
    }
    else if(var=="instlumi_second"){
      graph_title = " Slope of gas gain dependence on "+var+" after pressure/instlumi corr : "+year+" : "+type_channels+" channels";
    }
    else if(var=="intlumi"){
      graph_title = " Slope of gas gain dependence on integrated luminosity (after pressure/instlumi cor.) : "+year+" : "+type_channels+" channels";

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
    std::map <TString, int> marker_type = { {"ME11a", 20}, {"ME11b", 24}, 
      {"ME12HV1", 20}, {"ME12HV2", 20}, {"ME12HV3", 24},
      {"ME13HV1", 20}, {"ME13HV2", 20}, {"ME13HV3", 24},
      {"ME21HV1", 20}, {"ME21HV2", 20}, {"ME21HV3", 24},
      {"ME22HV1", 20},{"ME22HV2", 20},{"ME22HV3", 20},{"ME22HV4", 20},{"ME22HV5", 24},
      {"ME31HV1", 20}, {"ME31HV2", 20}, {"ME31HV3", 24},
      {"ME32HV1", 20},{"ME32HV2", 20},{"ME32HV3", 20},{"ME32HV4", 20},{"ME32HV5", 24},
      {"ME41HV1", 20}, {"ME41HV2", 20}, {"ME41HV3", 24},
      {"ME42HV1", 20},{"ME42HV2", 20},{"ME42HV3", 20},{"ME42HV4", 20},{"ME42HV5", 24}
    };

    //int marker_colour[32] = {8, 8, 2, 2, 2, 2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    //int marker_colour[32] = {8, 8, 2, 2, 2, 2,2,2,  4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    //int marker_colour[32] = {8, 2, 2, 2, 2, 2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2, 4,4,4, 2,2,2,2,2 };
    TLegend *legend_1 = new TLegend(0.7,0.7,0.9,0.9);
    for (Int_t j=0; j<chamber_name.size(); j++) {
    std::cout<<" declaring gr x arrayr "<<std::endl;
        TMarker *m = new TMarker(gr_xarray[j], gr_yarray[j], 20);
        m->SetMarkerColor(marker_colour[chamber_name[j]]);
        m->SetMarkerStyle(marker_type[chamber_name[j]]);
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

	 canv1->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/mean_slope_values_"+var+"_"+year+"_fit.pdf"); 
  }

void GasGain :: initialise_variables(TString reference_name, TString normalising_name, TString type_name){
//TString reference = "ME13HV3";
//TString normalising_term = "normalised_ME13HV3";

 reference = reference_name;
 normalising_term = normalising_name;
 type_of_channels = type_name;
//TString type_of_channels = "all";
//TString type_of_channels = "plus";
//TString type_of_channels = "minus";
//TString type_of_channels="lower_minus";
//TString type_of_channels="lower_plus";
//TString type_of_channels="upper_minus";
//TString type_of_channels="upper_plus";
//TString type_of_channels="odd_layers";
//TString type_of_channels="even_layers";
//TString type_of_channels="odd_chambers";
//TString type_of_channels="even_chambers";
  } //end of initialisation


void  reading_mean_graph_allyears(TString reference,TString normalising_term, TString type_of_channels ){

  GasGain *obj = new GasGain();
  obj->initialise_variables(reference, normalising_term, type_of_channels);
  obj->reading_gas_gain_cumulative(type_of_channels);
  delete obj;
  
//  GasGain *obj_pressure = new GasGain();
//  obj_pressure->initialise_variables(reference, normalising_term, type_of_channels);
//  obj_pressure->reading_pressure_instlumi_cumulative("pressure", type_of_channels);
//  delete obj_pressure;
//
//  GasGain *obj_pressure_second = new GasGain();
//  obj_pressure_second->initialise_variables(reference, normalising_term, type_of_channels);
//  obj_pressure_second->reading_pressure_instlumi_cumulative("pressure_second", type_of_channels);
//  delete obj_pressure_second;
//  GasGain *obj_instlumi = new GasGain();
//  obj_instlumi->initialise_variables(reference, normalising_term, type_of_channels);
//  obj_instlumi->reading_pressure_instlumi_cumulative("instlumi", type_of_channels);
//  delete obj_instlumi;
//
//  GasGain *obj_instlumi_second = new GasGain();
//  obj_instlumi_second->initialise_variables(reference, normalising_term, type_of_channels);
//  obj_instlumi_second->reading_pressure_instlumi_cumulative("instlumi_second", type_of_channels);
//  delete obj_instlumi_second;

}

void GasGain :: reading_pressure_instlumi_cumulative(TString thevar , TString type_channels){
    std::vector<float> mean_2016_value_list;
    std::vector<float> mean_2016_error_value_list;
    std::vector<float> mean_2017_value_list;
    std::vector<float> mean_2017_error_value_list;
    std::vector<float> mean_2018_value_list;
    std::vector<float> mean_2018_error_value_list;
  
    std::vector<TString> chambers = {"ME11a", "ME11b", 
    "ME12HV1" ,
    "ME12HV2" , "ME12HV3",
    "ME13HV1" ,"ME13HV2" , "ME13HV3",
    "ME21HV1", "ME21HV2", "ME21HV3", 
    "ME22HV1", 
    "ME22HV2" ,"ME22HV3","ME22HV4", 
    "ME22HV5",
    "ME31HV1", "ME31HV2" ,"ME31HV3",
    "ME32HV1", "ME32HV2" ,"ME32HV3", "ME32HV4", 
    "ME32HV5",
    "ME41HV1", "ME41HV2" ,"ME41HV3",
    "ME42HV1", "ME42HV2" ,"ME42HV3", "ME42HV4",
    "ME42HV5"
  };

    TString type_string; 
    if(type_channels=="all") type_string = "all channels";
    else if(type_channels=="plus") type_string = "+ endcap";
    else if(type_channels=="minus") type_string = "- endcap";
    else if(type_channels=="lower_minus") type_string = "lower side - endcap";
    else if(type_channels=="lower_plus") type_string = "lower side + endcap";
    else if(type_channels=="upper_plus") type_string = "upper side + endcap";
    else if(type_channels=="upper_minus") type_string = "upper side - endcap";
    else if(type_channels=="odd_chambers") type_string = "odd chambers";
    else if(type_channels=="even_chambers") type_string = "even chambers";
    else if(type_channels=="even_layers") type_string = "even layers";
    else if(type_channels=="odd_layers") type_string = "odd layers";
  
    // reading ME13HV3 histogram as reference histogram
    std::cout<<" will start looking into different chambers"<<std::endl;
    std::map<TString, TGraphAsymmErrors*> final_normalised_histograms;
  
   // iterating through all chambers 
   // Through all years  
    for (auto& chamber : chambers) {
      std::vector<TGraphAsymmErrors*> yearly_normalised_graphs;
      TGraphAsymmErrors* final_normalised_histogram;
      for (auto& year : years) {
          TDirectoryFile *dir;
          TGraphAsymmErrors* h_var;  

          TString input_file = year+"_old_full/dataset_output_"+chamber+"_"+year+".root";
    		  std::cout<<" input file "<<input_file<<std::endl;
				  TFile * file = TFile::Open(input_file, "READ");
    		
          dir = (TDirectoryFile*) file->Get(thevar);
     		  h_var = (TGraphAsymmErrors*) dir->Get("dataset_trimmed_"+chamber+"_allgoodchannelsvs_"+thevar+"_"+type_channels);

          // cosmetic
          h_var->GetYaxis()->SetTitle("Trimmed Mean Charge");
          if(thevar=="pressure") h_var->SetTitle("Gas gain dependence on pressure : "+chamber+"("+year+") : "+type_string);
          else if(thevar=="pressure_second") h_var->SetTitle("Gas gain dependence on pressure after pressure corr : "+chamber+"("+year+" : "+type_string);
          else if(thevar=="instlumi")  h_var->SetTitle("Gas gain dependence on instlumi : "+chamber+"("+year+") : "+type_string);
          else if(thevar=="instlumi_second")  h_var->SetTitle("Gas gain dependence on instlumi after pressure/instlumi corr : "+chamber+"("+year+" : "+type_string);
          TString saving_name;
          if(thevar=="pressure" || thevar=="pressure_second") 
           saving_name ="AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/pressure/"+year+"/plots_all/"+thevar+"_dependence_"+chamber+"_"+year+".pdf" ;
          else if(thevar=="instlumi" || thevar=="instlumi_second") 
           saving_name ="AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/instlumi/"+year+"/plots_all/"+thevar+"_dependence_"+chamber+"_"+year+".pdf" ;
          h_var->SetMarkerStyle(20);
		      h_var->SetMarkerSize(0.5);
          // end of cosmeetic
          TCanvas *c = new TCanvas();
          c->cd();
          h_var->Draw("AP");
          c->SaveAs(saving_name);

          TF1 *fa1 = (TF1*) h_var->GetFunction("expFit");
          double slope = fa1->GetParameter(1);
          double slope_error = fa1->GetParError(1);
		    	if(year=="2016") {
		       mean_2016_value_list.push_back(slope);
		       mean_2016_error_value_list.push_back(slope_error);
          } 
		    	else if(year=="2017"){
		       mean_2017_value_list.push_back(slope);
		       mean_2017_error_value_list.push_back(slope_error);
          }
		      else if(year=="2018"){
		       mean_2018_value_list.push_back(slope);
		       mean_2018_error_value_list.push_back(slope_error);
          }
       }
   } // end of all chambers
   if (std::find(years.begin(), years.end(), "2016") != years.end()) {
    this->draw_mean_plot(mean_2016_value_list, mean_2016_error_value_list, thevar, chambers, "2016", type_channels);
   }
   if (std::find(years.begin(), years.end(), "2017") != years.end()) {
    this->draw_mean_plot(mean_2017_value_list, mean_2017_error_value_list, thevar ,chambers, "2017", type_channels);
   }
   if (std::find(years.begin(), years.end(), "2018") != years.end()) {
    this->draw_mean_plot(mean_2018_value_list, mean_2018_error_value_list, thevar ,chambers, "2018", type_channels);
   }
}
void GasGain:: reading_gas_gain_cumulative(TString type_channels){
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
    "ME22HV2" ,"ME22HV3", "ME22HV4", 
    "ME22HV5",
    "ME31HV1", "ME31HV2" ,"ME31HV3",
    "ME32HV1", "ME32HV2" ,"ME32HV3", "ME32HV4",
    "ME32HV5",
    "ME41HV1", "ME41HV2" ,"ME41HV3",
    "ME42HV1", "ME42HV2" ,"ME42HV3", "ME42HV4",
    "ME42HV5"
  };


   std::map<TString, TString> ref_chambers;
  
  if(reference=="ME13HV3"){
   ref_chambers = {
     {"ME11a", "ME13HV3"}, {"ME11b", "ME13HV3"},  
     {"ME12HV1", "ME13HV3"}, {"ME12HV2", "ME13HV3"}, {"ME12HV3", "ME13HV3"}, 
     {"ME13HV1", "ME13HV3"}, {"ME13HV2", "ME13HV3"}, {"ME13HV3", "ME13HV3"},
     {"ME21HV1", "ME13HV3"}, {"ME21HV2", "ME13HV3"}, {"ME21HV3", "ME13HV3"},
     {"ME31HV1", "ME13HV3"}, {"ME31HV2", "ME13HV3"}, {"ME31HV3", "ME13HV3"},
     {"ME41HV1", "ME13HV3"}, {"ME41HV2", "ME13HV3"}, {"ME41HV3", "ME13HV3"},
     {"ME22HV1", "ME13HV3"}, {"ME22HV2", "ME13HV3"}, {"ME22HV3", "ME13HV3"},{"ME22HV4", "ME13HV3"},{"ME22HV5", "ME13HV3"},
     {"ME32HV1", "ME13HV3"}, {"ME32HV2", "ME13HV3"}, {"ME32HV3", "ME13HV3"},{"ME32HV4", "ME13HV3"},{"ME32HV5", "ME13HV3"},
     {"ME42HV1", "ME13HV3"}, {"ME42HV2", "ME13HV3"}, {"ME42HV3", "ME13HV3"},{"ME42HV4", "ME13HV3"},{"ME42HV5", "ME13HV3"}
   }; // end of map
  } // end of if
  else {
   ref_chambers = {
     {"ME11a", "ME11b"}, {"ME11b", "ME11b"},  
     {"ME12HV1", "ME12HV3"}, {"ME12HV2", "ME12HV3"}, {"ME12HV3", "ME12HV3"}, 
     {"ME13HV1", "ME13HV3"}, {"ME13HV2", "ME13HV3"}, {"ME13HV3", "ME13HV3"},
     {"ME21HV1", "ME21HV3"}, {"ME21HV2", "ME21HV3"}, {"ME21HV3", "ME21HV3"},
     {"ME31HV1", "ME31HV3"}, {"ME31HV2", "ME31HV3"}, {"ME31HV3", "ME31HV3"},
     {"ME41HV1", "ME41HV3"}, {"ME41HV2", "ME41HV3"}, {"ME41HV3", "ME41HV3"},
     {"ME22HV1", "ME22HV5"}, {"ME22HV2", "ME22HV5"}, {"ME22HV3", "ME22HV5"},{"ME22HV4", "ME22HV5"},{"ME22HV5", "ME22HV5"},
     {"ME32HV1", "ME32HV5"}, {"ME32HV2", "ME32HV5"}, {"ME32HV3", "ME32HV5"},{"ME32HV4", "ME32HV5"},{"ME32HV5", "ME32HV5"},
     {"ME42HV1", "ME42HV5"}, {"ME42HV2", "ME42HV5"}, {"ME42HV3", "ME42HV5"},{"ME42HV4", "ME42HV5"},{"ME42HV5", "ME42HV5"}
   }; // end of map

  } // end of else
  //end of deciding reference
 
    // reading ME13HV3 histogram as reference histogram
    std::cout<<" will start looking into different chambers"<<std::endl;


    std::map<TString, TGraphAsymmErrors*> final_normalised_histograms;

   for (auto& chamber : chambers) {
      std::vector<TGraphAsymmErrors*> yearly_normalised_graphs;
      std::vector<TGraphAsymmErrors*> yearly_normalised_graphs_orig;
      std::vector<TGraphAsymmErrors*> yearly_normalised_graphs_ref;
      TGraphAsymmErrors* final_combined_histogram_orig;
      TGraphAsymmErrors* final_combined_histogram_ref;
      TString reference_chamber = ref_chambers[chamber];
      for (auto& year : years) {

			    TString input_file_ref = year+"_old_full/dataset_output_"+reference_chamber+"_"+year+".root";
    			TFile * file_ref = TFile::Open(input_file_ref, "READ");
          TDirectoryFile *dir_intlumi_ref, *dir_intlumi_initial ; 
          TGraphAsymmErrors*  h_reference; 
          TGraphAsymmErrors* h_intlumi_initial;  
          
            dir_intlumi_ref = (TDirectoryFile*) file_ref->Get("intlumi_final");
	  		    h_reference = (TGraphAsymmErrors*) dir_intlumi_ref->Get("dataset_trimmed_"+reference_chamber+"_allgoodchannelsvs_intlumi_final_"+type_channels);

          TString input_file = year+"_old_full/dataset_output_"+chamber+"_"+year+".root";
    		  std::cout<<" input file "<<input_file<<std::endl;
				  TFile * file = TFile::Open(input_file, "READ");

          // reading intlumi initial and then plot it
          dir_intlumi_initial = (TDirectoryFile*) file->Get("intlumi_initial");
     		  h_intlumi_initial = (TGraphAsymmErrors*) dir_intlumi_initial->Get("dataset_trimmed_"+chamber+"_allgoodchannelsvs_intlumi_initial_"+type_channels);

          TCanvas *c1 = new TCanvas();
          c1->cd();
          h_intlumi_initial->GetXaxis()->SetTitle("Integrated luminosity (fb^{-1})");
          h_intlumi_initial->GetYaxis()->SetTitle("Trimmed mean Charge");
          h_intlumi_initial->SetTitle("Gas gain dependence on integrated luminosity before any correction");
          h_intlumi_initial->SetMarkerStyle(20);
          h_intlumi_initial->SetMarkerSize(0.5);

          h_intlumi_initial->Draw("AP");
          c1->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/"+year+"/plots_all/gas_gain_"+chamber+"_intlumi_initial.pdf");
          TDirectoryFile *dir_intlumi_final; 
          TGraphAsymmErrors* h_intlumi_final;
    	
          dir_intlumi_final = (TDirectoryFile*) file->Get("intlumi_final");
          h_intlumi_final = (TGraphAsymmErrors*) dir_intlumi_final->Get("dataset_trimmed_"+chamber+"_allgoodchannelsvs_intlumi_final_"+type_channels);
		  
          std::pair< TGraphAsymmErrors *  , TGraphAsymmErrors * > normalised_histograms; 
          // Normalising gas gain wrt intlumi of final chambers with reference chamber
          normalised_histograms = this->normalising(h_intlumi_final, chamber, year,  h_reference, reference_chamber, type_channels);  
           TGraphAsymmErrors * normalised_hist = normalised_histograms.first ;
           TGraphAsymmErrors * normalised_ratio_hist = normalised_histograms.second;
          yearly_normalised_graphs.push_back(normalised_hist);
          yearly_normalised_graphs_orig.push_back(h_intlumi_final);
          yearly_normalised_graphs_ref.push_back(h_reference);
          std::pair<float, float> mean_intlumi_pair = this->normalised_fit(normalised_ratio_hist, chamber, reference_chamber, year, type_channels);
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
 	 // final_combined_histogram = 	this->combining_different_years(yearly_normalised_graphs, chamber);
 	 // final_normalised_histogram = 	this->fitting_different_years(final_combined_histogram, chamber);
 	  final_combined_histogram_orig = 	this->combining_different_years(yearly_normalised_graphs_orig, chamber, type_channels);
 	  final_combined_histogram_ref = 	this->combining_different_years(yearly_normalised_graphs_ref, chamber, type_channels);
    std::pair<TGraphAsymmErrors *, TGraphAsymmErrors*> final_normalised_histogram;
    final_normalised_histogram = this->normalising(final_combined_histogram_orig,  chamber, "run2", final_combined_histogram_ref, reference_chamber, type_channels);  
    TGraphAsymmErrors * final_normalised_ratio_histogram = final_normalised_histogram.second;
    final_normalised_ratio_histogram->SetMarkerStyle(20);
    final_normalised_ratio_histogram->SetMarkerSize(0.5);
  
    TCanvas *c = new TCanvas();
    c->cd();
    final_normalised_ratio_histogram->Draw("AP");
    final_normalised_ratio_histogram->SetTitle(" Gas gain after pressure/instlumi corrections (normalised wrt "+reference_chamber+") : "+chamber);
    final_normalised_ratio_histogram->GetXaxis()->SetTitle("integrated luminosity (fb^{-1})");
    final_normalised_ratio_histogram->GetYaxis()->SetTitle("Gas gain");
    final_normalised_ratio_histogram->GetYaxis()->SetRangeUser(0.95,1.05);

    c->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/all/plots_all/final_normalised_"+chamber+"_"+reference_chamber+".pdf");
    double bin_low_edge, bin_up_edge;
    bool flag_first = false;
    int n = final_normalised_ratio_histogram->GetN();
    for (int i = 0; i < n; i++) {
         double x, y;
         final_normalised_ratio_histogram->GetPoint(i, x, y);
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
    final_normalised_ratio_histogram->Fit(fa1,"R");
    TCanvas *c1 = new TCanvas();
    c1->cd();
    gStyle->SetOptFit(1111);
    final_normalised_ratio_histogram->Draw("AP");
    c1->SaveAs("AllResults/results_gas_gain_"+type_channels+"_"+normalising_term+"/intlumi/all/fits/final_normalised_"+chamber+"_fit_"+reference_chamber+".pdf");

    mean_intlumi_run2_value_list.push_back(fa1->GetParameter(0));
    mean_intlumi_run2_error_value_list.push_back(fa1->GetParError(0));

  } // end of all chambers

   if (std::find(years.begin(), years.end(), "2016") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2016_value_list, mean_intlumi_2016_error_value_list, "intlumi", chambers, "2016", type_channels);
   }
   if (std::find(years.begin(), years.end(), "2017") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2017_value_list, mean_intlumi_2017_error_value_list, "intlumi", chambers, "2017", type_channels);
   }
   if (std::find(years.begin(), years.end(), "2018") != years.end()) {
    this->draw_mean_plot(mean_intlumi_2018_value_list, mean_intlumi_2018_error_value_list, "intlumi", chambers, "2018", type_channels);
   }
 
    this->draw_mean_plot(mean_intlumi_run2_value_list, mean_intlumi_run2_error_value_list, "intlumi", chambers, "run2", type_channels);
}
