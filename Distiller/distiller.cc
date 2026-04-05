/**
 * ==============================================================================
 * muEDM Analysis Distiller
 * ==============================================================================
 *
 * This module is dedicated to the high-level performance analysis of the muEDM 
 * tracking system. It processes the paired TTrees ("sim" and "rec") generated 
 * by the CHeTFitter using ROOT's RDataFrame.
 * 
 * Features:
 * - Efficiency & Resolution (Pulls)
 * - ROC Curves for Pattern Recognition
 * - Multiplicity & RAW Data Estimators
 * - Systematic Checks
 *
 * ==============================================================================
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <map>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TChain.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGraph.h>

// CHeT Includes
#include <CHeT/CHeTGlobalSettings.hh>

std::vector<double> GetTrueParams(double mc_x0, double mc_y0, double mc_z0, double mc_ux, double mc_uy, double mc_uz) {
    double trueLoc_x0 = mc_x0, trueLoc_y0 = mc_y0, trueLoc_z0 = mc_z0;
    CHeT::Config::ApplyInverseTransformation(trueLoc_x0, trueLoc_y0, trueLoc_z0);

    double trueLoc_ux = mc_ux, trueLoc_uy = mc_uy, trueLoc_uz = mc_uz;
    CHeT::Config::ApplyInverseRotation(trueLoc_ux, trueLoc_uy, trueLoc_uz);

    double t_to_y0 = -trueLoc_y0 / trueLoc_uy;
    double true_x0 = trueLoc_x0 + trueLoc_ux * t_to_y0;
    double true_z0 = trueLoc_z0 + trueLoc_uz * t_to_y0;
    double true_sx = trueLoc_ux / trueLoc_uy;
    double true_sz = trueLoc_uz / trueLoc_uy;

    return {true_x0, true_z0, true_sx, true_sz};
}

void DoEfficiencyAndResolution(ROOT::RDataFrame& df) {
    // 1. Definition of variables
    auto df_def = df
        // Momentum
        .Define("trk_p", "sqrt(trk_ux*trk_ux + trk_uy*trk_uy + trk_uz*trk_uz)")
        // Theta EDM formula: pz=cos(theta), pr=sin(theta)*(x*cos(phi)+y*sin(phi))/r
        .Define("trk_r", "sqrt(trk_x0*trk_x0 + trk_y0*trk_y0)")
        .Define("trk_phi", "atan2(trk_uy, trk_ux)")  // Note: coordinate meaning may vary depending on how mom is saved
        // Difference
        .Define("true_pars", GetTrueParams, {"trk_x0", "trk_y0", "trk_z0", "trk_ux", "trk_uy", "trk_uz"})
        .Define("diff_x", "rec.rec_x0 - true_pars[0]")
        .Define("diff_z", "rec.rec_z0 - true_pars[1]")
        .Define("diff_sx", "rec.rec_sx - true_pars[2]")
        .Define("diff_sz", "rec.rec_sz - true_pars[3]");

    auto df_acc = df_def.Filter("rec.rec_acceptance == true", "Events in acceptance");
    auto df_conv = df_acc.Filter("rec.rec_converged == true", "Converged fit events");

    // 2. Booking Histograms
    auto h_diff_x = df_conv.Histo1D({"h_diff_x", "Pull X; #Delta X [mm]; Events", 100, -5, 5}, "diff_x");
    auto h_diff_z = df_conv.Histo1D({"h_diff_z", "Pull Z; #Delta Z [mm]; Events", 100, -5, 5}, "diff_z");
    auto h_diff_sx = df_conv.Histo1D({"h_diff_sx", "Pull sX; #Delta sX; Events", 100, -0.1, 0.1}, "diff_sx");
    auto h_diff_sz = df_conv.Histo1D({"h_diff_sz", "Pull sZ; #Delta sZ; Events", 100, -0.1, 0.1}, "diff_sz");

    // We can also book TEfficiencies directly with DataFrame!
    auto fill_eff = [](TEfficiency& eff, double p, bool conv) { eff.Fill(conv, p); };
    TEfficiency eff("eff_p", "Efficiency vs Momentum; Momentum [MeV/c]; Efficiency", 20, 0, 70);
    
    // Instead of directly using Book which can be complex for TEfficiency without a custom helper,
    // we use a TH1D approach for passed/total or we manually loop/reduce.
    // For now let's just make two TH1Ds and divide them to get an efficiency plot using TEfficiency.
    
    auto h_passed = df_conv.Histo1D({"h_passed", "Passed; Momentum [MeV/c]; Events", 20, 0, 70}, "trk_p");
    auto h_total = df_acc.Histo1D({"h_total", "Total; Momentum [MeV/c]; Events", 20, 0, 70}, "trk_p");

    // 3. Draw Results
    TCanvas c("c_resolutions", "Resolutions", 800, 600);
    c.Divide(2, 2);
    c.cd(1); h_diff_x->Draw();
    c.cd(2); h_diff_z->Draw();
    c.cd(3); h_diff_sx->Draw();
    c.cd(4); h_diff_sz->Draw();
    c.SaveAs("Resolutions.pdf");

    TCanvas c2("c_eff", "Efficiency", 800, 600);
    if(TEfficiency::CheckConsistency(*h_passed, *h_total)) {
        TEfficiency* pEff = new TEfficiency(*h_passed, *h_total);
        pEff->SetTitle("Efficiency vs Momentum; Momentum [MeV/c]; Efficiency");
        pEff->Draw("AP");
    }
    c2.SaveAs("Efficiency.pdf");

    std::cout << ">>> Plotting finished. Created Resolutions.pdf and Efficiency.pdf" << std::endl;
}


void PlotHoughROC(const std::vector<std::string>& file_list, const std::vector<double>& thresholds) {
    if (file_list.size() != thresholds.size()) {
        std::cerr << "Mismatch between files and thresholds count!" << std::endl;
        return;
    }

    std::vector<double> eff_avg;
    std::vector<double> pur_avg;

    for (size_t i = 0; i < file_list.size(); ++i) {
        TChain chain_sim("sim");
        TChain chain_rec("rec");
        chain_sim.Add(file_list[i].c_str());
        chain_rec.Add(file_list[i].c_str());
        chain_sim.AddFriend(&chain_rec, "rec");

        ROOT::RDataFrame df(chain_sim);
        
        // Custom function to evaluate single event efficiency and purity
        // We compare the true hits (from sim) with the reconstructed hough hits (from rec)
        // Since we are working with vectors of hits, let's calculate intersection sizes
        
        auto df_eval = df.Define("tp", [](const std::vector<int>& sim_hits, const std::vector<int>& rec_hough) {
            int tp = 0;
            for (int h : rec_hough) {
                if (std::find(sim_hits.begin(), sim_hits.end(), h) != sim_hits.end()) tp++;
            }
            return tp;
        }, {"hit_ids", "rec.rec_hough2d_idx"})
        .Define("true_total", [](const std::vector<int>& sim_hits) { return (int)sim_hits.size(); }, {"hit_ids"})
        .Define("reco_total", [](const std::vector<int>& rec_hough) { return (int)rec_hough.size(); }, {"rec.rec_hough2d_idx"})
        .Define("eff", "true_total > 0 ? (double)tp / true_total : 0.0")
        .Define("pur", "reco_total > 0 ? (double)tp / reco_total : 0.0");

        double mean_eff = df_eval.Mean("eff").GetValue();
        double mean_pur = df_eval.Mean("pur").GetValue();

        eff_avg.push_back(mean_eff);
        pur_avg.push_back(mean_pur);
        
        std::cout << "Threshold: " << thresholds[i] << " -> Eff: " << mean_eff << ", Pur: " << mean_pur << std::endl;
    }

    TCanvas c("c_roc", "ROC Curve", 800, 600);
    TGraph* gROC = new TGraph(eff_avg.size(), pur_avg.data(), eff_avg.data());
    gROC->SetTitle("Hough Transform ROC; Purity; Efficiency");
    gROC->SetMarkerStyle(20);
    gROC->SetLineColor(kBlue);
    gROC->SetMarkerColor(kBlue);
    gROC->Draw("APL");

    // Add labels
    for (size_t i = 0; i < thresholds.size(); ++i) {
        auto *t = new TLatex(pur_avg[i], eff_avg[i], Form("%.2f", thresholds[i]));
        t->SetTextSize(0.03);
        t->Draw();
    }

    c.SaveAs("Hough_ROC.pdf");
    std::cout << ">>> Saved Hough_ROC.pdf" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> [analysis_type]" << std::endl;
        std::cerr << "       analysis_type: perf (default), raw, mult, est, align" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    std::string analysis_type = (argc > 2) ? argv[2] : "perf";

    ROOT::EnableImplicitMT();

    if (analysis_type == "perf") {
        TChain chain_sim("sim");
        TChain chain_rec("rec");

        chain_sim.Add(filename.c_str());
        chain_rec.Add(filename.c_str());

        // Join trees
        chain_sim.AddFriend(&chain_rec, "rec");
        ROOT::RDataFrame df(chain_sim);

        auto report = df.Report();
        report->Print();

        DoEfficiencyAndResolution(df);
    } 
    else if (analysis_type == "raw") {
        // Plot raw detector correlations (needs raw data runID extraction from filename or argument)
        // PlotRawCorrelation(0, 0, 1000, 0, 1000); 
    }
    else if (analysis_type == "mult") {
        // PlotMultiplicity(0, 0, 1000, 0, 1000);
    }
    else if (analysis_type == "est") {
        // PlotEstimators(0, 0, 1000, 0, 1000);
    }
        else if (analysis_type == "roc") {
        // Plot ROC curve using multiple input files
        std::vector<std::string> files;
        std::vector<double> thresholds;
        
        // Esempio: se da riga di comando passi "roc file1=0.2 file2=0.4 ..."
        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            size_t pos = arg.find("=");
            if (pos != std::string::npos) {
                files.push_back(arg.substr(0, pos));
                thresholds.push_back(std::stod(arg.substr(pos + 1)));
            }
        }
        
        if (files.empty()) {
            std::cerr << "Please provide files and thresholds! Example: ./Distiller file.root roc file1.root=0.2 file2.root=0.4" << std::endl;
        } else {
            PlotHoughROC(files, thresholds);
        }
    }
    else if (analysis_type == "align") {
        // Placeholder for Systematic Calculator
        // SolveSystematics();
        // FindAlignment(0, 0, 0, 0);
    }
    else {
        std::cerr << "Unknown analysis type: " << analysis_type << std::endl;
        return 1;
    }

    return 0;
}
