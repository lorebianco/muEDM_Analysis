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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT Includes
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TChain.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TProfile.h>

// CHeT Includes
#include <CHeT/CHeTGlobalSettings.hh>

void DoEfficiencyAndResolution(ROOT::RDF::RNode df, bool isMichel)
{
    // 1. Definition of variables
    auto df_base
        = df
              // Momentum
              .Define("trk_p", "sqrt(mc_px*mc_px + mc_py*mc_py + mc_pz*mc_pz)")
              // Theta EDM formula: pz=cos(theta), pr=sin(theta)*(x*cos(phi)+y*sin(phi))/r
              .Define("trk_r", "sqrt(mc_x*mc_x + mc_y*mc_y)")
              .Define("trk_phi", "atan2(mc_py, mc_px)"); // Note: coordinate meaning may vary
                                                         // depending on how mom is saved

    ROOT::RDF::RNode df_def = df_base;
    if(!isMichel)
    {
        df_def = df_base
                     // Difference
                     .Define("diff_x", "rec.rec_x0 - trk_x0")
                     .Define("diff_z", "rec.rec_z0 - trk_z0")
                     .Define("diff_sx", "rec.rec_sx - trk_sx")
                     .Define("diff_sz", "rec.rec_sz - trk_sz");
    }
    else
    {
        df_def = df_base
                     // Difference
                     .Define("diff_px", "rec_extrap_px * 1000.0 - mc_px")
                     .Define("diff_py", "rec_extrap_py * 1000.0 - mc_py")
                     .Define("diff_pz", "rec_extrap_pz * 1000.0 - mc_pz")
                     .Define("diff_x", "rec_extrap_x * 10.0 - mc_x")
                     .Define("diff_y", "rec_extrap_y * 10.0 - mc_y")
                     .Define("diff_z", "rec_extrap_z * 10.0 - mc_z");
    }

    auto df_acc = df_def.Filter("rec_acceptance == true", "Events in acceptance");
    auto df_conv = df_acc.Filter("rec_converged == true", "Converged fit events");

    // 2. Booking Histograms
    ROOT::RDF::RResultPtr<::TH1D> h_diff_x, h_diff_y, h_diff_z, h_diff_sx, h_diff_sz;
    ROOT::RDF::RResultPtr<::TH1D> h_diff_px, h_diff_py, h_diff_pz;

    if(!isMichel)
    {
        h_diff_x = df_conv.Histo1D(
            { "h_diff_x", "Pull X; #Delta X [mm]; Events", 100, -5, 5 }, "diff_x");
        h_diff_z = df_conv.Histo1D(
            { "h_diff_z", "Pull Z; #Delta Z [mm]; Events", 100, -5, 5 }, "diff_z");
        h_diff_sx = df_conv.Histo1D(
            { "h_diff_sx", "Pull sX; #Delta sX; Events", 100, -0.1, 0.1 }, "diff_sx");
        h_diff_sz = df_conv.Histo1D(
            { "h_diff_sz", "Pull sZ; #Delta sZ; Events", 100, -0.1, 0.1 }, "diff_sz");
    }
    else
    {
        h_diff_px = df_conv.Histo1D(
            { "h_diff_px", "Pull Px; #Delta Px [MeV/c]; Events", 100, -50, 50 }, "diff_px");
        h_diff_py = df_conv.Histo1D(
            { "h_diff_py", "Pull Py; #Delta Py [MeV/c]; Events", 100, -50, 50 }, "diff_py");
        h_diff_pz = df_conv.Histo1D(
            { "h_diff_pz", "Pull Pz; #Delta Pz [MeV/c]; Events", 100, -50, 50 }, "diff_pz");
        h_diff_x = df_conv.Histo1D(
            { "h_diff_x", "Pull X; #Delta X [mm]; Events", 100, -50, 50 }, "diff_x");
        h_diff_y = df_conv.Histo1D(
            { "h_diff_y", "Pull Y; #Delta Y [mm]; Events", 100, -50, 50 }, "diff_y");
        h_diff_z = df_conv.Histo1D(
            { "h_diff_z", "Pull Z; #Delta Z [mm]; Events", 100, -50, 50 }, "diff_z");
    }

    // We can also book TEfficiencies directly with DataFrame!
    // auto fill_eff = [](TEfficiency &eff, double p, bool conv) { eff.Fill(conv, p); };
    TEfficiency eff("eff_p", "Efficiency vs Momentum; Momentum [MeV/c]; Efficiency", 20, 0, 70);

    // Instead of directly using Book which can be complex for TEfficiency without a custom helper,
    // we use a TH1D approach for passed/total or we manually loop/reduce.
    // For now let's just make two TH1Ds and divide them to get an efficiency plot using
    // TEfficiency.

    auto h_passed
        = df_conv.Histo1D({ "h_passed", "Passed; Momentum [MeV/c]; Events", 20, 0, 70 }, "trk_p");
    auto h_total
        = df_acc.Histo1D({ "h_total", "Total; Momentum [MeV/c]; Events", 20, 0, 70 }, "trk_p");

    // 3. Draw Results
    TCanvas c("c_resolutions", "Resolutions", 800, 600);
    c.Divide(2, 2);
    if(!isMichel)
    {
        c.cd(1);
        h_diff_x->Draw();
        c.cd(2);
        h_diff_z->Draw();
        c.cd(3);
        h_diff_sx->Draw();
        c.cd(4);
        h_diff_sz->Draw();
    }
    else
    {
        c.Clear();
        c.Divide(3, 2);
        c.cd(1);
        h_diff_px->Draw();
        c.cd(2);
        h_diff_py->Draw();
        c.cd(3);
        h_diff_pz->Draw();
        c.cd(4);
        h_diff_x->Draw();
        c.cd(5);
        h_diff_y->Draw();
        c.cd(6);
        h_diff_z->Draw();
    }
    c.SaveAs("Resolutions.pdf");

    TCanvas c2("c_eff", "Efficiency", 800, 600);
    if(TEfficiency::CheckConsistency(*h_passed, *h_total))
    {
        TEfficiency *pEff = new TEfficiency(*h_passed, *h_total);
        pEff->SetTitle("Efficiency vs Momentum; Momentum [MeV/c]; Efficiency");
        pEff->Draw("AP");
    }
    c2.SaveAs("Efficiency.pdf");

    std::cout << ">>> Plotting finished. Created Resolutions.pdf and Efficiency.pdf" << std::endl;
}

void PlotHoughROC(const std::vector<std::string> &file_list, const std::vector<double> &thresholds)
{
    if(file_list.size() != thresholds.size())
    {
        std::cerr << "Mismatch between files and thresholds count!" << std::endl;
        return;
    }

    std::vector<double> eff_avg;
    std::vector<double> pur_avg;

    for(size_t i = 0; i < file_list.size(); ++i)
    {
        TChain chain_sim("sim");
        TChain chain_rec("rec");
        chain_sim.Add(file_list[i].c_str());
        chain_rec.Add(file_list[i].c_str());
        Long64_t rec_entries = chain_rec.GetEntries();
        chain_sim.AddFriend(&chain_rec, "rec");

        ROOT::RDataFrame df_full(chain_sim);
        auto df = df_full.Filter(
            "EventID < " + std::to_string(rec_entries), "Events processed by Fitter");

        // Custom function to evaluate single event efficiency and purity
        // We compare the true hits (from sim) with the reconstructed hough hits (from rec)
        // Since we are working with vectors of hits, let's calculate intersection sizes

        auto df_eval
            = df.Define("tp",
                    [](const std::vector<int> &sim_hits, const std::vector<int> &rec_hough)
                    {
                        int tp = 0;
                        for(int h : rec_hough)
                        {
                            if(std::find(sim_hits.begin(), sim_hits.end(), h) != sim_hits.end())
                                tp++;
                        }
                        return tp;
                    },
                    { "hit_ids", "rec.rec_hough2d_idx" })
                  .Define("true_total",
                      [](const std::vector<int> &sim_hits) { return (int)sim_hits.size(); },
                      { "hit_ids" })
                  .Define("reco_total",
                      [](const std::vector<int> &rec_hough) { return (int)rec_hough.size(); },
                      { "rec.rec_hough2d_idx" })
                  .Define("eff", "true_total > 0 ? (double)tp / true_total : 0.0")
                  .Define("pur", "reco_total > 0 ? (double)tp / reco_total : 0.0");

        double mean_eff = df_eval.Mean("eff").GetValue();
        double mean_pur = df_eval.Mean("pur").GetValue();

        eff_avg.push_back(mean_eff);
        pur_avg.push_back(mean_pur);

        std::cout << "Threshold: " << thresholds[i] << " -> Eff: " << mean_eff
                  << ", Pur: " << mean_pur << std::endl;
    }

    TCanvas c("c_roc", "ROC Curve", 800, 600);
    TGraph *gROC = new TGraph(eff_avg.size(), pur_avg.data(), eff_avg.data());
    gROC->SetTitle("Hough Transform ROC; Purity; Efficiency");
    gROC->SetMarkerStyle(20);
    gROC->SetLineColor(kBlue);
    gROC->SetMarkerColor(kBlue);
    gROC->Draw("APL");

    // Add labels
    for(size_t i = 0; i < thresholds.size(); ++i)
    {
        auto *t = new TLatex(pur_avg[i], eff_avg[i], Form("%.2f", thresholds[i]));
        t->SetTextSize(0.03);
        t->Draw();
    }

    c.SaveAs("Hough_ROC.pdf");
    std::cout << ">>> Saved Hough_ROC.pdf" << std::endl;
}

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> [analysis_type]" << std::endl;
        std::cerr << "       analysis_type: perf (default), raw, mult, est, align" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string analysis_type = (argc > 2) ? argv[2] : "perf";

    ROOT::EnableImplicitMT();

    if(analysis_type == "perf")
    {
        TChain chain_sim("sim");
        TChain chain_rec("rec");

        chain_sim.Add(filename.c_str());
        chain_rec.Add(filename.c_str());

        Long64_t rec_entries = chain_rec.GetEntries();
        std::cout << ">>> Found " << rec_entries << " events processed by Fitter." << std::endl;

        // Join trees
        chain_sim.AddFriend(&chain_rec, "rec");
        ROOT::RDataFrame df_full(chain_sim);

        // Use Filter to bound the processing to Fitter's output, avoiding Range and IMT crash
        auto df = df_full.Filter(
            "EventID < " + std::to_string(rec_entries), "Events processed by Fitter");

        bool isMichel = chain_rec.GetBranch("rec_R") != nullptr;

        auto report = df.Report();
        report->Print();

        DoEfficiencyAndResolution(df, isMichel);
    }
    else if(analysis_type == "raw")
    {
        // Plot raw detector correlations (needs raw data runID extraction from filename or
        // argument) PlotRawCorrelation(0, 0, 1000, 0, 1000);
    }
    else if(analysis_type == "mult")
    {
        // PlotMultiplicity(0, 0, 1000, 0, 1000);
    }
    else if(analysis_type == "est")
    {
        // PlotEstimators(0, 0, 1000, 0, 1000);
    }
    else if(analysis_type == "roc")
    {
        // Plot ROC curve using multiple input files
        std::vector<std::string> files;
        std::vector<double> thresholds;

        // Esempio: se da riga di comando passi "roc file1=0.2 file2=0.4 ..."
        for(int i = 3; i < argc; ++i)
        {
            std::string arg = argv[i];
            size_t pos = arg.find("=");
            if(pos != std::string::npos)
            {
                files.push_back(arg.substr(0, pos));
                thresholds.push_back(std::stod(arg.substr(pos + 1)));
            }
        }

        if(files.empty())
        {
            std::cerr << "Please provide files and thresholds! Example: ./Distiller file.root roc "
                         "file1.root=0.2 file2.root=0.4"
                      << std::endl;
        }
        else
        {
            PlotHoughROC(files, thresholds);
        }
    }
    else if(analysis_type == "align")
    {
        // Placeholder for Systematic Calculator
        // SolveSystematics();
        // FindAlignment(0, 0, 0, 0);
    }
    else
    {
        std::cerr << "Unknown analysis type: " << analysis_type << std::endl;
        return 1;
    }

    return 0;
}
