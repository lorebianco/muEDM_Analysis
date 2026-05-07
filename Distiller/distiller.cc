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

#include <RtypesCore.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

// ROOT Includes
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TSystem.h>

#include "lbrootstyle.hh"

// CHeT Includes
#include <CHeT/CHeTGlobalSettings.hh>

#include "lbrootstyle.hh"

void SetCHeTOptimalSelfAligment()
{
    // Config::SetOffsetExp(29. * TMath::DegToRad()); // 29 optimal
    CHeT::Config::SetDeltaI(1, CHeT::Config::GetDeltaI(1) + 4.5 * TMath::DegToRad()); // 4.5 optimal
}

struct EffStats
{
    long passed = 0;
    long total = 0;
};
using EffMap = std::map<int, EffStats>;

std::mutex eff_mutex;

std::pair<int, int> AccumulateEfficiency(double x0, double z0, double sx, double sz,
    const std::vector<int> &hit_ids, EffMap &bundleMap, EffMap &layerMap, EffMap &cylMap)
{
    int track_passed = 0;
    int track_total = 0;
    const double L_SAFE = CHeT::Config::L_HALF - 5.0;
    const int N_NEIGHBORS = 1;
    auto cylinders = CHeT::Config::GetCylinders();
    int global_bundle_offset = 0;
    int cyl_idx = 0;

    for(const auto &cyl : cylinders)
    {
        const CHeT::Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int lay_local_idx = 0; lay_local_idx < 2; ++lay_local_idx)
        {
            const auto &L = *layers[lay_local_idx];
            double R = L.radius;
            int lay_global_idx = cyl_idx * 2 + lay_local_idx;
            int current_layer_start = global_bundle_offset;
            double A = sx * sx + 1.0;
            double B = 2.0 * x0 * sx;
            double C = x0 * x0 - R * R;
            double delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                double y_sol[2] = { (-B + sqrt(delta)) / (2 * A), (-B - sqrt(delta)) / (2 * A) };
                for(double y : y_sol)
                {
                    double z = z0 + sz * y;
                    if(abs(z) <= L_SAFE)
                    {
                        double x = x0 + sx * y;
                        double phi_hit = atan2(y, x);
                        int best_bundle_id = -1;
                        double min_dist_norm = 1e9;
                        double bundle_width_rad = 2.0 * M_PI / L.nBundles;
                        for(int b = 0; b < L.nBundles; ++b)
                        {
                            int gid = current_layer_start + b;
                            auto prop = CHeT::Config::GetFiberProp(gid);
                            double alpha
                                = (z + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
                            double phi_fib = prop.phi0 + prop.dir * alpha * M_PI;
                            double dphi = abs(CHeT::Config::wrap0_2pi(phi_hit)
                                - CHeT::Config::wrap0_2pi(phi_fib));
                            if(dphi > M_PI)
                                dphi = 2.0 * M_PI - dphi;
                            double dist_norm = dphi / bundle_width_rad;
                            if(dist_norm < min_dist_norm)
                            {
                                min_dist_norm = dist_norm;
                                best_bundle_id = gid;
                            }
                        }

                        const double TOL_BUNDLE = 1.5;
                        if(best_bundle_id >= 0 && min_dist_norm < TOL_BUNDLE)
                        {
                            int hit_bundle_id = -1;
                            for(int h : hit_ids)
                            {
                                if(h == best_bundle_id)
                                {
                                    hit_bundle_id = best_bundle_id;
                                    break;
                                }
                            }
                            if(hit_bundle_id == -1)
                            {
                                auto best_prop = CHeT::Config::GetFiberProp(best_bundle_id);
                                for(int h : hit_ids)
                                {
                                    auto hit_prop = CHeT::Config::GetFiberProp(h);
                                    if(hit_prop.cylinderId != best_prop.cylinderId
                                        || hit_prop.layerId != best_prop.layerId)
                                        continue;
                                    int diff = abs(h - best_bundle_id);
                                    if((diff <= N_NEIGHBORS) || (diff >= L.nBundles - N_NEIGHBORS))
                                    {
                                        hit_bundle_id = h;
                                        break;
                                    }
                                }
                            }
                            int id_to_fill = (hit_bundle_id != -1) ? hit_bundle_id : best_bundle_id;
                            {
                                std::lock_guard<std::mutex> lock(eff_mutex);
                                bundleMap[id_to_fill].total++;
                                layerMap[lay_global_idx].total++;
                                cylMap[cyl_idx].total++;
                                track_total++;
                                if(hit_bundle_id != -1)
                                {
                                    bundleMap[id_to_fill].passed++;
                                    layerMap[lay_global_idx].passed++;
                                    cylMap[cyl_idx].passed++;
                                    track_passed++;
                                }
                            }
                        }
                    }
                }
            }
            global_bundle_offset += L.nBundles;
        }
        cyl_idx++;
    }
    return { track_passed, track_total };
}

std::pair<double, double> GetCPErrors(long passed, long total)
{
    if(total == 0)
        return { 0.0, 0.0 };
    double eff = (double)passed / total;
    double cl = 0.6827; // 1 sigma
    double lower = TEfficiency::ClopperPearson(total, passed, cl, false);
    double upper = TEfficiency::ClopperPearson(total, passed, cl, true);
    return { eff - lower, upper - eff };
}

void PlotEfficiencyResults(int runID, const EffMap &bMap, const EffMap &lMap, const EffMap &cMap,
    const std::map<int, int> &occMap)
{
    auto cylinders = CHeT::Config::GetCylinders();
    int num_cyls = cylinders.size();
    int num_lays = num_cyls * 2;

    std::vector<int> col_cyl;
    std::vector<std::string> lab_cyl;
    std::vector<int> col_lay;
    std::vector<std::string> lab_lay;

    for(const auto &cyl : cylinders)
    {
        col_cyl.push_back(cyl.inner.color);
        lab_cyl.push_back(Form("Cylinder %d", cyl.id));

        col_lay.push_back(cyl.inner.color);
        lab_lay.push_back(Form("Cyl %d - In", cyl.id));

        col_lay.push_back(cyl.outer.color);
        lab_lay.push_back(Form("Cyl %d - Out", cyl.id));
    }

    // ---------------------------------------------------------
    // 0. NEW: Noise / Occupancy Canvas
    // ---------------------------------------------------------
    TCanvas *cNoise = new TCanvas("cNoise", "Noise and Occupancy Map", 800, 800);
    cNoise->cd();
    cNoise->SetRightMargin(0.15);

    TH2Poly *h_noise_poly = new TH2Poly();
    h_noise_poly->SetName("h_noise_map");
    h_noise_poly->SetTitle("Unassigned Hits Map (z=0);X [mm];Y [mm]");

    int current_bundle_id_noise = 0;
    double max_R_noise = 0;
    for(const auto &cyl : cylinders)
    {
        if(cyl.outer.radius > max_R_noise)
            max_R_noise = cyl.outer.radius;
    }

    for(const auto &cyl : cylinders)
    {
        const CHeT::Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int lay_idx = 0; lay_idx < 2; ++lay_idx)
        {
            const auto &L = *layers[lay_idx];
            double dPhi = 2.0 * M_PI / L.nBundles;

            double t_draw = std::max(L.thickness, max_R_noise * 0.02);
            double R_cyl = (cyl.inner.radius + cyl.outer.radius) / 2.0;
            double R_draw = R_cyl;
            if(lay_idx == 0)
                R_draw -= t_draw * 0.6;
            else
                R_draw += t_draw * 0.6;

            for(int b = 0; b < L.nBundles; ++b)
            {
                auto prop = CHeT::Config::GetFiberProp(current_bundle_id_noise);
                double phi_center = prop.phi0 + prop.dir * M_PI / 2.0;

                double phi_gap = dPhi * 0.02;
                double phi1 = phi_center - dPhi / 2.0 + phi_gap;
                double phi2 = phi_center + dPhi / 2.0 - phi_gap;

                double r1 = R_draw - t_draw / 2.0;
                double r2 = R_draw + t_draw / 2.0;

                double px[4] = { r1 * cos(phi1), r2 * cos(phi1), r2 * cos(phi2), r1 * cos(phi2) };
                double py[4] = { r1 * sin(phi1), r2 * sin(phi1), r2 * sin(phi2), r1 * sin(phi2) };

                int binID = h_noise_poly->AddBin(4, px, py);

                int total = (occMap.count(current_bundle_id_noise))
                    ? occMap.at(current_bundle_id_noise)
                    : 0;
                int matched = (bMap.count(current_bundle_id_noise))
                    ? bMap.at(current_bundle_id_noise).passed
                    : 0;
                int noise = total - matched;
                if(noise < 0)
                    noise = 0;

                h_noise_poly->SetBinContent(binID, noise);

                current_bundle_id_noise++;
            }
        }
    }

    double plot_range_noise = max_R_noise + 5.0;
    TH1 *frame_noise = cNoise->DrawFrame(
        -plot_range_noise, -plot_range_noise, plot_range_noise, plot_range_noise);
    frame_noise->SetTitle("Unassigned Hits Map (z=0);X [mm];Y [mm]");
    h_noise_poly->SetStats(0);
    gStyle->SetPalette(kBird);
    h_noise_poly->Draw("COLZ L SAME");
    cNoise->Modified();
    cNoise->Update();
    cNoise->SaveAs("Efficiency_Noise.pdf");
    cNoise->SaveAs("Efficiency_Noise.png");

    // ---------------------------------------------------------
    // Canvas Efficiencies
    // ---------------------------------------------------------
    TCanvas *cEff = new TCanvas("cEff", Form("Efficiencies Run %d", runID), 900, 1200);
    cEff->Divide(1, 3);

    // ---------------------------------------------------------
    // 1. Cylinder Efficiency
    // ---------------------------------------------------------
    cEff->cd(1);
    gPad->SetGridy();

    double min_cyl_eff = 1.0;
    for(const auto &it : cMap)
    {
        if(it.second.total > 0)
        {
            double eff = (double)it.second.passed / it.second.total;
            auto errs = GetCPErrors(it.second.passed, it.second.total);
            if(eff - errs.first < min_cyl_eff)
                min_cyl_eff = eff - errs.first;
        }
    }
    double min_cyl_draw = std::max(0.0, min_cyl_eff - 0.05);

    TH1D *h_frame_cyl = new TH1D(
        "h_frame_cyl", "Cylinder Efficiency;;Efficiency", num_cyls, -0.5, num_cyls - 0.5);
    h_frame_cyl->SetMinimum(min_cyl_draw);
    h_frame_cyl->SetMaximum(1.05);
    h_frame_cyl->SetStats(0);
    for(int i = 0; i < num_cyls; ++i)
        h_frame_cyl->GetXaxis()->SetBinLabel(i + 1, lab_cyl[i].c_str());
    h_frame_cyl->GetXaxis()->SetLabelSize(0.07);
    h_frame_cyl->Draw("AXIS");

    TLine *lineC = new TLine(-0.5, 1.0, num_cyls - 0.5, 1.0);
    lineC->SetLineStyle(2);
    lineC->SetLineColor(kGray);
    lineC->Draw();

    std::vector<TH1D *> h_cyls(num_cyls);
    std::vector<TGraphAsymmErrors *> g_err_cyl(num_cyls);

    // --- Global Efficiency Calculation ---
    long g_passed = 0;
    long g_total = 0;

    for(const auto &it : cMap)
    {
        g_passed += it.second.passed;
        g_total += it.second.total;
    }

    if(g_total > 0)
    {
        double g_eff = (double)g_passed / g_total;

        TLatex tG;
        tG.SetNDC();
        tG.SetTextSize(0.08);
        tG.SetTextColor(kBlack);
        tG.DrawLatex(0.50, 0.75, Form("Global Efficiency: %.2f%%", g_eff * 100));

        printf("\n--- GLOBAL EFFICIENCY ---\n %.2f%% (%ld/%ld)\n", g_eff * 100, g_passed, g_total);
    }

    printf("\n--- CYLINDER EFFICIENCY ---\n");

    for(int i = 0; i < num_cyls; ++i)
    {
        h_cyls[i] = new TH1D(Form("h_cyl_%d", i), "", num_cyls, -0.5, num_cyls - 0.5);
        g_err_cyl[i] = new TGraphAsymmErrors();

        auto it = cMap.find(i);
        if(it != cMap.end() && it->second.total > 0)
        {
            auto errs = GetCPErrors(it->second.passed, it->second.total);
            double eff = (double)it->second.passed / it->second.total;

            h_cyls[i]->SetBinContent(i + 1, eff);

            g_err_cyl[i]->SetPoint(0, (double)i, eff);
            g_err_cyl[i]->SetPointError(0, 0.5, 0.5, errs.first, errs.second);

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_cyl[i].c_str(), eff * 100,
                errs.second * 100, errs.first * 100, it->second.passed, it->second.total);
        }

        h_cyls[i]->SetFillColorAlpha(col_cyl[i], 0.35);
        h_cyls[i]->SetLineColor(col_cyl[i]);
        h_cyls[i]->SetLineWidth(2);
        h_cyls[i]->SetStats(0);

        g_err_cyl[i]->SetLineColor(col_cyl[i]);
        g_err_cyl[i]->SetLineWidth(2);
        g_err_cyl[i]->SetMarkerStyle(0);

        h_cyls[i]->Draw("SAME HIST BAR");
        g_err_cyl[i]->Draw("P E1 SAME");
    }
    h_frame_cyl->Draw("AXIS SAME");

    // ---------------------------------------------------------
    // 2. Layer Efficiency
    // ---------------------------------------------------------
    cEff->cd(2);
    gPad->SetGridy();

    double min_lay_eff = 1.0;
    for(const auto &it : lMap)
    {
        if(it.second.total > 0)
        {
            double eff = (double)it.second.passed / it.second.total;
            auto errs = GetCPErrors(it.second.passed, it.second.total);
            if(eff - errs.first < min_lay_eff)
                min_lay_eff = eff - errs.first;
        }
    }
    double min_lay_draw = std::max(0.0, min_lay_eff - 0.05);

    TH1D *h_frame_lay = new TH1D(
        "h_frame_lay", "Sub-Layer Efficiency;;Efficiency", num_lays, -0.5, num_lays - 0.5);
    h_frame_lay->SetMinimum(min_lay_draw);
    h_frame_lay->SetMaximum(1.05);
    h_frame_lay->SetStats(0);
    for(int i = 0; i < num_lays; ++i)
        h_frame_lay->GetXaxis()->SetBinLabel(i + 1, lab_lay[i].c_str());
    h_frame_lay->GetXaxis()->SetLabelSize(0.07);
    h_frame_lay->Draw("AXIS");

    TLine *lineL = new TLine(-0.5, 1.0, num_lays - 0.5, 1.0);
    lineL->SetLineStyle(2);
    lineL->SetLineColor(kGray);
    lineL->Draw();

    std::vector<TH1D *> h_layers(num_lays);
    std::vector<TGraphAsymmErrors *> g_err_lay(num_lays);

    printf("\n--- LAYER EFFICIENCY ---\n");

    for(int i = 0; i < num_lays; ++i)
    {
        h_layers[i] = new TH1D(Form("h_lay_%d", i), "", num_lays, -0.5, num_lays - 0.5);
        g_err_lay[i] = new TGraphAsymmErrors();

        auto it = lMap.find(i);
        if(it != lMap.end() && it->second.total > 0)
        {
            auto errs = GetCPErrors(it->second.passed, it->second.total);
            double eff = (double)it->second.passed / it->second.total;

            h_layers[i]->SetBinContent(i + 1, eff);
            g_err_lay[i]->SetPoint(0, (double)i, eff);
            g_err_lay[i]->SetPointError(0, 0.5, 0.5, errs.first, errs.second);

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_lay[i].c_str(), eff * 100,
                errs.second * 100, errs.first * 100, it->second.passed, it->second.total);
        }

        h_layers[i]->SetFillColorAlpha(col_lay[i], 0.35);
        h_layers[i]->SetLineColor(col_lay[i]);
        h_layers[i]->SetLineWidth(2);
        h_layers[i]->SetStats(0);

        g_err_lay[i]->SetLineColor(col_lay[i]);
        g_err_lay[i]->SetLineWidth(2);
        g_err_lay[i]->SetMarkerStyle(0);

        h_layers[i]->Draw("SAME HISTO BAR");
        g_err_lay[i]->Draw("P E1 SAME");
    }
    h_frame_lay->Draw("AXIS SAME");

    // ---------------------------------------------------------
    // 3. Bundle Efficiency
    // ---------------------------------------------------------
    cEff->cd(3);
    gPad->SetGridy();

    TGraphAsymmErrors *g_bund = new TGraphAsymmErrors();
    int pt = 0;

    for(auto const &[id, cnt] : bMap)
    {
        if(cnt.total > 0)
        {
            double eff = (double)cnt.passed / cnt.total;
            auto errs = GetCPErrors(cnt.passed, cnt.total); // {low, high}

            g_bund->SetPoint(pt, (double)id, eff);
            g_bund->SetPointError(pt, 0, 0, errs.first, errs.second);
            // printf("Bundle %d eff %.2f total %ld\n", id, eff, cnt.total);
            pt++;
        }
    }

    double min_bund_eff = 1.0;
    for(const auto &it : bMap)
    {
        if(it.second.total > 0)
        {
            double eff = (double)it.second.passed / it.second.total;
            auto errs = GetCPErrors(it.second.passed, it.second.total);
            if(eff - errs.first < min_bund_eff)
                min_bund_eff = eff - errs.first;
        }
    }
    double min_bund_draw = std::max(0.0, min_bund_eff - 0.05);

    int max_id = -1;
    for(const auto &cyl : cylinders)
        max_id += cyl.inner.nBundles + cyl.outer.nBundles;

    TH1D *h_frame_bund = new TH1D("h_frame_bund", "Bundle Efficiency;Global Bundle ID;Efficiency",
        max_id + 1, -0.5, max_id + 0.5);
    h_frame_bund->SetMinimum(min_bund_draw);
    h_frame_bund->SetMaximum(1.05);
    h_frame_bund->SetStats(0);
    h_frame_bund->Draw("AXIS");

    // Ridisegno separatori anche qui
    int current_offset = 0;
    double yMin = h_frame_bund->GetMinimum();
    double yMax = h_frame_bund->GetMaximum();
    for(const auto &cyl : cylinders)
    {
        int nIn = cyl.inner.nBundles;
        double x_in = current_offset + nIn - 0.5;
        TLine *lIn = new TLine(x_in, yMin, x_in, yMax);
        lIn->SetLineStyle(2);
        lIn->SetLineColor(cyl.inner.color);
        lIn->SetLineWidth(2);
        lIn->Draw();
        current_offset += nIn;

        int nOut = cyl.outer.nBundles;
        double x_out = current_offset + nOut - 0.5;
        TLine *lOut = new TLine(x_out, yMin, x_out, yMax);
        lOut->SetLineStyle(2);
        lOut->SetLineColor(cyl.outer.color);
        lOut->SetLineWidth(2);
        lOut->Draw();
        current_offset += nOut;
    }

    g_bund->SetMarkerStyle(20);
    g_bund->SetMarkerSize(0.6);
    g_bund->Draw("P E1 SAME");

    TLine *line1 = new TLine(-0.5, 1.0, max_id + 0.5, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray);
    line1->Draw();
    h_frame_bund->Draw("AXIS SAME");

    // ---------------------------------------------------------
    // 4. 2D Detector Map (Efficiency at z=0)
    // ---------------------------------------------------------
    TCanvas *cMap2D = new TCanvas("cMap2D", "Detector Efficiency Map", 800, 800);
    cMap2D->cd();
    cMap2D->SetRightMargin(0.15);

    TH2Poly *h_poly = new TH2Poly();
    h_poly->SetName("h_eff_map");
    h_poly->SetTitle("Detector Efficiency Map (z=0);X [mm];Y [mm]");

    int current_bundle_id = 0;
    double max_R = 0;
    for(const auto &cyl : cylinders)
    {
        if(cyl.outer.radius > max_R)
            max_R = cyl.outer.radius;
    }

    for(const auto &cyl : cylinders)
    {
        const CHeT::Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int lay_idx = 0; lay_idx < 2; ++lay_idx)
        {
            const auto &L = *layers[lay_idx];
            double dPhi = 2.0 * M_PI / L.nBundles;

            // Artificially space out inner/outer layers to allow larger thickness
            // without overlap
            double t_draw = std::max(L.thickness, max_R * 0.02);
            double R_cyl = (cyl.inner.radius + cyl.outer.radius) / 2.0;
            double R_draw = R_cyl;
            if(lay_idx == 0)
                R_draw -= t_draw * 0.6;
            else
                R_draw += t_draw * 0.6;

            for(int b = 0; b < L.nBundles; ++b)
            {
                auto prop = CHeT::Config::GetFiberProp(current_bundle_id);
                // The phi angle at z=0 is phi0 + dir * (M_PI / 2)
                double phi_center = prop.phi0 + prop.dir * M_PI / 2.0;

                // Add a small gap between bundles azimuthally for better visualization
                double phi_gap = dPhi * 0.02;
                double phi1 = phi_center - dPhi / 2.0 + phi_gap;
                double phi2 = phi_center + dPhi / 2.0 - phi_gap;

                double r1 = R_draw - t_draw / 2.0;
                double r2 = R_draw + t_draw / 2.0;

                double px[4] = { r1 * cos(phi1), r2 * cos(phi1), r2 * cos(phi2), r1 * cos(phi2) };
                double py[4] = { r1 * sin(phi1), r2 * sin(phi1), r2 * sin(phi2), r1 * sin(phi2) };

                int binID = h_poly->AddBin(4, px, py);

                auto it = bMap.find(current_bundle_id);
                if(it != bMap.end() && it->second.total > 0)
                {
                    double eff = (double)it->second.passed / it->second.total;
                    h_poly->SetBinContent(binID, eff * 100.0);
                }

                current_bundle_id++;
            }
        }
    }

    double plot_range = max_R + 5.0;
    TH1 *frame = cMap2D->DrawFrame(-plot_range, -plot_range, plot_range, plot_range);
    frame->SetTitle("Detector Efficiency Map (z=0);X [mm];Y [mm]");
    h_poly->SetMinimum(20.0);
    h_poly->SetMaximum(100.0);
    h_poly->SetStats(0);

    // Colormap "Semaforo" (Rosso -> Giallo -> Verde)
    const Int_t NRGBs = 3;
    const Int_t NCont = 16;
    // Punti di controllo: 0% (rosso), 70% (giallo), 100% (verde scuro per
    // visibilità)
    Double_t stops[NRGBs] = { 0.00, 0.70, 1.00 };
    Double_t red[NRGBs] = { 0.80, 1.00, 0.00 };
    Double_t green[NRGBs] = { 0.00, 0.80, 0.60 };
    Double_t blue[NRGBs] = { 0.00, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    h_poly->Draw("COLZ L SAME");

    cEff->SaveAs("Efficiency_Canvases.pdf");
    cEff->SaveAs("Efficiency_Canvases.png");
    cMap2D->SaveAs("Efficiency_Map2D.pdf");
    cMap2D->SaveAs("Efficiency_Map2D.png");

    printf("------------------------\n");
}

auto fitResidui = [](auto h, bool useDoubleGauss, double limit = 0.0)
{
    if(!h || h->GetEntries() < 50)
        return;

    if(limit > 0.0)
        h->GetXaxis()->SetRangeUser(-limit, limit);

    if(!useDoubleGauss)
    {
        h->Fit("gaus", "QL");
        if(h->GetFunction("gaus"))
            h->GetFunction("gaus")->SetLineColor(kRed);
    }
    else
    {
        auto f2 = std::make_unique<TF1>("f2",
            "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
            h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

        double mean = h->GetMean();
        double rms = h->GetRMS();
        double peak = h->GetMaximum();

        f2->SetParNames("A_Core", "Mean", "Sigma_Core", "A_Tail", "Sigma_Tail");
        f2->SetParameter(0, peak);
        f2->SetParameter(1, mean);
        f2->SetParameter(2, rms * 0.5);
        f2->SetParameter(3, peak * 0.01);
        f2->SetParameter(4, rms * 2.0);

        f2->SetParLimits(0, peak * 0.75, peak * 2);
        f2->SetParLimits(2, 0, rms);
        f2->SetParLimits(3, 0, peak * 0.25);
        f2->SetParLimits(4, rms, rms * 10);

        f2->SetLineColor(kRed);
        f2->SetNpx(1000);

        auto fitptr = h->Fit(f2.get(), "SLMQ");

        if(fitptr->IsValid())
        {
            auto fCore = std::make_unique<TF1>(
                "fCore", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            fCore->SetParameters(f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2));
            fCore->SetLineColor(kGreen + 2);
            fCore->SetLineStyle(2);
            fCore->DrawCopy("same");

            auto fTail = std::make_unique<TF1>(
                "fTail", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            fTail->SetParameters(f2->GetParameter(3), f2->GetParameter(1), f2->GetParameter(4));
            fTail->SetLineColor(kMagenta);
            fTail->SetLineStyle(2);
            fTail->DrawCopy("same");

            f2->DrawCopy("same");
        }
        peak = h->GetMaximum();

        // Dynamic Y-axis range
        h->SetMaximum(peak * 1.15);
    }
};

void DoResolutions(ROOT::RDF::RNode df, bool isMichel)
{
    // 1. Definition of variables
    auto df_base = df
                       // Momentum
                       .Define("trk_p", "sqrt(mc_px*mc_px + mc_py*mc_py + mc_pz*mc_pz)")
                       // Theta EDM formula: pz=cos(theta),
                       // pr=sin(theta)*(x*cos(phi)+y*sin(phi))/r
                       .Define("trk_r", "sqrt(mc_x*mc_x + mc_y*mc_y)")
                       .Define("trk_phi",
                           "atan2(mc_py, mc_px)"); // Note: coordinate meaning may vary
                                                   // depending on how mom is saved

    ROOT::RDF::RNode df_def = df_base;
    if(!isMichel)
    {
        df_def = df_base
                     // Difference
                     .Define("diff_x", "rec.rec_x0 - trk_x0")
                     .Define("diff_z", "rec.rec_z0 - trk_z0")
                     .Define("diff_sx", "rec.rec_sx - trk_sx")
                     .Define("diff_sz", "rec.rec_sz - trk_sz")
                     .Define("diff_zi",
                         [](const std::vector<double> &rec_zi, const std::vector<double> &mc_z)
                         {
                             std::vector<double> diffs;
                             size_t n = std::min(rec_zi.size(), mc_z.size());
                             for(size_t i = 0; i < n; ++i)
                                 diffs.push_back(rec_zi[i] - mc_z[i]);
                             return diffs;
                         },
                         { "rec.rec_zi", "mc_hits_z" });
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
    ROOT::RDF::RResultPtr<::TH1D> h_diff_x, h_diff_y, h_diff_z, h_diff_sx, h_diff_sz, h_diff_zi;
    ROOT::RDF::RResultPtr<::TH1D> h_diff_px, h_diff_py, h_diff_pz;

    if(!isMichel)
    {
        h_diff_x = df_conv.Histo1D(
            { "h_diff_x", "Pull X; #Delta X [mm]; Events", 100, -1, 1 }, "diff_x");
        h_diff_z = df_conv.Histo1D(
            { "h_diff_z", "Pull Z; #Delta Z [mm]; Events", 100, -5, 5 }, "diff_z");
        h_diff_sx = df_conv.Histo1D(
            { "h_diff_sx", "Pull sX; #Delta sX; Events", 100, -0.025, 0.025 }, "diff_sx");
        h_diff_sz = df_conv.Histo1D(
            { "h_diff_sz", "Pull sZ; #Delta sZ; Events", 100, -0.1, 0.1 }, "diff_sz");
        h_diff_zi = df_conv.Histo1D(
            { "h_diff_zi", "Pull z_i; #Delta z_i [mm]; Hits", 100, -5, 5 }, "diff_zi");
    }
    else
    {
        h_diff_px = df_conv.Histo1D(
            { "h_diff_px", "Pull Px; #Delta Px [MeV/c]; Events", 200, -50, 50 }, "diff_px");
        h_diff_py = df_conv.Histo1D(
            { "h_diff_py", "Pull Py; #Delta Py [MeV/c]; Events", 200, -50, 50 }, "diff_py");
        h_diff_pz = df_conv.Histo1D(
            { "h_diff_pz", "Pull Pz; #Delta Pz [MeV/c]; Events", 200, -50, 50 }, "diff_pz");
        h_diff_x = df_conv.Histo1D(
            { "h_diff_x", "Pull X; #Delta X [mm]; Events", 200, -50, 50 }, "diff_x");
        h_diff_y = df_conv.Histo1D(
            { "h_diff_y", "Pull Y; #Delta Y [mm]; Events", 200, -50, 50 }, "diff_y");
        h_diff_z = df_conv.Histo1D(
            { "h_diff_z", "Pull Z; #Delta Z [mm]; Events", 200, -50, 50 }, "diff_z");
    }

    // We can also book TEfficiencies directly with DataFrame!
    // auto fill_eff = [](TEfficiency &eff, double p, bool conv) { eff.Fill(conv,
    // p); };
    TEfficiency eff("eff_p", "Efficiency vs Momentum; Momentum [MeV/c]; Efficiency", 20, 0, 70);

    // Instead of directly using Book which can be complex for TEfficiency without
    // a custom helper, we use a TH1D approach for passed/total or we manually
    // loop/reduce. For now let's just make two TH1Ds and divide them to get an
    // efficiency plot using TEfficiency.

    auto h_passed
        = df_conv.Histo1D({ "h_passed", "Passed; Momentum [MeV/c]; Events", 20, 0, 70 }, "trk_p");
    auto h_total
        = df_acc.Histo1D({ "h_total", "Total; Momentum [MeV/c]; Events", 20, 0, 70 }, "trk_p");

    // 3. Draw Results
    TCanvas c("c_resolutions", "Resolutions", 1200, 800);
    if(!isMichel)
    {
        c.Divide(2, 2);

        double lim_pos = std::max({ std::abs(h_diff_x->GetMean()) + 5 * h_diff_x->GetRMS(),
            std::abs(h_diff_z->GetMean()) + 5 * h_diff_z->GetRMS() });
        double lim_ang = std::max({ std::abs(h_diff_sx->GetMean()) + 5 * h_diff_sx->GetRMS(),
            std::abs(h_diff_sz->GetMean()) + 5 * h_diff_sz->GetRMS() });

        c.cd(1);
        fitResidui(h_diff_x.GetPtr(), true, lim_pos);
        c.cd(2);
        fitResidui(h_diff_z.GetPtr(), true, lim_pos);
        c.cd(3);
        fitResidui(h_diff_sx.GetPtr(), true, lim_ang);
        c.cd(4);
        fitResidui(h_diff_sz.GetPtr(), true, lim_ang);

        c.SaveAs("Resolutions.pdf");
        c.SaveAs("Resolutions.png");

        TCanvas c_zi("c_zi", "Z_i Resolutions", 800, 600);
        c_zi.cd();
        double lim_zi = std::abs(h_diff_zi->GetMean()) + 5 * h_diff_zi->GetRMS();
        fitResidui(h_diff_zi.GetPtr(), true, lim_zi);
        c_zi.SaveAs("Resolutions_zi.pdf");
        c_zi.SaveAs("Resolutions_zi.png");
    }
    else
    {
        c.Divide(3, 2);

        double lim_p = std::max({ std::abs(h_diff_px->GetMean()) + 5 * h_diff_px->GetRMS(),
            std::abs(h_diff_py->GetMean()) + 5 * h_diff_py->GetRMS(),
            std::abs(h_diff_pz->GetMean()) + 5 * h_diff_pz->GetRMS() });
        double lim_pos = std::max({ std::abs(h_diff_x->GetMean()) + 5 * h_diff_x->GetRMS(),
            std::abs(h_diff_y->GetMean()) + 5 * h_diff_y->GetRMS(),
            std::abs(h_diff_z->GetMean()) + 5 * h_diff_z->GetRMS() });

        c.cd(1);
        fitResidui(h_diff_px.GetPtr(), true, lim_p);
        c.cd(2);
        fitResidui(h_diff_py.GetPtr(), true, lim_p);
        c.cd(3);
        fitResidui(h_diff_pz.GetPtr(), true, lim_p);
        c.cd(4);
        fitResidui(h_diff_x.GetPtr(), true, lim_pos);
        c.cd(5);
        fitResidui(h_diff_y.GetPtr(), true, lim_pos);
        c.cd(6);
        fitResidui(h_diff_z.GetPtr(), true, lim_pos);

        c.SaveAs("Resolutions.pdf");
        c.SaveAs("Resolutions.png");
    }

    TCanvas c2("c_eff", "Efficiency", 800, 600);
    if(TEfficiency::CheckConsistency(*h_passed, *h_total))
    {
        TEfficiency *pEff = new TEfficiency(*h_passed, *h_total);
        pEff->SetTitle("Efficiency vs Momentum; Momentum [MeV/c]; Efficiency");
        pEff->Draw("AP");
    }
    c2.SaveAs("Efficiency.pdf");
    c2.SaveAs("Efficiency.png");

    std::cout << ">>> Plotting finished. Created Resolutions.pdf and Efficiency.pdf" << std::endl;
}

void PlotEfficiency(ROOT::RDF::RNode df, bool hasSim)
{
    std::cout << ">>> Plotting efficiency... (placeholder for full implementation)" << std::endl;
    EffMap bundleStats, layerStats, cylStats;
    std::map<int, int> occupancyMap;

    TProfile *p_eff_x0
        = new TProfile("p_eff_x0", "Track Efficiency vs x0;x0 [mm];Efficiency", 50, -0, 0);
    TProfile *p_eff_z0
        = new TProfile("p_eff_z0", "Track Efficiency vs z0;z0 [mm];Efficiency", 50, -0, 0);
    TProfile *p_eff_sx
        = new TProfile("p_eff_sx", "Track Efficiency vs sx;sx;Efficiency", 50, -0, 0);
    TProfile *p_eff_sz
        = new TProfile("p_eff_sz", "Track Efficiency vs sz;sz;Efficiency", 50, -0, 0);

    TProfile2D *p2_eff_x0_sx = new TProfile2D(
        "p2_eff_x0_sx", "Track Efficiency vs x0 and sx;x0 [mm];sx", 30, -80, 80, 30, -0.6, 0.6);
    TProfile2D *p2_eff_z0_sz = new TProfile2D(
        "p2_eff_z0_sz", "Track Efficiency vs z0 and sz;z0 [mm];sz", 30, -160, 160, 30, -1.0, 1.0);

    auto df_acc = df.Filter("rec_acceptance == true", "Events in acceptance");
    auto df_conv = df_acc.Filter("rec_converged == true", "Converged fit events");

    auto colNames = df.GetColumnNames();
    std::cout << "Available columns: ";
    for(const auto &c : colNames)
        std::cout << c << " ";
    std::cout << std::endl;

    std::vector<std::string> cols;
    if(hasSim)
    {
        std::string hits_col = "rec.rec_hits";
        if(std::find(colNames.begin(), colNames.end(), "hit_ids") != colNames.end())
        {
            hits_col = "hit_ids";
        }
        else if(std::find(colNames.begin(), colNames.end(), "hits") != colNames.end())
        {
            hits_col = "hits";
        }
        cols = { "rec.rec_x0", "rec.rec_z0", "rec.rec_sx", "rec.rec_sz", hits_col };
    }
    else
    {
        std::string hits_col = "rec_hits";
        auto it = std::find(colNames.begin(), colNames.end(), "rec_hits_idx");
        if(it != colNames.end())
        {
            hits_col = "rec_hits_idx";
        }
        cols = { "rec_x0", "rec_z0", "rec_sx", "rec_sz", hits_col };
    }

    df_conv.Foreach(
        [&](double x0, double z0, double sx, double sz, const std::vector<int> &hits)
        {
            {
                std::lock_guard<std::mutex> lock(eff_mutex);
                for(int h : hits)
                    occupancyMap[h]++;
            }
            auto [passed, total]
                = AccumulateEfficiency(x0, z0, sx, sz, hits, bundleStats, layerStats, cylStats);
            if(total > 0)
            {
                double eff = (double)passed / total;
                std::lock_guard<std::mutex> lock(eff_mutex);
                p_eff_x0->Fill(x0, eff);
                p_eff_z0->Fill(z0, eff);
                p_eff_sx->Fill(sx, eff);
                p_eff_sz->Fill(sz, eff);
                p2_eff_x0_sx->Fill(x0, sx, eff);
                p2_eff_z0_sz->Fill(z0, sz, eff);
            }
        },
        cols);

    TCanvas *cEffParams = new TCanvas("cEffParams", "Track Efficiency vs Parameters", 1200, 800);
    cEffParams->Divide(2, 2);

    cEffParams->cd(1);
    p_eff_x0->Draw();

    cEffParams->cd(2);
    p_eff_z0->Draw();

    cEffParams->cd(3);
    p_eff_sx->Draw();

    cEffParams->cd(4);
    p_eff_sz->Draw();

    cEffParams->SaveAs("Efficiency_vs_TrackParams.pdf");
    cEffParams->SaveAs("Efficiency_vs_TrackParams.png");

    TCanvas *cEffParams2D
        = new TCanvas("cEffParams2D", "Track Efficiency vs Parameters 2D", 1200, 600);
    cEffParams2D->Divide(2, 1);

    cEffParams2D->cd(1);
    gPad->SetRightMargin(0.15);
    p2_eff_x0_sx->Draw("COLZ");

    cEffParams2D->cd(2);
    gPad->SetRightMargin(0.15);
    p2_eff_z0_sz->Draw("COLZ");

    cEffParams2D->SaveAs("Efficiency_vs_TrackParams2D.pdf");
    cEffParams2D->SaveAs("Efficiency_vs_TrackParams2D.png");

    PlotEfficiencyResults(-1, bundleStats, layerStats, cylStats, occupancyMap);
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
        // We compare the true hits (from sim) with the reconstructed hough hits
        // (from rec) Since we are working with vectors of hits, let's calculate
        // intersection sizes

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
                    { "hit_ids", "rec.rec_hits" })
                  .Define("true_total",
                      [](const std::vector<int> &sim_hits) { return (int)sim_hits.size(); },
                      { "hit_ids" })
                  .Define("reco_total",
                      [](const std::vector<int> &rec_hough) { return (int)rec_hough.size(); },
                      { "rec.rec_hits" })
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
    std::string analysis_type = "perf";
    if(argc > 2 && std::string(argv[2]).find("--") != 0)
        analysis_type = argv[2];

    // Load custom style
    SetLBStyle();

    ROOT::EnableImplicitMT();
    Bool_t real = true;
    if(real)
    {
        CHeT::Config::SetActiveCylinders({ 0, 1 });
        SetCHeTOptimalSelfAligment();
    }

    if(analysis_type == "perf")
    {
        bool doEff = false;
        bool doRes = false;
        bool specified = false;

        for(int i = 2; i < argc; ++i)
        {
            std::string arg = argv[i];
            if(arg == "--eff")
            {
                doEff = true;
                specified = true;
            }
            if(arg == "--res")
            {
                doRes = true;
                specified = true;
            }
        }

        if(!specified)
        {
            doEff = true;
            doRes = true;
        }

        TFile *f = TFile::Open(filename.c_str(), "READ");
        bool hasSim = f && f->Get("sim") != nullptr;
        if(f)
        {
            TTree *metaTree = hasSim ? (TTree *)f->Get("sim") : (TTree *)f->Get("rec");
            if(metaTree && metaTree->GetUserInfo())
            {
                std::vector<int> active_cyls;
                bool found_active_meta = false;
                for(int i = 0; i < 6; ++i)
                {
                    auto param = (TParameter<int> *)metaTree->GetUserInfo()->FindObject(
                        Form("ActiveCyl_%d", i));
                    if(param)
                    {
                        found_active_meta = true;
                        if(param->GetVal() > 0)
                            active_cyls.push_back(i);
                    }
                }
                if(found_active_meta)
                {
                    CHeT::Config::SetActiveCylinders(active_cyls);
                    std::cout << ">>> Inherited Active Cylinders from metadata!" << std::endl;
                }
            }
            f->Close();
        }

        TChain chain_rec("rec");
        chain_rec.Add(filename.c_str());

        Long64_t rec_entries = chain_rec.GetEntries();
        std::cout << ">>> Found " << rec_entries << " events processed by Fitter." << std::endl;

        ROOT::RDataFrame *df_ptr = nullptr;
        TChain chain_sim("sim");

        if(hasSim)
        {
            chain_sim.Add(filename.c_str());
            chain_sim.AddFriend(&chain_rec, "rec");
            df_ptr = new ROOT::RDataFrame(chain_sim);
        }
        else
        {
            df_ptr = new ROOT::RDataFrame(chain_rec);
        }

        // Use Filter to bound the processing to Fitter's output, avoiding Range and
        // IMT crash
        auto df = df_ptr->Filter(
            "EventID < " + std::to_string(rec_entries), "Events processed by Fitter");

        bool isMichel = chain_rec.GetBranch("rec_R") != nullptr;

        auto report = df.Report();
        report->Print();

        if(hasSim && doRes)
        {
            DoResolutions(df, isMichel);
        }
        if(doEff)
        {
            PlotEfficiency(df, hasSim);
        }

        delete df_ptr;
    }
    else if(analysis_type == "raw")
    {
        // Plot raw detector correlations (needs raw data runID extraction from
        // filename or argument) PlotRawCorrelation(0, 0, 1000, 0, 1000);
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
            std::cerr << "Please provide files and thresholds! Example: ./Distiller "
                         "file.root roc "
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
