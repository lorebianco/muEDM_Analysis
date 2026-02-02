/**
 * UNIFIED ANALYSIS & SIMULATION FRAMEWORK
 *
 * Contains:
 * 1. Data Structures (Tracks, Hits, Fits)
 * 2. Math & Fitting Core (Minuit2, Hough)
 * 3. Simulation Engine (ToyMC, Digitization)
 * 4. Data Processing Engine (RDataFrame helpers)
 * 5. Plotting & High-Level Analysis Runners for Cosmics
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// ROOT Includes
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Types.h>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>

// Custom Includes
#include "CHeT/CHeTGlobalSettings.hh"
#include "CHeT/CHeTReader.hh"
#include "CHeT/CHeTVisualizer.hh"

using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;
using namespace CHeT;

// Alias per comodità RDataFrame
using RVecUI = RVec<unsigned int>;
using RVecUS = RVec<unsigned short>;
using RVecUC = RVec<unsigned char>;

// Debug
Bool_t DEBUG_TOY = false;

// ==============================================================================
// 1. DATA STRUCTURES
// ==============================================================================

// Struttura per il MC: Verità Monte Carlo
struct CosmicTrack
{
    Double_t x0, y0, z0, ux, uy, uz;
    bool active;
};

// Struttura per il risultato del Fit (Parametrizzazione: x = x0 + sx*y, z = z0
// + sz*y)
struct RecoTrack
{
    double x0; // Coordinata x a y=0
    double z0; // Coordinata z a y=0
    double sx; // Slope x (dx/dy)
    double sz; // Slope z (dz/dy)
    double chi2; // Chi2 ridotto
    bool converged; // Status del fit
};

// Struttura output completa del fit: parametri + punti grafici
struct FitOutput
{
    RecoTrack track;
    std::vector<Vis::VisPoint3D> fittedPoints;
};

// Struttura dati passata al funtore di minimizzazione
struct FitData
{
    vector<Config::FiberProp> props;
};

// Struttura per Hough Transform
struct HoughResult
{
    Double_t rho;
    Double_t theta;
    Double_t weight;
};

// Strutture per calcolo Efficienza
struct EffStats
{
    long passed = 0;
    long total = 0;
};
using EffMap = map<int, EffStats>;

// ==============================================================================
// 2. HELPER FUNCTIONS (Style & Summary)
// ==============================================================================

void SetGlobalStyle()
{
    static bool isStyleSet = false;
    if(isStyleSet)
        return;
    gStyle->SetPalette(60);
    TColor::InvertPalette();
    isStyleSet = true;
}

void PrintSummary(long long total, long long passed, Double_t t0, Double_t t1,
    UInt_t tot0, UInt_t tot1)
{
    double efficiency = (total > 0) ? (100.0 * passed / total) : 0;
    cout << "\n" << string(60, '=') << endl;
    cout << "  EVENT SELECTION SUMMARY" << endl;
    cout << string(60, '-') << endl;
    printf("  ToA Range       : [%.1f, %.1f] ns\n", t0, t1);
    printf("  ToT Range       : [%u, %u] LSB\n", tot0, tot1);
    cout << string(60, '-') << endl;
    printf("  Total Events    : %lld\n", total);
    printf("  Passed Filter   : %lld\n", passed);
    printf("  Efficiency      : %.2f %%\n", efficiency);
    cout << string(60, '=') << "\n" << endl;
}

pair<double, double> GetCPErrors(long passed, long total)
{
    if(total == 0)
        return { 0.0, 0.0 };
    double eff = (double)passed / total;
    double cl = 0.6827; // 1 sigma
    double lower = TEfficiency::ClopperPearson(total, passed, cl, false);
    double upper = TEfficiency::ClopperPearson(total, passed, cl, true);
    return { eff - lower, upper - eff };
}

// --- DEBUG ---

void DebugEfficiencyCalculation(const RecoTrack &tr, const vector<int> &hit_ids)
{
    printf("\n");
    printf("************************************************************\n");
    printf("               DEBUG EFFICIENCY LOGIC                       \n");
    printf("************************************************************\n");
    printf("Track Params: x0=%.2f, z0=%.2f, sx=%.4f, sz=%.4f\n", tr.x0, tr.z0,
        tr.sx, tr.sz);

    const double L_SAFE = Config::L_HALF - 5.0;
    const double TOL_BUNDLE = 1.5;

    auto cylinders = Config::GetCylinders();
    int global_bundle_offset = 0;
    int cyl_idx = 0;

    for(const auto &cyl : cylinders)
    {
        printf("\n--- Cylinder %d ---\n", cyl_idx);
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };

        for(int lay_local_idx = 0; lay_local_idx < 2; ++lay_local_idx)
        {
            const auto &L = *layers[lay_local_idx];
            double R = L.radius;
            int lay_global_idx = cyl_idx * 2 + lay_local_idx;

            printf("  > Layer %s [Global Index: %d] (R=%.1f mm, Global "
                   "Offset=%d)\n",
                (lay_local_idx == 0 ? "Inner" : "Outer"), lay_global_idx, R,
                global_bundle_offset);

            // Intersezione Retta-Cerchio
            double A = tr.sx * tr.sx + 1.0;
            double B = 2.0 * tr.x0 * tr.sx;
            double C = tr.x0 * tr.x0 - R * R;
            double delta = B * B - 4 * A * C;

            if(delta < 0)
            {
                printf("    [MISS] No geometric intersection (delta < 0)\n");
            }
            else
            {
                double y_sol[2] = { (-B + sqrt(delta)) / (2 * A),
                    (-B - sqrt(delta)) / (2 * A) };

                for(int sol_i = 0; sol_i < 2; ++sol_i)
                {
                    double y = y_sol[sol_i];
                    double z = tr.z0 + tr.sz * y;
                    double x = tr.x0 + tr.sx * y;

                    printf(
                        "    * Intersection %d at (x=%.1f, y=%.1f, z=%.1f): ",
                        sol_i, x, y, z);

                    if(abs(z) > L_SAFE)
                    {
                        printf("[SKIP] Z out of bounds (safe limit %.1f)\n",
                            L_SAFE);
                        continue;
                    }
                    printf("[OK] Inside volume.\n");

                    // Logica Bundle
                    double phi_hit = atan2(y, x);
                    int best_bundle_id = -1;
                    double min_dist_norm = 1e9;
                    double bundle_width_rad = 2.0 * M_PI / L.nBundles;

                    // Cerchiamo il bundle più vicino
                    for(int b = 0; b < L.nBundles; ++b)
                    {
                        int gid = global_bundle_offset + b;
                        auto prop = Config::GetFiberProp(gid);
                        double alpha
                            = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                        double phi_fib = prop.phi0 + prop.dir * alpha * M_PI;

                        double dphi = abs(Config::wrap0_2pi(phi_hit)
                            - Config::wrap0_2pi(phi_fib));
                        if(dphi > M_PI)
                            dphi = 2.0 * M_PI - dphi;

                        double dist_norm = dphi / bundle_width_rad;
                        if(dist_norm < min_dist_norm)
                        {
                            min_dist_norm = dist_norm;
                            best_bundle_id = gid;
                        }
                    }

                    printf("      -> Best Bundle: ID %d (Dist Norm: %.3f)\n",
                        best_bundle_id, min_dist_norm);

                    if(best_bundle_id >= 0 && min_dist_norm < TOL_BUNDLE)
                    {
                        printf("      -> [CHECKING HITS] Searching for ID %d "
                               "or neighbours...\n",
                            best_bundle_id);

                        bool hit_bundle = false;
                        int found_id = -1;

                        for(int h : hit_ids)
                        {
                            int diff = abs(h - best_bundle_id);
                            bool is_neighbour
                                = (diff <= 1) || (diff == L.nBundles - 1);

                            // Controllo extra: deve essere dello stesso layer
                            if(is_neighbour
                                && Config::GetFiberProp(h).layerId
                                    == Config::GetFiberProp(best_bundle_id)
                                           .layerId)
                            {
                                hit_bundle = true;
                                found_id = h;
                                break;
                            }
                        }

                        if(hit_bundle)
                        {
                            printf("         [SUCCESS] FOUND MATCH! Hit ID: %d "
                                   "(Delta: %d)\n",
                                found_id, found_id - best_bundle_id);
                        }
                        else
                        {
                            printf("         \033[31m[FAIL]\033[0m EXPECTED "
                                   "HIT NOT FOUND.\n");
                            printf("         (Hits in this event: ");
                            for(int h : hit_ids)
                                printf("%d ", h);
                            printf(")\n");
                        }
                    }
                    else
                    {
                        printf("      -> [SKIP] Too far from any bundle center "
                               "(Gap?).\n");
                    }
                }
            }
            global_bundle_offset += L.nBundles;
        }
        cyl_idx++;
    }
    printf("************************************************************\n\n");
}

// ==============================================================================
// 3. MATH & FITTING ENGINE (Minuit2 & Hough)
// ==============================================================================

// Funzione Chi2 / LogLikelihood
double Track3DNeg2LogL(const double *par, const FitData &data, bool usePrior)
{
    const double x0 = par[0];
    const double sx = par[1];
    const double z0 = par[2];
    const double sz = par[3];
    const double sigma2 = 0.3; // Risoluzione fibra (0.2mm / sqrt(12))^2 ~ 0.3
    const double sigmaL2 = 1e-6; // Penalty forte per hit fuori lunghezza fisica

    double n2ll = 0.0;
    if(usePrior)
        n2ll = 2.0 * log(1.0 + sx * sx + sz * sz);

    for(size_t i = 0; i < data.props.size(); ++i)
    {
        // Parametro minimizzato: dove la particella ha colpito lungo la fibra
        // (coordinata locale z)
        const double zi = par[4 + i];
        const auto &p = data.props[i];

        // Posizione geometrica sulla fibra data z_locale (spirale)
        const double alpha = (zi + Config::L_HALF) / (2.0 * Config::L_HALF);
        const double phi_fib = p.phi0 + p.dir * alpha * M_PI;

        const double x_f = p.r * cos(phi_fib);
        const double y_f = p.r * sin(phi_fib);

        // Posizione attesa dalla traccia (a quel y_f)
        const double x_trk = x0 + sx * y_f;
        const double z_trk = z0 + sz * y_f;

        // Residui (distanza nel piano XZ locale tangente + distanza lungo la
        // fibra) Nota: approssimazione, trattiamo residui su coordinate
        // cartesiane
        n2ll += ((x_f - x_trk) * (x_f - x_trk) + (zi - z_trk) * (zi - z_trk))
            / sigma2;

        // Penalty se il punto fittato esce dalla fibra
        if(abs(zi) > Config::L_HALF)
            n2ll += (pow(abs(zi) - Config::L_HALF, 2)) / sigmaL2;
    }
    return n2ll;
}

// Driver del Fit
FitOutput Do3DFit(const vector<int> &hit_ids, bool usePrior = false)
{
    FitData data;
    for(int id : hit_ids)
        data.props.push_back(Config::GetFiberProp(id));

    int n_hits = data.props.size();
    if(n_hits < 3)
        return { { 0, 0, 0, 0, -1.0, false }, {} };

    ROOT::Math::Minimizer *min
        = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetPrintLevel(0); // 0 = Silent

    auto chi2_func = [&](const double *par)
    { return Track3DNeg2LogL(par, data, usePrior); };
    ROOT::Math::Functor f(chi2_func, 4 + n_hits);
    min->SetFunction(f);

    // Variabili traccia: x0, sx, z0, sz
    min->SetVariable(0, "x0", 0.0, 0.1);
    min->SetVariable(1, "sx", 0.0, 0.01);
    min->SetVariable(2, "z0", 0.0, 0.1);
    min->SetVariable(3, "sz", 0.0, 0.01);

    // Variabili hit: z_locale per ogni fibra
    for(int i = 0; i < n_hits; ++i)
        min->SetVariable(
            4 + i, Form("z_hit_%d", i), 0, 5.0); // Step iniziale 5mm

    bool conv = min->Minimize();
    const double *res = min->X();
    double chi2
        = (n_hits * 2 > 4) ? min->MinValue() / (double)(n_hits * 2 - 4) : 0.0;

    FitOutput output;
    output.track = { res[0], res[2], res[1], res[3], chi2, conv };

    if(conv)
    {
        for(int i = 0; i < n_hits; ++i)
        {
            double z_fit = res[4 + i];
            const auto &p = data.props[i];

            double alpha = (z_fit + Config::L_HALF) / (2.0 * Config::L_HALF);
            double phi = p.phi0 + p.dir * alpha * M_PI;
            double x_fit = p.r * cos(phi);
            double y_fit = p.r * sin(phi);

            // Visualizzazione: punti neri cerchiati
            output.fittedPoints.emplace_back(
                x_fit, y_fit, z_fit, kBlack, 24, 1.0);
        }
    }
    return output;
}

// Trasformata di Hough (XY polare)
vector<HoughResult> DoHoughTransform(
    const vector<Int_t> &hit_ids, int nCandidates = 2)
{
    Int_t nBinsTheta = 80;
    Double_t thetaMin = 0.0, thetaMax = TMath::Pi();
    Int_t nBinsRho = 40;
    Double_t rhoMin = -30.0, rhoMax = 30.0;

    gROOT->Delete("h_xy_pol;1");
    TH2D *h_pol = new TH2D("h_xy_pol", "Hough Polare XY;#theta [rad];#rho [mm]",
        nBinsTheta, thetaMin, thetaMax, nBinsRho, rhoMin, rhoMax);

    auto inters = Config::FindIntersections(hit_ids);
    for(auto &p : inters)
    {
        for(int it = 1; it <= nBinsTheta; ++it)
        {
            double theta = h_pol->GetXaxis()->GetBinCenter(it);
            double rho = p.x * cos(theta) + p.y * sin(theta);
            h_pol->Fill(theta, rho);
        }
    }

    // Ricerca Picchi
    TH2D *h_work = (TH2D *)h_pol->Clone("h_work");
    vector<HoughResult> candidates;
    Double_t topWeight = 0.0;

    for(int n = 0; n < nCandidates; ++n)
    {
        Int_t bMax = h_work->GetMaximumBin();
        Int_t bx, by, bz;
        h_work->GetBinXYZ(bMax, bx, by, bz);
        double maxWeight = h_work->GetBinContent(bMax);

        if(maxWeight < 1.5)
            break;
        if(n == 0)
            topWeight = maxWeight;
        else if(maxWeight < (topWeight - 0.01))
            break;

        double sw = 0, st = 0, sr = 0;
        int win = 2;
        for(int ix = bx - win; ix <= bx + win; ++ix)
        {
            for(int iy = by - win; iy <= by + win; ++iy)
            {
                if(ix < 1 || ix > nBinsTheta || iy < 1 || iy > nBinsRho)
                    continue;
                double w = h_work->GetBinContent(ix, iy);
                sw += w;
                st += w * h_work->GetXaxis()->GetBinCenter(ix);
                sr += w * h_work->GetYaxis()->GetBinCenter(iy);
                h_work->SetBinContent(ix, iy, 0);
            }
        }
        if(sw > 0)
            candidates.push_back({ sr / sw, st / sw, maxWeight });
    }

    // Scommentare per debug grafico
    /*
    TCanvas *c = (TCanvas *)gROOT->FindObject("c_hough");
    if(!c)
        c = new TCanvas("c_hough", "Hough Polare Space", 600, 400);
    c->cd();
    h_pol->SetStats(0);
    h_pol->Draw("COLZ");
    */

    delete h_work;
    return candidates;
}

// ==============================================================================
// 4. SIMULATION ENGINE (Generator & Digitizer)
// ==============================================================================

CosmicTrack GenerateCosmic()
{
    static TRandom3 rnd(0);
    CosmicTrack track;
    track.active = true;

    Double_t halfX_up = 90., halfZ_up = 90., yUp = 386;
    Double_t halfX_down = 90., halfZ_down = 90., yDown = -476;
    Bool_t accepted = false;

    while(!accepted)
    {
        track.y0 = yUp;
        track.x0 = rnd.Uniform(-halfX_up, halfX_up);
        track.z0 = rnd.Uniform(-halfZ_up, halfZ_up);

        Double_t cosTheta, cosThetaSq_test, phi;
        do
        {
            cosTheta = -rnd.Uniform(0., 1.); // Verso il basso
            cosThetaSq_test = rnd.Uniform(0., 1.);
            phi = rnd.Uniform(0, 2. * M_PI);
        } while(cosThetaSq_test > (cosTheta * cosTheta));

        Double_t sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        track.uz = sinTheta * cos(phi);
        track.ux = sinTheta * sin(phi);
        track.uy = cosTheta;

        Double_t t = (yDown - track.y0) / track.uy;
        Double_t x_proj = track.x0 + track.ux * t;
        Double_t z_proj = track.z0 + track.uz * t;

        if(abs(x_proj) <= halfX_down && abs(z_proj) <= halfZ_down)
            accepted = true;
    }
    return track;
}

vector<Int_t> FindHitBundles(CosmicTrack tr, Double_t efficiency = 1.)
{
    vector<Int_t> hits;
    auto cylinders = Config::GetCylinders();
    int current_global_offset = 0;

    for(const auto &cyl : cylinders)
    {
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int i_lay = 0; i_lay < 2; ++i_lay)
        {
            const auto &lay = *layers[i_lay];
            Double_t R = lay.radius;

            // Intersezione Retta XY - Cerchio
            Double_t A = tr.ux * tr.ux + tr.uy * tr.uy;
            Double_t B = 2.0 * (tr.x0 * tr.ux + tr.y0 * tr.uy);
            Double_t C = tr.x0 * tr.x0 + tr.y0 * tr.y0 - R * R;
            Double_t delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                Double_t t_sol[2] = { (-B + sqrt(delta)) / (2 * A),
                    (-B - sqrt(delta)) / (2 * A) };
                for(int i = 0; i < 2; ++i)
                {
                    Double_t t = t_sol[i];
                    Double_t zi = tr.z0 + tr.uz * t;
                    // Check Z geometrico
                    if(abs(zi) <= Config::L_HALF)
                    {
                        Double_t xi = tr.x0 + tr.ux * t;
                        Double_t yi = tr.y0 + tr.uy * t;
                        Double_t phi_track = atan2(yi, xi);

                        // Check Bundles
                        for(int b = 0; b < lay.nBundles; ++b)
                        {
                            int b_id = current_global_offset + b;
                            Config::FiberProp p = Config::GetFiberProp(b_id);

                            Double_t alpha = (zi + Config::L_HALF)
                                / (2.0 * Config::L_HALF);
                            Double_t phi_f = p.phi0 + p.dir * alpha * M_PI;

                            Double_t dphi = abs(Config::wrap0_2pi(phi_track)
                                - Config::wrap0_2pi(phi_f));
                            if(dphi > M_PI)
                                dphi = 2.0 * M_PI - dphi;

                            // Condition di hit
                            if(dphi < (M_PI / lay.nBundles)
                                && gRandom->Rndm() <= efficiency)
                                hits.push_back(b_id);
                        }
                    }
                }
            }
            current_global_offset += lay.nBundles;
        }
    }
    sort(hits.begin(), hits.end());
    hits.erase(unique(hits.begin(), hits.end()), hits.end());
    return hits;
}

bool IsGeometricIntersection(const CosmicTrack &tr)
{
    auto cylinders = Config::GetCylinders();
    Double_t A = tr.ux * tr.ux + tr.uy * tr.uy;
    Double_t B = 2.0 * (tr.x0 * tr.ux + tr.y0 * tr.uy);

    for(const auto &cyl : cylinders)
    {
        double radii[2] = { cyl.inner.radius, cyl.outer.radius };
        for(double R : radii)
        {
            Double_t C = tr.x0 * tr.x0 + tr.y0 * tr.y0 - R * R;
            Double_t delta = B * B - 4 * A * C;
            if(delta >= 0)
            {
                Double_t t_sol[2] = { (-B + sqrt(delta)) / (2 * A),
                    (-B - sqrt(delta)) / (2 * A) };
                for(int i = 0; i < 2; ++i)
                    if(abs(tr.z0 + tr.uz * t_sol[i]) <= Config::L_HALF)
                        return true;
            }
        }
    }
    return false;
}

// ==============================================================================
// 5. DATA PROCESSING HELPERS (RDataFrame)
// ==============================================================================

void AccumulateEfficiency(const RecoTrack &tr, const vector<int> &hit_ids,
    EffMap &bundleMap, EffMap &layerMap, EffMap &cylMap)
{
    const double L_SAFE = Config::L_HALF - 5.0;
    const double TOL_BUNDLE = 1.5;
    const int N_NEIGHBORS = 1;

    auto cylinders = Config::GetCylinders();
    int global_bundle_offset = 0;
    int cyl_idx = 0;

    for(const auto &cyl : cylinders)
    {
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };

        for(int lay_local_idx = 0; lay_local_idx < 2; ++lay_local_idx)
        {
            const auto &L = *layers[lay_local_idx];
            double R = L.radius;
            int lay_global_idx = cyl_idx * 2 + lay_local_idx;

            // Salviamo l'inizio di questo layer per calcoli se necessario,
            // ma useremo principalmente GetFiberProp per sicurezza.
            int current_layer_start = global_bundle_offset;

            // Intersezione Retta-Cerchio
            double A = tr.sx * tr.sx + 1.0;
            double B = 2.0 * tr.x0 * tr.sx;
            double C = tr.x0 * tr.x0 - R * R;
            double delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                double y_sol[2] = { (-B + sqrt(delta)) / (2 * A),
                    (-B - sqrt(delta)) / (2 * A) };
                for(double y : y_sol)
                {
                    double z = tr.z0 + tr.sz * y;
                    if(abs(z) <= L_SAFE)
                    {
                        double x = tr.x0 + tr.sx * y;
                        double phi_hit = atan2(y, x);

                        int best_bundle_id = -1; // ID Globale
                        double min_dist_norm = 1e9;
                        double bundle_width_rad = 2.0 * M_PI / L.nBundles;

                        // --- TROVA IL BUNDLE TEORICO MIGLIORE ---
                        for(int b = 0; b < L.nBundles; ++b)
                        {
                            int gid = current_layer_start + b;

                            // Nota: GetFiberProp qui è sicura perché gid è
                            // sequenziale
                            auto prop = Config::GetFiberProp(gid);

                            double alpha
                                = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                            double phi_fib
                                = prop.phi0 + prop.dir * alpha * M_PI;

                            double dphi = abs(Config::wrap0_2pi(phi_hit)
                                - Config::wrap0_2pi(phi_fib));
                            if(dphi > M_PI)
                                dphi = 2.0 * M_PI - dphi;

                            double dist_norm = dphi / bundle_width_rad;
                            if(dist_norm < min_dist_norm)
                            {
                                min_dist_norm = dist_norm;
                                best_bundle_id = gid;
                            }
                        }

                        if(best_bundle_id >= 0 && min_dist_norm < TOL_BUNDLE)
                        {
                            int hit_bundle_id = -1;

                            // 1. Priorita' al match esatto
                            for(int h : hit_ids)
                            {
                                if(h == best_bundle_id)
                                {
                                    hit_bundle_id = best_bundle_id;
                                    break;
                                }
                            }

                            // 2. Cerca vicini (Solo se match esatto fallisce)
                            if(hit_bundle_id == -1)
                            {
                                // Recuperiamo le proprietà del 'best' per
                                // sapere chi siamo
                                auto best_prop
                                    = Config::GetFiberProp(best_bundle_id);

                                for(int h : hit_ids)
                                {
                                    auto hit_prop = Config::GetFiberProp(h);

                                    // --- FILTRO CRUCIALE ---
                                    // 1. Deve essere lo stesso CILINDRO
                                    if(hit_prop.cylinderId
                                        != best_prop.cylinderId)
                                        continue;

                                    // 2. Deve essere lo stesso LAYER
                                    // (Inner/Outer)
                                    if(hit_prop.layerId != best_prop.layerId)
                                        continue;

                                    // --- CALCOLO DISTANZA ---
                                    // Ora che siamo sicuri di essere sulla
                                    // stessa superficie logica, la differenza
                                    // tra ID globali corrisponde alla distanza
                                    // in indici.
                                    int diff = abs(h - best_bundle_id);

                                    // Logica circolare
                                    bool is_neighbour = (diff <= N_NEIGHBORS)
                                        || (diff >= L.nBundles - N_NEIGHBORS);

                                    if(is_neighbour)
                                    {
                                        hit_bundle_id = h;
                                        break;
                                    }
                                }
                            }

                            // --- AGGIORNAMENTO STATISTICHE ---
                            // Se abbiamo trovato un hit (esatto o vicino),
                            // usiamo quello. Altrimenti usiamo il best match
                            // teorico (che sarà contato come MISS).
                            int id_to_fill = (hit_bundle_id != -1)
                                ? hit_bundle_id
                                : best_bundle_id;

                            bundleMap[id_to_fill].total++;
                            layerMap[lay_global_idx].total++;
                            cylMap[cyl_idx].total++;

                            if(hit_bundle_id != -1)
                            {
                                bundleMap[id_to_fill].passed++;
                                layerMap[lay_global_idx].passed++;
                                cylMap[cyl_idx].passed++;
                            }
                        }
                    }
                }
            }
            // Avanziamo l'offset globale
            global_bundle_offset += L.nBundles;
        }
        cyl_idx++;
    }
}

// ==============================================================================
// 6. PLOTTING FUNCTIONS (Data)
// ==============================================================================

void DrawEvent(Int_t eventID, Int_t runID, Double_t toaMin, Double_t toaMax,
    UInt_t totMin, UInt_t totMax)
{
    if(ROOT::IsImplicitMTEnabled())
        ROOT::DisableImplicitMT();

    SetGlobalStyle();
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    if(gSystem->AccessPathName(filename.c_str()))
    {
        cerr << "[Error] File not found: " << filename << endl;
        return;
    }

    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    if(eventID == -1)
    {
        gRandom->SetSeed(0);
        long long nEntries = reader.GetRaw().Count().GetValue();
        eventID = gRandom->Integer(nEntries);
    }

    reader.SetSingleEvent(eventID);
    auto df = reader.GetEstimators();

    // The reader provides "All_Bundle", which contains the global IDs of
    // selected hits
    auto hitsResult = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");

    if(hitsResult->empty())
    {
        cerr << "[Warning] Event " << eventID << " seems empty or out of range."
             << endl;
        return;
    }

    // Convert RVec<int> to vector<int>
    const auto &rvec = hitsResult->at(0);
    vector<int> hit_ids(rvec.begin(), rvec.end());

    auto intersections = Config::FindIntersections(hit_ids);

    cout << "Displaying Run " << runID << " - Event " << eventID << endl;
    cout << "Selected Hits : " << hit_ids.size() << endl;
    cout << "Intersections : " << intersections.size() << endl;

    // Call Visualizer
    Vis::Draw2D(hit_ids, intersections);
    Vis::Draw3D(hit_ids);
}

void PlotRawCorrelation(
    Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    if(!ROOT::IsImplicitMTEnabled())
        ROOT::EnableImplicitMT(16);

    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    Data::Reader reader(filename);
    // Explicitly set cuts even though we define our own filtered vars below
    // (Reader defines FDxx_corrToA)
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df = reader.GetEstimators();
    auto rawCount = reader.GetRaw().Count();

    auto df_all_hits
        = df.Define("All_ToA_raw",
                [](const RVecD &t0, const RVecD &t1, const RVecD &t2,
                    const RVecD &t3) {
                    return Concatenate(
                        t0, Concatenate(t1, Concatenate(t2, t3)));
                },
                { "FD00_corrToA", "FD01_corrToA", "FD02_corrToA",
                    "FD03_corrToA" })
              .Define("All_ToT_raw",
                  [](const RVecUS &t0, const RVecUS &t1, const RVecUS &t2,
                      const RVecUS &t3) {
                      return RVecD(Concatenate(
                          t0, Concatenate(t1, Concatenate(t2, t3))));
                  },
                  { "FD00_ToT", "FD01_ToT", "FD02_ToT", "FD03_ToT" })
              .Define("All_ToA",
                  [=](const RVecD &toa, const RVecD &tot)
                  {
                      return toa[(toa > toaMin) && (toa < toaMax)
                          && (tot > totMin) && (tot < totMax)];
                  },
                  { "All_ToA_raw", "All_ToT_raw" })
              .Define("All_ToT",
                  [=](const RVecD &toa, const RVecD &tot)
                  {
                      return tot[(toa > toaMin) && (toa < toaMax)
                          && (tot > totMin) && (tot < totMax)];
                  },
                  { "All_ToA_raw", "All_ToT_raw" });

    auto nPassed = df_all_hits.Filter("All_ToA.size() > 0").Count();
    auto h_raw_corr = df_all_hits.Histo2D(
        { "h_raw_corr",
            Form("Raw Correlation Run %d;ToA [ns];ToT [LSB]", runID),
            static_cast<int>((toaMax - toaMin) / 4), toaMin, toaMax,
            static_cast<int>((totMax - totMin) / 4), (double)totMin,
            (double)totMax },
        "All_ToA", "All_ToT");

    TCanvas *cRaw = new TCanvas("cRaw", "Correlation Analysis", 800, 600);
    cRaw->SetLogz();
    h_raw_corr->DrawCopy("COLZ");
    PrintSummary(rawCount.GetValue(), nPassed.GetValue(), toaMin, toaMax,
        totMin, totMax);
}

void PlotMultiplicity(
    Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();

    // --- Setup Dati ---
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df_full = reader.GetEstimators();

    long long nRaw = reader.GetRaw().Count().GetValue();
    auto df = df_full.Filter("nHits_Total > 0");
    auto nPassed = df.Count();

    // ------------------------------------------
    // 1. Calcolo Dinamico dei Limiti Bundle ID
    // ------------------------------------------
    auto cylinders = Config::GetCylinders();
    double binEdges[5] = { 0, 0, 0, 0, 0 };
    int nBins[4] = { 0, 0, 0, 0 };
    double currentOffset = 0;

    if(cylinders.size() > 0)
    {
        nBins[0] = cylinders[0].inner.nBundles;
        binEdges[0] = currentOffset;
        currentOffset += nBins[0];
        binEdges[1] = currentOffset;
        nBins[1] = cylinders[0].outer.nBundles;
        currentOffset += nBins[1];
        binEdges[2] = currentOffset;
    }
    if(cylinders.size() > 1)
    {
        nBins[2] = cylinders[1].inner.nBundles;
        currentOffset += nBins[2];
        binEdges[3] = currentOffset;
        nBins[3] = cylinders[1].outer.nBundles;
        currentOffset += nBins[3];
        binEdges[4] = currentOffset;
    }

    // ------------------------------------------
    // 2. Definizione Istogrammi
    // ------------------------------------------
    auto h_total = df.Histo1D(
        { "h_total", "Total Multiplicity;n Hits;Events", 201, -0.5, 200.5 },
        "nHits_Total");
    auto h_c0 = df.Histo1D(
        { "h_c0", "Cylinders Multiplicity;n Hits;Events", 121, -0.5, 120.5 },
        "nHits_Cyl0");
    auto h_c1 = df.Histo1D({ "h_c1", "", 121, -0.5, 120.5 }, "nHits_Cyl1");
    auto h_l1 = df.Histo1D(
        { "h_l1", "Layers Multiplicity;n Hits;Events", 81, -0.5, 80.5 },
        "nHits_Cyl0_Inner");
    auto h_l2 = df.Histo1D({ "h_l2", "", 81, -0.5, 80.5 }, "nHits_Cyl0_Outer");
    auto h_l3 = df.Histo1D({ "h_l3", "", 81, -0.5, 80.5 }, "nHits_Cyl1_Inner");
    auto h_l4 = df.Histo1D({ "h_l4", "", 81, -0.5, 80.5 }, "nHits_Cyl1_Outer");

    // Occupancy
    auto df_occ = df.Define("Bundles_Cyl0_Lay0",
                        "All_Bundle[All_Cyl == 0 && All_Lay == 0]")
                      .Define("Bundles_Cyl0_Lay1",
                          "All_Bundle[All_Cyl == 0 && All_Lay == 1]")
                      .Define("Bundles_Cyl1_Lay0",
                          "All_Bundle[All_Cyl == 1 && All_Lay == 0]")
                      .Define("Bundles_Cyl1_Lay1",
                          "All_Bundle[All_Cyl == 1 && All_Lay == 1]");

    auto h_occ_c0l0 = df_occ.Histo1D(
        { "h_occ_c0l0", "Occupancy Cyl 0 - Inner;Global Bundle ID;Hits",
            nBins[0], binEdges[0], binEdges[1] },
        "Bundles_Cyl0_Lay0");
    auto h_occ_c0l1 = df_occ.Histo1D(
        { "h_occ_c0l1", "Occupancy Cyl 0 - Outer;Global Bundle ID;Hits",
            nBins[1], binEdges[1], binEdges[2] },
        "Bundles_Cyl0_Lay1");
    auto h_occ_c1l0 = df_occ.Histo1D(
        { "h_occ_c1l0", "Occupancy Cyl 1 - Inner;Global Bundle ID;Hits",
            nBins[2], binEdges[2], binEdges[3] },
        "Bundles_Cyl1_Lay0");
    auto h_occ_c1l1 = df_occ.Histo1D(
        { "h_occ_c1l1", "Occupancy Cyl 1 - Outer;Global Bundle ID;Hits",
            nBins[3], binEdges[3], binEdges[4] },
        "Bundles_Cyl1_Lay1");

    // ------------------------------------------
    // 3. Disegno
    // ------------------------------------------

    // Canvas 1: Multiplicity (Stile Classico)
    TCanvas *cMult
        = new TCanvas("cMult", Form("Multiplicity - Run %d", runID), 1200, 800);
    cMult->Divide(1, 2);

    cMult->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    TH1D *cp_total = (TH1D *)h_total->DrawCopy();
    cp_total->SetFillColorAlpha(kGray, 0.5);
    cp_total->SetLineColor(kBlack);
    cp_total->SetLineWidth(2);

    cMult->cd(2)->Divide(2, 1);
    cMult->cd(2)->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    TH1D *cp_c0 = (TH1D *)h_c0->DrawCopy();
    cp_c0->SetFillColorAlpha(kRed, 0.35);
    cp_c0->SetLineColor(kRed);
    cp_c0->SetLineWidth(2);
    TH1D *cp_c1 = (TH1D *)h_c1->DrawCopy("SAME");
    cp_c1->SetFillColorAlpha(kBlue, 0.35);
    cp_c1->SetLineColor(kBlue);
    cp_c1->SetLineWidth(2);

    auto legCyl = new TLegend(0.58, 0.45, 0.87, 0.54);
    legCyl->AddEntry(cp_c0, "Cylinder 0", "lf");
    legCyl->AddEntry(cp_c1, "Cylinder 1", "lf");
    legCyl->SetBorderSize(0);
    legCyl->Draw();

    cMult->cd(2)->cd(2);
    gPad->SetLogy();
    gPad->SetGridy();
    TH1D *cp_l1 = (TH1D *)h_l1->DrawCopy();
    cp_l1->SetFillColorAlpha(kRed, 0.3);
    cp_l1->SetLineColor(kRed);
    TH1D *cp_l2 = (TH1D *)h_l2->DrawCopy("SAME");
    cp_l2->SetFillColorAlpha(kOrange + 1, 0.3);
    cp_l2->SetLineColor(kOrange + 1);
    TH1D *cp_l3 = (TH1D *)h_l3->DrawCopy("SAME");
    cp_l3->SetFillColorAlpha(kBlue, 0.3);
    cp_l3->SetLineColor(kBlue);
    TH1D *cp_l4 = (TH1D *)h_l4->DrawCopy("SAME");
    cp_l4->SetFillColorAlpha(kCyan + 1, 0.3);
    cp_l4->SetLineColor(kCyan + 1);

    auto legLay = new TLegend(0.58, 0.4, 0.87, 0.58);
    legLay->AddEntry(cp_l1, "Cyl 0 - Inner", "lf");
    legLay->AddEntry(cp_l2, "Cyl 0 - Outer", "lf");
    legLay->AddEntry(cp_l3, "Cyl 1 - Inner", "lf");
    legLay->AddEntry(cp_l4, "Cyl 1 - Outer", "lf");
    legLay->SetBorderSize(0);
    legLay->Draw();

    // --- Canvas 2: Occupancy (Stile a Barre Separate) ---
    TCanvas *cOcc = new TCanvas(
        "cOcc", Form("Bundle Occupancy - Run %d", runID), 1200, 800);
    cOcc->Divide(2, 2);

    auto DrawOcc = [](TH1D *h, int col)
    {
        h->SetFillColor(col);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        h->DrawCopy("HIST BAR1");
    };

    cOcc->cd(1);
    gPad->SetGrid();
    DrawOcc((TH1D *)h_occ_c0l0.GetPtr(), kRed);
    cOcc->cd(2);
    gPad->SetGrid();
    DrawOcc((TH1D *)h_occ_c0l1.GetPtr(), kOrange + 1);
    cOcc->cd(3);
    gPad->SetGrid();
    DrawOcc((TH1D *)h_occ_c1l0.GetPtr(), kBlue);
    cOcc->cd(4);
    gPad->SetGrid();
    DrawOcc((TH1D *)h_occ_c1l1.GetPtr(), kCyan + 1);

    PrintSummary(nRaw, nPassed.GetValue(), toaMin, toaMax, totMin, totMax);
}

void PlotEstimators(
    Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    long long nRaw = reader.GetRaw().Count().GetValue();
    auto df = reader.GetEstimators().Filter("SumToT > 0");
    auto nPassed = df.Count();

    auto h_toa = df.Histo1D(
        { "h_toa", "First ToA;First ToA [ns];Events", 150, toaMin, toaMax },
        "FirstToA");
    auto h_tot = df.Histo1D(
        { "h_tot", "Sum of ToT;Sum ToT [LSB];Events", 100, 0, 5000 }, "SumToT");
    auto h_corr = df.Histo2D(
        { "h_corr",
            "Correlation First ToA vs Sum ToT;First ToA [ns];Sum ToT [LSB]",
            100, toaMin, toaMax, 100, 0, 5000 },
        "FirstToA", "SumToT");

    TCanvas *c
        = new TCanvas("c", Form("Estimators - Run %d", runID), 1200, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_toa->SetLineColor(kBlue + 1);
    h_toa->DrawCopy();
    c->cd(2);
    h_tot->SetLineColor(kRed + 1);
    h_tot->DrawCopy();
    c->cd(3);
    h_corr->DrawCopy("COLZ");

    PrintSummary(nRaw, nPassed.GetValue(), toaMin, toaMax, totMin, totMax);
}

void PlotMultiplicityCorrelations(
    Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    long long nRaw = reader.GetRaw().Count().GetValue();
    auto df = reader.GetEstimators().Filter("nHits_Total > 0");
    auto nPassed = df.Count();

    auto h_hits_toa = df.Histo2D(
        { "h_hits_toa", "Multiplicity vs First ToA;n Hits;First ToA [ns]", 101,
            -0.5, 100.5, 150, toaMin, toaMax },
        "nHits_Total", "FirstToA");
    auto h_hits_tot = df.Histo2D(
        { "h_hits_tot", "Multiplicity vs Sum ToT;n Hits;Sum ToT [LSB]", 101,
            -0.5, 100.5, 200, 0, 10000 },
        "nHits_Total", "SumToT");

    TCanvas *cMultCorr = new TCanvas("cMultCorr", "Correlations", 1400, 600);
    cMultCorr->Divide(2, 1);
    cMultCorr->cd(1);
    gPad->SetLogz();
    gPad->SetGrid();
    h_hits_toa->DrawCopy("COLZ");
    cMultCorr->cd(2);
    gPad->SetLogz();
    gPad->SetGrid();
    h_hits_tot->DrawCopy("COLZ");

    PrintSummary(nRaw, nPassed.GetValue(), toaMin, toaMax, totMin, totMax);
}

void PlotEfficiencyResults(int runID, const EffMap &bMap, const EffMap &lMap,
    const EffMap &cMap, const map<int, int> &occMap)
{
    // --- Configurazione Stili ---
    const int col_lay[] = { kRed, kOrange + 1, kBlue, kCyan + 1 };
    const char *lab_lay[]
        = { "Cyl 0 - In", "Cyl 0 - Out", "Cyl 1 - In", "Cyl 1 - Out" };

    const int col_cyl[] = { kRed, kBlue };
    const char *lab_cyl[] = { "Cylinder 0", "Cylinder 1" };

    // ---------------------------------------------------------
    // 0. NEW: Noise / Occupancy Canvas
    // ---------------------------------------------------------
    TCanvas *cNoise = new TCanvas("cNoise", "Noise and Occupancy", 1200, 600);
    cNoise->cd();
    gPad->SetGridy();

    // Calcolo Max ID per assi
    int max_id = 0;
    if(!bMap.empty())
        max_id = max(max_id, bMap.rbegin()->first);
    if(!occMap.empty())
        max_id = max(max_id, occMap.rbegin()->first);

    // Istogrammi per Stack
    TH1D *h_all = new TH1D("h_all",
        "Hit Occupancy (Signal vs Noise);Global Bundle ID;Counts", max_id + 1,
        -0.5, max_id + 0.5);
    TH1D *h_match = new TH1D(
        "h_match", "Matched", max_id + 1, -0.5, max_id + 0.5); // Verde
    TH1D *h_noise = new TH1D(
        "h_noise", "Unassigned", max_id + 1, -0.5, max_id + 0.5); // Rosso

    for(int i = 0; i <= max_id; ++i)
    {
        // Totale hits fisici
        int total = (occMap.count(i)) ? occMap.at(i) : 0;

        // Hits associati a tracce (usiamo 'passed' come proxy per matched)
        int matched = (bMap.count(i)) ? bMap.at(i).passed : 0;

        // Rumore = Totale - Matched (clamp a 0 per sicurezza)
        int noise = total - matched;
        if(noise < 0)
            noise = 0;

        h_all->SetBinContent(i + 1, total); // Serve per la scala
        h_match->SetBinContent(i + 1, matched);
        h_noise->SetBinContent(i + 1, noise);
    }

    // Disegno Frame
    h_all->SetStats(0);
    h_all->SetMaximum(h_all->GetMaximum() * 1.15); // +15% margine
    h_all->SetLineColor(0); // Invisibile
    h_all->Draw("AXIS");

    // Stack
    THStack *hs = new THStack("hs", "");
    h_match->SetFillColorAlpha(kGreen - 3, 0.35);
    h_match->SetLineColor(kGreen - 2);
    hs->Add(h_match);
    h_noise->SetFillColorAlpha(kRed - 3, 0.35);
    h_noise->SetLineColor(kRed - 2);
    hs->Add(h_noise);
    hs->Draw("SAME HIST");

    // Legenda
    TLegend *leg = new TLegend(0.78, 0.835, 0.88, 0.91);
    leg->AddEntry(h_noise, "Unassigned", "f");
    leg->AddEntry(h_match, "Track Signal", "f");
    leg->SetBorderSize(0);
    leg->Draw();

    // Separatori Verticali (Copia della logica esistente)
    auto cylinders = Config::GetCylinders();
    int current_offset = 0;
    double yMin = h_all->GetMinimum();
    double yMax = h_all->GetMaximum();

    for(const auto &cyl : cylinders)
    {
        // Inner
        int nIn = cyl.inner.nBundles;
        double x_line_in = current_offset + nIn - 0.5;
        TLine *lIn = new TLine(x_line_in, yMin, x_line_in, yMax);
        lIn->SetLineStyle(2);
        lIn->SetLineColor(kGray + 1);
        lIn->Draw();

        current_offset += nIn;

        // Outer
        int nOut = cyl.outer.nBundles;
        double x_line_out = current_offset + nOut - 0.5;
        TLine *lOut = new TLine(x_line_out, yMin, x_line_out, yMax);
        lOut->SetLineStyle(1);
        lOut->SetLineColor(kBlack);
        lOut->Draw();

        current_offset += nOut;
    }
    h_all->Draw("AXIS SAME"); // Ridisegna assi sopra

    // ---------------------------------------------------------
    // Canvas Originale Efficiencies
    // ---------------------------------------------------------
    TCanvas *cEff
        = new TCanvas("cEff", Form("Efficiencies Run %d", runID), 900, 1200);
    cEff->Divide(1, 3);

    // ---------------------------------------------------------
    // 1. Cylinder Efficiency
    // ---------------------------------------------------------
    cEff->cd(1);
    gPad->SetGridy();

    TH1D *h_frame_cyl = new TH1D(
        "h_frame_cyl", "Cylinder Efficiency;;Efficiency", 2, -0.5, 1.5);
    h_frame_cyl->SetMinimum(0.5);
    h_frame_cyl->SetMaximum(1.05);
    h_frame_cyl->SetStats(0);
    for(int i = 0; i < 2; ++i)
        h_frame_cyl->GetXaxis()->SetBinLabel(i + 1, lab_cyl[i]);
    h_frame_cyl->GetXaxis()->SetLabelSize(0.07);
    h_frame_cyl->Draw("AXIS");

    TLine *lineC = new TLine(-0.5, 1.0, 1.5, 1.0);
    lineC->SetLineStyle(2);
    lineC->SetLineColor(kGray);
    lineC->Draw();

    TH1D *h_cyls[2];
    TGraphAsymmErrors *g_err_cyl[2];

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
        tG.DrawLatex(
            0.60, 0.75, Form("Global Efficiency: %.2f%%", g_eff * 100));

        printf("\n--- GLOBAL EFFICIENCY ---\n %.2f%% (%ld/%ld)\n", g_eff * 100,
            g_passed, g_total);
    }

    printf("\n--- CYLINDER EFFICIENCY ---\n");

    for(int i = 0; i < 2; ++i)
    {
        h_cyls[i] = new TH1D(Form("h_cyl_%d", i), "", 2, -0.5, 1.5);
        g_err_cyl[i] = new TGraphAsymmErrors();

        auto it = cMap.find(i);
        if(it != cMap.end() && it->second.total > 0)
        {
            auto errs = GetCPErrors(it->second.passed, it->second.total);
            double eff = (double)it->second.passed / it->second.total;

            h_cyls[i]->SetBinContent(i + 1, eff);

            g_err_cyl[i]->SetPoint(0, (double)i, eff);
            g_err_cyl[i]->SetPointError(0, 0.5, 0.5, errs.first, errs.second);

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_cyl[i],
                eff * 100, errs.second * 100, errs.first * 100,
                it->second.passed, it->second.total);
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

    TH1D *h_frame_lay = new TH1D(
        "h_frame_lay", "Sub-Layer Efficiency;;Efficiency", 4, -0.5, 3.5);
    h_frame_lay->SetMinimum(0.5);
    h_frame_lay->SetMaximum(1.05);
    h_frame_lay->SetStats(0);
    for(int i = 0; i < 4; ++i)
        h_frame_lay->GetXaxis()->SetBinLabel(i + 1, lab_lay[i]);
    h_frame_lay->GetXaxis()->SetLabelSize(0.07);
    h_frame_lay->Draw("AXIS");

    TLine *lineL = new TLine(-0.5, 1.0, 3.5, 1.0);
    lineL->SetLineStyle(2);
    lineL->SetLineColor(kGray);
    lineL->Draw();

    TH1D *h_layers[4];
    TGraphAsymmErrors *g_err_lay[4];

    printf("\n--- LAYER EFFICIENCY ---\n");

    for(int i = 0; i < 4; ++i)
    {
        h_layers[i] = new TH1D(Form("h_lay_%d", i), "", 4, -0.5, 3.5);
        g_err_lay[i] = new TGraphAsymmErrors();

        auto it = lMap.find(i);
        if(it != lMap.end() && it->second.total > 0)
        {
            double eff = (double)it->second.passed / it->second.total;
            auto errs = GetCPErrors(it->second.passed, it->second.total);

            h_layers[i]->SetBinContent(i + 1, eff);

            g_err_lay[i]->SetPoint(0, (double)i, eff);
            g_err_lay[i]->SetPointError(0, 0.5, 0.5, errs.first, errs.second);

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_lay[i],
                eff * 100, errs.second * 100, errs.first * 100,
                it->second.passed, it->second.total);
        }

        h_layers[i]->SetFillColorAlpha(col_lay[i], 0.3);
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
        if(cnt.total > 5)
        {
            double eff = (double)cnt.passed / cnt.total;
            auto errs = GetCPErrors(cnt.passed, cnt.total); // {low, high}

            g_bund->SetPoint(pt, (double)id, eff);
            g_bund->SetPointError(pt, 0, 0, errs.first, errs.second);
            pt++;
        }
    }

    TH1D *h_frame_bund = new TH1D("h_frame_bund",
        "Bundle Efficiency;Global Bundle ID;Efficiency", max_id + 1, -0.5,
        max_id + 0.5);
    h_frame_bund->SetMinimum(0.0);
    h_frame_bund->SetMaximum(1.15);
    h_frame_bund->SetStats(0);
    h_frame_bund->Draw("AXIS");

    // Ridisegno separatori anche qui
    current_offset = 0;
    yMin = h_frame_bund->GetMinimum();
    yMax = h_frame_bund->GetMaximum();
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

    printf("------------------------\n");
}

// ==============================================================================
// 7. HIGH-LEVEL ANALYSIS RUNNERS FOR COSMICS
// ==============================================================================

// --- REAL DATA ANALYSIS ---

void FitCosmicEvent(Int_t eventID, Int_t runID, Double_t toaMin,
    Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    if(ROOT::IsImplicitMTEnabled())
        ROOT::DisableImplicitMT();

    // --- 1. Preparazione Dati (Lettura evento singolo) ---
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";
    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    // Usa Range per processare SOLO l'evento richiesto (molto veloce)
    if(eventID == -1)
    {
        gRandom->SetSeed(0);
        eventID = gRandom->Integer(reader.GetRaw().Count().GetValue());
    }

    reader.SetSingleEvent(eventID);
    auto df = reader.GetEstimators();

    // Take dei dati
    auto hits_ptr = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");
    if(hits_ptr->empty() || hits_ptr->at(0).empty())
    {
        cout << "[WARNING] Evento " << eventID << " vuoto o non trovato."
             << endl;
        return;
    }
    const auto &rvec = hits_ptr->at(0);
    vector<int> hit_ids(rvec.begin(), rvec.end());

    // Rimuovi duplicati e ordina
    sort(hit_ids.begin(), hit_ids.end());
    hit_ids.erase(unique(hit_ids.begin(), hit_ids.end()), hit_ids.end());

    // --- 2. Esecuzione del Fit ---
    cout << "\n=== Analysis Run " << runID << " Event " << eventID
         << " ===" << endl;
    cout << "Hits found: " << hit_ids.size() << endl;

    FitOutput fitOut = Do3DFit(hit_ids, false);
    RecoTrack tr = fitOut.track;

    if(tr.converged)
    {
        cout << "Fit CONVERGED!" << endl;
        printf("Position (y=0): x0=%.2f mm, z0=%.2f mm\n", tr.x0, tr.z0);
        printf("Slopes        : sx=%.4f, sz=%.4f\n", tr.sx, tr.sz);
        printf("Chi2/ndf      : %.3f\n", tr.chi2);
    }
    else
    {
        cout << "Fit FAILED or not enough hits." << endl;
    }

    // --- 3. Visualizzazione ---
    // Prepara la traccia fittata per il visualizer
    vector<Vis::VisLineTrack> visTracks;

    if(tr.converged)
        visTracks.emplace_back(tr.x0, 0.0, tr.z0, // Punto di start (a y=0)
            tr.sx, 1.0, tr.sz, // Vettore direzione (dx/dy, dy/dy, dz/dy)
            kBlack, 3 // Colore e spessore
        );

    // Disegna (Hit + Traccia Fit)
    Vis::Draw2D(hit_ids, Config::FindIntersections(hit_ids), visTracks);
    Vis::Draw3D(hit_ids, visTracks, fitOut.fittedPoints);
}

void AnalyzeCosmicRun(Int_t runID, Double_t toaMin, Double_t toaMax,
    UInt_t totMin, UInt_t totMax, Double_t minPValue = 0.01,
    Bool_t doEfficiency = true)
{
    // Disabilitiamo MT implicito per evitare conflitti con Minuit nel loop
    // manuale
    if(!ROOT::IsImplicitMTEnabled())
        ROOT::EnableImplicitMT();

    SetGlobalStyle();
    cout << "Starting Full Run Analysis for Run " << runID << "..." << endl;

    EffMap bundleStats, layerStats, cylStats;
    map<int, int> occupancyMap;

    // --- 1. Definizione Istogrammi ---
    TH1D *h_chi2 = new TH1D(
        "h_chi2", "Normalized #chi^{2};#chi^{2} / ndf;Events", 100, 0, 20);
    TH1D *h_prob = new TH1D("h_prob", "Prob(#chi^{2});Prob;Events", 50, 0, 1);
    TH1D *h_nhits
        = new TH1D("h_nhits", "Hits per Track;N Hits;Events", 15, -0.5, 50.5);

    TH1D *h_x0
        = new TH1D("h_x0", "Reconstructed X0;x_{0} [mm];Events", 80, -40, 40);
    TH1D *h_z0 = new TH1D(
        "h_z0", "Reconstructed Z0;z_{0} [mm];Events", 100, -150, 150);
    TH1D *h_sx
        = new TH1D("h_sx", "Slope X (sx);sx = dx/dy;Events", 100, -1.2, 1.2);
    TH1D *h_sz = new TH1D("h_sz", "Slope Z (sz);sz = dz/dy;Events", 100, -2, 2);
    TH1D *h_phi = new TH1D("h_phi", "Azimuthal Angle #phi; #phi [rad];Events",
        100, -TMath::Pi(), TMath::Pi());
    TH1D *h_cosTheta = new TH1D("h_cosTheta",
        "Polar Angle cos(#theta); cos(#theta); Events", 100, -1.0, 0.0);

    // --- 2. Estrazione Dati con RDataFrame ---
    string filename = "/home/lorenzo/muEDM_Project/Data/RootData/run00"
        + to_string(runID) + ".root";

    Data::Reader reader(filename);
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df = reader.GetEstimators();

    auto all_hits_ptr = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");
    const auto &all_hits = *all_hits_ptr;

    // --- 3. Loop di Fitting ---
    int nEvents = all_hits.size();
    cout << "Events to analyze: " << nEvents << endl;
    int nFitted = 0;

    for(int i = 0; i < nEvents; ++i)
    {
        const auto &rvec = all_hits[i];
        vector<int> hits(rvec.begin(), rvec.end());
        sort(hits.begin(), hits.end());
        hits.erase(unique(hits.begin(), hits.end()), hits.end());

        // Minimo 3 hits per tentare il fit
        if(hits.size() < 3)
            continue;

        if(i % 1000 == 0)
            printf("\rProcessing: %.1f%%", 100.0 * i / nEvents);

        // Fit geometrico (senza prior)
        FitOutput out = Do3DFit(hits, false);

        if(out.track.converged)
        {
            double chi2_red = out.track.chi2;
            double ndf = hits.size() * 2 - 4;
            double prob = TMath::Prob(chi2_red * ndf, (int)ndf);

            if(prob < minPValue)
                continue;

            nFitted++;

            // --- TRASFORMAZIONE COORDINATE (Analoga al ToyMC) ---
            double sx = out.track.sx;
            double sz = out.track.sz;

            // 1. Calcolo versori assumendo cosmici discendenti (uy < 0)
            double norm = sqrt(sx * sx + 1.0 + sz * sz);
            double uy = -1.0 / norm; // Forziamo verso il basso
            double ux = sx * uy; // Nota: sx = ux/uy -> ux = sx*uy
            double uz = sz * uy;

            // 2. Calcolo Angoli (Logica Z-X-Y: atan2(ux, uz))
            double phi = atan2(ux, uz);
            double cosTheta = uy;

            h_chi2->Fill(chi2_red);
            h_prob->Fill(prob);
            h_nhits->Fill(hits.size());

            h_x0->Fill(out.track.x0);
            h_z0->Fill(out.track.z0);

            h_sx->Fill(out.track.sx);
            h_sz->Fill(out.track.sz);
            h_phi->Fill(phi);
            h_cosTheta->Fill(cosTheta);

            // --- Calcolo efficienza su richiesta ---
            if(doEfficiency)
            {
                for(int h : hits)
                    occupancyMap[h]++;

                AccumulateEfficiency(
                    out.track, hits, bundleStats, layerStats, cylStats);
            }
        }
    }
    printf("\nAnalysis Complete. Fitted %d / %d events (%.1f%%)\n", nFitted,
        nEvents, 100.0 * nFitted / nEvents);

    // --- 4. Plotting ---
    TCanvas *cAn1
        = new TCanvas("cAn1", "Run Analysis Results - GOF", 800, 1200);
    cAn1->Divide(1, 3);

    cAn1->cd(1);
    gPad->SetLogy();
    h_chi2->SetLineColor(kBlue + 1);
    h_chi2->SetLineWidth(2);
    h_chi2->Draw();

    cAn1->cd(2);
    h_prob->SetLineColor(kMagenta + 1);
    h_prob->SetLineWidth(2);
    h_prob->SetMinimum(0);
    h_prob->DrawCopy();

    cAn1->cd(3);
    h_nhits->SetFillColor(kGray);
    h_nhits->DrawCopy();

    TCanvas *cAn2
        = new TCanvas("cAn2", "Run Analysis Results - Pars", 1200, 800);
    cAn2->Divide(2, 2);
    cAn2->cd(1);
    h_x0->SetLineColor(kRed + 1);
    h_x0->DrawCopy();

    cAn2->cd(2);
    h_z0->SetLineColor(kBlue + 1);
    h_z0->DrawCopy();

    cAn2->cd(3);
    h_sx->SetLineColor(kGreen + 2);
    h_sx->DrawCopy();
    // h_phi->SetLineColor(kGreen + 2);
    // h_phi->DrawCopy();

    cAn2->cd(4);
    h_sz->SetLineColor(kOrange + 1);
    h_sz->DrawCopy();
    // h_cosTheta->SetLineColor(kOrange + 1);
    // h_cosTheta->DrawCopy();

    if(doEfficiency)
        PlotEfficiencyResults(
            runID, bundleStats, layerStats, cylStats, occupancyMap);
}

// --- MONTE CARLO SIMULATION ---

void RunSingleHough()
{
    CosmicTrack tr = GenerateCosmic();
    vector<Int_t> hit_ids = FindHitBundles(tr);
    auto inters = Config::FindIntersections(hit_ids);
    vector<HoughResult> results = DoHoughTransform(hit_ids, 2);

    printf("\n--- RISULTATI HOUGH POLARE XY ---\n");
    for(size_t i = 0; i < results.size(); ++i)
    {
        printf("Candidato %lu: Rho=%.2f mm, Theta=%.2f rad (Peso: %.1f)\n",
            i + 1, results[i].rho, results[i].theta, results[i].weight);
    }

    // --- PREPARAZIONE VISUALIZZAZIONE ---
    vector<Vis::VisLineTrack> visTracks;

    // 1. Traccia Vera (Gialla)
    visTracks.emplace_back(
        tr.x0, tr.y0, tr.z0, tr.ux, tr.uy, tr.uz, kYellow, 3);

    // 2. Candidati Hough (convertiamo Rho/Theta in Linee 2D)
    int colors[2] = { kViolet - 9, kGreen + 2 };
    for(size_t i = 0; i < results.size(); ++i)
    {
        double r = results[i].rho;
        double t = results[i].theta;

        // Conversione Polare -> Punto + Direzione
        // Punto sulla retta più vicino all'origine: P(r*cos(t), r*sin(t))
        // Direzione perpendicolare alla normale: u(-sin(t), cos(t))
        double x0 = r * cos(t);
        double y0 = r * sin(t);
        double ux = -sin(t);
        double uy = cos(t);

        // Aggiungiamo come traccia puramente XY (z0=0, uz=0)
        // Style 2 = Dashed
        visTracks.emplace_back(x0, y0, 0, ux, uy, 0, colors[i % 2], 2, 2);
    }

    // Disegno
    Vis::Draw2D(hit_ids, inters, visTracks);
}

void RunSingleMC(double efficiency = 1.0)
{
    CosmicTrack trTrue = GenerateCosmic();
    vector<Int_t> hit_ids = FindHitBundles(trTrue, efficiency);

    if(hit_ids.size() < 3)
        return;

    FitOutput fitRes = Do3DFit(hit_ids, false);
    RecoTrack trFit = fitRes.track;

    // Stampa confronto
    if(trFit.converged)
    {
        // --- 1. Calcoli Lineari (Posizione e Slopes) ---
        double t_to_y0 = -trTrue.y0 / trTrue.uy;
        double true_x0 = trTrue.x0 + trTrue.ux * t_to_y0;
        double true_z0 = trTrue.z0 + trTrue.uz * t_to_y0;
        double true_sx = trTrue.ux / trTrue.uy;
        double true_sz = trTrue.uz / trTrue.uy;

        // --- 2. Calcoli Angolari (Phi e CosTheta) ---

        // A. True Angles
        double true_cosTheta = trTrue.uy;
        double true_phi = atan2(trTrue.ux, trTrue.uz);

        // B. Fitted Angles
        // Ricostruiamo il vettore dal fit (sx, 1, sz) e lo normalizziamo
        double norm = sqrt(trFit.sx * trFit.sx + 1.0 + trFit.sz * trFit.sz);

        // Forziamo uy negativo (cosmico discendente)
        double fit_uy = -1.0 / norm;
        double fit_ux = trFit.sx * fit_uy; // sx = ux/uy -> ux = sx*uy
        double fit_uz = trFit.sz * fit_uy; // sz = uz/uy -> uz = sz*uy

        double fit_cosTheta = fit_uy;
        double fit_phi = atan2(fit_ux, fit_uz);

        // --- 3. Stampa Tabella ---
        printf("\n");
        printf("==============================================================="
               "\n");
        printf("                   FIT RESULTS vs MC TRUTH                     "
               "\n");
        printf("==============================================================="
               "\n");
        printf(" %-10s | %12s | %12s | %12s \n", "Param", "True", "Fitted",
            "Residual");
        printf("---------------------------------------------------------------"
               "\n");
        // Parametri Spaziali
        printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "x0 [mm]", true_x0,
            trFit.x0, trFit.x0 - true_x0);
        printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "z0 [mm]", true_z0,
            trFit.z0, trFit.z0 - true_z0);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sx (dx/dy)", true_sx,
            trFit.sx, trFit.sx - true_sx);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sz (dz/dy)", true_sz,
            trFit.sz, trFit.sz - true_sz);
        printf("---------------------------------------------------------------"
               "\n");
        // Parametri Angolari
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "cos(theta)",
            true_cosTheta, fit_cosTheta, fit_cosTheta - true_cosTheta);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "phi [rad]", true_phi,
            fit_phi, fit_phi - true_phi);
        printf("---------------------------------------------------------------"
               "\n");
        printf(" Chi2/ndf  : %.4f\n", trFit.chi2);
        printf("==============================================================="
               "\n");

        DebugEfficiencyCalculation(trFit, hit_ids);
    }
    else
    {
        printf("\n[Fit Failed] Minuit did not converge.\n");
    }

    // --- PREPARAZIONE VISUALIZZAZIONE ---
    vector<Vis::VisLineTrack> visTracks;

    // 1. Traccia Vera
    visTracks.emplace_back(trTrue.x0, trTrue.y0, trTrue.z0, trTrue.ux,
        trTrue.uy, trTrue.uz, kYellow, 3);

    // 2. Traccia Fit
    if(trFit.converged)
        visTracks.emplace_back(trFit.x0, 0.0, trFit.z0, // Punto a y=0
            trFit.sx, 1.0, trFit.sz, // Vettore direzionale non normalizzato
            kBlack, 2, 7 // Tratteggiato
        );

    Vis::Draw2D(hit_ids, Config::FindIntersections(hit_ids), visTracks);
    Vis::Draw3D(hit_ids, visTracks, fitRes.fittedPoints);
}

void RunToyMC(int nEvents = 10000, double efficiency = 1.0,
    Double_t minPValue = 0.01, Bool_t doEfficiency = true,
    bool doDoubleFit = true)
{
    // --- 1. Definizione Istogrammi ---

    // a. Istogrammi Residui (Bias)
    TH1D *h_res_x0 = new TH1D("h_res_x0",
        "X0 Position Bias (Fit - True); #DeltaX0 [mm]; Counts", 100,
        -1 / efficiency, 1 / efficiency);
    TH1D *h_res_z0 = new TH1D("h_res_z0",
        "Z0 Position Bias (Fit - True); #DeltaZ0 [mm]; Counts", 100,
        -5 / efficiency, 5 / efficiency);
    TH1D *h_res_sx = new TH1D("h_res_sx",
        "X Slope Bias (sx_{fit} - sx_{true}); #Deltasx [rad]; Counts", 100,
        -0.1 / efficiency, 0.1 / efficiency);
    TH1D *h_res_sz = new TH1D("h_res_sz",
        "Z Slope Bias (sz_{fit} - sz_{true}); #Deltasz [rad]; Counts", 100,
        -0.5 / efficiency, 0.5 / efficiency);
    TH1D *h_res_phi = new TH1D("h_res_phi",
        "Bias Phi (Fit - True); #Delta #phi [rad]; Counts", 100, -3, 3);
    TH1D *h_res_cosTheta = new TH1D("h_res_cosTheta",
        "Bias CosTheta (Fit - True); #Delta cos#theta; Counts", 100, -0.02,
        0.02);

    // b. Istogrammi Variabili Fittate (Distribution)
    TH1D *h_fit_x0 = new TH1D(
        "h_fit_x0", "Reconstructed X0;x_{0} [mm];Events", 80, -40, 40);
    TH1D *h_fit_z0 = new TH1D(
        "h_fit_z0", "Reconstructed Z0;z_{0} [mm];Events", 100, -150, 150);
    TH1D *h_fit_sx = new TH1D("h_fit_sx",
        "Reconstructed Slope X (sx);sx = dx/dy;Events", 100, -1.2, 1.2);
    TH1D *h_fit_sz = new TH1D(
        "h_fit_sz", "Reconstructed Slope Z (sz);sz = dz/dy;Events", 100, -2, 2);
    TH1D *h_fit_phi = new TH1D("h_fit_phi", "Reco #phi; #phi [rad]; Events",
        100, -TMath::Pi(), TMath::Pi());
    TH1D *h_fit_cosTheta = new TH1D("h_fit_cosTheta",
        "Reco cos(#theta); cos(#theta); Events", 100, -1.0, 1.0);

    // c. Istogrammi Variabili True (Distribution)
    TH1D *h_true_x0 = new TH1D("h_true_x0",
        "True X0 Distribution;x_{0}^{true} [mm];Events", 80, -40, 40);
    TH1D *h_true_z0 = new TH1D("h_true_z0",
        "True Z0 Distribution;z_{0}^{true} [mm];Events", 100, -150, 150);
    TH1D *h_true_sx = new TH1D(
        "h_true_sx", "True Slope X (sx);sx^{true};Events", 100, -1.2, 1.2);
    TH1D *h_true_sz = new TH1D(
        "h_true_sz", "True Slope Z (sz);sz^{true};Events", 100, -2, 2);
    TH1D *h_true_phi = new TH1D("h_true_phi", "True #phi; #phi [rad];Events",
        100, -TMath::Pi(), TMath::Pi());
    TH1D *h_true_cosTheta = new TH1D("h_true_cosTheta",
        "True cos(#theta); cos(#theta); Events", 100, -1.0, 1.0);

    // Debug
    TH2D *hDebug = new TH2D("", "", 100, -3, 3, 100, -3, 3);
    TCanvas *cDbg = nullptr;

    // --- Contatori ---
    int n_triggers = 0;
    int n_geo_accepted = 0;
    int n_reconstructible = 0;
    int n_good_fit = 0;
    EffMap bundleStats, layerStats, cylStats;
    map<int, int> occupancyMap;

    printf("Starting Toy MC with %d triggered events...\n", nEvents);

    for(int i = 0; i < nEvents; ++i)
    {
        if(i % 100 == 0)
            printf("\rProcessing event %d/%d", i, nEvents);

        // 1. Generazione
        CosmicTrack trueTr = GenerateCosmic();
        n_triggers++;

        // 2. Accettanza Geometrica
        bool crosses_cylinder = IsGeometricIntersection(trueTr);
        if(crosses_cylinder)
            n_geo_accepted++;

        // 3. Digitizzazione
        vector<Int_t> hit_ids = FindHitBundles(trueTr, efficiency);

        // Debug
        if(DEBUG_TOY && efficiency > 0.99 && hit_ids.size() != 8
            && hit_ids.size() != 4 && hit_ids.size() != 0)
        {
            // DEBUG: Se siamo al 100% di efficienza, ci aspettiamo 0, 4 o 8
            // hits. Se ne troviamo un numero diverso, c'è un problema
            // geometrico o di logica.
            cout << "\n[DEBUG] Anomalous Event found! Hits: " << hit_ids.size()
                 << endl;

            if(!cDbg)
                cDbg = new TCanvas("cDbg", "Debug Anomalous Event", 1000, 800);
            else
                cDbg->Clear();

            cDbg->cd();

            // Visualizza la traccia vera in verde
            Vis::VisLineTrack vTr(trueTr.x0, trueTr.y0, trueTr.z0, trueTr.ux,
                trueTr.uy, trueTr.uz, kGreen + 2, 3);

            // Disegna tutto
            Vis::Draw3D(hit_ids, { vTr }, {}, true);

            cDbg->Update();
            gSystem->ProcessEvents();

            cout << "Press ENTER to continue processing..." << endl;
            if(cin.peek() == '\n')
                cin.ignore(); // Flush se necessario
            cin.get(); // Attesa input
        }

        if(hit_ids.size() < 3)
            continue;
        n_reconstructible++;

        // 4. Fitting
        FitOutput fitRes = Do3DFit(hit_ids, false);
        RecoTrack fitTr = fitRes.track;

        if(!fitTr.converged)
            continue;

        double chi2_red = fitTr.chi2;
        double ndf = hit_ids.size() * 2 - 4;
        double prob = TMath::Prob(chi2_red * ndf, (int)ndf);

        if(prob < minPValue)
            continue;

        n_good_fit++;

        // 5. Analisi e Confronto

        // A. Calcolo variabili TRUE nel formalismo del FIT (Slopes)
        // x = x0 + sx*y  => sx = ux/uy
        // z = z0 + sz*y  => sz = uz/uy
        // x(y=0) = x_gen + ux * t_0  dove t_0 porta a y=0 => t_0 = -y_gen/uy
        double t_to_y0 = -trueTr.y0 / trueTr.uy;
        double true_x0_at_y0 = trueTr.x0 + trueTr.ux * t_to_y0;
        double true_z0_at_y0 = trueTr.z0 + trueTr.uz * t_to_y0;
        double true_sx = trueTr.ux / trueTr.uy;
        double true_sz = trueTr.uz / trueTr.uy;
        double true_phi = atan2(trueTr.ux, trueTr.uz);
        double true_cosTheta = trueTr.uy;

        // B. Calcolo variabili FIT nel formalismo FISICO (Angoli)
        // Vettore direttore fit non normalizzato: V = (sx, 1, sz)
        double norm = sqrt(fitTr.sx * fitTr.sx + 1.0 + fitTr.sz * fitTr.sz);
        double fit_uy = -1.0 / norm;
        double fit_ux = fitTr.sx * fit_uy;
        double fit_uz = fitTr.sz * fit_uy;
        double fit_phi = atan2(fit_ux, fit_uz);
        double fit_cosTheta = fit_uy;

        // Riempimento Residui (Fit - True)
        h_res_x0->Fill(fitTr.x0 - true_x0_at_y0);
        h_res_z0->Fill(fitTr.z0 - true_z0_at_y0);
        h_res_sx->Fill(fitTr.sx - true_sx);
        h_res_sz->Fill(fitTr.sz - true_sz);

        double dPhi = fit_phi - true_phi;
        // Gestione attraversamento +/- Pi
        if(dPhi > M_PI)
            dPhi -= 2.0 * M_PI;
        if(dPhi < -M_PI)
            dPhi += 2.0 * M_PI;
        h_res_phi->Fill(dPhi);
        h_res_cosTheta->Fill(fit_cosTheta - true_cosTheta);

        // Riempimento Variabili Fittate
        h_fit_x0->Fill(fitTr.x0);
        h_fit_z0->Fill(fitTr.z0);
        h_fit_sx->Fill(fitTr.sx);
        h_fit_sz->Fill(fitTr.sz);
        h_fit_phi->Fill(fit_phi);
        h_fit_cosTheta->Fill(fit_cosTheta);

        // Riempimento Variabili True
        h_true_x0->Fill(true_x0_at_y0);
        h_true_z0->Fill(true_z0_at_y0);
        h_true_sx->Fill(true_sx);
        h_true_sz->Fill(true_sz);
        h_true_phi->Fill(true_phi);
        h_true_cosTheta->Fill(true_cosTheta);

        hDebug->Fill(true_phi, fit_phi);

        // Efficiency
        if(doEfficiency)
        {
            for(int h : hit_ids)
                occupancyMap[h]++;

            AccumulateEfficiency(
                fitTr, hit_ids, bundleStats, layerStats, cylStats);
        }
    }
    printf("\n\n");

    // --- Calcolo Statistiche Finali ---
    double acc_geo = 100.0 * (double)n_geo_accepted / (double)n_triggers;
    double eff_rec = 0.0;
    if(n_geo_accepted > 0)
        eff_rec = 100.0 * (double)n_reconstructible / (double)n_geo_accepted;

    printf("\n=== ACCEPTANCE RESULTS ===\n");
    printf("Total Triggers (Scintillators):       %d\n", n_triggers);
    printf("CHeT Geometric Acceptance:            %d (%.2f%%)\n",
        n_geo_accepted, acc_geo);
    printf("Events with Hits >= 3:                %d (%.2f%% of Geom. Acc.)\n",
        n_reconstructible, eff_rec);
    printf("Fit Converged:                        %d\n", n_good_fit);
    printf("=================================\n");

    // --- Plotting Bias (Canvas 1) ---
    TCanvas *c_res = new TCanvas("c_res", "Fit Residuals MC", 1200, 800);
    c_res->Divide(2, 2);

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);

    auto fitResidui = [](TH1D *h, bool useDoubleGauss)
    {
        if(!h || h->GetEntries() < 50)
            return;

        if(!useDoubleGauss)
        {
            h->Fit("gaus", "QL");
            if(h->GetFunction("gaus"))
                h->GetFunction("gaus")->SetLineColor(kRed);
        }
        else
        {
            TF1 *f2 = new TF1("f2",
                "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
                h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

            double mean = h->GetMean();
            double rms = h->GetRMS();
            double peak = h->GetMaximum();

            f2->SetParNames(
                "A_Core", "Mean", "Sigma_Core", "A_Tail", "Sigma_Tail");
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

            auto fitptr = h->Fit(f2, "SLMQ");

            if(fitptr->IsValid())
            {
                TF1 *fCore = new TF1("fCore", "gaus", h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
                fCore->SetParameters(f2->GetParameter(0), f2->GetParameter(1),
                    f2->GetParameter(2));
                fCore->SetLineColor(kGreen + 2);
                fCore->SetLineStyle(2);
                fCore->Draw("same");

                TF1 *fTail = new TF1("fTail", "gaus", h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
                fTail->SetParameters(f2->GetParameter(3), f2->GetParameter(1),
                    f2->GetParameter(4));
                fTail->SetLineColor(kMagenta);
                fTail->SetLineStyle(2);
                fTail->Draw("same");

                f2->Draw("same");
            }
        }
    };

    c_res->cd(1);
    h_res_x0->Draw();
    fitResidui(h_res_x0, doDoubleFit);
    c_res->cd(2);
    h_res_z0->Draw();
    fitResidui(h_res_z0, doDoubleFit);
    c_res->cd(3);
    h_res_sx->Draw();
    fitResidui(h_res_sx, doDoubleFit);
    // h_res_phi->Draw();
    // fitResidui(h_res_phi, doDoubleFit);
    c_res->cd(4);
    h_res_sz->Draw();
    fitResidui(h_res_sz, doDoubleFit);
    // h_res_cosTheta->Draw();
    // fitResidui(h_res_cosTheta, doDoubleFit);

    // --- Plotting Parametri Fittati (Canvas 2) ---
    TCanvas *c_pars
        = new TCanvas("c_pars", "Fitted Parameters Distribution", 1200, 800);
    c_pars->Divide(2, 2);

    c_pars->cd(1);
    h_fit_x0->SetLineColor(kRed + 1);
    h_fit_x0->Draw();

    c_pars->cd(2);
    h_fit_z0->SetLineColor(kBlue + 1);
    h_fit_z0->Draw();

    c_pars->cd(3);
    h_fit_sx->SetLineColor(kGreen + 2);
    h_fit_sx->Draw();
    // h_fit_phi->SetLineColor(kGreen + 2);
    // h_fit_phi->Draw();

    c_pars->cd(4);
    h_fit_sz->SetLineColor(kOrange + 1);
    h_fit_sz->Draw();
    // h_fit_cosTheta->SetLineColor(kOrange + 1);
    // h_fit_cosTheta->Draw();

    // --- Plotting Parametri True (Canvas 3) ---
    TCanvas *c_true = new TCanvas(
        "c_true", "True Parameters Distribution (Reconstructed)", 1200, 800);
    c_true->Divide(2, 2);

    // Stile: Azure + Fill per distinguerli
    int trueColor = kAzure - 3;

    c_true->cd(1);
    h_true_x0->SetLineColor(trueColor);
    h_true_x0->SetFillColorAlpha(trueColor, 0.3);
    h_true_x0->Draw();

    c_true->cd(2);
    h_true_z0->SetLineColor(trueColor);
    h_true_z0->SetFillColorAlpha(trueColor, 0.3);
    h_true_z0->Draw();

    c_true->cd(3);
    h_true_sx->SetLineColor(trueColor);
    h_true_sx->SetFillColorAlpha(trueColor, 0.3);
    h_true_sx->Draw();
    // h_true_phi->SetLineColor(trueColor);
    // h_true_phi->SetFillColorAlpha(trueColor, 0.3);
    // h_true_phi->Draw();

    c_true->cd(4);
    h_true_sz->SetLineColor(trueColor);
    h_true_sz->SetFillColorAlpha(trueColor, 0.3);
    h_true_sz->Draw();
    // h_true_cosTheta->SetLineColor(trueColor);
    // h_true_cosTheta->SetFillColorAlpha(trueColor, 0.3);
    // h_true_cosTheta->Draw();

    // new TCanvas();
    // hDebug->Draw();

    if(doEfficiency)
        PlotEfficiencyResults(
            -1, bundleStats, layerStats, cylStats, occupancyMap);

    if(DEBUG_TOY && cDbg)
        delete cDbg;

    return;
}
