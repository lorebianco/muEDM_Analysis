/**
 * UNIFIED ANALYSIS & SIMULATION SOFTWARE
 *
 * This file contains the complete pipeline for the muEDM analysis:
 * 1. Data Structures
 * 2. Math & Fitting Core (Minuit2, Hough)
 * 3. Simulation Engine (ToyMC, Digitization)
 * 4. Data Processing Engine (RDataFrame helpers)
 * 5. Plotting & High-Level Analysis Runners
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <sys/select.h>
#include <type_traits>
#include <unistd.h>
#include <vector>

// ROOT Includes
#include <Rtypes.h>
#include <RtypesCore.h>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Types.h>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TSpectrum2.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>

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
using RVecI = RVec<int>;
using RVecD = RVec<double>;

// Debug Flags
constexpr bool DEBUG_TOY = false;
constexpr bool DEBUG_RUN = false;

// ==============================================================================
// 1. DATA STRUCTURES
// ==============================================================================

// Struttura per il MC: Verità Monte Carlo
struct ToyEvent
{
    // Verità (Truth)
    double mc_x, mc_y, mc_z;
    double trk_x0, trk_y0, trk_z0;
    double trk_ux, trk_uy, trk_uz;
    double mc_E, trk_R, trk_cx, trk_cy; // Per Michel
    double trk_tmin, trk_tmax;

    // Dati (Digitization)
    vector<int> hit_ids;
};

struct CosmicTrack
{
    Double_t x0, y0, z0, ux, uy, uz;
    bool active;
};

// Struttura per il risultato del Fit (Parametrizzazione: x = x0 + sx*y, z = z0 + sz*y)
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
    vector<Vis::VisPoint3D> fittedPoints;
};

// Struttura dati passata al funtore di minimizzazione
struct FitData
{
    vector<Config::FiberProp> props;
};

// Struttura per Hough Transform
struct LinearHoughResult
{
    Double_t rho;
    Double_t theta;
    Double_t weight;
};

// Struttura per Hough Circolare (Piano XY)
struct CircularHoughResult
{
    Double_t xc;
    Double_t yc;
    Double_t R;
    Double_t weight;
};

struct Point3D
{
    double x, y, z;
};

struct ZHoughResult
{
    double z0;
    double dz_ds;
    double weight;
    vector<Point3D> track_hits;
    double t_min;
    double t_max;
};

// Strutture per calcolo Efficienza
struct EffStats
{
    long passed = 0;
    long total = 0;
};
using EffMap = map<int, EffStats>;

// Struttura per Michel Track Generation
struct MichelTrack
{
    double E_kin;
    double x0, y0; // Punto di emissione sul piano target
    double cx, cy, radius; // Parametri 2D Elica
    double z0, dz_dt; // Parametri Longitudinali Elica
    double t_min, t_max; // Limiti angolari di percorrenza
    double theta_rad, phi_dir; // Angoli cinematici al vertice
    bool hits_boundary; // Uscita dal limite radiale o dalle basi Z?
};

// Struttura per Systematics
struct SysParams
{
    double x0, z0, sx, sz;

    SysParams operator-(const SysParams &other) const
    {
        return { x0 - other.x0, z0 - other.z0, sx - other.sx, sz - other.sz };
    }
    SysParams &operator+=(const SysParams &other)
    {
        x0 += other.x0;
        z0 += other.z0;
        sx += other.sx;
        sz += other.sz;
        return *this;
    }
    SysParams operator/(double d) const
    {
        return { x0 / d, z0 / d, sx / d, sz / d };
    }
};

enum class AnalysisMode
{
    ALL,
    RANDOM,
    SINGLE
};

// ==============================================================================
// 2. GLOBAL UTILS & STYLE
// ==============================================================================

void SetOptimalConfig()
{
    // Config::SetOffsetExp(29. * TMath::DegToRad()); // 29 optimal
    Config::SetDelta2(Config::GetDelta2() + 4.5 * TMath::DegToRad()); // 4.5 optimal
}

void SetGlobalStyle()
{
    static bool isStyleSet = false;
    if(isStyleSet)
        return;
    gStyle->SetPalette(60);
    TColor::InvertPalette();
    isStyleSet = true;
}

void PrintSummary(
    long long total, long long passed, Double_t t0, Double_t t1, UInt_t tot0, UInt_t tot1)
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

// ==============================================================================
// 3. MATH & FITTING ENGINE
// ==============================================================================

// Da lavorare
double Track3DNeg2LogL_DOCA(const double *par, const FitData &data, bool usePrior)
{
    // Parametri della traccia
    const double x0 = par[0];
    const double sx = par[1];
    const double z0 = par[2];
    const double sz = par[3];

    // Vettore direzione della traccia: v = (sx, 1, sz)
    // Non serve normalizzarlo subito, lo faremo nel calcolo della distanza
    const double vx = sx;
    const double vy = 1.0;
    const double vz = sz;
    const double v2 = vx * vx + vy * vy + vz * vz; // Norma quadra del vettore direzione

    const double sigma2 = 0.3;
    const double sigmaL2 = 1e-6;

    double n2ll = 0.0;
    if(usePrior)
        n2ll = 2.0 * log(1.0 + sx * sx + sz * sz);

    for(size_t i = 0; i < data.props.size(); ++i)
    {
        // 1. Punto sulla fibra in base alla coordinata locale zi (parametro del
        // fit)
        const double zi_loc = par[4 + i];
        const auto &p = data.props[i];

        const double alpha = (zi_loc + Config::L_HALF) / (2.0 * Config::L_HALF);
        const double phi_f = p.phi0 + p.dir * alpha * M_PI;

        const double xf = p.r * cos(phi_f);
        const double yf = p.r * sin(phi_f);
        const double zf = zi_loc;

        // 2. Calcolo della DOCA al quadrato tra il punto P=(xf, yf, zf)
        //    e la retta Passante per P0=(x0, 0, z0) con direzione V=(sx, 1, sz)

        // Vettore W = P_fibra - P0_traccia
        const double wx = xf - x0;
        const double wy = yf - 0.0;
        const double wz = zf - z0;

        // Prodotto scalare (W . V)
        const double w_dot_v = wx * vx + wy * vy + wz * vz;

        // Distanza al quadrato: DOCA^2 = |W|^2 - (W.V)^2 / |V|^2
        const double w2 = wx * wx + wy * wy + wz * wz;
        const double doca2 = w2 - (w_dot_v * w_dot_v) / v2;

        // Aggiungiamo al log-likelihood (assumendo doca2 segua una gaussiana)
        n2ll += doca2 / sigma2;

        // 3. Penalty per sforamento lunghezza fisica della fibra
        if(abs(zi_loc) > Config::L_HALF)
        {
            double diff = abs(zi_loc) - Config::L_HALF;
            n2ll += (diff * diff) / sigmaL2;
        }
    }
    return n2ll;
}

double Track3DNeg2LogL(const double *par, const FitData &data, bool usePrior)
{
    const double x0 = par[0];
    const double sx = par[1];
    const double z0 = par[2];
    const double sz = par[3];
    const double sigma2 = 0.3; // Risoluzione bundle (2mm / sqrt(12))^2 ~ 0.3
    const double sigmaL2 = 1e-6; // Penalty forte per hit fuori lunghezza fisica

    double n2ll = 0.0;
    if(usePrior)
        n2ll = 2.0 * log(1.0 + sx * sx + sz * sz);

    for(size_t i = 0; i < data.props.size(); ++i)
    {
        const double zi = par[4 + i];
        const auto &p = data.props[i];

        const double alpha = (zi + Config::L_HALF) / (2.0 * Config::L_HALF);
        const double phi_fib = p.phi0 + p.dir * alpha * M_PI;

        const double x_f = p.r * cos(phi_fib);
        const double y_f = p.r * sin(phi_fib);

        const double x_trk = x0 + sx * y_f;
        const double z_trk = z0 + sz * y_f;

        n2ll += ((x_f - x_trk) * (x_f - x_trk) + (zi - z_trk) * (zi - z_trk)) / sigma2;

        if(abs(zi) > Config::L_HALF)
            n2ll += (pow(abs(zi) - Config::L_HALF, 2)) / sigmaL2;
    }
    return n2ll;
}

FitOutput Do3DFit(const vector<int> &hit_ids, bool usePrior = false)
{
    FitData data;
    for(int id : hit_ids)
        data.props.push_back(Config::GetFiberProp(id));

    int n_hits = data.props.size();
    if(n_hits < 3)
        return { { 0, 0, 0, 0, -1.0, false }, {} };

    unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min->SetPrintLevel(0);

    auto chi2_func = [&](const double *par) { return Track3DNeg2LogL(par, data, usePrior); };
    ROOT::Math::Functor f(chi2_func, 4 + n_hits);
    min->SetFunction(f);

    min->SetVariable(0, "x0", 0.0, 0.1);
    min->SetVariable(1, "sx", 0.0, 0.01);
    min->SetVariable(2, "z0", 0.0, 0.1);
    min->SetVariable(3, "sz", 0.0, 0.01);

    for(int i = 0; i < n_hits; ++i)
        min->SetVariable(4 + i, Form("z_hit_%d", i), 0, 5.0);

    bool conv = min->Minimize();
    const double *res = min->X();
    double chi2 = (n_hits * 2 > 4) ? min->MinValue() / (double)(n_hits * 2 - 4) : 0.0;

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
            output.fittedPoints.emplace_back(x_fit, y_fit, z_fit, kBlack, 24, 1.0, true);
        }
    }
    return output;
}

vector<LinearHoughResult> DoSimpleHoughTransform(
    const vector<Int_t> &hit_ids, int nCandidates = 2, bool drawGraphs = true)
{
    Int_t nBinsTheta = 80;
    Double_t thetaMin = 0.0, thetaMax = TMath::Pi();
    Int_t nBinsRho = 100;
    Double_t rhoMin = -100.0, rhoMax = 100.0;

    auto h_pol = make_unique<TH2D>("h_xy_pol", "Hough Polare XY;#theta [rad];#rho [mm]", nBinsTheta,
        thetaMin, thetaMax, nBinsRho, rhoMin, rhoMax);
    h_pol->SetDirectory(nullptr);

    auto inters = Config::FindIntersections(hit_ids);
    for(auto &p : inters)
    {
        for(int it = 1; it <= nBinsTheta; ++it)
        {
            double theta = h_pol->GetXaxis()->GetBinCenter(it);
            double rho = p.x_loc * cos(theta) + p.y_loc * sin(theta);
            h_pol->Fill(theta, rho);
        }
    }

    h_pol->Smooth(1);
    TSpectrum2 spec(nCandidates);
    TString opt = "noMarkov";
    if(!drawGraphs)
        opt += " goff";
    Int_t nFound = spec.Search(h_pol.get(), 2, opt, 0.1);

    Double_t *xPos = spec.GetPositionX();
    Double_t *yPos = spec.GetPositionY();

    vector<LinearHoughResult> candidates;
    for(int i = 0; i < nFound; ++i)
    {
        int bx = h_pol->GetXaxis()->FindBin(xPos[i]);
        int by = h_pol->GetYaxis()->FindBin(yPos[i]);
        candidates.push_back({ yPos[i], xPos[i], h_pol->GetBinContent(bx, by) });
    }
    sort(candidates.begin(), candidates.end(),
        [](const LinearHoughResult &a, const LinearHoughResult &b) { return a.weight > b.weight; });
    return candidates;
}

vector<CircularHoughResult> DoCircularHoughTransform(const vector<Int_t> &hit_ids,
    int nCandidates = 2,
    int combinatorial_threshold = 1000, // Nuovo: soglia per cambiare strategia
    int n_random_triplets = 10000, // Nuovo: numero di campioni per il metodo random
    bool drawGraphs = true)
{
    auto inters = Config::FindIntersections(hit_ids);
    if(inters.size() < 3)
        return {};

    double maxDim = 150.0;
    auto h_center = make_unique<TH2D>("h_circ_centers",
        "Hough Circle Centers (XY);X_c [mm];Y_c [mm]", 150, -maxDim, maxDim, 150, -maxDim, maxDim);
    h_center->SetDirectory(nullptr);

    // 1. POPOLAMENTO DELLO SPAZIO DEI CENTRI (con logica ibrida)
    int n = inters.size();

    // --- DECISIONE DELLA STRATEGIA ---
    if(n < combinatorial_threshold)
    {
        // STRATEGIA 1: COMBINATORIA ESATTA (veloce per pochi hit)
        if(drawGraphs)
            cout << "[Hough Info] Using exact combinatiorial logic (" << n << " hits < "
                 << combinatorial_threshold << ")" << endl;
        for(int i = 0; i < n - 2; ++i)
        {
            for(int j = i + 1; j < n - 1; ++j)
            {
                for(int k = j + 1; k < n; ++k)
                {
                    auto &p1 = inters[i];
                    auto &p2 = inters[j];
                    auto &p3 = inters[k];
                    double D = 2
                        * (p1.x_loc * (p2.y_loc - p3.y_loc) + p2.x_loc * (p3.y_loc - p1.y_loc)
                            + p3.x_loc * (p1.y_loc - p2.y_loc));
                    if(abs(D) < 1e-5)
                        continue;
                    double xc
                        = ((p1.x_loc * p1.x_loc + p1.y_loc * p1.y_loc) * (p2.y_loc - p3.y_loc)
                              + (p2.x_loc * p2.x_loc + p2.y_loc * p2.y_loc) * (p3.y_loc - p1.y_loc)
                              + (p3.x_loc * p3.x_loc + p3.y_loc * p3.y_loc) * (p1.y_loc - p2.y_loc))
                        / D;
                    double yc
                        = ((p1.x_loc * p1.x_loc + p1.y_loc * p1.y_loc) * (p3.x_loc - p2.x_loc)
                              + (p2.x_loc * p2.x_loc + p2.y_loc * p2.y_loc) * (p1.x_loc - p3.x_loc)
                              + (p3.x_loc * p3.x_loc + p3.y_loc * p3.y_loc) * (p2.x_loc - p1.x_loc))
                        / D;
                    if(abs(xc) < maxDim && abs(yc) < maxDim)
                    {
                        double R_triplet = sqrt(pow(p1.x_loc - xc, 2) + pow(p1.y_loc - yc, 2));
                        double d_origin = sqrt(xc * xc + yc * yc);
                        double R_target = 30.0, tol = 3.0;
                        bool intersects_target = (d_origin <= (R_triplet + R_target + tol))
                            && (d_origin >= (abs(R_triplet - R_target) - tol));
                        if(intersects_target)
                            h_center->Fill(xc, yc);
                    }
                }
            }
        }
    }
    else
    {
        // STRATEGIA 2: HOUGH RANDOMIZZATA (costo fisso per molti hit)
        if(drawGraphs)
            cout << "[Hough Info] Using random strategy (" << n
                 << " hits >= " << combinatorial_threshold << ")" << endl;
        gRandom->SetSeed(0); // Rende i risultati riproducibili
        for(int i = 0; i < n_random_triplets; ++i)
        {
            // Estrai 3 indici casuali DIVERSI
            int i1, i2, i3;
            do
            {
                i1 = gRandom->Integer(n);
                i2 = gRandom->Integer(n);
                i3 = gRandom->Integer(n);
            } while(i1 == i2 || i1 == i3 || i2 == i3);

            auto &p1 = inters[i1];
            auto &p2 = inters[i2];
            auto &p3 = inters[i3];
            double D = 2
                * (p1.x_loc * (p2.y_loc - p3.y_loc) + p2.x_loc * (p3.y_loc - p1.y_loc)
                    + p3.x_loc * (p1.y_loc - p2.y_loc));
            if(abs(D) < 1e-5)
                continue;
            double xc = ((p1.x_loc * p1.x_loc + p1.y_loc * p1.y_loc) * (p2.y_loc - p3.y_loc)
                            + (p2.x_loc * p2.x_loc + p2.y_loc * p2.y_loc) * (p3.y_loc - p1.y_loc)
                            + (p3.x_loc * p3.x_loc + p3.y_loc * p3.y_loc) * (p1.y_loc - p2.y_loc))
                / D;
            double yc = ((p1.x_loc * p1.x_loc + p1.y_loc * p1.y_loc) * (p3.x_loc - p2.x_loc)
                            + (p2.x_loc * p2.x_loc + p2.y_loc * p2.y_loc) * (p1.x_loc - p3.x_loc)
                            + (p3.x_loc * p3.x_loc + p3.y_loc * p3.y_loc) * (p2.x_loc - p1.x_loc))
                / D;
            if(abs(xc) < maxDim && abs(yc) < maxDim)
            {
                double R_triplet = sqrt(pow(p1.x_loc - xc, 2) + pow(p1.y_loc - yc, 2));
                double d_origin = sqrt(xc * xc + yc * yc);
                double R_target = 30.0, tol = 3.0;
                bool intersects_target = (d_origin <= (R_triplet + R_target + tol))
                    && (d_origin >= (abs(R_triplet - R_target) - tol));
                if(intersects_target)
                    h_center->Fill(xc, yc);
            }
        }
    }

    // 2. RICERCA E VISUALIZZAZIONE DEI PICCHI (logica invariata)
    vector<CircularHoughResult> candidates;
    TH2D *h_work = (TH2D *)h_center->Clone("h_work");
    int win = 2;
    for(int c = 0; c < nCandidates; ++c)
    {
        int bx, by, bz;
        h_work->GetMaximumBin(bx, by, bz);
        // double maxWeight = h_work->GetBinContent(bx, by);
        // if(maxWeight < 2.0)
        //     break;

        double sumW = 0, sumX = 0, sumY = 0;
        for(int i = bx - win; i <= bx + win; ++i)
        {
            for(int j = by - win; j <= by + win; ++j)
            {
                double w = h_work->GetBinContent(i, j);
                if(w > 0)
                {
                    sumW += w;
                    sumX += w * h_work->GetXaxis()->GetBinCenter(i);
                    sumY += w * h_work->GetYaxis()->GetBinCenter(j);
                    h_work->SetBinContent(i, j, 0);
                }
            }
        }

        if(sumW > 0)
        {
            double xc_peak = sumX / sumW;
            double yc_peak = sumY / sumW;

            TH1D h_radius("h_radius", "", 100, 0.0, 150.0);
            for(auto &p : inters)
            {
                double r = sqrt(pow(p.x_loc - xc_peak, 2) + pow(p.y_loc - yc_peak, 2));
                h_radius.Fill(r);
            }
            int best_r_bin = h_radius.GetMaximumBin();
            double R_est = h_radius.GetBinCenter(best_r_bin);

            double sumR = 0;
            int countR = 0;
            double bin_width = h_radius.GetBinWidth(1);
            for(auto &p : inters)
            {
                double r = sqrt(pow(p.x_loc - xc_peak, 2) + pow(p.y_loc - yc_peak, 2));
                if(abs(r - R_est) <= bin_width)
                {
                    sumR += r;
                    countR++;
                }
            }
            if(countR > 0)
                R_est = sumR / countR;

            candidates.push_back({ xc_peak, yc_peak, R_est, sumW });

            if(drawGraphs)
            {
                TString canvasName = Form("c_hough_radius_cand%d", c);
                TString histoName = Form("h_draw_radius_cand%d", c);
                TCanvas *c_radius = (TCanvas *)gROOT->FindObject(canvasName);
                if(!c_radius)
                    c_radius = new TCanvas(
                        canvasName, Form("Radius Hough Space (Candidate %d)", c), 500, 400);
                else
                    c_radius->Clear();
                c_radius->cd();
                c_radius->SetGridy();

                TH1D *h_draw_radius = (TH1D *)h_radius.Clone(histoName);
                h_draw_radius->SetDirectory(nullptr);
                h_draw_radius->SetTitle(
                    Form("Radius Distribution (Candidate %d);Radius [mm];Intersection Votes", c));
                h_draw_radius->SetFillColor(kAzure - 9);
                h_draw_radius->SetLineColor(kAzure + 1);
                h_draw_radius->DrawCopy("HIST BAR");

                double y_max = h_draw_radius->GetMaximum();
                h_draw_radius->SetMaximum(y_max);
                TLine *line = new TLine(R_est, 0, R_est, y_max);
                line->SetLineColor(kRed);
                line->SetLineWidth(2);
                line->SetLineStyle(2);
                line->Draw();
                c_radius->Update();
            }
        }
    }
    delete h_work;

    if(drawGraphs)
    {
        TCanvas *c_circ = (TCanvas *)gROOT->FindObject("c_hough_circ");
        if(!c_circ)
            c_circ = new TCanvas("c_hough_circ", "Circular Hough Space", 600, 500);
        else
            c_circ->Clear();
        c_circ->cd();

        h_center->SetStats(0);
        TH2D *h_draw = (TH2D *)h_center->DrawCopy("COLZ");
        if(h_draw)
            h_draw->SetDirectory(nullptr);

        for(auto &cand : candidates)
        {
            TMarker *m = new TMarker(cand.xc, cand.yc, 24);
            m->SetMarkerColor(kRed);
            m->SetMarkerSize(2.0);
            m->Draw();
        }
        c_circ->Update();
    }

    return candidates;
}

vector<ZHoughResult> DoZHoughTransform(const vector<Int_t> &hit_ids, double xc, double yc,
    double R_reco, int nCandidates = 1, bool drawGraphs = true, double tol_Z_fit = 15.0)
{
    struct Point3DLocal
    {
        double x, y, z, phi_raw;
    };
    vector<Point3DLocal> valid_points;
    double tolR = 5.;
    auto inters = Config::FindIntersections(hit_ids);
    for(auto &p : inters)
    {
        double r_point = sqrt(pow(p.x_loc - xc, 2) + pow(p.y_loc - yc, 2));
        if(abs(r_point - R_reco) <= tolR)
            valid_points.push_back(
                { p.x_loc, p.y_loc, p.z_loc, atan2(p.y_loc - yc, p.x_loc - xc) });
    }
    if(valid_points.empty())
        return {};

    double min_dz_ds = -3.0, max_dz_ds = 3.0;
    double min_z0 = -400.0, max_z0 = 400.0;
    auto h_z
        = make_unique<TH2D>("h_z_tiled", "Hough Z", 150, min_dz_ds, max_dz_ds, 200, min_z0, max_z0);
    h_z->SetDirectory(nullptr);

    int nSteps = 150;
    for(auto &p : valid_points)
    {
        for(int i = 0; i <= nSteps; ++i)
        {
            double dz_ds = min_dz_ds + i * (max_dz_ds - min_dz_ds) / nSteps;
            for(int k = -1; k <= 1; ++k)
            {
                double phi_k = p.phi_raw + 2.0 * M_PI * k;
                double s_k = R_reco * phi_k;
                double z0 = p.z - dz_ds * s_k;
                bool valid_cross = false;
                if(abs(dz_ds) > 1e-5)
                {
                    double phi_cross = -z0 / (dz_ds * R_reco);
                    if(phi_cross >= -M_PI && phi_cross <= M_PI)
                        valid_cross = true;
                }
                else if(k == 0)
                    valid_cross = true;

                if(valid_cross && z0 >= min_z0 && z0 <= max_z0)
                    h_z->Fill(dz_ds, z0);
            }
        }
    }

    vector<ZHoughResult> candidates;
    vector<pair<double, double>> raw_peaks;
    TH2D *h_work = (TH2D *)h_z->Clone("h_work_z");
    int win = 2;
    for(int c = 0; c < nCandidates; ++c)
    {
        int bx, by, bz;
        h_work->GetMaximumBin(bx, by, bz);
        if(h_work->GetBinContent(bx, by) < 2.0)
            break;
        double sumW = 0, sumX = 0, sumY = 0;
        for(int i = bx - win; i <= bx + win; ++i)
        {
            for(int j = by - win; j <= by + win; ++j)
            {
                double w = h_work->GetBinContent(i, j);
                if(w > 0)
                {
                    sumW += w;
                    sumX += w * h_work->GetXaxis()->GetBinCenter(i);
                    sumY += w * h_work->GetYaxis()->GetBinCenter(j);
                    h_work->SetBinContent(i, j, 0);
                }
            }
        }
        if(sumW > 0)
        {
            double dz_ds_peak = sumX / sumW;
            double z0_peak = sumY / sumW;
            raw_peaks.push_back({ dz_ds_peak, z0_peak });

            // --- ASSEGNAZIONE DELLE HIT ---
            // Associamo la hit alla Z attesa migliore, permettendo all'elica di evolversi (t_test
            // espanso)
            vector<Point3D> track_hits;
            // double tol_Z_fit = 15.0; // mm
            double min_t = 1e9, max_t = -1e9;

            for(auto &p : valid_points)
            {
                double min_dist = 1e9;
                double best_unrolled_t = 0;

                for(int k = -2; k <= 2; ++k)
                {
                    double t_test = p.phi_raw + 2.0 * M_PI * k;
                    double expected_z = z0_peak + dz_ds_peak * R_reco * t_test;
                    double dist = abs(p.z - expected_z);

                    if(dist < min_dist)
                    {
                        min_dist = dist;
                        best_unrolled_t = t_test;
                    }
                }

                if(min_dist <= tol_Z_fit)
                {
                    track_hits.push_back({ p.x, p.y, p.z });
                    if(best_unrolled_t < min_t)
                        min_t = best_unrolled_t;
                    if(best_unrolled_t > max_t)
                        max_t = best_unrolled_t;
                }
            }

            // Aggiungiamo il punto di incrocio del vertice a Z=0 nei bound visivi per garantire
            // che le TLine coprano anche l'origine e si veda che passa correttamente dentro il
            // range.
            double theta_v = 0.0;
            if(abs(dz_ds_peak) > 1e-5)
                theta_v = -z0_peak / (dz_ds_peak * R_reco);

            min_t = std::min(min_t, theta_v);
            max_t = std::max(max_t, theta_v);

            if(track_hits.empty())
            {
                min_t = -M_PI;
                max_t = M_PI;
            }

            candidates.push_back({ z0_peak, dz_ds_peak, sumW, track_hits, min_t, max_t });
        }
    }
    delete h_work;

    // 4. DISEGNO DEI GRAFICI
    if(drawGraphs)
    {
        // --- 4a. Canvas Hough Space ---
        TCanvas *c_z = (TCanvas *)gROOT->FindObject("c_hough_z");
        if(!c_z)
            c_z = new TCanvas("c_hough_z", "S-Z Linear Hough Space", 600, 500);
        else
            c_z->Clear();

        c_z->cd();
        h_z->SetStats(0);
        TH2D *h_draw = (TH2D *)h_z->DrawCopy("COLZ");
        if(h_draw)
            h_draw->SetDirectory(nullptr);

        for(auto &rp : raw_peaks)
        {
            TMarker *m = new TMarker(rp.first, rp.second, 24);
            m->SetMarkerColor(kRed);
            m->SetMarkerSize(2.0);
            m->Draw();
        }
        c_z->Update();

        // --- 4b. NUOVO PLOT: Hit mappate in Z-Phi (RAW + Tiled) ---
        TCanvas *c_phi_z_raw = (TCanvas *)gROOT->FindObject("c_phi_z_raw");
        if(!c_phi_z_raw)
            c_phi_z_raw = new TCanvas("c_phi_z_raw", "Raw Local Z-Phi", 800, 500);
        else
            c_phi_z_raw->Clear();

        c_phi_z_raw->cd();
        c_phi_z_raw->SetGrid();

        auto g_raw = new TGraph(); // Punti veri [-pi, pi]
        g_raw->SetMarkerStyle(20);
        g_raw->SetMarkerSize(1.2);
        g_raw->SetMarkerColor(kAzure + 1);

        auto g_tiled = new TGraph(); // Proiezioni visive (k = -1, 1)
        g_tiled->SetMarkerStyle(24);
        g_tiled->SetMarkerSize(1.2);
        g_tiled->SetMarkerColor(kAzure - 3);

        int pt_raw = 0;
        int pt_tiled = 0;
        for(auto &p : valid_points)
        {
            g_raw->SetPoint(pt_raw++, p.z, p.phi_raw);

            // Aggiungiamo k=-1 e k=1 per la visualizzazione continuata
            g_tiled->SetPoint(pt_tiled++, p.z, p.phi_raw - 2.0 * M_PI);
            g_tiled->SetPoint(pt_tiled++, p.z, p.phi_raw + 2.0 * M_PI);
        }

        if(g_raw->GetN() > 0)
        {
            TMultiGraph *mg = new TMultiGraph();
            mg->SetTitle(
                Form("Raw Hits w.r.t center (%.1f, %.1f) [Tiled];Z [mm];Local #phi [rad]", xc, yc));
            mg->Add(g_tiled, "P");
            mg->Add(g_raw, "P");
            mg->Draw("A");

            mg->GetXaxis()->SetLimits(-Config::L_HALF, Config::L_HALF);
            mg->GetYaxis()->SetRangeUser(
                -M_PI * 1.5, M_PI * 1.5); // Focus visivo leggermente espanso

            // --- Disegno delle Rette Fittate (Principale + Copie fantasma) ---
            for(size_t c = 0; c < candidates.size(); ++c)
            {
                auto &cand = candidates[c];

                // Disegniamo la retta "Canonica" (k=0) e le sue copie visive per far vedere
                // da dove passava la linearità.
                for(int k = -1; k <= 1; ++k)
                {
                    double shift_phi = 2.0 * M_PI * k;

                    // USIAMO i membri corretti t_min e t_max della struct
                    double phi_start = cand.t_min + shift_phi;
                    double phi_end = cand.t_max + shift_phi;

                    double z_start = cand.z0 + cand.dz_ds * R_reco * cand.t_min;
                    double z_end = cand.z0 + cand.dz_ds * R_reco * cand.t_max;

                    TLine *line = new TLine(z_start, phi_start, z_end, phi_end);

                    // La retta vera (k=0) è Rossa. Le sue proiezioni (k!=0) sono arancioni
                    // tratteggiate.
                    if(k == 0 && c == 0)
                    {
                        line->SetLineColor(kRed);
                        line->SetLineWidth(2);
                        line->SetLineStyle(1);
                    }
                    else
                    {
                        line->SetLineColor(kOrange + 1);
                        line->SetLineWidth(2);
                        line->SetLineStyle(2);
                    }
                    line->Draw("SAME");
                }
            }
        }
        c_phi_z_raw->Update();
    }

    return candidates;
}

// ==============================================================================
// 4. SIMULATION ENGINE
// ==============================================================================

CosmicTrack GenerateCosmic(
    double xC_up = 0., double zC_up = 0., double xC_down = 0., double zC_down = 0.)
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
        track.x0 = rnd.Uniform(xC_up - halfX_up, xC_up + halfX_up);
        track.z0 = rnd.Uniform(zC_up - halfZ_up, zC_up + halfZ_up);

        Double_t cosTheta, cosThetaSq_test, phi;
        do
        {
            cosTheta = -rnd.Uniform(0., 1.);
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

        if(abs(x_proj - xC_down) <= halfX_down && abs(z_proj - zC_down) <= halfZ_down)
            accepted = true;
    }
    return track;
}

vector<Int_t> FindHitBundles(const CosmicTrack &track_in, double efficiency = 1.0)
{
    CosmicTrack track = track_in;
    Config::ApplyInverseTransformation(track.x0, track.y0, track.z0);
    Config::ApplyInverseRotation(track.ux, track.uy, track.uz);

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
            Double_t A = track.ux * track.ux + track.uy * track.uy;
            Double_t B = 2.0 * (track.x0 * track.ux + track.y0 * track.uy);
            Double_t C = track.x0 * track.x0 + track.y0 * track.y0 - R * R;
            Double_t delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                Double_t t_sol[2] = { (-B + sqrt(delta)) / (2 * A), (-B - sqrt(delta)) / (2 * A) };
                for(int i = 0; i < 2; ++i)
                {
                    Double_t t = t_sol[i];
                    Double_t zi = track.z0 + track.uz * t;
                    if(abs(zi) <= Config::L_HALF)
                    {
                        Double_t xi = track.x0 + track.ux * t;
                        Double_t yi = track.y0 + track.uy * t;
                        Double_t phi_track = atan2(yi, xi);
                        for(int b = 0; b < lay.nBundles; ++b)
                        {
                            int b_id = current_global_offset + b;
                            Config::FiberProp p = Config::GetFiberProp(b_id);
                            Double_t alpha = (zi + Config::L_HALF) / (2.0 * Config::L_HALF);
                            Double_t phi_f = p.phi0 + p.dir * alpha * M_PI;
                            Double_t dphi
                                = abs(Config::wrap0_2pi(phi_track) - Config::wrap0_2pi(phi_f));
                            if(dphi > M_PI)
                                dphi = 2.0 * M_PI - dphi;
                            if(dphi < (M_PI / lay.nBundles) && gRandom->Rndm() <= efficiency)
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

bool IsGeometricIntersection(const CosmicTrack &tr_in)
{
    CosmicTrack tr = tr_in;
    Config::ApplyInverseTransformation(tr.x0, tr.y0, tr.z0);
    Config::ApplyInverseRotation(tr.ux, tr.uy, tr.uz);

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
                Double_t t_sol[2] = { (-B + sqrt(delta)) / (2 * A), (-B - sqrt(delta)) / (2 * A) };
                for(int i = 0; i < 2; ++i)
                    if(abs(tr.z0 + tr.uz * t_sol[i]) <= Config::L_HALF)
                        return true;
            }
        }
    }
    return false;
}

double GenerateMichelEnergy()
{
    while(true)
    {
        double x = gRandom->Uniform(0.0, 1.0);
        double y = gRandom->Uniform(0.0, 1.5);
        if(y <= (2.0 * x * x * (3.0 - 2.0 * x)))
            return x;
    }
}

MichelTrack GenerateMichelTrack(bool manualMode = false, double E_MeV = 50.0, double gamma = 0,
    double theta_deg = 45.0, double phi_deg = 0.0)
{
    const double E_max = 52.8;
    const double B_tesla = 2.89;
    const double m_e = 0.511;
    const double R_start = 30.0;
    const double R_limit = 199.0 / 2.0;
    MichelTrack tr;

    if(manualMode)
    {
        tr.E_kin = E_MeV;
        tr.theta_rad = theta_deg * M_PI / 180.0;
        tr.phi_dir = phi_deg * M_PI / 180.0;
    }
    else
    {
        double x = GenerateMichelEnergy();
        tr.E_kin = x * E_max;
        gamma = gRandom->Uniform(0.0, 2.0 * M_PI);
        double cosTheta = gRandom->Uniform(-0.99, 0.99);
        tr.theta_rad = std::acos(cosTheta);
        tr.phi_dir = gRandom->Uniform(0.0, 2.0 * M_PI);
    }
    tr.x0 = R_start * std::cos(gamma);
    tr.y0 = R_start * std::sin(gamma);
    double p = std::sqrt(tr.E_kin * tr.E_kin + 2.0 * m_e * tr.E_kin);
    double pt = p * std::sin(tr.theta_rad);
    tr.radius = pt / (0.29979 * std::abs(B_tesla));

    double omega0 = (B_tesla > 0) ? -1.0 : 1.0;
    tr.t_min = tr.phi_dir - omega0 * (M_PI / 2.0);
    while(tr.t_min > M_PI)
        tr.t_min -= 2.0 * M_PI;
    while(tr.t_min < -M_PI)
        tr.t_min += 2.0 * M_PI;

    tr.cx = tr.x0 - tr.radius * std::cos(tr.t_min);
    tr.cy = tr.y0 - tr.radius * std::sin(tr.t_min);
    double rho_c = std::sqrt(tr.cx * tr.cx + tr.cy * tr.cy);
    double phi_c = std::atan2(tr.cy, tr.cx);

    tr.dz_dt = omega0 * tr.radius / std::tan(tr.theta_rad);
    tr.z0 = -tr.dz_dt * tr.t_min;
    double delta_t_z = (std::abs(tr.dz_dt) > 1e-9) ? (Config::L_HALF / std::abs(tr.dz_dt)) : 1e9;
    double t_end_z = tr.t_min + omega0 * delta_t_z;
    tr.t_max = t_end_z;
    tr.hits_boundary = false;

    double arg
        = (R_limit * R_limit - rho_c * rho_c - tr.radius * tr.radius) / (2.0 * tr.radius * rho_c);
    if(std::abs(arg) <= 1.0)
    {
        double dt_acos = std::acos(arg);
        double closest_dt = 1e9;
        bool found_radial = false;
        for(double ang : { phi_c + dt_acos, phi_c - dt_acos })
        {
            double diff = ang - tr.t_min;
            if(omega0 < 0)
            {
                while(diff > 0)
                    diff -= 2.0 * M_PI;
                while(diff < -2.0 * M_PI)
                    diff += 2.0 * M_PI;
            }
            else
            {
                while(diff < 0)
                    diff += 2.0 * M_PI;
                while(diff > 2.0 * M_PI)
                    diff += 2.0 * M_PI;
            }
            if(std::abs(diff) < std::abs(closest_dt))
            {
                closest_dt = diff;
                found_radial = true;
            }
        }
        if(found_radial)
        {
            double t_radial = tr.t_min + closest_dt;
            if(std::abs(t_radial - tr.t_min) < std::abs(t_end_z - tr.t_min))
            {
                tr.t_max = t_radial;
                tr.hits_boundary = true;
            }
        }
    }
    return tr;
}

vector<Int_t> FindMichelHits(const MichelTrack &tr, double efficiency = 1.0)
{
    std::vector<int> hits;
    auto cylinders = Config::GetCylinders();
    int current_global_offset = 0;
    double rho_c = std::sqrt(tr.cx * tr.cx + tr.cy * tr.cy);
    double phi_c = std::atan2(tr.cy, tr.cx);
    double t_start_search = std::min(tr.t_min, tr.t_max);
    double t_end_search = std::max(tr.t_min, tr.t_max);

    for(const auto &cyl : cylinders)
    {
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int l = 0; l < 2; ++l)
        {
            const auto &lay = *layers[l];
            double arg_layer = (lay.radius * lay.radius - rho_c * rho_c - tr.radius * tr.radius)
                / (2.0 * tr.radius * rho_c);
            if(std::abs(arg_layer) <= 1.0)
            {
                double dt_layer = std::acos(arg_layer);
                double t_base[2] = { phi_c + dt_layer, phi_c - dt_layer };
                for(int k = -25; k <= 25; ++k)
                {
                    for(double t_b : t_base)
                    {
                        double t = t_b + 2.0 * M_PI * k;
                        if(t >= t_start_search && t <= t_end_search)
                        {
                            double z = tr.z0 + tr.dz_dt * t;
                            if(std::abs(z) <= Config::L_HALF)
                            {
                                double x = tr.cx + tr.radius * std::cos(t);
                                double y = tr.cy + tr.radius * std::sin(t);
                                double phi_hit = std::atan2(y, x);
                                for(int b = 0; b < lay.nBundles; ++b)
                                {
                                    int b_id = current_global_offset + b;
                                    auto prop = Config::GetFiberProp(b_id);
                                    double alpha = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                                    double phi_f = prop.phi0 + prop.dir * alpha * M_PI;
                                    double dphi = std::abs(
                                        Config::wrap0_2pi(phi_hit) - Config::wrap0_2pi(phi_f));
                                    if(dphi > M_PI)
                                        dphi = 2.0 * M_PI - dphi;
                                    if(dphi < (M_PI / lay.nBundles)
                                        && gRandom->Rndm() <= efficiency)
                                        hits.push_back(b_id);
                                }
                            }
                        }
                    }
                }
            }
            current_global_offset += lay.nBundles;
        }
    }
    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    return hits;
}

// ==============================================================================
// 5. DATA ANALYSIS HELPERS
// ==============================================================================

void AccumulateEfficiency(const RecoTrack &tr, const vector<int> &hit_ids, EffMap &bundleMap,
    EffMap &layerMap, EffMap &cylMap)
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
            int current_layer_start = global_bundle_offset;
            double A = tr.sx * tr.sx + 1.0;
            double B = 2.0 * tr.x0 * tr.sx;
            double C = tr.x0 * tr.x0 - R * R;
            double delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                double y_sol[2] = { (-B + sqrt(delta)) / (2 * A), (-B - sqrt(delta)) / (2 * A) };
                for(double y : y_sol)
                {
                    double z = tr.z0 + tr.sz * y;
                    if(abs(z) <= L_SAFE)
                    {
                        double x = tr.x0 + tr.sx * y;
                        double phi_hit = atan2(y, x);
                        int best_bundle_id = -1;
                        double min_dist_norm = 1e9;
                        double bundle_width_rad = 2.0 * M_PI / L.nBundles;
                        for(int b = 0; b < L.nBundles; ++b)
                        {
                            int gid = current_layer_start + b;
                            auto prop = Config::GetFiberProp(gid);
                            double alpha = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                            double phi_fib = prop.phi0 + prop.dir * alpha * M_PI;
                            double dphi
                                = abs(Config::wrap0_2pi(phi_hit) - Config::wrap0_2pi(phi_fib));
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
                                auto best_prop = Config::GetFiberProp(best_bundle_id);
                                for(int h : hit_ids)
                                {
                                    auto hit_prop = Config::GetFiberProp(h);
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
            global_bundle_offset += L.nBundles;
        }
        cyl_idx++;
    }
}

void PlotEfficiencyResults(int runID, const EffMap &bMap, const EffMap &lMap, const EffMap &cMap,
    const map<int, int> &occMap)
{
    // --- Configurazione Stili ---
    const int col_lay[] = { kRed, kOrange + 1, kBlue, kCyan + 1 };
    const char *lab_lay[] = { "Cyl 0 - In", "Cyl 0 - Out", "Cyl 1 - In", "Cyl 1 - Out" };

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
    TH1D *h_all = new TH1D("h_all", "Hit Occupancy (Signal vs Noise);Global Bundle ID;Counts",
        max_id + 1, -0.5, max_id + 0.5);
    TH1D *h_match = new TH1D("h_match", "Matched", max_id + 1, -0.5, max_id + 0.5); // Verde
    TH1D *h_noise = new TH1D("h_noise", "Unassigned", max_id + 1, -0.5,
        max_id + 0.5); // Rosso

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
    // Canvas Efficiencies
    // ---------------------------------------------------------
    TCanvas *cEff = new TCanvas("cEff", Form("Efficiencies Run %d", runID), 900, 1200);
    cEff->Divide(1, 3);

    // ---------------------------------------------------------
    // 1. Cylinder Efficiency
    // ---------------------------------------------------------
    cEff->cd(1);
    gPad->SetGridy();

    TH1D *h_frame_cyl = new TH1D("h_frame_cyl", "Cylinder Efficiency;;Efficiency", 2, -0.5, 1.5);
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
        tG.DrawLatex(0.60, 0.75, Form("Global Efficiency: %.2f%%", g_eff * 100));

        printf("\n--- GLOBAL EFFICIENCY ---\n %.2f%% (%ld/%ld)\n", g_eff * 100, g_passed, g_total);
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

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_cyl[i], eff * 100,
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

    TH1D *h_frame_lay = new TH1D("h_frame_lay", "Sub-Layer Efficiency;;Efficiency", 4, -0.5, 3.5);
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

            printf("%-15s: %.2f +%.2f/-%.2f %% (%ld/%ld)\n", lab_lay[i], eff * 100,
                errs.second * 100, errs.first * 100, it->second.passed, it->second.total);
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

    TH1D *h_frame_bund = new TH1D("h_frame_bund", "Bundle Efficiency;Global Bundle ID;Efficiency",
        max_id + 1, -0.5, max_id + 0.5);
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

void DebugEfficiencyCalculation(const RecoTrack &tr, const vector<int> &hit_ids)
{
    printf("\n");
    printf("************************************************************\n");
    printf("               DEBUG EFFICIENCY LOGIC                       \n");
    printf("************************************************************\n");
    printf("Track Params: x0=%.2f, z0=%.2f, sx=%.4f, sz=%.4f\n", tr.x0, tr.z0, tr.sx, tr.sz);

    const double L_SAFE = Config::L_HALF; // - 5.0;
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
                (lay_local_idx == 0 ? "Inner" : "Outer"), lay_global_idx, R, global_bundle_offset);

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
                double y_sol[2] = { (-B + sqrt(delta)) / (2 * A), (-B - sqrt(delta)) / (2 * A) };

                for(int sol_i = 0; sol_i < 2; ++sol_i)
                {
                    double y = y_sol[sol_i];
                    double z = tr.z0 + tr.sz * y;
                    double x = tr.x0 + tr.sx * y;

                    printf("    * Intersection %d at (x=%.1f, y=%.1f, z=%.1f): ", sol_i, x, y, z);

                    if(abs(z) > L_SAFE)
                    {
                        printf("[SKIP] Z out of bounds (safe limit %.1f)\n", L_SAFE);
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
                        double alpha = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                        double phi_fib = prop.phi0 + prop.dir * alpha * M_PI;

                        double dphi = abs(Config::wrap0_2pi(phi_hit) - Config::wrap0_2pi(phi_fib));
                        if(dphi > M_PI)
                            dphi = 2.0 * M_PI - dphi;

                        double dist_norm = dphi / bundle_width_rad;
                        if(dist_norm < min_dist_norm)
                        {
                            min_dist_norm = dist_norm;
                            best_bundle_id = gid;
                        }
                    }

                    printf("      -> Best Bundle: ID %d (Dist Norm: %.3f)\n", best_bundle_id,
                        min_dist_norm);

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
                            bool is_neighbour = (diff <= 1) || (diff == L.nBundles - 1);

                            // Controllo extra: deve essere dello stesso layer
                            if(is_neighbour
                                && Config::GetFiberProp(h).layerId
                                    == Config::GetFiberProp(best_bundle_id).layerId)
                            {
                                hit_bundle = true;
                                found_id = h;
                                break;
                            }
                        }

                        if(hit_bundle)
                        {
                            printf("         \033[1;32m[SUCCESS]\033[0m FOUND "
                                   "MATCH! Hit ID: %d "
                                   "(Delta: %d)\n",
                                found_id, found_id - best_bundle_id);
                        }
                        else
                        {
                            printf("         \033[1;31m[FAIL]\033[0m EXPECTED "
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

void InspectEventImpl(bool condition, const vector<int> &hit_ids,
    const vector<Vis::VisLineTrack> &tracks, const string &msg)
{
    static bool stop_debug = false;
    if(stop_debug || !condition)
        return;
    cout << "\n[DEBUG] " << msg << " | Hits: " << hit_ids.size() << endl;
    Vis::Draw3D(hit_ids, tracks, {}, true);
    cout << "Press ENTER for next event, 'q' to stop debugging..." << endl;
    while(true)
    {
        gSystem->ProcessEvents();
        fd_set readfds;
        FD_ZERO(&readfds);
        FD_SET(STDIN_FILENO, &readfds);
        struct timeval timeout;
        timeout.tv_sec = 0;
        timeout.tv_usec = 50000;
        int ready = select(STDIN_FILENO + 1, &readfds, nullptr, nullptr, &timeout);
        if(ready > 0)
        {
            char c = cin.get();
            if(c == 'q' || c == 'Q')
            {
                stop_debug = true;
                break;
            }
            if(c == '\n')
                break;
        }
    }
}

template <typename... Args>
void InspectEvent(bool condition, const vector<int> &hit_ids, Args &&...args)
{
    if(!condition)
        return;
    vector<Vis::VisLineTrack> tracks;
    string msg = "Event Inspection";
    auto process = [&](auto &&arg)
    {
        using T = decay_t<decltype(arg)>;
        if constexpr(is_same_v<T, CosmicTrack>)
            tracks.emplace_back(arg.x0, arg.y0, arg.z0, arg.ux, arg.uy, arg.uz, kYellow + 1, 3);
        else if constexpr(is_same_v<T, RecoTrack>)
        {
            if(arg.converged)
                tracks.emplace_back(arg.x0, 0.0, arg.z0, arg.sx, 1.0, arg.sz, kRed, 2, 7, true);
            else
                msg += " [Fit Failed]";
        }
        else if constexpr(is_convertible_v<T, string>)
            msg = string(arg);
    };
    (process(std::forward<Args>(args)), ...);
    InspectEventImpl(true, hit_ids, tracks, msg);
}

// ==============================================================================
// 6. SYSTEMATICS CORE
// ==============================================================================

struct SystematicCalculator
{
    static constexpr double M[6][4]
        = { { -0.0710, -0.0144, 0.7249, -0.0148 }, { 0.0006, -0.9636, -0.0305, 0.0410 },
              { 7.2e-05, 0.0014, -0.0049, -1.0240 }, { 0.0003, -0.0003, 1.0115, 0.0071 },
              { 0.0041, -0.0005, 1.0589, -0.0213 }, { 0.0019, 0.0501, 0.0197, 1.1175 } };

    static void Predict(
        double shiftX, double shiftZ, double rotX, double rotZ, double shearX, double shearZ)
    {
        double S[6] = { shiftX, shiftZ, rotX, rotZ, shearX, shearZ };
        double P[4] = { 0, 0, 0, 0 };
        for(int j = 0; j < 4; ++j)
            for(int i = 0; i < 6; ++i)
                P[j] += M[i][j] * S[i];
        printf("\n=== SYSTEMATIC PREDICTION ===\n");
        printf("Predicted: dx0=%.4f, dz0=%.4f, dsx=%.4f, dsz=%.4f\n", P[0], P[1], P[2], P[3]);
    }

    static void CalculateAlignment(double dx0, double dz0, double dsx, double dsz)
    {
        TMatrixD A(4, 6);
        for(int obs = 0; obs < 4; ++obs)
            for(int sys = 0; sys < 6; ++sys)
                A[obs][sys] = M[sys][obs];
        TVectorD b(4);
        b[0] = dx0;
        b[1] = dz0;
        b[2] = dsx;
        b[3] = dsz;
        TMatrixD AT(TMatrixD::kTransposed, A);
        TDecompSVD svd(AT);
        if(!svd.Decompose())
        {
            printf("SVD Failed.\n");
            return;
        }
        TMatrixD AT_pinv = svd.Invert();
        TMatrixD A_pinv(TMatrixD::kTransposed, AT_pinv);
        TVectorD x = A_pinv * b;
        printf("\n=== ALIGNMENT ESTIMATION ===\n");
        printf("Est: sX=%.4f, sZ=%.4f, rX=%.4f, rZ=%.4f, shX=%.4f, shZ=%.4f\n", x[0], x[1], x[2],
            x[3], x[4], x[5]);
    }
};

// ==============================================================================
// 7. HIGH LEVEL RUNNERS - SINGLE EVENT / VISUAL
// ==============================================================================

void RunSimpleHough()
{
    CosmicTrack tr = GenerateCosmic();
    vector<Int_t> hit_ids = FindHitBundles(tr);
    vector<LinearHoughResult> results = DoSimpleHoughTransform(hit_ids, 2);
    vector<Vis::VisLineTrack> visTracks;
    visTracks.emplace_back(tr.x0, tr.y0, tr.z0, tr.ux, tr.uy, tr.uz, kYellow, 3);
    int colors[2] = { kViolet - 9, kGreen + 2 };
    for(size_t i = 0; i < results.size(); ++i)
    {
        double r = results[i].rho;
        double t = results[i].theta;
        double x0 = r * cos(t);
        double y0 = r * sin(t);
        double ux = -sin(t);
        double uy = cos(t);
        visTracks.emplace_back(x0, y0, 0, ux, uy, 0, colors[i % 2], 2, 2, true);
    }
    Vis::Draw2D(hit_ids, visTracks);
}

void SimulateCosmicEvent(double efficiency = 1.0)
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
        // Convertiamo la traccia VERA (Globale) nel sistema LOCALE per il
        // confronto
        CosmicTrack trTrueLoc = trTrue;
        Config::ApplyInverseTransformation(trTrueLoc.x0, trTrueLoc.y0, trTrueLoc.z0);
        Config::ApplyInverseRotation(trTrueLoc.ux, trTrueLoc.uy, trTrueLoc.uz);

        double t_to_y0 = -trTrueLoc.y0 / trTrueLoc.uy;
        double true_x0 = trTrueLoc.x0 + trTrueLoc.ux * t_to_y0;
        double true_z0 = trTrueLoc.z0 + trTrueLoc.uz * t_to_y0;
        double true_sx = trTrueLoc.ux / trTrueLoc.uy;
        double true_sz = trTrueLoc.uz / trTrueLoc.uy;

        // --- 2. Calcoli Angolari (Phi e CosTheta) ---

        // A. True Angles
        double true_cosTheta = trTrueLoc.uy;
        double true_phi = atan2(trTrueLoc.ux, trTrueLoc.uz);

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
        printf(" %-10s | %12s | %12s | %12s \n", "Param", "True", "Fitted", "Residual");
        printf("---------------------------------------------------------------"
               "\n");
        // Parametri Spaziali
        printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "x0 [mm]", true_x0, trFit.x0,
            trFit.x0 - true_x0);
        printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "z0 [mm]", true_z0, trFit.z0,
            trFit.z0 - true_z0);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sx (dx/dy)", true_sx, trFit.sx,
            trFit.sx - true_sx);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sz (dz/dy)", true_sz, trFit.sz,
            trFit.sz - true_sz);
        printf("---------------------------------------------------------------"
               "\n");
        // Parametri Angolari
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "cos(theta)", true_cosTheta, fit_cosTheta,
            fit_cosTheta - true_cosTheta);
        printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "phi [rad]", true_phi, fit_phi,
            fit_phi - true_phi);
        printf("---------------------------------------------------------------"
               "\n");
        printf(" Chi2/ndf  : %.4f\n", trFit.chi2);
        printf("==============================================================="
               "\n");

        if(DEBUG_TOY)
            DebugEfficiencyCalculation(trFit, hit_ids);
    }
    else
    {
        printf("\n[Fit Failed] Minuit did not converge.\n");
    }

    // --- PREPARAZIONE VISUALIZZAZIONE ---
    vector<Vis::VisLineTrack> visTracks;

    // 1. Traccia Vera
    visTracks.emplace_back(
        trTrue.x0, trTrue.y0, trTrue.z0, trTrue.ux, trTrue.uy, trTrue.uz, kYellow, 3);

    // 2. Traccia Fit
    if(trFit.converged)
        visTracks.emplace_back(trFit.x0, 0.0, trFit.z0, // Punto a y=0
            trFit.sx, 1.0,
            trFit.sz, // Vettore direzionale non normalizzato
            kBlack, 2, 7, true // Tratteggiato, isLocal
        );

    Vis::Draw2D(hit_ids, visTracks);

    Vis::Draw3D(hit_ids, visTracks, fitRes.fittedPoints);
}

void SimulateMichelEvent(bool manualMode = false, double E_MeV = 50.0, double gamma = 0,
    double theta_deg = 45.0, double phi_deg = 0.0)
{
    gRandom->SetSeed(0);

    // 1. Genera Verità Monte Carlo (MC Truth)
    MichelTrack tr = GenerateMichelTrack(manualMode, E_MeV, gamma, theta_deg, phi_deg);
    vector<Int_t> hits = FindMichelHits(tr, 1.0);

    if(hits.empty())
    {
        cout << "[WARNING] Zero hits generated!" << endl;
        return;
    }

    // 2. Parametri di Ricerca
    int nCands2D = 1; // Cerchiamo fino a 2 cerchi
    int nCandsZ = 1; // Per ogni cerchio, cerchiamo fino a 2 soluzioni Z (tiling)

    // Esecuzione Hough 2D
    auto circ_cands = DoCircularHoughTransform(hits, nCands2D, 1000, 10000, true);

    // --- HEADER OUTPUT BELLO ---
    cout << "\n" << string(80, '=') << endl;
    cout << "   MICHEL EVENT RECONSTRUCTION SUMMARY (MULTI-CANDIDATE)" << endl;
    cout << string(80, '=') << endl;
    printf(" MC TRUTH | E: %5.2f MeV | Theta: %5.1f° | Phi: %5.1f° | R: %6.2f mm\n", tr.E_kin,
        tr.theta_rad * 180 / M_PI, tr.phi_dir * 180 / M_PI, tr.radius);
    printf(" MC TRUTH | Center: (%6.2f, %6.2f) | z0: %6.2f | dz/ds: %6.3f\n", tr.cx, tr.cy, tr.z0,
        tr.dz_dt / tr.radius);
    cout << string(80, '-') << endl;
    printf(" %-15s | %-10s | %-10s | %-10s | %-10s | %-8s\n", "CANDIDATE ID", "XC [mm]", "YC [mm]",
        "R [mm]", "Z0 [mm]", "dz/ds");
    cout << string(80, '-') << endl;

    std::vector<Vis::VisHelixTrack> visTracks;
    std::vector<Vis::VisPoint2D> extraPoints;

    // Traccia Vera (Gialla, spessa)
    visTracks.emplace_back(
        tr.cx, tr.cy, tr.radius, tr.z0, tr.dz_dt, tr.t_min, tr.t_max, kYellow, 6);
    extraPoints.emplace_back(tr.x0, tr.y0, kRed, 20, 1.2); // Punto di nascita

    int colorCycle[] = { kRed + 1, kAzure + 1, kSpring - 1, kOrange + 1 };
    int totalRecoCands = 0;

    // 3. Loop sui Candidati 2D (Cerchi)
    for(size_t i = 0; i < circ_cands.size(); ++i)
    {
        auto &c2D = circ_cands[i];

        // Per ogni cerchio, cerchiamo le soluzioni Z
        auto z_cands = DoZHoughTransform(hits, c2D.xc, c2D.yc, c2D.R, nCandsZ, true);

        if(z_cands.empty())
        {
            printf(" [2D-Cand %zu]     | %10.2f | %10.2f | %10.2f | %-21s\n", i, c2D.xc, c2D.yc,
                c2D.R, "NO Z-MATCH");
            continue;
        }

        // 4. Loop sui Candidati Z per il cerchio corrente
        for(size_t j = 0; j < z_cands.size(); ++j)
        {
            auto &cZ = z_cands[j];
            totalRecoCands++;

            // Stampa riga della tabella
            string id_str = Form("C%zu-Z%zu (W:%.0f)", i, j, cZ.weight);
            printf(" %-15s | %10.2f | %10.2f | %10.2f | %10.2f | %8.3f\n", id_str.c_str(), c2D.xc,
                c2D.yc, c2D.R, cZ.z0, cZ.dz_ds);

            // Calcolo residui rispetto al vero per questo candidato specifico
            double dR = c2D.R - tr.radius;
            double dZ0 = cZ.z0 - tr.z0;
            printf(" %-15s | %10s | %10s | %10.2f | %10.2f | %8.3f (RESIDUALS)\n", "", "", "", dR,
                dZ0, cZ.dz_ds - (tr.dz_dt / tr.radius));

            // Aggiunta al visualizer con colori diversi
            double reco_dz_dt = cZ.dz_ds * c2D.R;

            Vis::VisHelixTrack reco_helix(
                c2D.xc, c2D.yc, c2D.R, cZ.z0, reco_dz_dt, cZ.t_min, cZ.t_max, colorCycle[j], 4, 2);
            visTracks.push_back(reco_helix);
        }
        cout << string(80, '.') << endl;
    }

    cout << string(80, '=') << endl;
    cout << " TOTAL HITS: " << hits.size() << " | CLUSTERS FOUND: " << totalRecoCands << endl;
    cout << string(80, '=') << "\n" << endl;

    // --- Visualizzazione ---
    Vis::Draw3D(hits, visTracks);
    Vis::Draw2D(hits, visTracks, extraPoints);

    if(gSystem)
        gSystem->ProcessEvents();
}

// ==============================================================================
// 8. HIGH LEVEL RUNNERS - BATCH / ANALYSIS
// ==============================================================================

void RunCosmicToyMC(int nEvents = 10000, double efficiency = 1.0, Double_t minPValue = 0.,
    Bool_t doEfficiency = true, bool doDoubleFit = true)
{
    // --- 1. Definizione Istogrammi ---
    TDirectory::AddDirectory(false);

    // a. Istogrammi Residui (Bias)
    auto h_res_x0
        = make_unique<TH1D>("h_res_x0", "X0 Position Bias (Fit - True); #DeltaX0 [mm]; Counts", 100,
            -1 / efficiency, 1 / efficiency);
    auto h_res_z0
        = make_unique<TH1D>("h_res_z0", "Z0 Position Bias (Fit - True); #DeltaZ0 [mm]; Counts", 100,
            -5 / efficiency, 5 / efficiency);
    auto h_res_sx = make_unique<TH1D>("h_res_sx",
        "X Slope Bias (sx_{fit} - sx_{true}); #Deltasx [rad]; Counts", 100, -0.1 / efficiency,
        0.1 / efficiency);
    auto h_res_sz = make_unique<TH1D>("h_res_sz",
        "Z Slope Bias (sz_{fit} - sz_{true}); #Deltasz [rad]; Counts", 100, -0.5 / efficiency,
        0.5 / efficiency);
    auto h_res_phi = make_unique<TH1D>(
        "h_res_phi", "Bias Phi (Fit - True); #Delta #phi [rad]; Counts", 100, -3, 3);
    auto h_res_cosTheta = make_unique<TH1D>(
        "h_res_cosTheta", "Bias CosTheta (Fit - True); #Delta cos#theta; Counts", 100, -0.02, 0.02);
    auto h_res_z_hit = make_unique<TH1D>(
        "h_res_z_hit", "Z Hit Residuals (Fit - True); #Delta z_{hit} [mm]; Counts", 100, -10, 10);

    // b. Istogrammi Variabili Fittate (Distribution)
    auto h_fit_x0
        = make_unique<TH1D>("h_fit_x0", "Reconstructed X0;x_{0} [mm];Events", 80, -40, 40);
    auto h_fit_z0
        = make_unique<TH1D>("h_fit_z0", "Reconstructed Z0;z_{0} [mm];Events", 100, -150, 150);
    auto h_fit_sx = make_unique<TH1D>(
        "h_fit_sx", "Reconstructed Slope X (sx);sx = dx/dy;Events", 100, -1.2, 1.2);
    auto h_fit_sz
        = make_unique<TH1D>("h_fit_sz", "Reconstructed Slope Z (sz);sz = dz/dy;Events", 100, -2, 2);
    auto h_fit_phi = make_unique<TH1D>(
        "h_fit_phi", "Reco #phi; #phi [rad]; Events", 100, -TMath::Pi(), TMath::Pi());
    auto h_fit_cosTheta = make_unique<TH1D>(
        "h_fit_cosTheta", "Reco cos(#theta); cos(#theta); Events", 100, -1.0, 1.0);

    // c. Istogrammi Variabili True (Distribution)
    auto h_true_x0 = make_unique<TH1D>(
        "h_true_x0", "True X0 Distribution;x_{0}^{true} [mm];Events", 80, -40, 40);
    auto h_true_z0 = make_unique<TH1D>(
        "h_true_z0", "True Z0 Distribution;z_{0}^{true} [mm];Events", 100, -150, 150);
    auto h_true_sx
        = make_unique<TH1D>("h_true_sx", "True Slope X (sx);sx^{true};Events", 100, -1.2, 1.2);
    auto h_true_sz
        = make_unique<TH1D>("h_true_sz", "True Slope Z (sz);sz^{true};Events", 100, -2, 2);
    auto h_true_phi = make_unique<TH1D>(
        "h_true_phi", "True #phi; #phi [rad];Events", 100, -TMath::Pi(), TMath::Pi());
    auto h_true_cosTheta = make_unique<TH1D>(
        "h_true_cosTheta", "True cos(#theta); cos(#theta); Events", 100, -1.0, 1.0);

    // Debug
    // TH2D *hDebug = new TH2D("", "", 100, -3, 3, 100, -3, 3);

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
        CosmicTrack trueTr = GenerateCosmic(); // -9.75, 10, -9.75, -50 //
                                               // -9.75, -19, -9.75, -19 //
                                               // -41.3, 14.5, 41.3, -14.5
        n_triggers++;

        // 2. Accettanza Geometrica
        // IsGeometricIntersection gestisce internamente la rotazione inversa
        bool crosses_cylinder = IsGeometricIntersection(trueTr);
        if(crosses_cylinder)
            n_geo_accepted++;

        // 3. Digitizzazione
        // FindHitBundles gestisce internamente la rotazione inversa, quindi
        // passiamo trueTr (Globale)
        vector<Int_t> hit_ids = FindHitBundles(trueTr, efficiency);

        if(hit_ids.size() < 4)
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

        // Debug
        if(DEBUG_TOY)
        {
            // bool condition = hit_ids.size() != 8 && hit_ids.size() != 4 && hit_ids.size() != 0;
            Bool_t condition2 = fitTr.x0 >= -17 && fitTr.x0 <= -16;
            InspectEvent(condition2, hit_ids, trueTr, fitTr);
        }

        // if(fitTr.x0 <= -21 || fitTr.x0 >= -20)
        //     continue;

        if(prob < minPValue)
            continue;

        n_good_fit++;

        // 5. Analisi e Confronto

        // A. Calcolo variabili TRUE nel formalismo del FIT (Slopes)
        // x = x0 + sx*y  => sx = ux/uy
        // z = z0 + sz*y  => sz = uz/uy
        // x(y=0) = x_gen + ux * t_0  dove t_0 porta a y=0 => t_0 = -y_gen/uy
        // USIAMO LA TRACCIA RUOTATA (LOCALE) PER IL CONFRONTO
        CosmicTrack trueTrLoc = trueTr;
        Config::ApplyInverseTransformation(trueTrLoc.x0, trueTrLoc.y0, trueTrLoc.z0);
        Config::ApplyInverseRotation(trueTrLoc.ux, trueTrLoc.uy, trueTrLoc.uz);

        double t_to_y0 = -trueTrLoc.y0 / trueTrLoc.uy;
        double true_x0_at_y0 = trueTrLoc.x0 + trueTrLoc.ux * t_to_y0;
        double true_z0_at_y0 = trueTrLoc.z0 + trueTrLoc.uz * t_to_y0;
        double true_sx = trueTrLoc.ux / trueTrLoc.uy;
        double true_sz = trueTrLoc.uz / trueTrLoc.uy;
        double true_phi = atan2(trueTrLoc.ux, trueTrLoc.uz);
        double true_cosTheta = trueTrLoc.uy;

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

        // Fill Z hit residuals
        for(size_t k = 0; k < hit_ids.size(); ++k)
        {
            int id = hit_ids[k];
            auto prop = Config::GetFiberProp(id);
            double R = prop.r;

            // Intersection trueTrLoc with cylinder R (Local Frame)
            double A = trueTrLoc.ux * trueTrLoc.ux + trueTrLoc.uy * trueTrLoc.uy;
            double B = 2 * (trueTrLoc.x0 * trueTrLoc.ux + trueTrLoc.y0 * trueTrLoc.uy);
            double C = trueTrLoc.x0 * trueTrLoc.x0 + trueTrLoc.y0 * trueTrLoc.y0 - R * R;
            double delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                double t1 = (-B - sqrt(delta)) / (2 * A);
                double t2 = (-B + sqrt(delta)) / (2 * A);
                double z1 = trueTrLoc.z0 + trueTrLoc.uz * t1;
                double z2 = trueTrLoc.z0 + trueTrLoc.uz * t2;

                // Select closest to fiber phi to identify correct intersection
                double x1 = trueTrLoc.x0 + trueTrLoc.ux * t1;
                double y1 = trueTrLoc.y0 + trueTrLoc.uy * t1;
                double phi1 = atan2(y1, x1);

                double x2 = trueTrLoc.x0 + trueTrLoc.ux * t2;
                double y2 = trueTrLoc.y0 + trueTrLoc.uy * t2;
                double phi2 = atan2(y2, x2);

                double alpha1 = (z1 + Config::L_HALF) / (2.0 * Config::L_HALF);
                double phi_f1 = prop.phi0 + prop.dir * alpha1 * M_PI;
                double dphi1 = abs(Config::wrap0_2pi(phi1) - Config::wrap0_2pi(phi_f1));
                if(dphi1 > M_PI)
                    dphi1 = 2.0 * M_PI - dphi1;

                double alpha2 = (z2 + Config::L_HALF) / (2.0 * Config::L_HALF);
                double phi_f2 = prop.phi0 + prop.dir * alpha2 * M_PI;
                double dphi2 = abs(Config::wrap0_2pi(phi2) - Config::wrap0_2pi(phi_f2));
                if(dphi2 > M_PI)
                    dphi2 = 2.0 * M_PI - dphi2;

                double z_true = (dphi1 < dphi2) ? z1 : z2;

                if(k < fitRes.fittedPoints.size())
                {
                    double z_fit = fitRes.fittedPoints[k].z;
                    h_res_z_hit->Fill(z_fit - z_true);
                }
            }
        }

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

        // hDebug->Fill(true_phi, fit_phi);

        // Efficiency
        if(doEfficiency)
        {
            for(int h : hit_ids)
                occupancyMap[h]++;

            AccumulateEfficiency(fitTr, hit_ids, bundleStats, layerStats, cylStats);
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
    printf("CHeT Geometric Acceptance:            %d (%.2f%%)\n", n_geo_accepted, acc_geo);
    printf("Events with Hits >= 3:                %d (%.2f%% of Geom. Acc.)\n", n_reconstructible,
        eff_rec);
    printf("Fit Converged:                        %d\n", n_good_fit);
    printf("=================================\n");

    // --- Plotting Bias (Canvas 1) ---
    TCanvas *c_res = (TCanvas *)gROOT->FindObject("c_res");
    if(!c_res)
        c_res = new TCanvas("c_res", "Fit Residuals MC", 1200, 800);
    else
        c_res->Clear();
    c_res->Divide(2, 2);

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);

    auto fitResidui = [](auto h, bool useDoubleGauss)
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
            auto f2 = make_unique<TF1>("f2",
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
                auto fCore = make_unique<TF1>(
                    "fCore", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                fCore->SetParameters(f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2));
                fCore->SetLineColor(kGreen + 2);
                fCore->SetLineStyle(2);
                fCore->DrawCopy("same");

                auto fTail = make_unique<TF1>(
                    "fTail", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                fTail->SetParameters(f2->GetParameter(3), f2->GetParameter(1), f2->GetParameter(4));
                fTail->SetLineColor(kMagenta);
                fTail->SetLineStyle(2);
                fTail->DrawCopy("same");

                f2->DrawCopy("same");
            }
        }
    };

    c_res->cd(1);
    fitResidui(h_res_x0->DrawCopy(), doDoubleFit);
    c_res->cd(2);
    fitResidui(h_res_z0->DrawCopy(), doDoubleFit);
    c_res->cd(3);
    fitResidui(h_res_sx->DrawCopy(), doDoubleFit);
    // h_res_phi->Draw();
    // fitResidui(h_res_phi.get(), doDoubleFit);
    c_res->cd(4);
    fitResidui(h_res_sz->DrawCopy(), doDoubleFit);
    // h_res_cosTheta->Draw();
    // fitResidui(h_res_cosTheta, doDoubleFit);

    TCanvas *c_zres = (TCanvas *)gROOT->FindObject("c_zres");
    if(!c_zres)
        c_zres = new TCanvas("c_zres", "Z Hit Residuals", 600, 400);
    else
        c_zres->Clear();
    c_zres->cd();
    fitResidui(h_res_z_hit->DrawCopy(), doDoubleFit);

    // --- Plotting Parametri Fittati (Canvas 2) ---
    TCanvas *c_pars = (TCanvas *)gROOT->FindObject("c_pars");
    if(!c_pars)
        c_pars = new TCanvas("c_pars", "Fitted Parameters Distribution", 1200, 800);
    else
        c_pars->Clear();
    c_pars->Divide(2, 2);

    c_pars->cd(1);
    h_fit_x0->SetLineColor(kRed + 1);
    h_fit_x0->DrawCopy();

    c_pars->cd(2);
    h_fit_z0->SetLineColor(kBlue + 1);
    h_fit_z0->DrawCopy();

    c_pars->cd(3);
    h_fit_sx->SetLineColor(kGreen + 2);
    h_fit_sx->DrawCopy();
    // h_fit_phi->SetLineColor(kGreen + 2);
    // h_fit_phi->Draw();

    c_pars->cd(4);
    h_fit_sz->SetLineColor(kOrange + 1);
    h_fit_sz->DrawCopy();
    // h_fit_cosTheta->SetLineColor(kOrange + 1);
    // h_fit_cosTheta->Draw();

    // --- Plotting Parametri True (Canvas 3) ---
    TCanvas *c_true = (TCanvas *)gROOT->FindObject("c_true");
    if(!c_true)
        c_true = new TCanvas("c_true", "True Parameters Distribution (Reconstructed)", 1200, 800);
    else
        c_true->Clear();
    c_true->Divide(2, 2);

    // Stile: Azure + Fill per distinguerli
    int trueColor = kAzure - 3;

    c_true->cd(1);
    h_true_x0->SetLineColor(trueColor);
    h_true_x0->SetFillColorAlpha(trueColor, 0.3);
    h_true_x0->DrawCopy();

    c_true->cd(2);
    h_true_z0->SetLineColor(trueColor);
    h_true_z0->SetFillColorAlpha(trueColor, 0.3);
    h_true_z0->DrawCopy();

    c_true->cd(3);
    h_true_sx->SetLineColor(trueColor);
    h_true_sx->SetFillColorAlpha(trueColor, 0.3);
    h_true_sx->DrawCopy();
    // h_true_phi->SetLineColor(trueColor);
    // h_true_phi->SetFillColorAlpha(trueColor, 0.3);
    // h_true_phi->Draw();

    c_true->cd(4);
    h_true_sz->SetLineColor(trueColor);
    h_true_sz->SetFillColorAlpha(trueColor, 0.3);
    h_true_sz->DrawCopy();
    // h_true_cosTheta->SetLineColor(trueColor);
    // h_true_cosTheta->SetFillColorAlpha(trueColor, 0.3);
    // h_true_cosTheta->Draw();

    // new TCanvas();
    // hDebug->Draw();

    if(doEfficiency)
        PlotEfficiencyResults(-1, bundleStats, layerStats, cylStats, occupancyMap);

    return;
}

// ==============================================================================
// X. HOUGH EFFICIENCY & PURITY EVALUATION
// ==============================================================================

// Aggiunta per Efficienza e Purezza: calcolo DOCA tra un'intersezione e la traccia simulata
double GetDistanceToHelix(const CHeT::Config::BundlesIntersection &inter, const MichelTrack &tr)
{
    // 1. Distanza dal cilindro teorico (2D)
    double r_hit = std::sqrt(std::pow(inter.x - tr.cx, 2) + std::pow(inter.y - tr.cy, 2));
    double dist_xy = std::abs(r_hit - tr.radius);

    // 2. Distanza longitudinale (Z) calcolata srotolando la fase
    double phi = std::atan2(inter.y - tr.cy, inter.x - tr.cx);
    double dist_z = 1e9;

    if(std::abs(tr.dz_dt) > 1e-6)
    {
        double t_z = (inter.z - tr.z0) / tr.dz_dt;
        double k_float = (t_z - phi) / (2.0 * M_PI);
        int k = std::round(k_float);
        double t_best = phi + 2.0 * M_PI * k;
        dist_z = std::abs(inter.z - (tr.z0 + tr.dz_dt * t_best));
    }
    else
    {
        dist_z = std::abs(inter.z - tr.z0);
    }

    return std::sqrt(dist_xy * dist_xy + dist_z * dist_z);
}

struct EfficiencyResult
{
    int true_intersections_total = 0;
    int true_positives = 0;
    int false_positives = 0;

    double GetEfficiency() const
    {
        return true_intersections_total > 0 ? (double)true_positives / true_intersections_total
                                            : 0.0;
    }
    double GetPurity() const
    {
        return (true_positives + false_positives) > 0
            ? (double)true_positives / (true_positives + false_positives)
            : 0.0;
    }
};

EfficiencyResult EvaluateHoughTrack(const vector<Point3D> &cand_hits_local,
    const vector<int> &hit_ids, const MichelTrack &tr, double tol_mm = 2.0)
{
    EfficiencyResult res;
    auto all_inters = Config::FindIntersections(hit_ids);

    vector<bool> is_true(all_inters.size(), false);
    for(size_t i = 0; i < all_inters.size(); ++i)
    {
        double dist = GetDistanceToHelix(all_inters[i], tr);
        if(dist <= tol_mm)
        {
            is_true[i] = true;
            res.true_intersections_total++;
        }
    }

    for(const auto &ch : cand_hits_local)
    {
        bool found = false;
        bool matched_as_true = false;
        for(size_t i = 0; i < all_inters.size(); ++i)
        {
            double d_match = std::sqrt(std::pow(ch.x - all_inters[i].x_loc, 2)
                + std::pow(ch.y - all_inters[i].y_loc, 2)
                + std::pow(ch.z - all_inters[i].z_loc, 2));
            if(d_match < 1e-4)
            {
                found = true;
                if(is_true[i])
                    matched_as_true = true;
                break;
            }
        }

        if(found)
        {
            if(matched_as_true)
                res.true_positives++;
            else
                res.false_positives++;
        }
    }

    return res;
}

void RunROCForMichelHough(int nEvents = 50, double efficiency = 1.0)
{
    TDirectory::AddDirectory(false);
    gRandom->SetSeed(0);

    printf("Starting ROC Evaluation for Michel Hough Transform with %d events...\n", nEvents);

    // Valuteremo la ROC variando tol_Z_fit
    std::vector<double> thresholds = { 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25 };
    // std::vector<double> thresholds = { 15.0, 50.0, 100, 200 };
    std::vector<double> eff_avg(thresholds.size(), 0.0);
    std::vector<double> pur_avg(thresholds.size(), 0.0);
    std::vector<int> counts(thresholds.size(), 0);

    for(int i = 0; i < nEvents; ++i)
    {
        if(i % 5 == 0)
            cout << "\rProcessing event " << i << "/" << nEvents << flush;

        MichelTrack tr = GenerateMichelTrack(false);
        vector<int> hits = FindMichelHits(tr, efficiency);
        if(hits.size() < 3)
            continue;

        auto cands2D = DoCircularHoughTransform(hits, 1., 1000, 10000, false);
        if(cands2D.empty())
            continue;
        auto cand2D = cands2D[0];

        for(size_t t_idx = 0; t_idx < thresholds.size(); ++t_idx)
        {
            double tolZ = thresholds[t_idx];
            auto candsZ = DoZHoughTransform(hits, cand2D.xc, cand2D.yc, cand2D.R, 1, false, tolZ);
            if(candsZ.empty())
                continue;
            auto candZ = candsZ[0];

            EfficiencyResult res
                = EvaluateHoughTrack(candZ.track_hits, hits, tr, 2.0); // 2.0 mm true distance

            if(res.true_intersections_total > 0)
            {
                eff_avg[t_idx] += res.GetEfficiency();
                pur_avg[t_idx] += res.GetPurity();
                counts[t_idx]++;
            }
        }
    }
    cout << endl;

    TGraph *gROC = new TGraph();
    for(size_t t_idx = 0; t_idx < thresholds.size(); ++t_idx)
    {
        if(counts[t_idx] > 0)
        {
            eff_avg[t_idx] /= counts[t_idx];
            pur_avg[t_idx] /= counts[t_idx];
            gROC->SetPoint(gROC->GetN(), eff_avg[t_idx], pur_avg[t_idx]);
            printf("Threshold tolZ = %5.1f mm | Eff = %5.3f | Pur = %5.3f\n", thresholds[t_idx],
                eff_avg[t_idx], pur_avg[t_idx]);
        }
    }

    TCanvas *cROC = new TCanvas("cROC", "Efficiency vs Purity (ROC)", 800, 600);
    cROC->SetGrid();
    gROC->SetTitle("Hough Selection ROC;Efficiency;Purity");
    gROC->SetMarkerStyle(20);
    gROC->SetMarkerSize(1.5);
    gROC->SetMarkerColor(kAzure + 1);
    gROC->SetLineColor(kAzure + 1);
    gROC->SetLineWidth(2);
    gROC->Draw("APL");
}

void RunMichelToyMC(int nEvents = 1000, double efficiency = 1.0, bool doDoubleFit = true)
{
    TDirectory::AddDirectory(false);
    gRandom->SetSeed(0); // Start reproducible run

    TStopwatch timer;
    timer.Start();

    // --- 1. Definizione Istogrammi dei Residui ---
    auto h_res_xc = make_unique<TH1D>(
        "h_res_xc", "X_{c} Bias (Reco - True); #Delta X_{c}[mm]; Counts", 100, -20, 20);
    auto h_res_yc = make_unique<TH1D>(
        "h_res_yc", "Y_{c} Bias (Reco - True); #Delta Y_{c} [mm]; Counts", 100, -20, 20);
    auto h_res_R = make_unique<TH1D>(
        "h_res_R", "Radius Bias (Reco - True); #Delta R [mm]; Counts", 100, -15, 15);
    auto h_res_z0 = make_unique<TH1D>(
        "h_res_z0", "Z_{0} Bias (Reco - True); #Delta Z_{0}[mm]; Counts", 100, -60, 60);
    auto h_res_dz = make_unique<TH1D>(
        "h_res_dz", "dz/ds Bias (Reco - True); #Delta (dz/ds) [mm/mm]; Counts", 100, -0.4, 0.4);

    int n_generated = 0;
    int n_reconstructible = 0;
    int n_reconstructed = 0;

    printf("Starting Michel Toy MC with %d triggered events...\n", nEvents);

    // --- 2. Loop Eventi ---
    for(int i = 0; i < nEvents; ++i)
    {
        // Usiamo printf perché va su stdout, mentre cout è bloccato
        if(i % 10 == 0)
            cout << Form("\rProcessing event %d/%d", i, nEvents) << flush;

        // A. Generazione e Hit
        MichelTrack tr = GenerateMichelTrack(false);
        n_generated++;

        vector<int> hits = FindMichelHits(tr, efficiency);
        if(hits.size() < 3)
            continue;
        n_reconstructible++;

        // B. Ricostruzione 2D
        auto cands2D = DoCircularHoughTransform(hits, 1., 1000, 10000, false);
        if(cands2D.empty())
            continue;
        auto cand2D = cands2D[0];

        // C. Ricostruzione 3D
        auto candsZ = DoZHoughTransform(hits, cand2D.xc, cand2D.yc, cand2D.R, 1, false);
        if(candsZ.empty())
            continue;
        auto candZ = candsZ[0];

        n_reconstructed++;

        // D. Calcolo Bias
        double true_dz_ds = tr.dz_dt / tr.radius;

        h_res_xc->Fill(cand2D.xc - tr.cx);
        h_res_yc->Fill(cand2D.yc - tr.cy);
        h_res_R->Fill(cand2D.R - tr.radius);
        h_res_z0->Fill(candZ.z0 - tr.z0);
        h_res_dz->Fill(candZ.dz_ds - true_dz_ds);
    }

    // --- 1. Fermo il Timer e calcolo i minuti ---
    timer.Stop();
    double realTimeSeconds = timer.RealTime(); // Tempo effettivo trascorso
    // double cpuTimeSeconds = timer.CpuTime(); // Tempo totale di calcolo CPU

    double realTimeMinutes = realTimeSeconds / 60.0;
    double eventsPerSec = (realTimeSeconds > 0) ? nEvents / realTimeSeconds : 0;

    // --- 2. Stampa dei risultati ---
    string lineDouble(40, '=');
    string lineSimple(40, '-');

    printf("\n\n%s\n", lineDouble.c_str());
    printf("   MICHEL HOUGH RECO RESULTS\n");
    printf("%s\n", lineSimple.c_str());
    printf("  Total Generated    : %d\n", n_generated);
    printf("  Geom Accepted      : %d (%.1f%%)\n", n_reconstructible,
        100.0 * n_reconstructible / n_generated);
    printf("  Hough Reco         : %d (%.1f%%)\n", n_reconstructed,
        100.0 * n_reconstructed / n_reconstructible);
    printf("%s\n", lineSimple.c_str());

    printf("  Simulation Duration: %.2f minutes\n", realTimeMinutes);
    printf("  Processing Speed   : %.1f events/sec\n", eventsPerSec);
    printf("%s\n\n", lineDouble.c_str());

    // --- 3. Visualizzazione Risultati ---
    TCanvas *c_res = (TCanvas *)gROOT->FindObject("c_res_michel");
    if(!c_res)
        c_res = new TCanvas("c_res_michel", "Michel Hough Residuals", 1400, 800);
    else
        c_res->Clear();
    c_res->Divide(3, 2);

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);

    // Lambda per il fit Gaussiano (uguale al tuo ToyMC cosmico)
    auto fitResidui = [](auto h, bool useDoubleGauss)
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
            auto f2 = make_unique<TF1>("f2",
                "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
                h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            double mean = h->GetMean();
            double rms = h->GetRMS();
            double peak = h->GetMaximum();
            f2->SetParameters(peak, mean, rms * 0.5, peak * 0.01, rms * 2.0);
            f2->SetParLimits(0, peak * 0.75, peak * 2);
            f2->SetParLimits(2, 0, rms);
            f2->SetParLimits(3, 0, peak * 0.25);
            f2->SetParLimits(4, rms, rms * 10);
            f2->SetLineColor(kRed);
            f2->SetNpx(1000);
            if(h->Fit(f2.get(), "SLMQ")->IsValid())
            {
                auto fC = make_unique<TF1>(
                    "fC", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                fC->SetParameters(f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2));
                fC->SetLineColor(kGreen + 2);
                fC->SetLineStyle(2);
                fC->DrawCopy("same");
                auto fT = make_unique<TF1>(
                    "fT", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                fT->SetParameters(f2->GetParameter(3), f2->GetParameter(1), f2->GetParameter(4));
                fT->SetLineColor(kMagenta);
                fT->SetLineStyle(2);
                fT->DrawCopy("same");
                f2->DrawCopy("same");
            }
        }
    };

    c_res->cd(1);
    fitResidui(h_res_xc->DrawCopy(), doDoubleFit);
    c_res->cd(2);
    fitResidui(h_res_yc->DrawCopy(), doDoubleFit);
    c_res->cd(3);
    fitResidui(h_res_R->DrawCopy(), doDoubleFit);
    c_res->cd(4);
    fitResidui(h_res_z0->DrawCopy(), doDoubleFit);
    c_res->cd(5);
    fitResidui(h_res_dz->DrawCopy(), doDoubleFit);

    // Ultimo Pad informativo
    c_res->cd(6);
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->AddText("Michel Track Seeding Stats");
    pt->AddText(Form("Events Generated: %d", n_generated));
    pt->AddText(Form("Acc. Geometric: %d", n_reconstructible));
    pt->AddText(Form("Pattern Reco: %d", n_reconstructed));
    pt->SetFillColor(kWhite);
    pt->SetTextFont(42);
    pt->SetTextSize(0.06);
    pt->Draw();

    c_res->Update();
}

// ==============================================================================
// 9. DATASET GENERATORS
// ==============================================================================

void GenerateCosmicDataset(string filename, int nEvents, double efficiency = 1.0)
{
    TFile f(filename.c_str(), "RECREATE");
    TTree *tree = new TTree("CosmicToy", "Toy Cosmic Events");
    ToyEvent ev;
    tree->Branch("trk_x0", &ev.trk_x0);
    tree->Branch("trk_y0", &ev.trk_y0);
    tree->Branch("trk_z0", &ev.trk_z0);
    tree->Branch("trk_ux", &ev.trk_ux);
    tree->Branch("trk_uy", &ev.trk_uy);
    tree->Branch("trk_uz", &ev.trk_uz);
    tree->Branch("hits", &ev.hit_ids);
    for(int i = 0; i < nEvents; ++i)
    {
        CosmicTrack tr = GenerateCosmic();
        ev.trk_x0 = tr.x0;
        ev.trk_y0 = tr.y0;
        ev.trk_z0 = tr.z0;
        ev.trk_ux = tr.ux;
        ev.trk_uy = tr.uy;
        ev.trk_uz = tr.uz;
        ev.hit_ids = FindHitBundles(tr, efficiency);
        if(!ev.hit_ids.empty())
            tree->Fill();
    }
    tree->Write();
    f.Close();
}

void GenerateMichelDataset(string filename, int nEvents, double efficiency = 1.0)
{
    TFile f(filename.c_str(), "RECREATE");
    TTree *tree = new TTree("MichelToy", "Toy Michel Events");
    ToyEvent ev;
    tree->Branch("mc_E", &ev.mc_E);
    tree->Branch("trk_R", &ev.trk_R);
    tree->Branch("trk_cx", &ev.trk_cx);
    tree->Branch("trk_cy", &ev.trk_cy);
    tree->Branch("trk_z0", &ev.trk_z0);
    tree->Branch("trk_uz", &ev.trk_uz);
    tree->Branch("trk_tmin", &ev.trk_tmin);
    tree->Branch("trk_tmax", &ev.trk_tmax);

    for(int i = 0; i < nEvents; ++i)
    {
        MichelTrack tr = GenerateMichelTrack(false);
        ev.mc_E = tr.E_kin;
        ev.trk_R = tr.radius;
        ev.trk_cx = tr.cx;
        ev.trk_cy = tr.cy;
        ev.trk_z0 = tr.z0;
        ev.trk_uz = tr.dz_dt;
        ev.trk_tmin = tr.t_min;
        ev.trk_tmax = tr.t_max;
        ev.hit_ids = FindMichelHits(tr, efficiency);
        if(!ev.hit_ids.empty())
            tree->Fill();
    }
    tree->Write();
    f.Close();
}

// ==============================================================================
// 10. SYSTEMATICS RUNNERS
// ==============================================================================

SysParams RunSystematicToy(int nEvents, double efficiency, double shearX = 0., double shearZ = 0.,
    bool useResiduals = false)
{
    double tx, ty, tz, rx, ry, rz;
    Config::GetTranslation(tx, ty, tz);
    Config::GetRotation(rx, ry, rz);
    SysParams sumObs = { 0, 0, 0, 0 };
    int nGood = 0;
    for(int i = 0; i < nEvents; ++i)
    {
        double xC_up = shearX / 2.0;
        double zC_up = shearZ / 2.0;
        CosmicTrack trTrue = GenerateCosmic(xC_up, zC_up, -xC_up, -zC_up);
        vector<int> hits = FindHitBundles(trTrue, efficiency);
        if(hits.size() < 4)
            continue;
        Config::SetTranslation(0, 0, 0);
        Config::SetRotation(0, 0, 0);
        FitOutput fitOut = Do3DFit(hits, false);
        Config::SetTranslation(tx, ty, tz);
        Config::SetRotation(rx, ry, rz);
        if(fitOut.track.converged)
        {
            sumObs.x0 += fitOut.track.x0;
            sumObs.z0 += fitOut.track.z0;
            sumObs.sx += fitOut.track.sx;
            sumObs.sz += fitOut.track.sz;
            if(useResiduals)
            {
                // To get residuals we need true track in local frame of fit (0,0,0)
                // but the track was generated with shear.
                // For simplified systematics, usually we just look at the shift of the mean
                // so we don't subtract truth here unless explicitly asked.
                // If asked, we need to project trTrue to y=0.
                CosmicTrack trTrueLoc = trTrue; // Assuming "perfect" detector for truth
                double t = -trTrueLoc.y0 / trTrueLoc.uy;
                sumObs.x0 -= (trTrueLoc.x0 + trTrueLoc.ux * t);
                sumObs.z0 -= (trTrueLoc.z0 + trTrueLoc.uz * t);
                sumObs.sx -= (trTrueLoc.ux / trTrueLoc.uy);
                sumObs.sz -= (trTrueLoc.uz / trTrueLoc.uy);
            }
            nGood++;
        }
    }
    if(nGood > 0)
        return sumObs / nGood;
    return { 0, 0, 0, 0 };
}

void StudySystematicSensitivityGrid(Double_t efficiency, Double_t rangePos = 2.0,
    double rangeRot = 0.005, Double_t rangeShear = 20., int nSteps = 5, int nEventsPerPoint = 5000,
    bool useResiduals = false)
{
    if(ROOT::IsImplicitMTEnabled())
        ROOT::DisableImplicitMT();
    const char *sysLabs[] = { "Shift X", "Shift Z", "Rot X", "Rot Z", "Shear X", "Shear Z" };
    const char *parLabs[] = { "x0 [mm]", "z0 [mm]", "sx [mrad]", "sz [mrad]" };
    double sensitivity[6][4];
    TCanvas *cGrid = new TCanvas("cSensGrid", "Systematic Sensitivity Grid", 1200, 1500);
    cGrid->Divide(4, 6, 0.001, 0.001);

    Config::SetTranslation(0, 0, 0);
    Config::SetRotation(0, 0, 0);
    SysParams baseline = RunSystematicToy(nEventsPerPoint, efficiency, 0, 0, useResiduals);

    for(int iSys = 0; iSys < 6; ++iSys)
    {
        double minVal, maxVal;
        if(iSys < 2)
        {
            minVal = -rangePos;
            maxVal = rangePos;
        }
        else if(iSys < 4)
        {
            minVal = -rangeRot;
            maxVal = rangeRot;
        }
        else
        {
            minVal = -rangeShear;
            maxVal = rangeShear;
        }
        double step = (maxVal - minVal) / (nSteps - 1);
        TGraph *graphs[4];
        for(int k = 0; k < 4; ++k)
            graphs[k] = new TGraph();

        printf("Scanning %s ... ", sysLabs[iSys]);
        fflush(stdout);

        for(int pt = 0; pt < nSteps; ++pt)
        {
            double val = minVal + pt * step;
            Config::SetTranslation(0, 0, 0);
            Config::SetRotation(0, 0, 0);
            double sX = 0., sZ = 0.;
            if(iSys == 0)
                Config::SetTranslation(val, 0, 0);
            if(iSys == 1)
                Config::SetTranslation(0, 0, val);
            if(iSys == 2)
                Config::SetRotation(val, 0, 0);
            if(iSys == 3)
                Config::SetRotation(0, 0, val);
            if(iSys == 4)
                sX = val;
            if(iSys == 5)
                sZ = val;

            SysParams res = RunSystematicToy(nEventsPerPoint, efficiency, sX, sZ, useResiduals);
            SysParams net = res - baseline;
            double plotVal = (iSys >= 2 && iSys < 4) ? val * 1000.0 : val;

            graphs[0]->SetPoint(pt, plotVal, net.x0);
            graphs[1]->SetPoint(pt, plotVal, net.z0);
            graphs[2]->SetPoint(pt, plotVal, net.sx * 1000.0);
            graphs[3]->SetPoint(pt, plotVal, net.sz * 1000.0);
        }
        printf("Done.\n");

        for(int iPar = 0; iPar < 4; ++iPar)
        {
            int padIdx = iSys * 4 + iPar + 1;
            cGrid->cd(padIdx);
            gPad->SetGrid();
            graphs[iPar]->Fit("pol1", "Q");
            TF1 *fit = graphs[iPar]->GetFunction("pol1");
            double slope = (fit) ? fit->GetParameter(1) : 0.0;
            sensitivity[iSys][iPar] = slope;
            graphs[iPar]->SetMarkerStyle(20);
            graphs[iPar]->Draw("AP");
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.08);
            lat.DrawLatex(0.2, 0.8, Form("m = %.3f", slope));
        }
    }

    TCanvas *cMat = new TCanvas("cSensMat", "Sensitivity Matrix Heatmap", 800, 600);
    cMat->cd();
    gStyle->SetPaintTextFormat("5.3f");
    gStyle->SetPalette(kTemperatureMap);
    TH2D *hMat = new TH2D("hSensMat", "Sensitivity Coefficients;Fitted Parameter;Systematic Source",
        4, 0, 4, 6, 0, 6);
    for(int r = 0; r < 6; ++r)
    {
        for(int c = 0; c < 4; ++c)
            hMat->SetBinContent(c + 1, 6 - r, sensitivity[r][c]);
        hMat->GetYaxis()->SetBinLabel(6 - r, sysLabs[r]);
    }
    for(int c = 0; c < 4; ++c)
        hMat->GetXaxis()->SetBinLabel(c + 1, parLabs[c]);
    hMat->Draw("COLZ TEXT");
}

void SolveSystematics(
    double sX = 0, double sZ = 0, double rX = 0, double rZ = 0, double shX = 0, double shZ = 0)
{
    SystematicCalculator::Predict(sX, sZ, rX, rZ, shX, shZ);
}

void FindAlignment(double dx0, double dz0, double dsx, double dsz)
{
    SystematicCalculator::CalculateAlignment(dx0, dz0, dsx * 1e3, dsz * 1e3);
}

// ==============================================================================
// 11. UNIFIED ANALYZER
// ==============================================================================

void AnalyzeToyDataset(
    string filename, bool isCosmic = true, AnalysisMode mode = AnalysisMode::ALL, int entryID = -1)
{
    bool doDoubleFit = true;

    TFile *f = TFile::Open(filename.c_str(), "READ");
    if(!f || f->IsZombie())
    {
        cerr << "[Error] Could not open file: " << filename << endl;
        return;
    }
    TTree *treeCHeT = (TTree *)f->Get("chet");
    TTree *treeSIM = (TTree *)f->Get("sim");
    if(!treeCHeT || !treeSIM)
    {
        cerr << "[Error] Required trees (Event, sim) not found." << endl;
        return;
    }

    ToyEvent ev;
    std::vector<int> *hits_ptr = &ev.hit_ids;
    treeCHeT->SetBranchAddress("All_Bundle", &hits_ptr);

    treeSIM->SetBranchAddress("mc_x", &ev.mc_x);
    treeSIM->SetBranchAddress("mc_y", &ev.mc_y);
    treeSIM->SetBranchAddress("mc_z", &ev.mc_z);

    if(isCosmic)
    {
        treeSIM->SetBranchAddress("trk_x0", &ev.trk_x0);
        treeSIM->SetBranchAddress("trk_y0", &ev.trk_y0);
        treeSIM->SetBranchAddress("trk_z0", &ev.trk_z0);
        treeSIM->SetBranchAddress("trk_ux", &ev.trk_ux);
        treeSIM->SetBranchAddress("trk_uy", &ev.trk_uy);
        treeSIM->SetBranchAddress("trk_uz", &ev.trk_uz);
    }
    else
    {
        treeSIM->SetBranchAddress("mc_E", &ev.mc_E);
        treeSIM->SetBranchAddress("trk_R", &ev.trk_R);
        treeSIM->SetBranchAddress("trk_cx", &ev.trk_cx);
        treeSIM->SetBranchAddress("trk_cy", &ev.trk_cy);
        treeSIM->SetBranchAddress("trk_z0", &ev.trk_z0);
        treeSIM->SetBranchAddress("trk_uz", &ev.trk_uz);
        treeSIM->SetBranchAddress("trk_tmin", &ev.trk_tmin);
        treeSIM->SetBranchAddress("trk_tmax", &ev.trk_tmax);
    }

    TDirectory::AddDirectory(false);
    // Histograms definition
    // Cosmics
    auto h_res_x0 = make_unique<TH1D>(
        "h_res_x0", "X0 Position Bias (Fit - True); #DeltaX0 [mm]; Counts", 100, 0, 0);
    auto h_res_z0 = make_unique<TH1D>(
        "h_res_z0", "Z0 Position Bias (Fit - True); #DeltaZ0 [mm]; Counts", 100, 0, 0);
    auto h_res_sx = make_unique<TH1D>(
        "h_res_sx", "X Slope Bias (sx_{fit} - sx_{true}); #Deltasx [rad]; Counts", 100, 0, 0);
    auto h_res_sz = make_unique<TH1D>(
        "h_res_sz", "Z Slope Bias (sz_{fit} - sz_{true}); #Deltasz [rad]; Counts", 100, 0, 0);
    auto h_fit_x0
        = make_unique<TH1D>("h_fit_x0", "Reconstructed X0;x_{0} [mm];Events", 80, -40, 40);
    auto h_fit_z0
        = make_unique<TH1D>("h_fit_z0", "Reconstructed Z0;z_{0} [mm];Events", 100, -150, 150);
    auto h_fit_sx = make_unique<TH1D>(
        "h_fit_sx", "Reconstructed Slope X (sx);sx = dx/dy;Events", 100, -1.2, 1.2);
    auto h_fit_sz
        = make_unique<TH1D>("h_fit_sz", "Reconstructed Slope Z (sz);sz = dz/dy;Events", 100, -2, 2);

    // Michel
    auto h_res_xc = make_unique<TH1D>("h_res_xc", "X_{c} Residual; #Delta X_{c}[mm]", 100, -20, 20);
    auto h_res_yc
        = make_unique<TH1D>("h_res_yc", "Y_{c} Residual; #Delta Y_{c} [mm]", 100, -20, 20);
    auto h_res_R = make_unique<TH1D>("h_res_R", "Radius Residual; #Delta R [mm]", 100, -15, 15);
    auto h_res_z = make_unique<TH1D>("h_res_z", "Z_{0} Residual; #Delta Z_{0}[mm]", 100, -60, 60);
    auto h_res_dz = make_unique<TH1D>("h_res_dz", "dz/ds Residual; #Delta (dz/ds)", 100, -0.4, 0.4);

    int nEntries = treeCHeT->GetEntries();
    vector<int> entriesToProcess;
    if(mode == AnalysisMode::ALL)
        for(int i = 0; i < nEntries; ++i)
            entriesToProcess.push_back(i);
    else if(mode == AnalysisMode::RANDOM)
    {
        gRandom->SetSeed(0);
        entriesToProcess.push_back(gRandom->Integer(nEntries));
    }
    else if(entryID >= 0)
        entriesToProcess.push_back(entryID);

    // --- 4. Loop di Analisi ---
    int nGood = 0;
    int current = 0;
    int total = entriesToProcess.size();

    for(int idx : entriesToProcess)
    {
        current++; // Incrementiamo la posizione attuale

        // Contatore elegante ogni 100 eventi o all'ultimo evento
        if(mode == AnalysisMode::ALL && (current % 100 == 0 || current == total))
        {
            double progress = (double)current / total * 100.0;
            printf("\r[Progress] Analyzing event %d / %d (%5.1f%%) ", current, total, progress);
            fflush(stdout);
        }

        treeCHeT->GetEntry(idx);
        treeSIM->GetEntry(idx);

        if(isCosmic)
        {
            if(ev.hit_ids.size() < 4)
                continue;

            FitOutput out = Do3DFit(ev.hit_ids, false);
            if(out.track.converged)
            {
                Config::ApplyInverseTransformation(ev.trk_x0, ev.trk_y0, ev.trk_z0);
                Config::ApplyInverseRotation(ev.trk_ux, ev.trk_uy, ev.trk_uz);

                double t = -ev.trk_y0 / ev.trk_uy;
                double true_x0 = ev.trk_x0 + ev.trk_ux * t;
                double true_z0 = ev.trk_z0 + ev.trk_uz * t;
                double true_sx = ev.trk_ux / ev.trk_uy;
                double true_sz = ev.trk_uz / ev.trk_uy;

                // Variables
                h_fit_x0->Fill(out.track.x0);
                h_fit_z0->Fill(out.track.z0);
                h_fit_sx->Fill(out.track.sx);
                h_fit_sz->Fill(out.track.sz);
                // Residuals
                h_res_x0->Fill(out.track.x0 - true_x0);
                h_res_z0->Fill(out.track.z0 - true_z0);
                h_res_sx->Fill(out.track.sx - true_sx);
                h_res_sz->Fill(out.track.sz - true_sz);

                nGood++;

                if(mode != AnalysisMode::ALL)
                {
                    vector<Vis::VisLineTrack> tracks;
                    tracks.emplace_back(ev.trk_x0, ev.trk_y0, ev.trk_z0, ev.trk_ux, ev.trk_uy,
                        ev.trk_uz, kYellow, 3);
                    tracks.emplace_back(out.track.x0, 0, out.track.z0, out.track.sx, 1,
                        out.track.sz, kRed, 2, 7, true);
                    Vis::Draw3D(ev.hit_ids, tracks, out.fittedPoints);
                }
            }
        }
        else
        {
            auto cands2D
                = DoCircularHoughTransform(ev.hit_ids, 1, 1000, 10000, (mode != AnalysisMode::ALL));
            if(!cands2D.empty())
            {
                auto cand2D = cands2D[0];
                auto candsZ = DoZHoughTransform(
                    ev.hit_ids, cand2D.xc, cand2D.yc, cand2D.R, 1, (mode != AnalysisMode::ALL));
                if(!candsZ.empty())
                {
                    auto candZ = candsZ[0];
                    h_res_xc->Fill(cand2D.xc - ev.trk_cx);
                    h_res_yc->Fill(cand2D.yc - ev.trk_cy);
                    h_res_R->Fill(cand2D.R - ev.trk_R);
                    h_res_z->Fill(candZ.z0 - ev.trk_z0);
                    double true_dz_ds = ev.trk_uz / ev.trk_R;
                    h_res_dz->Fill(candZ.dz_ds - true_dz_ds);
                    nGood++;

                    if(mode != AnalysisMode::ALL)
                    {
                        vector<Vis::VisHelixTrack> tracks;
                        tracks.emplace_back(ev.trk_cx, ev.trk_cy, ev.trk_R, ev.trk_z0, ev.trk_uz,
                            ev.trk_tmin, ev.trk_tmax, kYellow, 4);
                        tracks.emplace_back(cand2D.xc, cand2D.yc, cand2D.R, candZ.z0,
                            candZ.dz_ds * cand2D.R, candZ.t_min, candZ.t_max, kRed, 2);
                        Vis::Draw3D(ev.hit_ids, tracks);
                    }
                }
            }
        }
    }

    if(mode == AnalysisMode::ALL)
    {
        printf("\n");
        printf("Reconstruction: %d / %zu events.\n", nGood, entriesToProcess.size());
        gStyle->SetOptFit(111);
        if(isCosmic)
        {
            TCanvas *c_res = (TCanvas *)gROOT->FindObject("c_res");
            if(!c_res)
                c_res = new TCanvas("c_res", "Fit Residuals MC", 1200, 800);
            else
                c_res->Clear();
            c_res->Divide(2, 2);

            gStyle->SetOptStat(1110);
            gStyle->SetOptFit(111);

            auto fitResidui = [](auto h, bool useDoubleGauss)
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
                    auto f2 = make_unique<TF1>("f2",
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
                        auto fCore = make_unique<TF1>(
                            "fCore", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                        fCore->SetParameters(
                            f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2));
                        fCore->SetLineColor(kGreen + 2);
                        fCore->SetLineStyle(2);
                        fCore->DrawCopy("same");

                        auto fTail = make_unique<TF1>(
                            "fTail", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
                        fTail->SetParameters(
                            f2->GetParameter(3), f2->GetParameter(1), f2->GetParameter(4));
                        fTail->SetLineColor(kMagenta);
                        fTail->SetLineStyle(2);
                        fTail->DrawCopy("same");

                        f2->DrawCopy("same");
                    }
                }
            };

            c_res->cd(1);
            fitResidui(h_res_x0->DrawCopy(), doDoubleFit);
            c_res->cd(2);
            fitResidui(h_res_z0->DrawCopy(), doDoubleFit);
            c_res->cd(3);
            fitResidui(h_res_sx->DrawCopy(), doDoubleFit);
            c_res->cd(4);
            fitResidui(h_res_sz->DrawCopy(), doDoubleFit);

            TCanvas *c_pars = (TCanvas *)gROOT->FindObject("c_pars");
            if(!c_pars)
                c_pars = new TCanvas("c_pars", "Fitted Parameters Distribution", 1200, 800);
            else
                c_pars->Clear();
            c_pars->Divide(2, 2);

            c_pars->cd(1);
            h_fit_x0->SetLineColor(kRed + 1);
            h_fit_x0->DrawCopy();
            c_pars->cd(2);
            h_fit_z0->SetLineColor(kBlue + 1);
            h_fit_z0->DrawCopy();
            c_pars->cd(3);
            h_fit_sx->SetLineColor(kGreen + 2);
            h_fit_sx->DrawCopy();
            c_pars->cd(4);
            h_fit_sz->SetLineColor(kOrange + 1);
            h_fit_sz->DrawCopy();
        }
        else
        {
            TCanvas *c_res = (TCanvas *)gROOT->FindObject("c_res");
            if(!c_res)
                c_res = new TCanvas("c_res", "Fit Residuals MC", 1200, 800);
            else
                c_res->Clear();

            c_res->Divide(3, 2);
            c_res->cd(1);
            h_res_xc->Fit("gaus", "Q");
            h_res_xc->DrawCopy();
            c_res->cd(2);
            h_res_yc->Fit("gaus", "Q");
            h_res_yc->DrawCopy();
            c_res->cd(3);
            h_res_R->Fit("gaus", "Q");
            h_res_R->DrawCopy();
            c_res->cd(4);
            h_res_z->Fit("gaus", "Q");
            h_res_z->DrawCopy();
            c_res->cd(5);
            h_res_dz->Fit("gaus", "Q");
            h_res_dz->DrawCopy();
        }
    }
}

// ==============================================================================
// --- REAL DATA ANALYSIS ---
// ==============================================================================

void DrawEvent(
    Int_t eventID, Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    if(ROOT::IsImplicitMTEnabled())
        ROOT::DisableImplicitMT();

    SetGlobalStyle();
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    if(gSystem->AccessPathName(filename.c_str()))
    {
        cerr << "[Error] File not found: " << filename << endl;
        return;
    }

    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    if(eventID == -1)
    {
        gRandom->SetSeed(0);
        long long nEntries = reader.GetRaw().Count().GetValue();
        eventID = gRandom->Integer(nEntries);
    }

    reader.SetSingleEntry(eventID);
    auto df = reader.GetCHeTTree();

    // The reader provides "All_Bundle", which contains the global IDs of
    // selected hits
    auto hitsResult = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");

    if(hitsResult->empty())
    {
        cerr << "[Warning] Event " << eventID << " seems empty or out of range." << endl;
        return;
    }

    // Convert RVec<int> to vector<int>
    const auto &rvec = hitsResult->at(0);
    vector<int> hit_ids(rvec.begin(), rvec.end());

    cout << "Displaying Run " << runID << " - Event " << eventID << endl;
    cout << "Selected Hits : " << hit_ids.size() << endl;
    cout << "Intersections : " << Config::FindIntersections(hit_ids).size() << endl;

    // Call Visualizer
    Vis::Draw2D(hit_ids);
    Vis::Draw3D(hit_ids);
}

void PlotRawCorrelation(Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    if(!ROOT::IsImplicitMTEnabled())
        ROOT::EnableImplicitMT(16);

    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    Data::Reader reader(filename, "auto");
    // Explicitly set cuts even though we define our own filtered vars below
    // (Reader defines FDxx_corrToA)
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df = reader.GetCHeTTree();
    auto rawCount = reader.GetRaw().Count();

    auto df_all_hits
        = df.Define("All_ToA_raw",
                [](const RVecD &t0, const RVecD &t1, const RVecD &t2, const RVecD &t3)
                { return Concatenate(t0, Concatenate(t1, Concatenate(t2, t3))); },
                { "FD00_corrToA", "FD01_corrToA", "FD02_corrToA", "FD03_corrToA" })
              .Define("All_ToT_raw",
                  [](const RVecUS &t0, const RVecUS &t1, const RVecUS &t2, const RVecUS &t3)
                  { return RVecD(Concatenate(t0, Concatenate(t1, Concatenate(t2, t3)))); },
                  { "FD00_ToT", "FD01_ToT", "FD02_ToT", "FD03_ToT" })
              .Define("All_ToA",
                  [=](const RVecD &toa, const RVecD &tot) {
                      return toa[(toa > toaMin) && (toa < toaMax) && (tot > totMin)
                          && (tot < totMax)];
                  },
                  { "All_ToA_raw", "All_ToT_raw" })
              .Define("All_ToT",
                  [=](const RVecD &toa, const RVecD &tot) {
                      return tot[(toa > toaMin) && (toa < toaMax) && (tot > totMin)
                          && (tot < totMax)];
                  },
                  { "All_ToA_raw", "All_ToT_raw" });

    auto nPassed = df_all_hits.Filter("All_ToA.size() > 0").Count();
    auto h_raw_corr = df_all_hits.Histo2D(
        { "h_raw_corr", Form("Raw Correlation Run %d;ToA [ns];ToT [LSB]", runID),
            static_cast<int>((toaMax - toaMin) / 4), toaMin, toaMax,
            static_cast<int>((totMax - totMin) / 4), (double)totMin, (double)totMax },
        "All_ToA", "All_ToT");

    TCanvas *cRaw = (TCanvas *)gROOT->FindObject("cRaw");
    if(!cRaw)
        cRaw = new TCanvas("cRaw", "Correlation Analysis", 800, 600);
    else
        cRaw->Clear();

    cRaw->cd();
    cRaw->SetLogz();
    h_raw_corr->DrawCopy("COLZ");
    PrintSummary(rawCount.GetValue(), nPassed.GetValue(), toaMin, toaMax, totMin, totMax);
}

void PlotMultiplicity(Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();

    // --- Setup Dati ---
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df_full = reader.GetCHeTTree();

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
        { "h_total", "Total Multiplicity;n Hits;Events", 201, -0.5, 200.5 }, "nHits_Total");
    auto h_c0 = df.Histo1D(
        { "h_c0", "Cylinders Multiplicity;n Hits;Events", 121, -0.5, 120.5 }, "nHits_Cyl0");
    auto h_c1 = df.Histo1D({ "h_c1", "", 121, -0.5, 120.5 }, "nHits_Cyl1");
    auto h_l1 = df.Histo1D(
        { "h_l1", "Layers Multiplicity;n Hits;Events", 81, -0.5, 80.5 }, "nHits_Cyl0_Inner");
    auto h_l2 = df.Histo1D({ "h_l2", "", 81, -0.5, 80.5 }, "nHits_Cyl0_Outer");
    auto h_l3 = df.Histo1D({ "h_l3", "", 81, -0.5, 80.5 }, "nHits_Cyl1_Inner");
    auto h_l4 = df.Histo1D({ "h_l4", "", 81, -0.5, 80.5 }, "nHits_Cyl1_Outer");

    // Occupancy
    auto df_occ = df.Define("Bundles_Cyl0_Lay0", "All_Bundle[All_Cyl == 0 && All_Lay == 0]")
                      .Define("Bundles_Cyl0_Lay1", "All_Bundle[All_Cyl == 0 && All_Lay == 1]")
                      .Define("Bundles_Cyl1_Lay0", "All_Bundle[All_Cyl == 1 && All_Lay == 0]")
                      .Define("Bundles_Cyl1_Lay1", "All_Bundle[All_Cyl == 1 && All_Lay == 1]");

    auto h_occ_c0l0
        = df_occ.Histo1D({ "h_occ_c0l0", "Occupancy Cyl 0 - Inner;Global Bundle ID;Hits", nBins[0],
                             binEdges[0], binEdges[1] },
            "Bundles_Cyl0_Lay0");
    auto h_occ_c0l1
        = df_occ.Histo1D({ "h_occ_c0l1", "Occupancy Cyl 0 - Outer;Global Bundle ID;Hits", nBins[1],
                             binEdges[1], binEdges[2] },
            "Bundles_Cyl0_Lay1");
    auto h_occ_c1l0
        = df_occ.Histo1D({ "h_occ_c1l0", "Occupancy Cyl 1 - Inner;Global Bundle ID;Hits", nBins[2],
                             binEdges[2], binEdges[3] },
            "Bundles_Cyl1_Lay0");
    auto h_occ_c1l1
        = df_occ.Histo1D({ "h_occ_c1l1", "Occupancy Cyl 1 - Outer;Global Bundle ID;Hits", nBins[3],
                             binEdges[3], binEdges[4] },
            "Bundles_Cyl1_Lay1");

    // ------------------------------------------
    // 3. Disegno
    // ------------------------------------------

    // Canvas 1: Multiplicity (Stile Classico)
    TCanvas *cMult = (TCanvas *)gROOT->FindObject("cMult");
    if(!cMult)
        cMult = new TCanvas("cMult", Form("Multiplicity - Run %d", runID), 1200, 800);
    else
        cMult->Clear();
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
    TCanvas *cOcc = (TCanvas *)gROOT->FindObject("cOcc");
    if(!cOcc)
        cOcc = new TCanvas("cOcc", Form("Bundle Occupancy - Run %d", runID), 1200, 800);
    else
        cOcc->Clear();
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

void PlotEstimators(Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    long long nRaw = reader.GetRaw().Count().GetValue();
    auto df = reader.GetCHeTTree().Filter("SumToT > 0");
    auto nPassed = df.Count();

    auto h_toa = df.Histo1D(
        { "h_toa", "First ToA;First ToA [ns];Events", 150, toaMin, toaMax }, "FirstToA");
    auto h_tot = df.Histo1D({ "h_tot", "Sum of ToT;Sum ToT [LSB];Events", 100, 0, 5000 }, "SumToT");
    auto h_corr
        = df.Histo2D({ "h_corr", "Correlation First ToA vs Sum ToT;First ToA [ns];Sum ToT [LSB]",
                         100, toaMin, toaMax, 100, 0, 5000 },
            "FirstToA", "SumToT");

    TCanvas *cEst = (TCanvas *)gROOT->FindObject("cEst");
    if(!cEst)
        cEst = new TCanvas("cEst", Form("Estimators - Run %d", runID), 1200, 800);
    else
        cEst->Clear();
    cEst->Divide(2, 2);
    cEst->cd(1);
    h_toa->SetLineColor(kBlue + 1);
    h_toa->DrawCopy();
    cEst->cd(2);
    h_tot->SetLineColor(kRed + 1);
    h_tot->DrawCopy();
    cEst->cd(3);
    h_corr->DrawCopy("COLZ");

    PrintSummary(nRaw, nPassed.GetValue(), toaMin, toaMax, totMin, totMax);
}

void PlotMultiplicityCorrelations(
    Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    SetGlobalStyle();
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    long long nRaw = reader.GetRaw().Count().GetValue();
    auto df = reader.GetCHeTTree().Filter("nHits_Total > 0");
    auto nPassed = df.Count();

    auto h_hits_toa = df.Histo2D({ "h_hits_toa", "Multiplicity vs First ToA;n Hits;First ToA [ns]",
                                     101, -0.5, 100.5, 150, toaMin, toaMax },
        "nHits_Total", "FirstToA");
    auto h_hits_tot = df.Histo2D({ "h_hits_tot", "Multiplicity vs Sum ToT;n Hits;Sum ToT [LSB]",
                                     101, -0.5, 100.5, 200, 0, 10000 },
        "nHits_Total", "SumToT");

    TCanvas *cMultCorr = (TCanvas *)gROOT->FindObject("cMultCorr");
    if(!cMultCorr)
        cMultCorr = new TCanvas("cMultCorr", "Correlations", 1400, 600);
    else
        cMultCorr->Clear();
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

void FitCosmicEvent(
    Int_t eventID, Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax)
{
    if(ROOT::IsImplicitMTEnabled())
        ROOT::DisableImplicitMT();

    // --- 1. Preparazione Dati (Lettura evento singolo) ---
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";
    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);

    // Usa Range per processare SOLO l'evento richiesto (molto veloce)
    if(eventID == -1)
    {
        gRandom->SetSeed(0);
        eventID = gRandom->Integer(reader.GetRaw().Count().GetValue());
    }

    reader.SetSingleEntry(eventID);
    auto df = reader.GetCHeTTree();

    // Take dei dati
    auto hits_ptr = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");
    if(hits_ptr->empty() || hits_ptr->at(0).empty())
    {
        cout << "[WARNING] Evento " << eventID << " vuoto o non trovato." << endl;
        return;
    }
    const auto &rvec = hits_ptr->at(0);
    vector<int> hit_ids(rvec.begin(), rvec.end());

    // Rimuovi duplicati e ordina
    sort(hit_ids.begin(), hit_ids.end());
    hit_ids.erase(unique(hit_ids.begin(), hit_ids.end()), hit_ids.end());

    // --- 2. Esecuzione del Fit ---
    cout << "\n=== Analysis Run " << runID << " Event " << eventID << " ===" << endl;
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
            tr.sx, 1.0,
            tr.sz, // Vettore direzione (dx/dy, dy/dy, dz/dy)
            kBlack, 3, 1,
            true // Colore, spessore, stile, isLocal
        );

    // Disegna (Hit + Traccia Fit)
    Vis::Draw2D(hit_ids, visTracks);

    Vis::Draw3D(hit_ids, visTracks, fitOut.fittedPoints);
}

void AnalyzeCosmicRun(Int_t runID, Double_t toaMin, Double_t toaMax, UInt_t totMin, UInt_t totMax,
    Double_t minPValue = 0., Bool_t doEfficiency = true)
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
    TDirectory::AddDirectory(false);
    auto h_chi2
        = make_unique<TH1D>("h_chi2", "Normalized #chi^{2};#chi^{2} / ndf;Events", 100, 0, 20);
    auto h_prob = make_unique<TH1D>("h_prob", "Prob(#chi^{2});Prob;Events", 50, 0, 1);
    auto h_nhits = make_unique<TH1D>("h_nhits", "Hits per Track;N Hits;Events", 15, -0.5, 50.5);

    auto h_x0 = make_unique<TH1D>("h_x0", "Reconstructed X0;x_{0} [mm];Events", 80, -40, 40);
    auto h_z0 = make_unique<TH1D>("h_z0", "Reconstructed Z0;z_{0} [mm];Events", 100, -150, 150);
    auto h_sx = make_unique<TH1D>("h_sx", "Slope X (sx);sx = dx/dy;Events", 100, -1.2, 1.2);
    auto h_sz = make_unique<TH1D>("h_sz", "Slope Z (sz);sz = dz/dy;Events", 100, -2, 2);
    auto h_phi = make_unique<TH1D>(
        "h_phi", "Azimuthal Angle #phi; #phi [rad];Events", 100, -TMath::Pi(), TMath::Pi());
    auto h_cosTheta = make_unique<TH1D>(
        "h_cosTheta", "Polar Angle cos(#theta); cos(#theta); Events", 100, -1.0, 0.0);

    // --- 2. Estrazione Dati con RDataFrame ---
    string filename
        = "/home/lorenzo/muEDM_Project/Data/RootData/run00" + to_string(runID) + ".root";

    // Geometry test
    // Config::SetOffsetExp(29. * TMath::DegToRad()); // 29 optimal
    // Config::SetDelta2(
    //     Config::GetDelta2() + 4.5 * TMath::DegToRad()); // 4.5 optimal
    // Config::SetDelta1(Config::GetDelta1() - 18.0 * (M_PI / 180.0));

    // Data reader
    Data::Reader reader(filename, "auto");
    reader.SetCuts(toaMin, toaMax, totMin, totMax);
    auto df = reader.GetCHeTTree();

    auto all_hits_ptr = df.Take<RVecI>("All_Bundle");
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
        if(hits.size() < 4)
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

            if(DEBUG_RUN)
            {
                // Bool_t condition = prob > 0 && prob < 0.01;
                Bool_t condition2 = out.track.x0 >= -17 && out.track.x0 <= -16;
                InspectEvent(condition2, hits, out.track);
            }

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

                AccumulateEfficiency(out.track, hits, bundleStats, layerStats, cylStats);
            }
        }
    }
    printf("\nAnalysis Complete. Fitted %d / %d events (%.1f%%)\n", nFitted, nEvents,
        100.0 * nFitted / nEvents);

    // --- 4. Plotting ---
    TCanvas *cAn1 = (TCanvas *)gROOT->FindObject("cAn1");
    if(!cAn1)
        cAn1 = new TCanvas("cAn1", "Run Analysis Results - GOF", 800, 1200);
    else
        cAn1->Clear();
    cAn1->Divide(1, 3);

    cAn1->cd(1);
    gPad->SetLogy();
    h_chi2->SetLineColor(kBlue + 1);
    h_chi2->SetLineWidth(2);
    h_chi2->DrawCopy();

    cAn1->cd(2);
    h_prob->SetLineColor(kMagenta + 1);
    h_prob->SetLineWidth(2);
    h_prob->SetMinimum(0);
    h_prob->DrawCopy();

    cAn1->cd(3);
    h_nhits->SetFillColor(kGray);
    h_nhits->DrawCopy();

    TCanvas *cAn2 = (TCanvas *)gROOT->FindObject("cAn2");
    if(!cAn2)
        cAn2 = new TCanvas("cAn2", "Run Analysis Results - Pars", 1200, 800);
    else
        cAn2->Clear();
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
        PlotEfficiencyResults(runID, bundleStats, layerStats, cylStats, occupancyMap);
}
