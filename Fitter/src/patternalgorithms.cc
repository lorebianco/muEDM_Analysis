#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMultiGraph.h>

#include "CHeT/CHeTGlobalSettings.hh"
#include "patternalgorithms.hh"

using namespace std;

pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SortVectorZ(
    vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Sort vectors according to z(third) value

    vector<Double_t> z_of_vectors;
    for(auto v : vectors)
        z_of_vectors.push_back(abs(v.at(2)));

    // initialize original index locations
    vector<size_t> idx(z_of_vectors.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using stable_sort instead of sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
        [&z_of_vectors](size_t i1, size_t i2) { return z_of_vectors[i1] < z_of_vectors[i2]; });

    vector<vector<Double_t>> copy_vector;
    for(auto i : idx)
        copy_vector.push_back(vectors.at(i));

    vector<Int_t> copy_cylinders;
    for(auto i : idx)
        copy_cylinders.push_back(cylinders.at(i));

    return { copy_vector, copy_cylinders };
}

pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::ShuffleVectorZ(
    vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Shuffle the vectors and cylinders randomly

    // Create a random engine and a distribution
    random_device rd;
    mt19937 g(rd());

    // Combine vectors and cylinders into a single vector of pairs
    vector<pair<vector<Double_t>, Int_t>> combined;
    for(size_t i = 0; i < vectors.size(); ++i)
    {
        combined.push_back({ vectors[i], cylinders[i] });
    }

    // Shuffle the combined vector randomly
    shuffle(combined.begin(), combined.end(), g);

    // Extract the shuffled vectors and cylinders back
    vector<vector<Double_t>> shuffled_vectors;
    vector<Int_t> shuffled_cylinders;
    for(const auto &p : combined)
    {
        shuffled_vectors.push_back(p.first);
        shuffled_cylinders.push_back(p.second);
    }

    return { shuffled_vectors, shuffled_cylinders };
}

Int_t PTTALG::CountTurns(const vector<vector<Double_t>> hitsCoordinates)
{
    if(hitsCoordinates.size() < 3)
        return 0; // Servono almeno 3 punti per trovare un massimo o minimo

    Int_t nTurns = 0;
    Bool_t foundMax = false, foundMin = false;

    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se abbiamo sia un massimo che un minimo -> un giro completato
        if(foundMax && foundMin)
        {
            nTurns++;
            foundMax = false;
            foundMin = false;
        }
    }

    return nTurns;
}

pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SelectTurn(
    Float_t turnID, const vector<vector<Double_t>> &hitsCoordinates, const vector<Int_t> &cylinders)
{
    vector<vector<Double_t>> turnHits;
    vector<Int_t> turnCylinders;

    if(hitsCoordinates.size() < 3)
        return { turnHits, turnCylinders }; // Troppi pochi punti per definire un giro

    Int_t totalTurns = CountTurns(hitsCoordinates);
    if(turnID > totalTurns)
        return { hitsCoordinates,
            cylinders }; // Se voglio più giri di quelli presenti, prendo tutta la traccia

    Double_t nHalfTurns = 0.0;
    turnHits.push_back(hitsCoordinates.front()); // Includi il primo punto
    turnCylinders.push_back(cylinders.front());

    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Identificazione di un massimo o minimo locale
        if((yCurr > yPrev && yCurr > yNext) || (yCurr < yPrev && yCurr < yNext))
        {
            nHalfTurns += 1.0; // Ora conto direttamente i mezzi giri
        }

        // Se il numero di **giri completi** supera `turnID`, interrompo
        if(nHalfTurns / 2.0 >= turnID)
        {
            break;
        }

        turnHits.push_back(hitsCoordinates[i]);
        turnCylinders.push_back(cylinders[i]);
    }

    // Includi sempre l'ultimo punto del semigiro
    turnHits.push_back(hitsCoordinates[turnHits.size()]);
    turnCylinders.push_back(cylinders[turnCylinders.size()]);

    return { turnHits, turnCylinders };
}

vector<Int_t> PTTALG::SplitTurns(const vector<vector<Double_t>> &hitsCoordinates)
{
    vector<Int_t> turnIndices;
    if(hitsCoordinates.size() < 3)
        return turnIndices;

    turnIndices.push_back(0); // Il primo indice è sempre 0

    Bool_t foundMax = false, foundMin = false;

    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se troviamo un massimo e poi un minimo (o viceversa), aggiungiamo l'indice
        if(foundMax && foundMin)
        {
            turnIndices.push_back(i);
            foundMax = false;
            foundMin = false;
        }
    }

    return turnIndices;
}

Int_t PTTALG::CountCylinders(const vector<Int_t> &cylinders)
{
    set<Int_t> uniqueValues(cylinders.begin(), cylinders.end());
    return uniqueValues.size();
}

pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SelectCylinders(vector<Int_t> targetCylinders,
    const vector<vector<Double_t>> &hitsCoordinates, const vector<Int_t> &cylinders)
{
    vector<vector<Double_t>> selectedHits;
    vector<Int_t> selectedCylinders;

    if(hitsCoordinates.size() != cylinders.size())
        return { selectedHits, selectedCylinders };

    for(size_t i = 0; i < hitsCoordinates.size(); i++)
    {
        Int_t currentID = cylinders[i];
        Bool_t isTarget = false;

        for(size_t j = 0; j < targetCylinders.size(); j++)
        {
            if(targetCylinders[j] == currentID)
            {
                isTarget = true;
                break;
            }
        }

        if(isTarget)
        {
            selectedHits.push_back(hitsCoordinates[i]);
            selectedCylinders.push_back(cylinders[i]);
        }
    }

    return { selectedHits, selectedCylinders };
}

pair<TVector3, TVector3> PTTALG::SmearSeed(TVector3 pos, TVector3 mom)
{
    TVector3 smearPos(pos);
    TVector3 smearMom(mom);

    smearPos.SetX(gRandom->Gaus(pos.X(), ANS::sigmaSeedPos));
    smearPos.SetY(gRandom->Gaus(pos.Y(), ANS::sigmaSeedPos));
    smearPos.SetZ(gRandom->Gaus(pos.Z(), ANS::sigmaSeedPos));

    smearMom.SetPhi(gRandom->Gaus(mom.Phi(), ANS::sigmaSeedMomDir));
    smearMom.SetTheta(gRandom->Gaus(mom.Theta(), ANS::sigmaSeedMomDir));
    smearMom.SetMag(gRandom->Gaus(mom.Mag(), ANS::sigmaSeedMomMag * mom.Mag()));

    return { smearPos, smearMom };
}

vector<CircularHoughResult> PTTALG::DoCircularHoughTransform(const vector<Int_t> &hit_ids,
    int nCandidates, int combinatorial_threshold, int n_random_triplets, bool drawGraphs)
{
    auto inters = CHeT::Config::FindIntersections(hit_ids);
    if(inters.size() < 3)
        return {};

    double maxDim = 150.0;
    auto h_center = make_unique<TH2D>("h_circ_centers",
        "Hough Circle Centers (XY);X_c [mm];Y_c [mm]", 150, -maxDim, maxDim, 150, -maxDim, maxDim);
    h_center->SetDirectory(nullptr);

    int n = inters.size();

    if(n < combinatorial_threshold)
    {
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
        if(drawGraphs)
            cout << "[Hough Info] Using random strategy (" << n
                 << " hits >= " << combinatorial_threshold << ")" << endl;
        gRandom->SetSeed(0);
        for(int i = 0; i < n_random_triplets; ++i)
        {
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

    vector<CircularHoughResult> candidates;
    TH2D *h_work = (TH2D *)h_center->Clone("h_work");
    int win = 2;
    for(int c = 0; c < nCandidates; ++c)
    {
        int bx, by, bz;
        h_work->GetMaximumBin(bx, by, bz);

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
            std::vector<int> track_hits_idx;
            for(size_t i = 0; i < inters.size(); ++i)
            {
                auto &p = inters[i];
                double r = sqrt(pow(p.x_loc - xc_peak, 2) + pow(p.y_loc - yc_peak, 2));
                if(abs(r - R_est) <= bin_width)
                {
                    sumR += r;
                    countR++;
                    track_hits_idx.push_back(i);
                }
            }
            if(countR > 0)
                R_est = sumR / countR;

            candidates.push_back({ xc_peak, yc_peak, R_est, sumW, track_hits_idx });

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

vector<ZHoughResult> PTTALG::DoZHoughTransform(const vector<Int_t> &hit_ids, double xc, double yc,
    double R_reco, int nCandidates, bool drawGraphs, double tol_Z_fit)
{
    struct Point3DLocal
    {
        double x, y, z, phi_raw;
        int orig_idx;
    };
    vector<Point3DLocal> valid_points;
    double tolR = 5.;
    auto inters = CHeT::Config::FindIntersections(hit_ids);
    for(size_t i = 0; i < inters.size(); ++i)
    {
        auto &p = inters[i];
        double r_point = sqrt(pow(p.x_loc - xc, 2) + pow(p.y_loc - yc, 2));
        if(abs(r_point - R_reco) <= tolR)
            valid_points.push_back(
                { p.x_loc, p.y_loc, p.z_loc, atan2(p.y_loc - yc, p.x_loc - xc), (int)i });
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

            vector<Point3D> track_hits;
            std::vector<int> track_hits_idx;
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
                    track_hits_idx.push_back(p.orig_idx);
                    if(best_unrolled_t < min_t)
                        min_t = best_unrolled_t;
                    if(best_unrolled_t > max_t)
                        max_t = best_unrolled_t;
                }
            }

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

            candidates.push_back(
                { z0_peak, dz_ds_peak, sumW, track_hits, track_hits_idx, min_t, max_t });
        }
    }
    delete h_work;

    if(drawGraphs)
    {
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

        TCanvas *c_phi_z_raw = (TCanvas *)gROOT->FindObject("c_phi_z_raw");
        if(!c_phi_z_raw)
            c_phi_z_raw = new TCanvas("c_phi_z_raw", "Raw Local Z-Phi", 800, 500);
        else
            c_phi_z_raw->Clear();

        c_phi_z_raw->cd();
        c_phi_z_raw->SetGrid();

        auto g_raw = new TGraph();
        g_raw->SetMarkerStyle(20);
        g_raw->SetMarkerSize(1.2);
        g_raw->SetMarkerColor(kAzure + 1);

        auto g_tiled = new TGraph();
        g_tiled->SetMarkerStyle(24);
        g_tiled->SetMarkerSize(1.2);
        g_tiled->SetMarkerColor(kAzure - 3);

        int pt_raw = 0;
        int pt_tiled = 0;
        for(auto &p : valid_points)
        {
            g_raw->SetPoint(pt_raw++, p.z, p.phi_raw);
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

            mg->GetXaxis()->SetLimits(-CHeT::Config::L_HALF, CHeT::Config::L_HALF);
            mg->GetYaxis()->SetRangeUser(-M_PI * 1.5, M_PI * 1.5);

            for(size_t c = 0; c < candidates.size(); ++c)
            {
                auto &cand = candidates[c];

                for(int k = -1; k <= 1; ++k)
                {
                    double shift_phi = 2.0 * M_PI * k;

                    double phi_start = cand.t_min + shift_phi;
                    double phi_end = cand.t_max + shift_phi;

                    double z_start = cand.z0 + cand.dz_ds * R_reco * cand.t_min;
                    double z_end = cand.z0 + cand.dz_ds * R_reco * cand.t_max;

                    TLine *line = new TLine(z_start, phi_start, z_end, phi_end);

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
