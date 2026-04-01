#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "ToyGenerator.hh"

void PrintUsage(const char *progName)
{
    std::cout << "Usage: " << progName << " <mode> <n_events> <output_file> [efficiency]\n"
              << "  mode:        cosmic | michel\n"
              << "  n_events:    number of events to generate\n"
              << "  output_file: output ROOT file (e.g., toy_data.root)\n"
              << "  efficiency:  (optional) detector efficiency [0.0 - 1.0], default 1.0\n";
}

int main(int argc, char **argv)
{
    if(argc < 4)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    int nEvents = std::stoi(argv[2]);
    std::string outFile = "../../data/input/toy/" + std::string(argv[3]);
    double efficiency = 1.0;
    if(argc >= 5)
    {
        efficiency = std::stod(argv[4]);
    }

    if(mode != "cosmic" && mode != "michel")
    {
        std::cerr << "Error: Unknown mode '" << mode << "'. Must be 'cosmic' or 'michel'.\n";
        return 1;
    }

    std::cout << "--- muEDM ToyMC Generator ---\n"
              << "Mode:       " << mode << "\n"
              << "Events:     " << nEvents << "\n"
              << "Output:     " << outFile << "\n"
              << "Efficiency: " << efficiency << "\n"
              << "-----------------------------\n";

    // Create the output file
    TFile *fOut = TFile::Open(outFile.c_str(), "RECREATE");
    if(!fOut || fOut->IsZombie())
    {
        std::cerr << "Error: Could not open output file " << outFile << "\n";
        return 1;
    }

    // =========================================================================
    // 1. CHeT Tree (Detector Data - High Level)
    // =========================================================================
    TTree *treeCHeT = new TTree("chet", "Single Event Detector Data");

    int eventID = 0;
    std::vector<int> all_bundle;

    treeCHeT->Branch("EventID", &eventID);
    treeCHeT->Branch("All_Bundle", &all_bundle);

    // =========================================================================
    // 2. SIM Tree (MonteCarlo Truth Data)
    // =========================================================================
    TTree *treeSIM = new TTree("sim", "MonteCarlo Truth Data");

    treeSIM->Branch("EventID", &eventID);

    std::vector<double> true_hit_x, true_hit_y, true_hit_z;
    treeSIM->Branch("mc_hits_x", &true_hit_x);
    treeSIM->Branch("mc_hits_y", &true_hit_y);
    treeSIM->Branch("mc_hits_z", &true_hit_z);

    int trackID = 0;
    int particleID = 0; // e.g., 11 for electron, 0 for cosmic
    treeSIM->Branch("TrackID", &trackID);
    treeSIM->Branch("ParticleID", &particleID);

    // Common Truth Variables (can be padded with 0 if not used in a specific mode)
    // Track Variables (Mathematical parametrization for fit comparison)
    double trk_x0 = 0, trk_y0 = 0, trk_z0 = 0;
    double trk_ux = 0, trk_uy = 0, trk_uz = 0;
    double trk_R = 0, trk_cx = 0, trk_cy = 0, trk_tmin = 0, trk_tmax = 0;

    // MC Truth Variables (Physical generation points and sanity checks)
    double mc_x = 0, mc_y = 0, mc_z = 0;
    double mc_E = 0;

    treeSIM->Branch("mc_x", &mc_x);
    treeSIM->Branch("mc_y", &mc_y);
    treeSIM->Branch("mc_z", &mc_z);

    if(mode == "cosmic")
    {
        particleID = 0; // Cosmic muon placeholder
        treeSIM->Branch("trk_x0", &trk_x0);
        treeSIM->Branch("trk_y0", &trk_y0);
        treeSIM->Branch("trk_z0", &trk_z0);
        treeSIM->Branch("trk_ux", &trk_ux);
        treeSIM->Branch("trk_uy", &trk_uy);
        treeSIM->Branch("trk_uz", &trk_uz);
    }
    else if(mode == "michel")
    {
        particleID = 11; // Michel Electron
        treeSIM->Branch("mc_E", &mc_E);
        treeSIM->Branch("trk_R", &trk_R);
        treeSIM->Branch("trk_cx", &trk_cx);
        treeSIM->Branch("trk_cy", &trk_cy);
        treeSIM->Branch("trk_z0", &trk_z0);
        treeSIM->Branch("trk_uz", &trk_uz);
        treeSIM->Branch("trk_tmin", &trk_tmin);
        treeSIM->Branch("trk_tmax", &trk_tmax);
    }

    std::cout << "Generating events...\n";
    for(eventID = 0; eventID < nEvents; ++eventID)
    {
        all_bundle.clear();
        true_hit_x.clear();
        true_hit_y.clear();
        true_hit_z.clear();
        trackID = eventID; // Assuming 1 track per event for now

        if(mode == "cosmic")
        {
            ToyMC::CosmicTrack tr = ToyMC::GenerateCosmic();
            mc_x = tr.x0;
            mc_y = tr.y0;
            mc_z = tr.z0;

            trk_x0 = tr.x0;
            trk_y0 = tr.y0;
            trk_z0 = tr.z0;
            trk_ux = tr.ux;
            trk_uy = tr.uy;
            trk_uz = tr.uz;
            ToyMC::HitResult res = ToyMC::FindCosmicHits(tr, efficiency);
            all_bundle = res.bundles;
            true_hit_x = res.x;
            true_hit_y = res.y;
            true_hit_z = res.z;
        }
        else if(mode == "michel")
        {
            ToyMC::MichelTrack tr = ToyMC::GenerateMichelTrack(false);
            mc_x = tr.x0;
            mc_y = tr.y0;
            mc_z = 0.0;
            mc_E = tr.E_kin;

            trk_R = tr.radius;
            trk_cx = tr.cx;
            trk_cy = tr.cy;
            trk_z0 = tr.z0;
            trk_uz = tr.dz_dt;
            trk_tmin = tr.t_min;
            trk_tmax = tr.t_max;
            ToyMC::HitResult res = ToyMC::FindMichelHits(tr, efficiency);
            all_bundle = res.bundles;
            true_hit_x = res.x;
            true_hit_y = res.y;
            true_hit_z = res.z;
        }

        // Save parallel entries to both trees
        treeCHeT->Fill();
        treeSIM->Fill();

        if((eventID + 1) % 1000 == 0)
        {
            std::cout << "Processed " << (eventID + 1) << " / " << nEvents << " events\r"
                      << std::flush;
        }
    }
    std::cout << "\nGeneration complete. Writing to disk...\n";

    fOut->cd();
    treeCHeT->Write();
    treeSIM->Write();
    fOut->Close();

    std::cout << "Saved successfully to " << outFile << "\n";

    return 0;
}
