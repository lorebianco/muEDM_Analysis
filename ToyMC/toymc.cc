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
    // 1. SEV Tree (Detector Data - High Level)
    // =========================================================================
    TTree *treeSEV = new TTree("Event", "Single Event Detector Data");

    int eventID = 0;
    std::vector<int> all_bundle;

    treeSEV->Branch("EventID", &eventID);
    treeSEV->Branch("All_Bundle", &all_bundle);

    // =========================================================================
    // 2. SIM Tree (MonteCarlo Truth Data)
    // =========================================================================
    TTree *treeSIM = new TTree("sim", "MonteCarlo Truth Data");

    treeSIM->Branch("EventID", &eventID);

    int trackID = 0;
    int particleID = 0; // e.g., 11 for electron, 0 for cosmic
    treeSIM->Branch("TrackID", &trackID);
    treeSIM->Branch("ParticleID", &particleID);

    // Common Truth Variables (can be padded with 0 if not used in a specific mode)
    double mc_x0 = 0, mc_y0 = 0, mc_z0 = 0;
    double mc_ux = 0, mc_uy = 0, mc_uz = 0;
    double mc_E = 0, mc_R = 0, mc_cx = 0, mc_cy = 0, mc_tmin = 0, mc_tmax = 0;

    if(mode == "cosmic")
    {
        particleID = 0; // Cosmic muon placeholder
        treeSIM->Branch("mc_x0", &mc_x0);
        treeSIM->Branch("mc_y0", &mc_y0);
        treeSIM->Branch("mc_z0", &mc_z0);
        treeSIM->Branch("mc_ux", &mc_ux);
        treeSIM->Branch("mc_uy", &mc_uy);
        treeSIM->Branch("mc_uz", &mc_uz);
    }
    else if(mode == "michel")
    {
        particleID = 11; // Michel Electron
        treeSIM->Branch("mc_E", &mc_E);
        treeSIM->Branch("mc_R", &mc_R);
        treeSIM->Branch("mc_cx", &mc_cx);
        treeSIM->Branch("mc_cy", &mc_cy);
        treeSIM->Branch("mc_z0", &mc_z0);
        treeSIM->Branch("mc_uz", &mc_uz);
        treeSIM->Branch("mc_tmin", &mc_tmin);
        treeSIM->Branch("mc_tmax", &mc_tmax);
        // Let's add x0, y0 emission point just in case
        treeSIM->Branch("mc_x0", &mc_x0);
        treeSIM->Branch("mc_y0", &mc_y0);
    }

    std::cout << "Generating events...\n";
    for(eventID = 0; eventID < nEvents; ++eventID)
    {
        all_bundle.clear();
        trackID = eventID; // Assuming 1 track per event for now

        if(mode == "cosmic")
        {
            ToyMC::CosmicTrack tr = ToyMC::GenerateCosmic();
            mc_x0 = tr.x0;
            mc_y0 = tr.y0;
            mc_z0 = tr.z0;
            mc_ux = tr.ux;
            mc_uy = tr.uy;
            mc_uz = tr.uz;
            all_bundle = ToyMC::FindCosmicHits(tr, efficiency);
        }
        else if(mode == "michel")
        {
            ToyMC::MichelTrack tr = ToyMC::GenerateMichelTrack(false);
            mc_E = tr.E_kin;
            mc_R = tr.radius;
            mc_cx = tr.cx;
            mc_cy = tr.cy;
            mc_z0 = tr.z0;
            mc_uz = tr.dz_dt;
            mc_tmin = tr.t_min;
            mc_tmax = tr.t_max;
            mc_x0 = tr.x0;
            mc_y0 = tr.y0;
            all_bundle = ToyMC::FindMichelHits(tr, efficiency);
        }

        // Save parallel entries to both trees
        treeSEV->Fill();
        treeSIM->Fill();

        if((eventID + 1) % 1000 == 0)
        {
            std::cout << "Processed " << (eventID + 1) << " / " << nEvents << " events\r"
                      << std::flush;
        }
    }
    std::cout << "\nGeneration complete. Writing to disk...\n";

    fOut->cd();
    treeSEV->Write();
    treeSIM->Write();
    fOut->Close();

    std::cout << "Saved successfully to " << outFile << "\n";

    return 0;
}
