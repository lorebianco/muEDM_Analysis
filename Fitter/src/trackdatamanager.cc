#include "trackdatamanager.hh"

using namespace std;

TrackDataManager::TrackDataManager()
{
    // >>>>> Input <<<<< //
    simChain = new TChain("sim");

    for(const auto &fileName : Config::get().inputDataFiles)
    {
        cout << ">>> Adding file: " << fileName << endl;
        simChain->Add(fileName.c_str());
    }

    // Set branch addresses
    simChain->SetBranchAddress("EventID", &EventID);

    simChain->SetBranchAddress("mc_x", &trueDecay_posX);
    simChain->SetBranchAddress("mc_y", &trueDecay_posY);
    simChain->SetBranchAddress("mc_z", &trueDecay_posZ);
    simChain->SetBranchAddress("mc_px", &trueDecay_momX);
    simChain->SetBranchAddress("mc_py", &trueDecay_momY);
    simChain->SetBranchAddress("mc_pz", &trueDecay_momZ);

    if(simChain->GetBranch("mc_E"))
    {
        simChain->SetBranchAddress("mc_E", &mc_E);
        simChain->SetBranchAddress("trk_R", &trk_R);
        simChain->SetBranchAddress("trk_cx", &trk_cx);
        simChain->SetBranchAddress("trk_cy", &trk_cy);
        simChain->SetBranchAddress("trk_z0", &trk_z0);
        simChain->SetBranchAddress("trk_uz", &trk_uz);
        simChain->SetBranchAddress("trk_tmin", &trk_tmin);
        simChain->SetBranchAddress("trk_tmax", &trk_tmax);
    }

    reader = std::make_unique<CHeT::Data::Reader>(
        Config::get().inputDataFiles[0], Config::get().inputTreeName);

    if(simChain->GetEntries() == 0)
    {
        reader->SetCuts(Config::get().cutToAMin, Config::get().cutToAMax, Config::get().cutToTMin,
            Config::get().cutToTMax);
    }

    auto df = reader->GetCHeTTree();
    auto hits_ptr = df.Take<ROOT::VecOps::RVec<int>>("All_Bundle");
    allEventsHits = *hits_ptr;

    // Set nEvents and processedEvents
    nEvents = simChain->GetEntries();
    if(nEvents == 0)
        nEvents = allEventsHits.size();

    if(Config::get().processSingle)
        processedEvents = 1;
    else
        processedEvents = (Config::get().rangeLoop.second - Config::get().rangeLoop.first) < nEvents
            ? (Config::get().rangeLoop.second - Config::get().rangeLoop.first)
            : nEvents;

    // Summary
    cout << "\n>>> There are " << nEvents << " events in the chain" << endl;

    if(Config::get().processSingle)
        cout << ">>> "
             << "1"
             << " event will be processed\n"
             << endl;
    else
        cout << ">>> " << processedEvents << " events will be processed\n\n";

    // >>>>> Output <<<<< //

    // Initialize output file and trees
    if(!Config::get().outputFile.empty())
    {
        outputFile = new TFile(Config::get().outputFile.c_str(), "RECREATE");

        // Try to clone input trees if they exist
        TTree *simTree = (TTree *)simChain->GetTree();
        if(simTree)
        {
            outputFile->cd();
            TTree *clonedSim = simTree->CloneTree();
            clonedSim->Write();
        }
    }
}

void TrackDataManager::InitRecoTree(bool isMichel)
{
    if(!outputFile)
        return;
    outputFile->cd();
    if(!recTree)
        recTree = new TTree("rec", "Reconstructed Data");

    recTree->Branch("EventID", &EventID, "EventID/I");
    recTree->Branch("rec_chi2", &rec_chi2, "rec_chi2/D");
    recTree->Branch("rec_acceptance", &rec_acceptance, "rec_acceptance/O");
    recTree->Branch("rec_converged", &rec_converged, "rec_converged/O");
    recTree->Branch("rec_hits", &rec_hits);
    recTree->Branch("rec_hough2d_idx", &rec_hough2d_idx);
    recTree->Branch("rec_houghz_idx", &rec_houghz_idx);

    if(!isMichel)
    {
        recTree->Branch("rec_x0", &rec_x0, "rec_x0/D");
        recTree->Branch("rec_z0", &rec_z0, "rec_z0/D");
        recTree->Branch("rec_sx", &rec_sx, "rec_sx/D");
        recTree->Branch("rec_sz", &rec_sz, "rec_sz/D");
    }
    else
    {
        recTree->Branch("rec_R", &rec_R, "rec_R/D");
        recTree->Branch("rec_cx", &rec_cx, "rec_cx/D");
        recTree->Branch("rec_cy", &rec_cy, "rec_cy/D");
        recTree->Branch("rec_z0", &rec_z0, "rec_z0/D");
        recTree->Branch("rec_dz_ds", &rec_dz_ds, "rec_dz_ds/D");
        recTree->Branch("rec_phi0", &rec_phi0, "rec_phi0/D");
        recTree->Branch("rec_t_min", &rec_t_min, "rec_t_min/D");
        recTree->Branch("rec_t_max", &rec_t_max, "rec_t_max/D");

        recTree->Branch("rec_extrap_x", &rec_extrap_x, "rec_extrap_x/D");
        recTree->Branch("rec_extrap_y", &rec_extrap_y, "rec_extrap_y/D");
        recTree->Branch("rec_extrap_z", &rec_extrap_z, "rec_extrap_z/D");
        recTree->Branch("rec_extrap_px", &rec_extrap_px, "rec_extrap_px/D");
        recTree->Branch("rec_extrap_py", &rec_extrap_py, "rec_extrap_py/D");
        recTree->Branch("rec_extrap_pz", &rec_extrap_pz, "rec_extrap_pz/D");

        recTree->Branch("rec_n_candidates_2d", &rec_n_candidates_2d, "rec_n_candidates_2d/I");
        recTree->Branch("rec_n_candidates_z", &rec_n_candidates_z, "rec_n_candidates_z/I");
    }
}

void TrackDataManager::SaveRecoTree()
{
    if(outputFile)
    {
        outputFile->cd();
        if(recTree)
        {
            recTree->Write();
        }
        outputFile->Close();
        delete outputFile;
        outputFile = nullptr;
    }
}

TrackDataManager::~TrackDataManager()
{
    SaveRecoTree();
    delete simChain;
}

void TrackDataManager::ProcessAndFilterEvent(Long64_t eventID)
{
    if(eventID >= 0)
        currentEventIndex = eventID;

    this->EventID = currentEventIndex;

    // Clean containers
    hitsCoordinates.clear();
    hitsCylinderID.clear();
    rec_hits.clear();
    rec_hough2d_idx.clear();
    rec_houghz_idx.clear();

    if(currentEventIndex >= (Long64_t)allEventsHits.size())
        return;

    const auto &rvec = allEventsHits[currentEventIndex];
    std::vector<int> hit_ids(rvec.begin(), rvec.end());

    std::sort(hit_ids.begin(), hit_ids.end());
    hit_ids.erase(std::unique(hit_ids.begin(), hit_ids.end()), hit_ids.end());

    rec_hits = hit_ids;

    SetTrueDecayData(currentEventIndex);

    // SAVE DATA
    // Calculate all possible 3D intersections from the fired bundles.
    // Spacepoint and Helix fitters will use these true 3D points.
    const Double_t mm2cm = 0.1;
    auto intersections = CHeT::Config::FindIntersections(hit_ids);

    for(const auto &inter : intersections)
    {
        std::vector<Double_t> point = { inter.x * mm2cm, inter.y * mm2cm, inter.z * mm2cm };
        hitsCoordinates.push_back(point);
        hitsCylinderID.push_back(inter.cylinderId);
    }

    currentEventIndex++;
}

Bool_t TrackDataManager::SetTrueDecayData(Int_t eventID)
{
    if(simChain->GetEntries() == 0 || eventID >= simChain->GetEntries())
        return false;

    simChain->GetEntry(eventID);

    // Set time of decay
    trueDecayTime = 0; // Not available in simple sim

    // Conversion mm to cm
    const Double_t mm2cm = 0.1;

    // Set TVector3 for decay pos and mom
    trueDecayPos.SetXYZ(trueDecay_posX * mm2cm, trueDecay_posY * mm2cm, trueDecay_posZ * mm2cm);

    trueDecayMom.SetXYZ(trueDecay_momX, trueDecay_momY, trueDecay_momZ);

    // Compute and set the EDM theta angle
    Double_t x = trueDecayPos.X();
    Double_t y = trueDecayPos.Y();
    // Protezione per divisione per zero se x=0 e y=0
    Double_t r = sqrt(x * x + y * y);

    if(r > 1e-9)
    {
        Double_t pz = cos(trueDecayMom.Theta());
        Double_t pr = sin(trueDecayMom.Theta())
            * (x * cos(trueDecayMom.Phi()) + y * sin(trueDecayMom.Phi())) / r;
        trueDecayThetaEDM = TMath::ATan2(pz, pr);
    }
    else
    {
        trueDecayThetaEDM = 0; // O un valore di default sensato
    }

    return true; // All good!
}
