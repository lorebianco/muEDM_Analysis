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
    simChain->SetBranchAddress("EventID", &truedecay_EventID);
    simChain->SetBranchAddress("trk_x0", &trueDecay_posX);
    simChain->SetBranchAddress("trk_y0", &trueDecay_posY);
    simChain->SetBranchAddress("trk_z0", &trueDecay_posZ);
    simChain->SetBranchAddress("trk_ux", &trueDecay_momX);
    simChain->SetBranchAddress("trk_uy", &trueDecay_momY);
    simChain->SetBranchAddress("trk_uz", &trueDecay_momZ);

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

    // Init histograms
    accPhi = new TEfficiency("accPhi", "Acceptance: Phi vs Momentum; Momentum [MeV/c];#phi [rad]",
        10, 0, 68.9, 10, -TMath::Pi(), TMath::Pi());
    accTheta = new TEfficiency("accTheta",
        "Acceptance: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20,
        -TMath::Pi(), TMath::Pi());
    effPhi = new TEfficiency("effPhi", "Efficiency: Phi vs Momentum; Momentum [MeV/c];#phi [rad]",
        10, 0, 68.9, 10, -TMath::Pi(), TMath::Pi());
    effTheta = new TEfficiency("effTheta",
        "Efficiency: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20,
        -TMath::Pi(), TMath::Pi());

    histTurns = new TH1I("histTurns", "Number of Turns;nTurns;Counts", 20, 0, 20);
    effTurns = new TEfficiency("effTurns", "nTurns Efficiency;nTurns;Efficiency", 10, 0, 10);
    histCylinders = new TH1I("histCylinders", "Number of Cylinders;nCylinders;Counts", 7, 0, 7);
    effCylinders
        = new TEfficiency("effCylinders", "nCylinders Efficiency;nCylinders;Efficiency", 7, 0, 7);
    histCylVMom = new TProfile("histCylVMom",
        "Number of Cylinders vs Momentum;Momentum [MeV/c];nCylinders", 20, 0, 68.9, 0, 7);
    histTurnsVMom = new TProfile(
        "histTurnsVMom", "Number of Turns vs Momentum;Momentum [MeV/c];nTurns", 20, 0, 68.9, 0, 20);

    graphMom
        = new TH2D("graphMom", "Linearity: Momentum; Momentum_{MC} [MeV/c]; Momentum_{fit} [MeV/c]",
            100, 0., 70., 100, 0., 100.);
    graphTheta = new TH2D("graphTheta", "Linearity: Theta; #theta_{MC} [rad]; #theta_{fit} [rad]",
        100, 0., 0., 100, 0., 0.);
    graphPhi = new TH2D(
        "graphPhi", "Linearity: Phi; #phi_{MC} [rad]; #phi_{fit} [rad]", 100, 0., 0., 100, 0., 0.);

    histDiffX = new TH1D("histDiffX", "Pull Plot: decay X position;Pulls;Counts", 50, -10, 10);
    histDiffY = new TH1D("histDiffY", "Pull Plot: decay Y position;Pulls;Counts", 50, -10, 10);
    histDiffZ = new TH1D("histDiffZ", "Pull Plot: decay Z position;Pulls;Counts", 50, -10, 10);
    histDiffMom = new TH1D("histDiffMom", "Pull Plot: Momentum;Pulls;Counts", 50, -10, 10);
    histDiffTheta = new TH1D("histDiffTheta", "Pull Plot: #theta;Pulls;Counts", 50, -10, 10);
    histDiffPhi = new TH1D("histDiffPhi", "Pull Plot: #phi;Pulls;Counts", 50, -10, 10);

    hist2MomRes = new TH2D("hist2MomRes",
        "Histo 2D: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]", 20, 0., 68.9, 50, 0,
        20);
    hist2ThetaRes
        = new TH2D("hist2ThetaRes", "Histo 2D: Theta resolution;#theta [rad];#sigma_{#theta} [rad]",
            40, -TMath::Pi(), TMath::Pi(), 50, 0., 0.5);
    hist2PhiRes = new TH2D("hist2PhiRes", "Histo 2D: Phi resolution;#phi [rad];#sigma_{#phi} [rad]",
        20, -TMath::Pi(), TMath::Pi(), 50, 0., 0.5);

    profMomRes = new TProfile("profMomRes",
        "Profile plot: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]", 20, 0., 68.9, 0,
        40.);
    profThetaRes = new TProfile("profThetaRes",
        "Profile plot: Theta resolution;#theta [rad];#sigma_{#theta} [rad]", 40, -TMath::Pi(),
        TMath::Pi(), 0., 0.5);
    profPhiRes
        = new TProfile("profPhiRes", "Profile plot: Phi resolution;#phi [rad];#sigma_{#phi} [rad]",
            20, -TMath::Pi(), TMath::Pi(), 0., 0.5);

    histFakeHits = new TH1I("histFakeHits", "Number of pre-fitted hits added", 10, 0, 0);

    histTime = new TH1F("histTime", "Time of fitted tracks; Time [ns]; Counts", 100, 0, 0);
}

void TrackDataManager::InitRecoTree(bool isMichel)
{
    if(!outputFile)
        return;
    outputFile->cd();
    if(!recTree)
        recTree = new TTree("rec", "Reconstructed Data");

    recTree->Branch("rec_chi2", &rec_chi2, "rec_chi2/D");
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
        recTree->Branch("is_cosmic_valid", &is_cosmic_valid, "is_cosmic_valid/O");
    }
    else
    {
        recTree->Branch("rec_R", &rec_R, "rec_R/D");
        recTree->Branch("rec_cx", &rec_cx, "rec_cx/D");
        recTree->Branch("rec_cy", &rec_cy, "rec_cy/D");
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
        recTree->Branch("is_michel_valid", &is_michel_valid, "is_michel_valid/O");
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
