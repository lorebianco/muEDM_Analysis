#include "trackdatamanager.hh"

using namespace std;

TrackDataManager::TrackDataManager()
{
    // >>>>> Input <<<<< //
    tracksChain = new TChain("hits_discretized");
    decayChain = new TChain("MuDecayOutput");

    for(const auto &fileName : Config::get().inputDataFiles)
    {
        cout << ">>> Adding file: " << fileName << endl;
        tracksChain->Add(fileName.c_str());
        decayChain->Add(fileName.c_str());
    }

    // Set branch addresses
    hits_EventID = nullptr, hits_TrackID = nullptr, hits_ParticleID = nullptr,
    hits_DetectorID = nullptr, hits_DetName = nullptr, hits_time = nullptr, hits_posX = nullptr,
    hits_posY = nullptr, hits_posZ = nullptr;

    decayChain->SetBranchAddress("EventID", &truedecay_EventID);
    decayChain->SetBranchAddress("muDecayPosX", &trueDecay_posX);
    decayChain->SetBranchAddress("muDecayPosY", &trueDecay_posY);
    decayChain->SetBranchAddress("muDecayPosZ", &trueDecay_posZ);
    decayChain->SetBranchAddress("posIniMomX", &trueDecay_momX);
    decayChain->SetBranchAddress("posIniMomY", &trueDecay_momY);
    decayChain->SetBranchAddress("posIniMomZ", &trueDecay_momZ);
    decayChain->SetBranchAddress("time", &trueDecay_time);

    tracksChain->SetBranchAddress("EventID", &hits_EventID);
    tracksChain->SetBranchAddress("TrackID", &hits_TrackID);
    tracksChain->SetBranchAddress("ParticleID", &hits_ParticleID);
    tracksChain->SetBranchAddress("DetectorID", &hits_DetectorID);
    tracksChain->SetBranchAddress("DetName", &hits_DetName);
    tracksChain->SetBranchAddress("det_x", &hits_posX);
    tracksChain->SetBranchAddress("det_y", &hits_posY);
    tracksChain->SetBranchAddress("det_z", &hits_posZ);
    tracksChain->SetBranchAddress("time", &hits_time);

    // Construct the map
    BuildDecayIndex();

    // Set nEvents and processedEvents
    nEvents = tracksChain->GetEntries();

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

TrackDataManager::~TrackDataManager()
{
    delete decayChain;
    delete tracksChain;
}

void TrackDataManager::BuildDecayIndex()
{
    cout << ">>> Building Decay Event Map..." << endl;

    decayChain->SetBranchStatus("EventID", 1);
    Int_t tempID;
    decayChain->SetBranchAddress("EventID", &tempID);

    Long64_t nEntries = decayChain->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++)
    {
        decayChain->GetEntry(i);
        decayIndexMap[tempID] = i; // Map ID -> Index i
    }

    decayChain->SetBranchStatus("*", 1);
    decayChain->SetBranchAddress("EventID", &truedecay_EventID);

    cout << ">>> Map built. Indexed " << decayIndexMap.size() << " decay events" << endl;
}

void TrackDataManager::ProcessAndFilterEvent()
{
    // Clean containers
    hitsCoordinates.clear();
    hitsCylinderID.clear();

    // Check pointers
    if(!hits_posX || hits_posX->empty())
        return;

    // Synchro TTrees
    Int_t currentEventID = hits_EventID->at(0);

    Bool_t truthFound = SetTrueDecayData(currentEventID);

    if(!truthFound)
        return;

    // Conversion: 1 mm = 0.1 cm
    const Double_t mm2cm = 0.1;

    // Copy in temporary vectors
    vector<string> local_DetName = *hits_DetName;

    // Group floats
    vector<vector<Float_t>> local_fltVecs = {
        *hits_posX, *hits_posY, *hits_posZ, // Index 0, 1, 2
    };

    // Group ints
    vector<vector<Int_t>> local_intVecs = {
        *hits_DetectorID, // Index 0
        *hits_TrackID, // Index 1
        *hits_ParticleID // Index 2
    };

    // FIRST FILTER: Detector Exclusions
    FilterVectors(local_DetName, local_fltVecs, local_intVecs);

    if(local_DetName.empty())
        return;

    // Find positron
    Int_t particleTarget = -11;
    Int_t targetTrackID = FindFirstTrack(local_intVecs[1], local_intVecs[2], particleTarget);

    // SECOND FILTER: Keep only positron track
    if(targetTrackID != -1)
    {
        FilterParticleAndTrack(
            local_DetName, local_fltVecs, local_intVecs, targetTrackID, particleTarget);
    }
    else
    {
        return;
    }

    // SAVE DATA
    size_t nHitsRemaining = local_fltVecs[0].size();
    hitsCoordinates.reserve(nHitsRemaining);
    hitsCylinderID.reserve(nHitsRemaining);

    for(size_t i = 0; i < nHitsRemaining; ++i)
    {
        vector<Double_t> point = { local_fltVecs[0][i] * mm2cm, local_fltVecs[1][i] * mm2cm,
            local_fltVecs[2][i] * mm2cm };
        hitsCoordinates.push_back(point);
        hitsCylinderID.push_back(local_intVecs[0][i]);
    }
}

// >>> Helper Functions for filter <<< //
void TrackDataManager::FilterVectors(
    vector<string> &strVec, vector<vector<float>> &fltVecs, vector<vector<int>> &intVecs)
{
    // Indices to keep
    vector<size_t> indicesToKeep;
    size_t size = strVec.size();

    for(size_t i = 0; i < size; ++i)
    {
        // If string is NOT in the set excludeStrings, keep index
        if(excludeStrings.find(strVec[i]) == excludeStrings.end())
            indicesToKeep.push_back(i);
    }

    // Reconstruct vectors (filter)
    // Strings
    vector<string> newStrVec;
    newStrVec.reserve(indicesToKeep.size());
    for(auto idx : indicesToKeep)
        newStrVec.push_back(strVec[idx]);
    strVec = move(newStrVec);

    // Float
    for(auto &vec : fltVecs)
    {
        vector<Float_t> newVec;
        newVec.reserve(indicesToKeep.size());
        for(auto idx : indicesToKeep)
            newVec.push_back(vec[idx]);
        vec = move(newVec);
    }

    // Int
    for(auto &vec : intVecs)
    {
        vector<Int_t> newVec;
        newVec.reserve(indicesToKeep.size());
        for(auto idx : indicesToKeep)
            newVec.push_back(vec[idx]);
        vec = move(newVec);
    }
}

Int_t TrackDataManager::FindFirstTrack(
    const vector<Int_t> &trackIDs, const vector<Int_t> &particleIDs, Int_t targetParticleID)
{
    Int_t minTrackID = 999999;
    Bool_t found = false;

    for(size_t i = 0; i < trackIDs.size(); ++i)
    {
        if(particleIDs[i] == targetParticleID)
        {
            if(trackIDs[i] < minTrackID)
            {
                minTrackID = trackIDs[i];
                found = true;
            }
        }
    }
    return found ? minTrackID : -1;
}

void TrackDataManager::FilterParticleAndTrack(vector<string> &strVec,
    vector<vector<float>> &fltVecs, vector<vector<Int_t>> &intVecs, Int_t trackID, Int_t particleID)
{
    // Assume: intVecs[1] = TrackID, intVecs[2] = ParticleID
    const auto &vTrackID = intVecs[1];
    const auto &vPartID = intVecs[2];

    vector<size_t> indicesToKeep;
    size_t size = vTrackID.size();

    for(size_t i = 0; i < size; ++i)
    {
        if(vTrackID[i] == trackID && vPartID[i] == particleID)
            indicesToKeep.push_back(i);
    }

    // Apply filter (copy-paste of previous logic)
    vector<string> newStrVec;
    newStrVec.reserve(indicesToKeep.size());
    for(auto idx : indicesToKeep)
        newStrVec.push_back(strVec[idx]);
    strVec = move(newStrVec);

    for(auto &vec : fltVecs)
    {
        vector<Float_t> newVec;
        newVec.reserve(indicesToKeep.size());
        for(auto idx : indicesToKeep)
            newVec.push_back(vec[idx]);
        vec = move(newVec);
    }

    for(auto &vec : intVecs)
    {
        vector<Int_t> newVec;
        newVec.reserve(indicesToKeep.size());
        for(auto idx : indicesToKeep)
            newVec.push_back(vec[idx]);
        vec = move(newVec);
    }
}

Bool_t TrackDataManager::SetTrueDecayData(Int_t eventID)
{
    // --- PARTE 1: Sincronizzazione ---

    // Cerchiamo l'ID nella mappa (che hai riempito nel costruttore)
    auto it = decayIndexMap.find(eventID);

    if(it == decayIndexMap.end())
    {
        // Evento di verità non trovato per questo ID!
        return false;
    }

    // Trovato: carichiamo i dati grezzi dal Tree
    Long64_t entryIndex = it->second;
    decayChain->GetEntry(entryIndex);

    // Ora trueDecay_posX, trueDecay_momX, ecc. contengono i dati giusti.
    // --- PARTE 2: Calcoli e Conversione (La tua vecchia logica) ---

    // Set time of decay
    trueDecayTime = trueDecay_time;

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
