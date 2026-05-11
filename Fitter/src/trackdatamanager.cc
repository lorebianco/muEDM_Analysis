#include <TNamed.h>
#include <TParameter.h>

#include "auxiliaryalgorithms.hh"
#include "config.hh"
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

    simChain->SetBranchAddress("mc_hits_x", &mc_hits_x);
    simChain->SetBranchAddress("mc_hits_y", &mc_hits_y);
    simChain->SetBranchAddress("mc_hits_z", &mc_hits_z);
    simChain->SetBranchAddress("mc_x", &trueDecay_posX);
    simChain->SetBranchAddress("mc_y", &trueDecay_posY);
    simChain->SetBranchAddress("mc_z", &trueDecay_posZ);
    simChain->SetBranchAddress("mc_px", &trueDecay_momX);
    simChain->SetBranchAddress("mc_py", &trueDecay_momY);
    simChain->SetBranchAddress("mc_pz", &trueDecay_momZ);
    simChain->SetBranchAddress("mc_E", &mc_E);

    if(simChain->GetBranch("trk_R"))
    {
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

    // --- EXTRACT GEOMETRY MISALIGNMENTS FROM USERINFO ---
    CHeT::Config::GetTranslation(
        Config::get().geom_tx, Config::get().geom_ty, Config::get().geom_tz);
    CHeT::Config::GetRotation(Config::get().geom_rx, Config::get().geom_ry, Config::get().geom_rz);

    simChain->LoadTree(0);
    if(TTree *firstTree = simChain->GetTree())
    {
        if(TList *userInfo = firstTree->GetUserInfo())
        {
            auto getParam = [&](const char *name, double defVal) -> double
            {
                auto *p = dynamic_cast<TParameter<double> *>(userInfo->FindObject(name));
                return p ? p->GetVal() : defVal;
            };

            double tx = getParam("Geom_Tx", 0.0);
            double ty = getParam("Geom_Ty", 0.0);
            double tz = getParam("Geom_Tz", 0.0);
            double rx = getParam("Geom_Rx", 0.0);
            double ry = getParam("Geom_Ry", 0.0);
            double rz = getParam("Geom_Rz", 0.0);
            double offset = getParam("Geom_OffsetExp", 0.0);
            std::vector<double> deltas(6, 0.0);
            for(int i = 0; i < 6; ++i)
            {
                deltas[i] = getParam(Form("Geom_Delta_%d", i), 0.0);
            }

            std::vector<int> active_cyls;
            bool found_active_meta = false;
            for(int i = 0; i < 6; ++i)
            {
                double val = getParam(Form("ActiveCyl_%d", i), -1.0);
                if(val >= 0.0)
                {
                    found_active_meta = true;
                    if(val > 0.5)
                        active_cyls.push_back(i);
                }
            }

            if(found_active_meta)
            {
                CHeT::Config::SetActiveCylinders(active_cyls);
                cout << ">>> Inherited Active Cylinders from input file!" << endl;
            }

            if(Config::get().useTrueMCGeom)
            {
                CHeT::Config::SetTranslation(tx, ty, tz);
                CHeT::Config::SetRotation(rx, ry, rz);
                CHeT::Config::SetOffsetExp(offset);
                CHeT::Config::SetDeltas(deltas);
                cout << ">>> Inherited TRUE MC Geometry Misalignments from input file!" << endl;
            }
            else if(Config::get().overrideGeom)
            {
                std::vector<int> act_cyls;
                for(int i = 0; i < 6; ++i)
                {
                    if(Config::get().active_cyls[i])
                        act_cyls.push_back(i);
                }
                CHeT::Config::SetActiveCylinders(act_cyls);

                CHeT::Config::SetTranslation(
                    Config::get().geom_tx, Config::get().geom_ty, Config::get().geom_tz);
                CHeT::Config::SetRotation(
                    Config::get().geom_rx, Config::get().geom_ry, Config::get().geom_rz);
                CHeT::Config::SetOffsetExp(Config::get().geom_offset);
                CHeT::Config::SetDeltas(Config::get().geom_deltas);
                cout << ">>> GUI Override Geometry applied!" << endl;
            }
            else
            {
                // Nominal
                // as implemented in the CHeTLibrary
                cout << ">>> Nominal Geometry applied for Reconstruction!" << endl;
            }
        }
    }

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

            // Clone tree actually keeps UserInfo.

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

    // Copia i metadati dal simTree (la geometria e le configurazioni)
    if(simChain)
    {
        if(TTree *st = simChain->GetTree())
        {
            if(TList *userInfo = st->GetUserInfo())
            {
                TList *recUserInfo = recTree->GetUserInfo();
                for(int i = 0; i < userInfo->GetSize(); ++i)
                {
                    recUserInfo->Add(userInfo->At(i)->Clone());
                }
            }
        }
    }

    if(!isMichel)
    {
        recTree->Branch("rec_x0", &rec_x0, "rec_x0/D");
        recTree->Branch("rec_z0", &rec_z0, "rec_z0/D");
        recTree->Branch("rec_sx", &rec_sx, "rec_sx/D");
        recTree->Branch("rec_sz", &rec_sz, "rec_sz/D");
        recTree->Branch("rec_cov_cosmic", "TMatrixDSym", &rec_cov_cosmic);
        recTree->Branch("rec_zi", &rec_zi);
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
        recTree->Branch("rec_cov_extrap", "TMatrixDSym", &rec_cov_extrap);

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
    rec_zi.clear();
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

    if(Config::get().useTrueMCHits && mc_hits_x && mc_hits_y && mc_hits_z)
    {
        size_t nHitsMC = mc_hits_x->size();
        size_t nBundleIDs = allEventsHits[currentEventIndex].size();

        struct TempHit
        {
            double x, y, z;
            int cyl;
            bool processed;
        };

        std::vector<TempHit> tempHits;
        for(size_t i = 0; i < nHitsMC; ++i)
        {
            int b_id = (i < nBundleIDs) ? allEventsHits[currentEventIndex][i] : -1;
            tempHits.push_back({ mc_hits_x->at(i) * mm2cm, mc_hits_y->at(i) * mm2cm,
                mc_hits_z->at(i) * mm2cm, CHeT::Config::GetFiberProp(b_id).cylinderId, false });
        }

        hitsCoordinates.clear();
        hitsCylinderID.clear();

        for(size_t i = 0; i < tempHits.size(); ++i)
        {
            if(tempHits[i].processed)
                continue;

            double sumX = tempHits[i].x;
            double sumY = tempHits[i].y;
            double sumZ = tempHits[i].z;
            int count = 1;
            tempHits[i].processed = true;

            // Cerca altre hit nello stesso cilindro e nella stessa posizione
            for(size_t j = i + 1; j < tempHits.size(); ++j)
            {
                if(tempHits[j].processed)
                    continue;

                if(tempHits[i].cyl == tempHits[j].cyl)
                {
                    double dist = std::sqrt(std::pow(tempHits[i].x - tempHits[j].x, 2)
                        + std::pow(tempHits[i].y - tempHits[j].y, 2)
                        + std::pow(tempHits[i].z - tempHits[j].z, 2));

                    // Se la distanza è minore di 1 mm (0.1 cm), le consideriamo la stessa "misura"
                    if(dist < 0.1)
                    {
                        sumX += tempHits[j].x;
                        sumY += tempHits[j].y;
                        sumZ += tempHits[j].z;
                        count++;
                        tempHits[j].processed = true;
                    }
                }
            }

            // Calcola la media
            std::vector<Double_t> finalPoint = { sumX / count, sumY / count, sumZ / count };
            int finalCyl = tempHits[i].cyl;

            if(Config::get().useSmearing)
                finalPoint = AUXALG::SmearMeasurement(finalCyl, finalPoint);

            hitsCoordinates.push_back(finalPoint);
            hitsCylinderID.push_back(finalCyl);
        }
    }
    else
    {
        auto intersections = CHeT::Config::FindIntersections(hit_ids);

        for(const auto &inter : intersections)
        {
            std::vector<Double_t> point = { inter.x * mm2cm, inter.y * mm2cm, inter.z * mm2cm };
            hitsCoordinates.push_back(point);
            hitsCylinderID.push_back(inter.cylinderId);
        }
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
