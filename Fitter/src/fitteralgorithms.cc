#include "analyticviewer.hh"
#include "fitteralgorithms.hh"

using namespace std;
using namespace ROOT;

constexpr Int_t DEBUG_LVL = 0;
constexpr Bool_t NORMALIZED_PULLS = false;
constexpr Double_t FAKEHIT_ANGLE_LIM = TMath::Pi() / 20;
constexpr Int_t pdg = -11; // PID for positron

double FITALG::Track3DNeg2LogL(const double *par, const FITALG::FitData &data, bool usePrior)
{
    const double x0 = par[0], sx = par[1], z0 = par[2], sz = par[3];
    const double sigma2 = 0.3;
    const double sigmaL2 = 1e-6;

    double n2ll = 0.0;
    if(usePrior)
        n2ll = 2.0 * log(1.0 + sx * sx + sz * sz);

    for(size_t i = 0; i < data.props.size(); ++i)
    {
        const double zi = par[4 + i];
        const auto &p = data.props[i];
        const double alpha = (zi + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
        const double phi_fib = p.phi0 + p.dir * alpha * M_PI;

        const double x_f = p.r * cos(phi_fib), y_f = p.r * sin(phi_fib);
        const double x_trk = x0 + sx * y_f, z_trk = z0 + sz * y_f;

        n2ll += ((x_f - x_trk) * (x_f - x_trk) + (zi - z_trk) * (zi - z_trk)) / sigma2;
        if(abs(zi) > CHeT::Config::L_HALF)
            n2ll += (pow(abs(zi) - CHeT::Config::L_HALF, 2)) / sigmaL2;
    }
    return n2ll;
}

FITALG::FitOutput FITALG::Do3DFit(const std::vector<int> &hit_ids, bool usePrior)
{
    FITALG::FitData data;
    for(int id : hit_ids)
        data.props.push_back(CHeT::Config::GetFiberProp(id));

    int n_hits = data.props.size();
    if(n_hits < 3)
        return { { 0, 0, 0, 0, -1.0, false }, {} };

    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min->SetPrintLevel(0);

    auto chi2_func
        = [&](const double *par) { return FITALG::Track3DNeg2LogL(par, data, usePrior); };
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

    FITALG::FitOutput output;
    output.track = { res[0], res[2], res[1], res[3], chi2, conv };

    if(conv)
    {
        for(int i = 0; i < n_hits; ++i)
        {
            double z_fit = res[4 + i];
            const auto &p = data.props[i];
            double alpha = (z_fit + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
            double phi = p.phi0 + p.dir * alpha * M_PI;
            double x_fit = p.r * cos(phi), y_fit = p.r * sin(phi);
            output.fittedPoints.emplace_back(x_fit, y_fit, z_fit, kBlack, 24, 1.0, true);
        }
    }
    return output;
}

void FITALG::CosmicFitter()
{
    cout << ">>> Running CosmicFitter..." << endl;

    TrackDataManager data;
    data.InitRecoTree(false);

    auto &config = Config::get();

    // 1. Istanzia il viewer personalizzato
    FITALG::AnalyticViewer *myViewer = new FITALG::AnalyticViewer();

    // Carica la geometria globale (GDML o ROOT file) tramite il viewer
    if(!config.geometryFile.empty())
    {
        myViewer->LoadGeometry(config.geometryFile.c_str());
    }

    Int_t nEvents = data.GetNEvents();
    bool hasSIM = (data.GetChain()->GetEntries() > 0);

    Long64_t maxEvents = (config.rangeLoop.second < nEvents) ? config.rangeLoop.second : nEvents;
    Long64_t startEvent = config.rangeLoop.first;

    gRandom->SetSeed(0);
    if(config.processSingle && config.event == -1)
        config.event = gRandom->Integer(nEvents);

    if(!config.quietMode && !config.processSingle)
    {
        cout << ">>> Processing multiple events. Results won't be printed, but they are stored in "
                "TEve if available."
             << endl;
    }

    for(Long64_t i = startEvent; i < maxEvents; ++i)
    {
        if(config.processSingle && i != config.event)
            continue;

        if(!config.quietMode && !config.processSingle && i % 1000 == 0)
        {
            cout << ">>> Processed " << i << " events..." << endl;
        }

        // Delego I/O a TrackDataManager
        data.ProcessAndFilterEvent(i);

        const auto &hit_ids = data.rec_hits;

        data.rec_acceptance = false;
        data.rec_converged = false;
        data.rec_x0 = 0;
        data.rec_z0 = 0;
        data.rec_sx = 0;
        data.rec_sz = 0;
        data.rec_chi2 = -1;
        data.rec_hough2d_idx.clear();
        data.rec_houghz_idx.clear();

        bool is_converged = false;
        FitOutput fitRes;
        RecoTrack trFit;

        if(hit_ids.size() >= 4)
        {
            data.rec_acceptance = true;
            fitRes = Do3DFit(hit_ids, false);
            trFit = fitRes.track;

            data.rec_x0 = trFit.x0;
            data.rec_z0 = trFit.z0;
            data.rec_sx = trFit.sx;
            data.rec_sz = trFit.sz;
            data.rec_chi2 = trFit.chi2;
            data.rec_converged = trFit.converged;

            is_converged = trFit.converged;
        }

        if(data.recTree)
            data.recTree->Fill();

        double mc_x0 = 0, mc_y0 = 0, mc_z0 = 0, mc_ux = 0, mc_uy = 0, mc_uz = 0;
        if(hasSIM)
        {
            // Convert da cm a mm se necessario (TrackDataManager memorizza in cm per default)
            mc_x0 = data.trueDecay_posX / 0.1;
            mc_y0 = data.trueDecay_posY / 0.1;
            mc_z0 = data.trueDecay_posZ / 0.1;
            mc_ux = data.trueDecay_momX;
            mc_uy = data.trueDecay_momY;
            mc_uz = data.trueDecay_momZ;
        }

        if(is_converged)
        {
            myViewer->AddEvent(i, fitRes);

            if(hasSIM)
            {
                double trueLoc_x0 = mc_x0, trueLoc_y0 = mc_y0, trueLoc_z0 = mc_z0;
                CHeT::Config::ApplyInverseTransformation(trueLoc_x0, trueLoc_y0, trueLoc_z0);

                double trueLoc_ux = mc_ux, trueLoc_uy = mc_uy, trueLoc_uz = mc_uz;
                CHeT::Config::ApplyInverseRotation(trueLoc_ux, trueLoc_uy, trueLoc_uz);

                double t_to_y0 = -trueLoc_y0 / trueLoc_uy;
                double true_x0 = trueLoc_x0 + trueLoc_ux * t_to_y0;
                double true_z0 = trueLoc_z0 + trueLoc_uz * t_to_y0;
                double true_sx = trueLoc_ux / trueLoc_uy;
                double true_sz = trueLoc_uz / trueLoc_uy;

                double true_cosTheta = trueLoc_uy;
                double true_phi = atan2(trueLoc_ux, trueLoc_uz);

                double norm = sqrt(trFit.sx * trFit.sx + 1.0 + trFit.sz * trFit.sz);
                double fit_uy = -1.0 / norm;
                double fit_ux = trFit.sx * fit_uy;
                double fit_uz = trFit.sz * fit_uy;

                double fit_cosTheta = fit_uy;
                double fit_phi = atan2(fit_ux, fit_uz);

                if(config.processSingle)
                {
                    printf("\nEvent %lld\n", i);
                    printf("===============================================================\n");
                    printf("                   FIT RESULTS vs MC TRUTH                     \n");
                    printf("===============================================================\n");
                    printf(" %-10s | %12s | %12s | %12s \n", "Param", "True", "Fitted", "Residual");
                    printf("---------------------------------------------------------------\n");
                    printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "x0 [mm]", true_x0, trFit.x0,
                        trFit.x0 - true_x0);
                    printf(" %-10s | %12.4f | %12.4f | %12.4f \n", "z0 [mm]", true_z0, trFit.z0,
                        trFit.z0 - true_z0);
                    printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sx (dx/dy)", true_sx, trFit.sx,
                        trFit.sx - true_sx);
                    printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "sz (dz/dy)", true_sz, trFit.sz,
                        trFit.sz - true_sz);
                    printf("---------------------------------------------------------------\n");
                    printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "cos(theta)", true_cosTheta,
                        fit_cosTheta, fit_cosTheta - true_cosTheta);
                    printf(" %-10s | %12.5f | %12.5f | %12.5f \n", "phi [rad]", true_phi, fit_phi,
                        fit_phi - true_phi);
                    printf("---------------------------------------------------------------\n");
                    printf(" Chi2/ndf  : %.4f\n", trFit.chi2);
                    printf("===============================================================\n");
                }
            }
            else
            {
                double norm = sqrt(trFit.sx * trFit.sx + 1.0 + trFit.sz * trFit.sz);
                double fit_uy = -1.0 / norm;
                double fit_ux = trFit.sx * fit_uy;
                double fit_uz = trFit.sz * fit_uy;
                double fit_cosTheta = fit_uy;
                double fit_phi = atan2(fit_ux, fit_uz);

                if(config.processSingle)
                {
                    printf("\nEvent %lld\n", i);
                    printf("===============================================================\n");
                    printf("                   FIT RESULTS                                 \n");
                    printf("===============================================================\n");
                    printf(" %-10s | %12.4f \n", "x0 [mm]", trFit.x0);
                    printf(" %-10s | %12.4f \n", "z0 [mm]", trFit.z0);
                    printf(" %-10s | %12.5f \n", "sx (dx/dy)", trFit.sx);
                    printf(" %-10s | %12.5f \n", "sz (dz/dy)", trFit.sz);
                    printf(" %-10s | %12.5f \n", "cos(theta)", fit_cosTheta);
                    printf(" %-10s | %12.5f \n", "phi [rad]", fit_phi);
                    printf(" Chi2/ndf  : %.4f\n", trFit.chi2);
                    printf("===============================================================\n");
                }
            }
        }
        else
        {
            if(config.processSingle)
                printf("\nEvent %lld: [Fit Failed] Minuit did not converge.\n", i);
        }

        if(!config.quietMode)
        {
            std::vector<CHeT::Vis::VisLineTrack> visTracks;
            if(hasSIM)
                visTracks.emplace_back(mc_x0, mc_y0, mc_z0, mc_ux, mc_uy, mc_uz, kYellow, 3);

            if(is_converged)
            {
                visTracks.emplace_back(
                    trFit.x0, 0.0, trFit.z0, trFit.sx, 1.0, trFit.sz, kBlack, 2, 7, true);
            }

            if(config.processSingle)
            {
                CHeT::Vis::Draw2D(hit_ids, visTracks);
                CHeT::Vis::Draw3D(hit_ids, visTracks, fitRes.fittedPoints);
            }
        }
    }

    if(!config.quietMode && gApplication)
    {
        myViewer->Show();

        // Necessario per tenere aperta l'interfaccia se lanciato come applicazione
        if(gApplication)
            gApplication->Run(true);
    }

    // Finally
    data.SaveRecoTree();
    gSystem->Exit(0);
    exit(0);
}

void FITALG::SpacepointFitter()
{
    // Load events
    TrackDataManager data;
    data.InitRecoTree(true);

    Int_t nEvents = data.GetNEvents();
    Int_t processedEvents = data.GetProcessedEvents();

    // Init geometry and magnetic field
    new TGeoManager("DetectorGeometry", "CHET geometry");
    TGeoManager::Import(Config::get().geometryFile.c_str());
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

    // B field
    try
    {
        muEDM::Fields::Configure(Config::get().bFieldType, // Enum (Multi, Single, Const)
            Config::get().bFieldPath, // String path
            Config::get().bScale, // Scale factor
            Config::get().bConstZ // Constant Bz value
        );

        // Call the setup
        muEDM::Fields::SetupBField();
    }
    catch(const exception &e)
    {
        cerr << ">>> Critical error during field maps loading: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Init event display
    genfit::EventDisplay *display = genfit::EventDisplay::getInstance();
    display->reset();
    cout << endl;

    // Init fitter (maxIterations, deltaPVal) (Possible values = 20, 1.E-3)
    genfit::AbsKalmanFitter *fitter = new genfit::KalmanFitterRefTrack(20, 1.E-3);
    fitter->setDebugLvl(DEBUG_LVL);

    // Create array of hits
    TClonesArray chetHitArray("genfit::mySpacepointDetectorHit");

    // Init the factory
    Int_t detId = 0;
    genfit::MeasurementFactory<genfit::AbsMeasurement> chetFactory;
    genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement>
        cylProducer(&chetHitArray);
    chetFactory.addProducer(detId, &cylProducer);

    // Create Track
    genfit::Track *fitTrack = nullptr;

    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(Config::get().processSingle && Config::get().event == -1)
        Config::get().event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0;
    Int_t inEfficiency = 0;
    Int_t nTurns, nCylinders;
    for(Long_t ev = Config::get().rangeLoop.first; ev < Config::get().rangeLoop.second; ev++)
    {
        if(ev >= nEvents)
            break;
        if(Config::get().processSingle)
            if(ev != Config::get().event)
                continue;

        // Clean up
        delete fitTrack;
        fitTrack = nullptr;
        chetHitArray.Clear();

        // Get event and filter
        data.ProcessAndFilterEvent(ev);

        Int_t nHits = data.hitsCoordinates.size();

        if(Config::get().processSingle)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Check momentum: RKTrackRep can handle if > 4 MeV
        // if(mom.Mag()*1E3 < 5)
        //{
        //    processedEvents--;
        //    continue;
        //}

        if(Config::get().processSingle)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", data.trueDecayPos.X(),
                data.trueDecayPos.Y(), data.trueDecayPos.Z())
                 << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV/c", data.trueDecayMom.X(),
                data.trueDecayMom.Y(), data.trueDecayMom.Z())
                 << endl;
            cout << Form(">>> (p, thetaEDM, phi) = (%f MeV/c, %f pi rad, %f pi rad)",
                data.trueDecayMom.Mag(), data.trueDecayThetaEDM / TMath::Pi(),
                data.trueDecayMom.Phi() / TMath::Pi())
                 << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(Config::get().processSingle)
                cout << ">>> Track is not in acceptance!" << endl;

            data.accTheta->Fill(false, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
            data.accPhi->Fill(false, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

            continue;
        }

        // Sort and if in single event mode draw hits
        // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = PTTALG::SortVectorZ(data.hitsCoordinates, data.hitsCylinderID);
        data.hitsCoordinates = sortedHits.first;
        data.hitsCylinderID = sortedHits.second;
        nTurns = PTTALG::CountTurns(data.hitsCoordinates);
        nCylinders = PTTALG::CountCylinders(data.hitsCylinderID);

        // Apply turn analysis
        if(Config::get().turnMode)
        {
            if(Config::get().processSingle)
            {
                cout << ">>> nTurns = " << nTurns << endl;
                cout << ">>> nCylinders = " << nCylinders << endl;
            }

            auto turnHits = PTTALG::SelectTurn(
                Config::get().turnID, data.hitsCoordinates, data.hitsCylinderID);
            data.hitsCoordinates = turnHits.first;
            data.hitsCylinderID = turnHits.second;

            //! Cylinders Filter... Will be removed
            auto selCylHits
                = PTTALG::SelectCylinders({ 0, 1 }, data.hitsCoordinates, data.hitsCylinderID);
            data.hitsCoordinates = selCylHits.first;
            data.hitsCylinderID = selCylHits.second;

            // Update nHits
            nHits = data.hitsCoordinates.size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(Config::get().processSingle)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                data.accTheta->Fill(false, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
                data.accPhi->Fill(false, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

                continue;
            }
        }

        // Track is in acceptance!
        inAcceptance++;
        data.accTheta->Fill(true, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
        data.accPhi->Fill(true, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

        // Fast detector simulation
        if(Config::get().useSmearing)
        {
            for(Int_t i = 0; i < nHits; i++)
                data.hitsCoordinates[i]
                    = AUXALG::SmearMeasurement(data.hitsCylinderID[i], data.hitsCoordinates[i]);
        }

        // Prefitter
        Double_t xC, yC, R, phi0, z0, tanLambda;
        vector<Double_t> s_cumulative;
        if(Config::get().usePrefitter)
        {
            auto helixPars = HelixPrefitter(data.hitsCoordinates, data.hitsCylinderID);

            xC = helixPars[0], yC = helixPars[1], R = helixPars[2], phi0 = helixPars[3],
            z0 = helixPars[4], tanLambda = helixPars[5];

            // Cumulative arc length computation
            Double_t previous_phi = phi0;
            Double_t previous_s = 0;

            for(const auto &point : data.hitsCoordinates)
            {
                Double_t dx = point[0] - xC;
                Double_t dy = point[1] - yC;
                Double_t phi = atan2(dy, dx);

                Double_t dphi = phi - previous_phi;
                if(dphi > M_PI)
                    dphi -= 2 * M_PI;
                if(dphi < -M_PI)
                    dphi += 2 * M_PI;

                Double_t s = previous_s + R * abs(dphi);

                s_cumulative.push_back(s);
                previous_phi = phi;
                previous_s = s;
            }
        }

        // Track candidate
        genfit::TrackCand trackCand;

        // Resolution of detectors
        const CHeT::Resolutions hitCov;

        // Fill the candidate
        vector<TVector3> measuredCoordinates;
        vector<TVector3> virtualCoordinates;
        Int_t hitID = 0;
        Int_t nAddedHits = 0;
        for(Int_t i = 0; i < nHits; i++)
        {
            TVector3 hitCoords;
            vector<Double_t> measuredCoords = data.hitsCoordinates[i];

            hitCoords[0] = measuredCoords[0];
            hitCoords[1] = measuredCoords[1];
            hitCoords[2] = measuredCoords[2];

            measuredCoordinates.push_back(hitCoords);

            new(chetHitArray[hitID]) genfit::mySpacepointDetectorHit(hitCoords,
                hitCov.GetMatrixCartesian(
                    data.hitsCylinderID[i], TMath::ATan2(hitCoords[1], hitCoords[0])));

            trackCand.addHit(detId, hitID, -1, i);
            hitID++;

            // LeTrick
            if(Config::get().usePrefitter && ((i + 1) < nHits)
                && ((s_cumulative[i + 1] - s_cumulative[i]) / R > FAKEHIT_ANGLE_LIM))
            {
                const Double_t s_middle = 0.5 * (s_cumulative[i] + s_cumulative[i + 1]);
                AddFakeHitFromHelix(trackCand, hitID, i + 0.5, s_middle, xC, yC, R, z0, phi0,
                    tanLambda, chetHitArray, virtualCoordinates, 1.);
                nAddedHits++;
                hitID++;
            }
        }

        // Plotting
        if(Config::get().processSingle)
        {
            vector<vector<Double_t>> plottedCoordsVec;
            for(const auto &vec : measuredCoordinates)
                plottedCoordsVec.push_back({ vec.X(), vec.Y(), vec.Z() });

            vector<vector<Double_t>> virtualCoordsVec;
            for(const auto &vec : virtualCoordinates)
                virtualCoordsVec.push_back({ vec.X(), vec.Y(), vec.Z() });

            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ", "canvHitsYZ", 900, 500);
            TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ");

            AUXALG::DrawXYView_hits(
                &data.trueDecayPos, plottedCoordsVec, canvHitsXY, virtualCoordsVec);
            AUXALG::DrawYZView_hits(
                &data.trueDecayPos, plottedCoordsVec, canvHitsYZ, virtualCoordsVec);
            AUXALG::DrawXYZView_hits(plottedCoordsVec, canvHitsXYZ, virtualCoordsVec);
        }

        // Start values for fit (would come from pattern recognition -> from helix prefitter)
        // Initial guess for cov
        const Double_t seedResolution = 0.1; // ?
        TMatrixDSym covSeed(6);
        for(Int_t i = 0; i < 3; i++)
            covSeed(i, i) = seedResolution * seedResolution;
        for(Int_t i = 3; i < 6; i++)
            covSeed(i, i) = pow(seedResolution / nHits / sqrt(3), 2);

        // Set start values
        TVector3 posSeed(data.trueDecayPos); // cm
        TVector3 momSeed(data.trueDecayMom * 1E-3); // GeV
        if(Config::get().pttrecMode)
            tie(posSeed, momSeed) = PTTALG::SmearSeed(data.trueDecayPos, data.trueDecayMom);

        if(Config::get().processSingle)
        {
            cout << endl;
            cout << Form(">>> Seed position = (%f, %f, %f) cm", posSeed[0], posSeed[1], posSeed[2])
                 << endl;
            cout << Form(">>> Seed momentum = (%f, %f, %f) MeV", momSeed[0] * 1E3, momSeed[1] * 1E3,
                momSeed[2] * 1E3)
                 << endl;
        }

        trackCand.setPosMomSeedAndPdgCode(posSeed, momSeed, pdg);
        trackCand.setCovSeed(covSeed);

        // Track rep
        genfit::AbsTrackRep *rep = new genfit::RKTrackRep(pdg);

        // Create track
        fitTrack = new genfit::Track(trackCand, chetFactory, rep);

        // Check
        fitTrack->checkConsistency();
        // fitTrack->Print();

        // Do the fit
        try
        {
            fitter->processTrack(fitTrack, true);
        }
        catch(genfit::Exception &e)
        {
            cerr << e.what();
            cerr << "Exception, next track" << endl;
            continue;
        }

        // Fit result
        // fitTrack->Print();
        Bool_t isFitConverged = fitTrack->getFitStatus(rep)->isFitConverged();

        if(isFitConverged)
        {
            inEfficiency++;

            auto [fRes, fSigma, fPulls]
                = AUXALG::GetResults(fitTrack, rep, data.trueDecayPos, data.trueDecayMom.Mag(),
                    data.trueDecayThetaEDM, data.trueDecayMom.Phi(), NORMALIZED_PULLS);
            if(fRes.size() == 0)
                continue;
            data.histDiffX->Fill(fPulls[0]);
            data.histDiffY->Fill(fPulls[1]);
            data.histDiffZ->Fill(fPulls[2]);
            data.histDiffMom->Fill(fPulls[3]);
            data.histDiffTheta->Fill(fPulls[4]);
            data.histDiffPhi->Fill(fPulls[5]);

            data.graphMom->Fill(data.trueDecayMom.Mag(), fRes[3]);
            data.graphTheta->Fill(data.trueDecayThetaEDM, fRes[4]);
            data.graphPhi->Fill(data.trueDecayMom.Phi(), fRes[5]);

            data.hist2MomRes->Fill(data.trueDecayMom.Mag(), fSigma[3]);
            data.hist2ThetaRes->Fill(data.trueDecayThetaEDM, fSigma[4]);
            data.hist2PhiRes->Fill(data.trueDecayMom.Phi(), fSigma[5]);

            data.profMomRes->Fill(data.trueDecayMom.Mag(), fSigma[3]);
            data.profThetaRes->Fill(data.trueDecayThetaEDM, fSigma[4]);
            data.profPhiRes->Fill(data.trueDecayMom.Phi(), fSigma[5]);

            data.histTime->Fill(data.trueDecayTime);

            if(Config::get().processSingle)
                cout << Form(">>> Fitted (p, thetaEDM, phi) = (%f MeV/c, %f pi rad, %f pi rad)",
                    fRes[3], fRes[4] / TMath::Pi(), fRes[5] / TMath::Pi())
                     << endl;
        }

        if(Config::get().processSingle)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n"
                 << endl;
        }

        // Track is in efficiency?
        data.effTheta->Fill(isFitConverged, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
        data.effPhi->Fill(isFitConverged, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

        // Store turns and cylinders data
        data.histTurns->Fill(nTurns);
        data.effTurns->Fill(isFitConverged, nTurns);
        data.histCylinders->Fill(nCylinders);
        data.effCylinders->Fill(isFitConverged, nCylinders);
        data.histTurnsVMom->Fill(data.trueDecayMom.Mag(), nTurns);
        data.histCylVMom->Fill(data.trueDecayMom.Mag(), nCylinders);

        // Store pre-fit data
        data.histFakeHits->Fill(nAddedHits);

        // Check
        fitTrack->checkConsistency();

        // Last steps
        display->addEvent(fitTrack);

        cout << "\r>>> Processed event number " << (ev - Config::get().rangeLoop.first) << flush;
    }
    cout << endl;

    // Print results and draw graphs
    if(!Config::get().processSingle)
    {
        // Recap
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t)(inAcceptance * 100) / processedEvents << endl;
        cout << ">>> Efficiency = " << (Float_t)(inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;

        // Graphs
        if(Config::get().quietMode)
            gROOT->SetBatch(true);

        TCanvas *canvEfficiency = new TCanvas("canvEfficiency");
        canvEfficiency->Divide(2, 2);
        canvEfficiency->cd(1);
        data.accTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(2);
        data.accPhi->Draw("COLZ TEXT");
        canvEfficiency->cd(3);
        data.effTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(4);
        data.effPhi->Draw("COLZ TEXT");

        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        data.histTurns->Draw();
        canvTurns->cd(2);
        data.effTurns->Draw("AP");

        TCanvas *canvCylinders = new TCanvas("canvCylinders");
        canvCylinders->Divide(2, 2);
        canvCylinders->cd(1);
        data.histCylinders->Draw();
        canvCylinders->cd(2);
        data.effCylinders->Draw("AP");
        canvCylinders->cd(3);
        data.histTurnsVMom->Draw();
        canvCylinders->cd(4);
        data.histCylVMom->Draw();

        TCanvas *canvLinearity = new TCanvas("canvLinearity");
        canvLinearity->Divide(3);
        canvLinearity->cd(1);
        data.graphMom->SetMarkerStyle(20);
        data.graphMom->Draw("SCAT");
        canvLinearity->cd(2);
        data.graphTheta->SetMarkerStyle(20);
        data.graphTheta->Draw("SCAT");
        canvLinearity->cd(3);
        data.graphPhi->SetMarkerStyle(20);
        data.graphPhi->Draw("SCAT");

        TCanvas *canvfPullsPos = new TCanvas("Pulls Position");
        canvfPullsPos->Divide(3);

        canvfPullsPos->cd(1);
        data.histDiffX->Draw();

        canvfPullsPos->cd(2);
        data.histDiffY->Draw();

        canvfPullsPos->cd(3);
        data.histDiffZ->Draw();

        TCanvas *canvfPullsMom = new TCanvas("Pulls Momentum");
        canvfPullsMom->Divide(3);

        canvfPullsMom->cd(1);
        data.histDiffMom->Draw();

        canvfPullsMom->cd(2);
        data.histDiffTheta->Draw();

        canvfPullsMom->cd(3);
        data.histDiffPhi->Draw();

        TCanvas *canvProfRes = new TCanvas("ProfResolutions");
        canvProfRes->Divide(3);
        canvProfRes->cd(1);
        data.profMomRes->Draw();
        canvProfRes->cd(2);
        data.profThetaRes->Draw();
        canvProfRes->cd(3);
        data.profPhiRes->Draw();

        TCanvas *canvResMom = new TCanvas("canvResMom");
        canvResMom->cd();
        data.hist2MomRes->Draw();

        TCanvas *canvResTheta = new TCanvas("canvResTheta");
        canvResTheta->cd();
        data.hist2ThetaRes->Draw();

        TCanvas *canvResPhi = new TCanvas("canvResPhi");
        canvResPhi->cd();
        data.hist2PhiRes->Draw();

        TCanvas *canvPrefit = new TCanvas("canvPrefit");
        canvPrefit->cd();
        data.histFakeHits->Draw();

        TCanvas *canvTime = new TCanvas("canvTime");
        canvTime->cd();
        data.histTime->Draw();

        if(Config::get().quietMode)
        {
            canvEfficiency->SaveAs("canvEfficiency.pdf");
            canvResMom->SaveAs("canvResMom.pdf");
            canvResTheta->SaveAs("canvResTheta.pdf");
            canvResPhi->SaveAs("canvResPhi.pdf");
        }
    }

    // Delete fitter
    delete fitter;

    // Open event display
    if(Config::get().quietMode)
        display->setOptions("X");
    else
        display->setOptions("ABDEFGHMPT");
    display->open();

    // Finally
    data.SaveRecoTree();
    gSystem->Exit(0);
    exit(0);
}

void FITALG::HelixFitter()
{
    // Load events
    TrackDataManager data;
    data.InitRecoTree(true);

    Int_t nEvents = data.GetNEvents();
    Int_t processedEvents = data.GetProcessedEvents();

    genfit::EventDisplay *display = genfit::EventDisplay::getInstance();

    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(Config::get().processSingle && Config::get().event == -1)
        Config::get().event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0;
    Int_t inEfficiency = 0;
    Int_t nTurns, nCylinders;
    for(Long_t ev = Config::get().rangeLoop.first; ev < Config::get().rangeLoop.second; ev++)
    {
        if(ev >= nEvents)
            break;
        if(Config::get().processSingle)
            if(ev != Config::get().event)
                continue;

        // Get event
        data.ProcessAndFilterEvent(ev);

        Int_t nHits = data.hitsCoordinates.size();

        if(Config::get().processSingle)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        if(Config::get().processSingle)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", data.trueDecayPos.X(),
                data.trueDecayPos.Y(), data.trueDecayPos.Z())
                 << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV/c", data.trueDecayMom.X(),
                data.trueDecayMom.Y(), data.trueDecayMom.Z())
                 << endl;
            cout << Form(">>> (p, thetaEDM, phi) = (%f MeV/c, %f pi rad, %f pi rad)",
                data.trueDecayMom.Mag(), data.trueDecayThetaEDM / TMath::Pi(),
                data.trueDecayMom.Phi() / TMath::Pi())
                 << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(Config::get().processSingle)
                cout << ">>> Track is not in acceptance!" << endl;

            data.accTheta->Fill(false, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
            data.accPhi->Fill(false, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
        // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = PTTALG::SortVectorZ(data.hitsCoordinates, data.hitsCylinderID);
        data.hitsCoordinates = sortedHits.first;
        data.hitsCylinderID = sortedHits.second;
        nTurns = PTTALG::CountTurns(data.hitsCoordinates);
        nCylinders = PTTALG::CountCylinders(data.hitsCylinderID);

        // Apply turn analysis
        if(Config::get().turnMode)
        {
            if(Config::get().processSingle)
            {
                cout << ">>> nTurns = " << nTurns << endl;
                cout << ">>> nCylinders = " << nCylinders << endl;
            }

            auto turnHits = PTTALG::SelectTurn(
                Config::get().turnID, data.hitsCoordinates, data.hitsCylinderID);
            data.hitsCoordinates = turnHits.first;
            data.hitsCylinderID = turnHits.second;

            nHits = data.hitsCoordinates.size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(Config::get().processSingle)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                inAcceptance--;

                data.accTheta->Fill(false, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
                data.accPhi->Fill(false, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

                continue;
            }
        }

        // Track is in acceptance!
        data.accTheta->Fill(true, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
        data.accPhi->Fill(true, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

        // Resolution of detectors
        const CHeT::Resolutions hitCov;

        // Fill the candidate
        vector<TVector3> measuredCoordinates;

        RVecD r_1, r_2;
        RVecD w;
        RVec<RVecD> V(2 * nHits, RVecD(2 * nHits));

        for(Int_t i = 0; i < nHits; i++)
        {
            TVector3 hitCoords;

            // Measurements
            vector<Double_t> measuredCoords;
            if(Config::get().useSmearing)
                measuredCoords
                    = AUXALG::SmearMeasurement(data.hitsCylinderID[i], data.hitsCoordinates[i]);
            else
                measuredCoords = data.hitsCoordinates[i];

            hitCoords[0] = measuredCoords.at(0);
            hitCoords[1] = measuredCoords.at(1);
            hitCoords[2] = measuredCoords.at(2);

            measuredCoordinates.push_back(hitCoords);

            Double_t uu = measuredCoords.at(0);
            Double_t vv = measuredCoords.at(1);
            Double_t phi_global = atan2(vv, uu);

            TMatrixDSym cov_xy
                = hitCov.GetMatrixCartesian(data.hitsCylinderID[i], phi_global).GetSub(0, 1, 0, 1);

            TMatrixDSym cov_rphiz = hitCov.GetMatrixCylindrical(data.hitsCylinderID[i]);

            // Fill m, V and w
            r_1.push_back(uu);
            r_2.push_back(vv);
            w.push_back(1. / cov_rphiz(1, 1));

            for(auto j = 0; j < 2; j++)
                for(auto k = 0; k < 2; k++)
                    V[nHits * j + i][nHits * k + i] = cov_xy(j, k);
        }

        // ... Centering ...
        // u -= Mean(u);
        // v -= Mean(v);
        RVecD m_c = Concatenate(r_1, r_2);

        // ... Scaling ...
        // Double_t b = 0.5;
        // auto q = Dot(m_c, m_c);
        // auto Q = sqrt(q/nHits);
        // auto m_cs = m_c * b/Q;

        TMatrixD V_11(nHits, nHits), V_12(nHits, nHits), V_21(nHits, nHits), V_22(nHits, nHits);

        for(auto i = 0; i < nHits; i++)
            for(auto j = 0; j < nHits; j++)
            {
                V_11(i, j) = V[i][j];
                V_12(i, j) = V[i][nHits + j];
                V_21(i, j) = V[nHits + i][j];
                V_22(i, j) = V[nHits + i][nHits + j];
            }

        // ... Mapping ...
        RVecD r_3 = r_1 * r_1 + r_2 * r_2;

        RVecD r = Concatenate(m_c, r_3);

        // C Matrix
        map<pair<Int_t, Int_t>, TMatrixD> C;
        for(auto i = 1; i <= 3; ++i)
            for(auto j = 1; j <= 3; ++j)
                C.insert({ { i, j }, TMatrixD(nHits, nHits) });

        C[{ 1, 1 }] = V_11;
        C[{ 1, 2 }] = V_12;
        C[{ 2, 1 }] = V_21;
        C[{ 2, 2 }] = V_22;

        // C_13, C_23
        TMatrixD C13(nHits, nHits), C23(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
            {
                C13(i, j) = 2 * V_11(i, j) * r_1[j] + 2 * V_12(i, j) * r_2[j];
                C23(i, j) = 2 * V_21(i, j) * r_1[j] + 2 * V_22(i, j) * r_2[j];
            }

        C[{ 1, 3 }] = C13;
        C[{ 2, 3 }] = C23;
        C[{ 3, 1 }] = TMatrixD(TMatrixD::kTransposed, C13);
        C[{ 3, 2 }] = TMatrixD(TMatrixD::kTransposed, C23);

        // C_33
        TMatrixD C33(nHits, nHits);
        for(auto i = 0; i < 2; ++i)
            for(auto j = 0; j < 2; ++j)
            {
                const TMatrixD &Vii = (i == 0 ? V_11 : V_22);
                const TMatrixD &Vij = (i == 0 && j == 0) ? V_11
                    : (i == 0 && j == 1)                 ? V_12
                    : (i == 1 && j == 0)                 ? V_21
                                                         : V_22;

                const RVecD &ri = (i == 0 ? r_1 : r_2);
                const RVecD &rj = (j == 0 ? r_1 : r_2);

                for(auto m = 0; m < nHits; ++m)
                    for(auto n = 0; n < nHits; ++n)
                        C33(m, n) += 2 * Vii(m, n) * Vij(m, n) + 4 * Vij(m, n) * ri[m] * rj[n];
            }

        C[{ 3, 3 }] = C33;

        // ... Center of gravity ...
        w /= Sum(w);
        TVectorD w_vec(w.size());
        for(auto i = 0; i < w.size(); ++i)
            w_vec(i) = w[i];

        TMatrixD r_mat(nHits, 3);
        for(auto j = 0; j < 3; ++j)
            for(auto i = 0; i < nHits; ++i)
                r_mat(i, j) = r[j * nHits + i];

        TMatrixD r_mat_tr(TMatrixD::kTransposed, r_mat);
        TVectorD r_0 = r_mat_tr * w_vec;

        // Var(r_0)
        TMatrixD C_0(3, 3);

        for(auto i = 1; i <= 3; ++i)
        {
            for(auto j = 1; j <= 3; ++j)
            {
                const TMatrixD &Cij = C[{ i, j }];
                C_0(i - 1, j - 1) = Cij.Similarity(w_vec);
            }
        }

        // ... Substract ...
        TMatrixD H(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
                H(i, j) = (i == j ? 1 : 0) - w_vec(j);

        // s and D Matrix
        TMatrixD s = H * r_mat;

        map<Int_t, TVectorD> s_map;
        for(auto i = 1; i <= 3; ++i)
            s_map.insert({ i, TVectorD(nHits) });

        for(auto i = 0; i < nHits; ++i)
        {
            s_map[1](i) = s(i, 0);
            s_map[2](i) = s(i, 1);
            s_map[3](i) = s(i, 2);
        }

        map<pair<Int_t, Int_t>, TMatrixD> D;
        for(auto i = 1; i <= 3; ++i)
            for(auto j = 1; j <= 3; ++j)
                D.insert({ { i, j }, TMatrixD(nHits, nHits) });

        for(auto i = 1; i <= 3; ++i)
        {
            for(auto j = 1; j <= 3; ++j)
            {
                TMatrixD &Dij = D[{ i, j }];
                const TMatrixD &Cij = C[{ i, j }];
                TMatrixD H_tr(TMatrixD::kTransposed, H);

                Dij = H * Cij * H_tr;
            }
        }

        // ... Computation of weighted sample covariance matrix 𝑨 ...
        map<Int_t, pair<Int_t, Int_t>> nu = { { 1, { 1, 1 } }, { 2, { 1, 2 } }, { 3, { 1, 3 } },
            { 4, { 2, 2 } }, { 5, { 2, 3 } }, { 6, { 3, 3 } } };

        map<Int_t, Double_t> A_alpha;
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            TVectorD w_sj(nHits);
            for(auto k = 0; k < w.size(); ++k)
                w_sj[k] = w_vec(k) * s_map[j](k); // w ⊙ s_j

            A_alpha[alpha] = s_map[i] * w_sj; // s_iᵀ ⋅ (w ⊙ s_j)
        }

        // Compute covariance matrix of A E_{alpha,beta} for alpha, beta = 1..6
        TMatrixDSym E(6);

        // Precompute W2 = w * w^T (outer product)
        TMatrixD W2(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
                W2(i, j) = w_vec(i) * w_vec(j);

        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            const auto &si = s_map[i];
            const auto &sj = s_map[j];

            for(auto beta = alpha; beta <= 6; ++beta)
            {
                auto [k, l] = nu[beta];
                const auto &sk = s_map[k];
                const auto &sl = s_map[l];

                const TMatrixD &Dik = D[{ i, k }];
                const TMatrixD &Dil = D[{ i, l }];
                const TMatrixD &Djk = D[{ j, k }];
                const TMatrixD &Djl = D[{ j, l }];

                // Hadamard products
                Double_t S_term = 0.;

                for(auto qi = 0; qi < nHits; ++qi)
                    for(auto qj = 0; qj < nHits; ++qj)
                        S_term += (Dik(qi, qj) * W2(qi, qj) * Djl(qi, qj)
                            + Dil(qi, qj) * W2(qi, qj) * Djk(qi, qj));

                TMatrixD temp_jl2(nHits, nHits), temp_jk2(nHits, nHits), temp_il2(nHits, nHits),
                    temp_ik2(nHits, nHits);
                for(auto qu = 0; qu < nHits; ++qu)
                    for(auto qv = 0; qv < nHits; ++qv)
                    {
                        temp_jl2(qu, qv) = Djl(qu, qv) * W2(qu, qv);
                        temp_jk2(qu, qv) = Djk(qu, qv) * W2(qu, qv);
                        temp_il2(qu, qv) = Dil(qu, qv) * W2(qu, qv);
                        temp_ik2(qu, qv) = Dik(qu, qv) * W2(qu, qv);
                    }

                Double_t scalar = si * (temp_jl2 * sk) + si * (temp_jk2 * sl) + sj * (temp_il2 * sk)
                    + sj * (temp_ik2 * sl);

                // Finally
                E(alpha - 1, beta - 1) = S_term + scalar;

                if(alpha != beta)
                    E(beta - 1, alpha - 1) = E(alpha - 1, beta - 1); // symmetry
            }
        }

        // E.Print();

        // ... Computation of n and c
        // A matrix
        TMatrixDSym A(3);
        A.Zero();
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            A(i - 1, j - 1) = A_alpha[alpha];
            A(j - 1, i - 1) = A_alpha[alpha];
        }

        // Diagonalizzazione
        TVectorD eigenVals(3);
        TMatrixD eigenVecs = A.EigenVectors(eigenVals);

        // Normale = autovettore con autovalore minimo
        Int_t minIdx = (eigenVals(0) < eigenVals(1)) ? ((eigenVals(0) < eigenVals(2)) ? 0 : 2)
                                                     : ((eigenVals(1) < eigenVals(2)) ? 1 : 2);

        TVectorD n(3);
        for(Int_t i = 0; i < 3; ++i)
            n[i] = eigenVecs(i, minIdx);

        Double_t c = -(n * r_0);

        // Compute the Jacobian
        TMatrixD J2(3, 6);
        const Double_t epsilon = 1e-2;

        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            // Copia A_alpha e perturba il solo alpha-esimo parametro
            auto A_plus = A_alpha;
            auto A_minus = A_alpha;

            A_plus[alpha] += epsilon;
            A_minus[alpha] -= epsilon;

            auto calc_n = [&](const map<Int_t, Double_t> &A_mod) -> TVectorD
            {
                TMatrixDSym A_mat(3);
                A_mat.Zero();
                for(auto a = 1; a <= 6; ++a)
                {
                    auto [i, j] = nu[a];
                    A_mat(i - 1, j - 1) = A_mod.at(a);
                    A_mat(j - 1, i - 1) = A_mod.at(a);
                }

                TVectorD evals(3);
                TMatrixD evecs = A_mat.EigenVectors(evals);

                Int_t minIdx = (evals(0) < evals(1)) ? ((evals(0) < evals(2)) ? 0 : 2)
                                                     : ((evals(1) < evals(2)) ? 1 : 2);

                TVectorD n(3);

                for(auto i = 0; i < 3; ++i)
                    n[i] = evecs(i, minIdx);

                return n;
            };

            TVectorD n_plus = calc_n(A_plus);
            TVectorD n_minus = calc_n(A_minus);

            for(auto i = 0; i < 3; ++i)
                J2(i, alpha - 1) = (n_plus[i] - n_minus[i]) / (2. * epsilon);
        }

        // TMatrixDSym C_n(3, 3);
        auto C_n = E.Similarity(J2); // J2 * E * J2ᵀ

        // Joint covariance matrix of n and c
        // Parte alta-sinistra: Cn
        TMatrixD Cnc(4, 4);
        for(auto i = 0; i < 3; ++i)
            for(auto j = 0; j < 3; ++j)
                Cnc(i, j) = C_n(i, j);

        // Parte in alto a destra: -Cn * r0
        TVectorD Cn_r0(3);
        Cn_r0 = C_n * r_0;

        for(auto i = 0; i < 3; ++i)
        {
            Cnc(i, 3) = -Cn_r0[i]; // colonna finale
            Cnc(3, i) = -Cn_r0[i]; // riga finale (simmetrico)
        }

        // Calcolo var[c]
        Double_t ncnrcr = C_0.Similarity(n) + C_n.Similarity(r_0);
        Double_t S_trace = 0.;
        for(auto i = 0; i < 3; ++i)
            for(auto j = 0; j < 3; ++j)
                S_trace += C_n(i, j) * C_0(i, j); // Hadamard product

        Double_t var_c = ncnrcr + S_trace;
        Cnc(3, 3) = var_c;

        // ... Circle parameters ...
        const Double_t xC = -n(0) / (2 * n(2));
        const Double_t yC = -n(1) / (2 * n(2));
        Double_t r2 = (1 - n(2) * n(2) - 4 * c * n(2)) / (4 * n(2) * n(2));
        const Double_t R = sqrt(r2);

        // Jacobian
        Double_t h = sqrt(1 - n(2) * n(2) - 4 * c * n(2));
        TMatrixD J3(3, 4);
        J3.Zero();
        J3(0, 0) = -1. / (2 * n(2));
        J3(0, 2) = n(0) / (2 * n(2) * n(2));
        J3(1, 1) = -1. / (2 * n(2));
        J3(1, 2) = n(1) / (2 * n(2) * n(2));
        J3(2, 2) = -h / (2 * n(2) * n(2)) - (4 * c + 2 * n(2)) / (4 * h * n(2));
        J3(2, 3) = -1. / h;

        // Covariance matrix
        TMatrixD J3_tr = TMatrixD(TMatrixD::kTransposed, J3);
        TMatrixD covCircle(3, 3);
        covCircle = J3 * Cnc * J3_tr;

        if(Config::get().processSingle)
            cout << Form("\nxC = %f +/- %f\nyC = %f +/- %f\nR = %f +/- %f\n\n", xC,
                sqrt(covCircle(0, 0)), yC, sqrt(covCircle(1, 1)), R, sqrt(covCircle(2, 2)));

        // --- Fit helix in Z vs arc length s ---
        vector<Double_t> s_values;
        vector<Double_t> z_values;

        // Choose the pivot
        const TVector3 pivot = measuredCoordinates[0];

        // Compute phi0 and dr
        const Double_t phi0 = atan2(pivot.Y() - yC, pivot.X() - xC);
        const Double_t cos_phi0 = cos(phi0);
        const Double_t sin_phi0 = sin(phi0);
        const Double_t dr = sqrt(pow(pivot.X() - xC, 2) + pow(pivot.Y() - yC, 2)) - R;

        // Compute arc lengths s
        Double_t previous_phi = phi0;
        Double_t previous_s = 0;

        for(const auto &point : measuredCoordinates)
        {
            Double_t dx = point.X() - xC;
            Double_t dy = point.Y() - yC;
            Double_t phi = atan2(dy, dx);

            // Angle unwrapping
            Double_t dphi = phi - previous_phi;
            if(dphi > TMath::Pi())
                dphi -= 2 * TMath::Pi();
            if(dphi < -TMath::Pi())
                dphi += 2 * TMath::Pi();
            // dphi *= -1;

            Double_t s = previous_s + R * abs(dphi);

            s_values.push_back(s);
            z_values.push_back(point.Z());

            previous_phi = phi;
            previous_s = s;
        }

        // --- Fit z vs s ---
        auto *graphZvsS = new TGraphErrors(s_values.size());
        for(size_t i = 0; i < s_values.size(); ++i)
        {
            Double_t s = s_values[i];
            Double_t z = z_values[i];

            // Get the hit
            const auto &hit = measuredCoordinates[i];
            const TMatrixDSym &matCov
                = hitCov.GetMatrixCartesian(data.hitsCylinderID[i], atan2(hit.Y(), hit.X()));

            // Uncertainties on Z
            Double_t sigmaZ = sqrt(matCov(2, 2));

            graphZvsS->SetPoint(i, s, z);
            graphZvsS->SetPointError(i, 0, sigmaZ);
        }

        TF1 *fitZvsS = new TF1("fitZvsS", "[0] + x*[1]", -10, 10);
        fitZvsS->SetParNames("z0", "tan(lambda)");
        fitZvsS->SetParameters(0, 1);
        TFitResultPtr fitlinePtr = nullptr;
        if(Config::get().processSingle)
            fitlinePtr = graphZvsS->Fit(fitZvsS, "S");
        else
            fitlinePtr = graphZvsS->Fit(fitZvsS, "SQ");

        // Get results
        Double_t z0 = fitZvsS->GetParameter(0);
        Double_t tanLambda = fitZvsS->GetParameter(1);

        TMatrixDSym covLine = fitlinePtr->GetCovarianceMatrix();

        // Helix Fit result
        Bool_t isFitConverged = fitlinePtr->IsValid();
        if(isFitConverged)
        {
            inEfficiency++;

            // Construct the full params vector
            TVectorD parHelix(5);
            parHelix(0) = xC;
            parHelix(1) = yC;
            parHelix(2) = R;
            parHelix(3) = z0;
            parHelix(4) = tanLambda;

            // Construct the full covariance matrix
            TMatrixDSym covHelix(5);
            for(Int_t i = 0; i < 3; ++i)
                for(Int_t j = 0; j < 3; ++j)
                    covHelix(i, j) = covCircle(i, j);
            for(Int_t i = 0; i < 2; ++i)
                for(Int_t j = 0; j < 2; ++j)
                    covHelix(i + 3, j + 3) = covLine(i, j);

            // Fit info and plots
            if(Config::get().processSingle)
            {
                // --- Output parameters ---
                cout << "\n\nHelix parameters:";
                parHelix.Print();
                cout << "Phi0 = " << phi0 << endl;

                // --- Covariance matrix ---
                covHelix.Print();

                // Convert hits
                vector<vector<Double_t>> plottedCoordsVec;
                for(const auto &vec : measuredCoordinates)
                    plottedCoordsVec.push_back({ vec.X(), vec.Y(), vec.Z() });

                // Canvas
                TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "XY View", 700, 700);
                TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ", "YZ View", 900, 500);
                TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ", "3D Helix Fit", 800, 600);
                TCanvas *canvZvsS = new TCanvas("canvZvsS", "z vs s", 700, 500);

                AUXALG::DrawXYView_hits(&data.trueDecayPos, plottedCoordsVec, canvHitsXY);
                AUXALG::DrawYZView_hits(&data.trueDecayPos, plottedCoordsVec, canvHitsYZ);
                AUXALG::DrawXYZView_hits(plottedCoordsVec, canvHitsXYZ);

                AUXALG::DrawXYView_arc(
                    xC, yC, R, plottedCoordsVec, canvHitsXY, nTurns, Config::get().turnID);
                AUXALG::DrawZvsSFit(graphZvsS, fitZvsS, canvZvsS);
                AUXALG::DrawXYZView_helixFromHits(
                    xC, yC, R, z0, phi0, tanLambda, plottedCoordsVec, canvHitsXYZ);
            }

            // Extrapolate status at decay vertex
            // z_vertex = z0 + s * tanLambda
            Double_t z_vertex = data.trueDecayPos.Z();
            Double_t s_at_vertex = (z_vertex - z0) / tanLambda;

            // phi(s) = phi0 + h * (s/R)
            Double_t delta_phi = s_at_vertex / R;
            Double_t phi_at_vertex = phi0 + muEDM::h * delta_phi;

            // Vertex decay
            Double_t x_at_vertex = xC + R * cos(phi_at_vertex);
            Double_t y_at_vertex = yC + R * sin(phi_at_vertex);
            // Double_t x_at_z0 = pivot.X() - dr * cos_phi0 + R * (cos(phi0 - phi_at_z0) -
            // cos_phi0); Double_t y_at_z0 = pivot.Y() - dr * sin_phi0 + R * (sin(phi0 - phi_at_z0)
            // - sin_phi0);
            // Vertex momentum
            const Double_t k = 2.99792458; // MeV/c * T * cm
            TVector3 pFitted(
                -muEDM::h * sin(phi_at_vertex), muEDM::h * cos(phi_at_vertex), tanLambda);
            pFitted *= k * 0.1 * abs(muEDM::Fields::constBz) * R;

            // Uncertainties
            TMatrixD Jac = ANS::ComputeHelixJacobian(
                pivot, xC, yC, R, z0, tanLambda, 0.1 * muEDM::Fields::constBz);
            TMatrixDSym covFittedState = covHelix.Similarity(Jac);
            covFittedState = ANS::CovFromCardinalToCylindricalMom(covFittedState, pFitted);

            // Compute angles with muEDM convention
            Double_t pFittedPhi = pFitted.Phi();
            Double_t pFitted_z = cos(pFitted.Theta());
            Double_t pFitted_r = sin(pFitted.Theta())
                * (x_at_vertex * cos(pFittedPhi) + y_at_vertex * sin(pFittedPhi))
                / sqrt(x_at_vertex * x_at_vertex + y_at_vertex * y_at_vertex);
            Double_t pFittedTheTheta = atan2(pFitted_z, pFitted_r);

            // Fitted final extrapolated state
            // x, y, z, p, TheTheta, phi
            TVectorD fittedState(6);
            fittedState(0) = x_at_vertex;
            fittedState(1) = y_at_vertex;
            fittedState(2) = 0.;
            fittedState(3) = pFitted.Mag();
            fittedState(4) = pFittedTheTheta;
            fittedState(5) = pFittedPhi;

            // Print results
            if(Config::get().processSingle)
            {
                cout << "\nFitted vertex position: (" << fittedState(0) << ", " << fittedState(1)
                     << ", 0) cm" << endl;
                cout << endl;
                cout << Form(
                    "Fidded p = (%.2f, %.2f, %.2f) MeV/c", pFitted.X(), pFitted.Y(), pFitted.Z())
                     << endl;
                cout << Form("(p, thetaEDM, phi) = (%.2f MeV/c, %.2f pi rad, %.2f pi rad)",
                    fittedState(3), fittedState(4) / TMath::Pi(), fittedState(5) / TMath::Pi());
                cout << endl;

                covFittedState.Print();
            }

            // Fill histos
            auto [fRes, fSigma, fPulls] = AUXALG::GetResults(fittedState, covFittedState,
                data.trueDecayPos, data.trueDecayMom.Mag(), data.trueDecayThetaEDM,
                data.trueDecayMom.Phi(), NORMALIZED_PULLS);
            if(fRes.size() == 0)
                continue;

            data.histDiffX->Fill(fPulls[0]);
            data.histDiffY->Fill(fPulls[1]);
            data.histDiffZ->Fill(fPulls[2]);
            data.histDiffMom->Fill(fPulls[3]);
            data.histDiffTheta->Fill(fPulls[4]);
            data.histDiffPhi->Fill(fPulls[5]);

            data.graphMom->Fill(data.trueDecayMom.Mag(), fRes[3]);
            data.graphTheta->Fill(data.trueDecayThetaEDM, fRes[4]);
            data.graphPhi->Fill(data.trueDecayMom.Phi(), fRes[5]);

            data.hist2MomRes->Fill(data.trueDecayMom.Mag(), fSigma[3]);
            data.hist2ThetaRes->Fill(data.trueDecayThetaEDM, fSigma[4]);
            data.hist2PhiRes->Fill(data.trueDecayMom.Phi(), fSigma[5]);

            data.profMomRes->Fill(data.trueDecayMom.Mag(), fSigma[3]);
            data.profThetaRes->Fill(data.trueDecayThetaEDM, fSigma[4]);
            data.profPhiRes->Fill(data.trueDecayMom.Phi(), fSigma[5]);
        }

        if(Config::get().processSingle)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n"
                 << endl;
        }

        // Track is in efficiency?
        data.effTheta->Fill(isFitConverged, data.trueDecayMom.Mag(), data.trueDecayThetaEDM);
        data.effPhi->Fill(isFitConverged, data.trueDecayMom.Mag(), data.trueDecayMom.Phi());

        // Store turns and cylinders data
        data.histTurns->Fill(nTurns);
        data.effTurns->Fill(isFitConverged, nTurns);
        data.histCylinders->Fill(nCylinders);
        data.effCylinders->Fill(isFitConverged, nCylinders);
        data.histTurnsVMom->Fill(data.trueDecayMom.Mag(), nTurns);
        data.histCylVMom->Fill(data.trueDecayMom.Mag(), nCylinders);

        cout << "\r>>> Processed event number " << (ev - Config::get().rangeLoop.first) << flush;
    }

    cout << endl;

    // Print results and draw graphs
    if(!Config::get().processSingle)
    {
        // Recap
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t)(inAcceptance * 100) / processedEvents << endl;
        cout << ">>> Efficiency = " << (Float_t)(inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;

        // Graphs
        if(Config::get().quietMode)
            gROOT->SetBatch(true);

        TCanvas *canvEfficiency = new TCanvas("canvEfficiency");
        canvEfficiency->Divide(2, 2);
        canvEfficiency->cd(1);
        data.accTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(2);
        data.accPhi->Draw("COLZ TEXT");
        canvEfficiency->cd(3);
        data.effTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(4);
        data.effPhi->Draw("COLZ TEXT");

        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        data.histTurns->Draw();
        canvTurns->cd(2);
        data.effTurns->Draw("AP");

        TCanvas *canvCylinders = new TCanvas("canvCylinders");
        canvCylinders->Divide(2, 2);
        canvCylinders->cd(1);
        data.histCylinders->Draw();
        canvCylinders->cd(2);
        data.effCylinders->Draw("AP");
        canvCylinders->cd(3);
        data.histTurnsVMom->Draw();
        canvCylinders->cd(4);
        data.histCylVMom->Draw();

        TCanvas *canvLinearity = new TCanvas("canvLinearity");
        canvLinearity->Divide(3);
        canvLinearity->cd(1);
        data.graphMom->SetMarkerStyle(20);
        data.graphMom->Draw("SCAT");
        canvLinearity->cd(2);
        data.graphTheta->SetMarkerStyle(20);
        data.graphTheta->Draw("SCAT");
        canvLinearity->cd(3);
        data.graphPhi->SetMarkerStyle(20);
        data.graphPhi->Draw("SCAT");

        TCanvas *canvfPullsPos = new TCanvas("Pulls Position");
        canvfPullsPos->Divide(3);

        canvfPullsPos->cd(1);
        data.histDiffX->Draw();

        canvfPullsPos->cd(2);
        data.histDiffY->Draw();

        canvfPullsPos->cd(3);
        data.histDiffZ->Draw();

        TCanvas *canvfPullsMom = new TCanvas("Pulls Momentum");
        canvfPullsMom->Divide(3);

        canvfPullsMom->cd(1);
        data.histDiffMom->Draw();

        canvfPullsMom->cd(2);
        data.histDiffTheta->Draw();

        canvfPullsMom->cd(3);
        data.histDiffPhi->Draw();

        TCanvas *canvProfRes = new TCanvas("ProfResolutions");
        canvProfRes->Divide(3);
        canvProfRes->cd(1);
        data.profMomRes->Draw();
        canvProfRes->cd(2);
        data.profThetaRes->Draw();
        canvProfRes->cd(3);
        data.profPhiRes->Draw();

        TCanvas *canvResMom = new TCanvas("canvResMom");
        canvResMom->cd();
        data.hist2MomRes->Draw();

        TCanvas *canvResTheta = new TCanvas("canvResTheta");
        canvResTheta->cd();
        data.hist2ThetaRes->Draw();

        TCanvas *canvResPhi = new TCanvas("canvResPhi");
        canvResPhi->cd();
        data.hist2PhiRes->Draw();

        if(Config::get().quietMode)
        {
            canvEfficiency->SaveAs("canvEfficiency.pdf");
            canvResMom->SaveAs("canvResMom.pdf");
            canvResTheta->SaveAs("canvResTheta.pdf");
            canvResPhi->SaveAs("canvResPhi.pdf");
        }
    }

    // Trick
    display->open();

    // Finally
    data.SaveRecoTree();
    exit(0);
}

array<Double_t, 6> FITALG::HelixPrefitter(
    const vector<vector<Double_t>> &hitsCoordinates, const vector<Int_t> &cylinders)
{
    // Return parameters: xC, yC, R, phi0, z0, tanLambda
    array<Double_t, 6> fParameters = { 0, 0, 0, 0, 0, 0 };

    // Resolution of detectors
    const CHeT::Resolutions hitCov;

    // Fill the candidate
    vector<TVector3> measuredCoordinates;
    Int_t nHits = hitsCoordinates.size();

    RVecD r_1, r_2;
    RVecD w;
    RVec<RVecD> V(2 * nHits, RVecD(2 * nHits));

    for(Int_t i = 0; i < nHits; i++)
    {
        TVector3 hitCoords;

        // Measurements
        vector<Double_t> measuredCoords = hitsCoordinates[i];

        hitCoords[0] = measuredCoords.at(0);
        hitCoords[1] = measuredCoords.at(1);
        hitCoords[2] = measuredCoords.at(2);

        measuredCoordinates.push_back(hitCoords);

        Double_t uu = measuredCoords.at(0);
        Double_t vv = measuredCoords.at(1);
        Double_t phi_global = atan2(vv, uu);

        TMatrixDSym cov_xy = hitCov.GetMatrixCartesian(cylinders[i], phi_global).GetSub(0, 1, 0, 1);

        TMatrixDSym cov_rphiz = hitCov.GetMatrixCylindrical(cylinders[i]);

        // Fill m, V and w
        r_1.push_back(uu);
        r_2.push_back(vv);
        w.push_back(1. / cov_rphiz(1, 1));

        for(auto j = 0; j < 2; j++)
            for(auto k = 0; k < 2; k++)
                V[nHits * j + i][nHits * k + i] = cov_xy(j, k);
    }

    // --- RIEMANN CIRCLE FIT ---
    RVecD m_c = Concatenate(r_1, r_2);

    TMatrixD V_11(nHits, nHits), V_12(nHits, nHits), V_21(nHits, nHits), V_22(nHits, nHits);

    for(auto i = 0; i < nHits; i++)
        for(auto j = 0; j < nHits; j++)
        {
            V_11(i, j) = V[i][j];
            V_12(i, j) = V[i][nHits + j];
            V_21(i, j) = V[nHits + i][j];
            V_22(i, j) = V[nHits + i][nHits + j];
        }

    // ... Mapping ...
    RVecD r_3 = r_1 * r_1 + r_2 * r_2;

    RVecD r = Concatenate(m_c, r_3);

    // C Matrix
    map<pair<Int_t, Int_t>, TMatrixD> C;
    for(auto i = 1; i <= 3; ++i)
        for(auto j = 1; j <= 3; ++j)
            C.insert({ { i, j }, TMatrixD(nHits, nHits) });

    C[{ 1, 1 }] = V_11;
    C[{ 1, 2 }] = V_12;
    C[{ 2, 1 }] = V_21;
    C[{ 2, 2 }] = V_22;

    // C_13, C_23
    TMatrixD C13(nHits, nHits), C23(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
        {
            C13(i, j) = 2 * V_11(i, j) * r_1[j] + 2 * V_12(i, j) * r_2[j];
            C23(i, j) = 2 * V_21(i, j) * r_1[j] + 2 * V_22(i, j) * r_2[j];
        }

    C[{ 1, 3 }] = C13;
    C[{ 2, 3 }] = C23;
    C[{ 3, 1 }] = TMatrixD(TMatrixD::kTransposed, C13);
    C[{ 3, 2 }] = TMatrixD(TMatrixD::kTransposed, C23);

    // C_33
    TMatrixD C33(nHits, nHits);
    for(auto i = 0; i < 2; ++i)
        for(auto j = 0; j < 2; ++j)
        {
            const TMatrixD &Vii = (i == 0 ? V_11 : V_22);
            const TMatrixD &Vij = (i == 0 && j == 0) ? V_11
                : (i == 0 && j == 1)                 ? V_12
                : (i == 1 && j == 0)                 ? V_21
                                                     : V_22;

            const RVecD &ri = (i == 0 ? r_1 : r_2);
            const RVecD &rj = (j == 0 ? r_1 : r_2);

            for(auto m = 0; m < nHits; ++m)
                for(auto n = 0; n < nHits; ++n)
                    C33(m, n) += 2 * Vii(m, n) * Vij(m, n) + 4 * Vij(m, n) * ri[m] * rj[n];
        }

    C[{ 3, 3 }] = C33;

    // ... Center of gravity ...
    w /= Sum(w);
    TVectorD w_vec(w.size());
    for(auto i = 0; i < w.size(); ++i)
        w_vec(i) = w[i];

    TMatrixD r_mat(nHits, 3);
    for(auto j = 0; j < 3; ++j)
        for(auto i = 0; i < nHits; ++i)
            r_mat(i, j) = r[j * nHits + i];

    TMatrixD r_mat_tr(TMatrixD::kTransposed, r_mat);
    TVectorD r_0 = r_mat_tr * w_vec;

    // Var(r_0)
    TMatrixD C_0(3, 3);

    for(auto i = 1; i <= 3; ++i)
    {
        for(auto j = 1; j <= 3; ++j)
        {
            const TMatrixD &Cij = C[{ i, j }];
            C_0(i - 1, j - 1) = Cij.Similarity(w_vec);
        }
    }

    // ... Substract ...
    TMatrixD H(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
            H(i, j) = (i == j ? 1 : 0) - w_vec(j);

    // s and D Matrix
    TMatrixD s = H * r_mat;

    map<Int_t, TVectorD> s_map;
    for(auto i = 1; i <= 3; ++i)
        s_map.insert({ i, TVectorD(nHits) });

    for(auto i = 0; i < nHits; ++i)
    {
        s_map[1](i) = s(i, 0);
        s_map[2](i) = s(i, 1);
        s_map[3](i) = s(i, 2);
    }

    map<pair<Int_t, Int_t>, TMatrixD> D;
    for(auto i = 1; i <= 3; ++i)
        for(auto j = 1; j <= 3; ++j)
            D.insert({ { i, j }, TMatrixD(nHits, nHits) });

    for(auto i = 1; i <= 3; ++i)
    {
        for(auto j = 1; j <= 3; ++j)
        {
            TMatrixD &Dij = D[{ i, j }];
            const TMatrixD &Cij = C[{ i, j }];
            TMatrixD H_tr(TMatrixD::kTransposed, H);

            Dij = H * Cij * H_tr;
        }
    }

    // ... Computation of weighted sample covariance matrix 𝑨 ...
    map<Int_t, pair<Int_t, Int_t>> nu = { { 1, { 1, 1 } }, { 2, { 1, 2 } }, { 3, { 1, 3 } },
        { 4, { 2, 2 } }, { 5, { 2, 3 } }, { 6, { 3, 3 } } };

    map<Int_t, Double_t> A_alpha;
    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        TVectorD w_sj(nHits);
        for(auto k = 0; k < w.size(); ++k)
            w_sj[k] = w_vec(k) * s_map[j](k); // w ⊙ s_j

        A_alpha[alpha] = s_map[i] * w_sj; // s_iᵀ ⋅ (w ⊙ s_j)
    }

    // Compute covariance matrix of A E_{alpha,beta} for alpha, beta = 1..6
    TMatrixDSym E(6);

    // Precompute W2 = w * w^T (outer product)
    TMatrixD W2(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
            W2(i, j) = w_vec(i) * w_vec(j);

    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        const auto &si = s_map[i];
        const auto &sj = s_map[j];

        for(auto beta = alpha; beta <= 6; ++beta)
        {
            auto [k, l] = nu[beta];
            const auto &sk = s_map[k];
            const auto &sl = s_map[l];

            const TMatrixD &Dik = D[{ i, k }];
            const TMatrixD &Dil = D[{ i, l }];
            const TMatrixD &Djk = D[{ j, k }];
            const TMatrixD &Djl = D[{ j, l }];

            // Hadamard products
            Double_t S_term = 0.;

            for(auto qi = 0; qi < nHits; ++qi)
                for(auto qj = 0; qj < nHits; ++qj)
                    S_term += (Dik(qi, qj) * W2(qi, qj) * Djl(qi, qj)
                        + Dil(qi, qj) * W2(qi, qj) * Djk(qi, qj));

            TMatrixD temp_jl2(nHits, nHits), temp_jk2(nHits, nHits), temp_il2(nHits, nHits),
                temp_ik2(nHits, nHits);
            for(auto qu = 0; qu < nHits; ++qu)
                for(auto qv = 0; qv < nHits; ++qv)
                {
                    temp_jl2(qu, qv) = Djl(qu, qv) * W2(qu, qv);
                    temp_jk2(qu, qv) = Djk(qu, qv) * W2(qu, qv);
                    temp_il2(qu, qv) = Dil(qu, qv) * W2(qu, qv);
                    temp_ik2(qu, qv) = Dik(qu, qv) * W2(qu, qv);
                }

            Double_t scalar = si * (temp_jl2 * sk) + si * (temp_jk2 * sl) + sj * (temp_il2 * sk)
                + sj * (temp_ik2 * sl);

            // Finally
            E(alpha - 1, beta - 1) = S_term + scalar;

            if(alpha != beta)
                E(beta - 1, alpha - 1) = E(alpha - 1, beta - 1); // symmetry
        }
    }

    // ... Computation of n and c
    // A matrix
    TMatrixDSym A(3);
    A.Zero();
    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        A(i - 1, j - 1) = A_alpha[alpha];
        A(j - 1, i - 1) = A_alpha[alpha];
    }

    // Diagonalizzazione
    TVectorD eigenVals(3);
    TMatrixD eigenVecs = A.EigenVectors(eigenVals);

    // Normale = autovettore con autovalore minimo
    Int_t minIdx = (eigenVals(0) < eigenVals(1)) ? ((eigenVals(0) < eigenVals(2)) ? 0 : 2)
                                                 : ((eigenVals(1) < eigenVals(2)) ? 1 : 2);

    TVectorD n(3);
    for(Int_t i = 0; i < 3; ++i)
        n[i] = eigenVecs(i, minIdx);

    Double_t c = -(n * r_0);

    // Compute the Jacobian
    TMatrixD J2(3, 6);
    const Double_t epsilon = 1e-2;

    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        // Copia A_alpha e perturba il solo alpha-esimo parametro
        auto A_plus = A_alpha;
        auto A_minus = A_alpha;

        A_plus[alpha] += epsilon;
        A_minus[alpha] -= epsilon;

        auto calc_n = [&](const map<Int_t, Double_t> &A_mod) -> TVectorD
        {
            TMatrixDSym A_mat(3);
            A_mat.Zero();
            for(auto a = 1; a <= 6; ++a)
            {
                auto [i, j] = nu[a];
                A_mat(i - 1, j - 1) = A_mod.at(a);
                A_mat(j - 1, i - 1) = A_mod.at(a);
            }

            TVectorD evals(3);
            TMatrixD evecs = A_mat.EigenVectors(evals);

            Int_t minIdx = (evals(0) < evals(1)) ? ((evals(0) < evals(2)) ? 0 : 2)
                                                 : ((evals(1) < evals(2)) ? 1 : 2);

            TVectorD n(3);

            for(auto i = 0; i < 3; ++i)
                n[i] = evecs(i, minIdx);

            return n;
        };

        TVectorD n_plus = calc_n(A_plus);
        TVectorD n_minus = calc_n(A_minus);

        for(auto i = 0; i < 3; ++i)
            J2(i, alpha - 1) = (n_plus[i] - n_minus[i]) / (2. * epsilon);
    }

    // TMatrixDSym C_n(3, 3);
    auto C_n = E.Similarity(J2); // J2 * E * J2ᵀ

    // Joint covariance matrix of n and c
    // Parte alta-sinistra: Cn
    TMatrixD Cnc(4, 4);
    for(auto i = 0; i < 3; ++i)
        for(auto j = 0; j < 3; ++j)
            Cnc(i, j) = C_n(i, j);

    // Parte in alto a destra: -Cn * r0
    TVectorD Cn_r0(3);
    Cn_r0 = C_n * r_0;

    for(auto i = 0; i < 3; ++i)
    {
        Cnc(i, 3) = -Cn_r0[i]; // colonna finale
        Cnc(3, i) = -Cn_r0[i]; // riga finale (simmetrico)
    }

    // Calcolo var[c]
    Double_t ncnrcr = C_0.Similarity(n) + C_n.Similarity(r_0);
    Double_t S_trace = 0.;
    for(auto i = 0; i < 3; ++i)
        for(auto j = 0; j < 3; ++j)
            S_trace += C_n(i, j) * C_0(i, j); // Hadamard product

    Double_t var_c = ncnrcr + S_trace;
    Cnc(3, 3) = var_c;

    // ... Circle parameters ...
    const Double_t xC = -n(0) / (2 * n(2));
    const Double_t yC = -n(1) / (2 * n(2));
    Double_t r2 = (1 - n(2) * n(2) - 4 * c * n(2)) / (4 * n(2) * n(2));
    const Double_t R = sqrt(r2);

    // Jacobian
    Double_t h = sqrt(1 - n(2) * n(2) - 4 * c * n(2));
    TMatrixD J3(3, 4);
    J3.Zero();
    J3(0, 0) = -1. / (2 * n(2));
    J3(0, 2) = n(0) / (2 * n(2) * n(2));
    J3(1, 1) = -1. / (2 * n(2));
    J3(1, 2) = n(1) / (2 * n(2) * n(2));
    J3(2, 2) = -h / (2 * n(2) * n(2)) - (4 * c + 2 * n(2)) / (4 * h * n(2));
    J3(2, 3) = -1. / h;

    // Covariance matrix
    TMatrixD J3_tr = TMatrixD(TMatrixD::kTransposed, J3);
    TMatrixD covCircle(3, 3);
    covCircle = J3 * Cnc * J3_tr;

    if(Config::get().processSingle)
        cout << Form("\n>>> Prefitter =\nxC = %f +/- %f\nyC = %f +/- %f\nR = %f +/- %f\n\n", xC,
            sqrt(covCircle(0, 0)), yC, sqrt(covCircle(1, 1)), R, sqrt(covCircle(2, 2)));

    // --- Z VS S FIT ---
    vector<Double_t> s_values;
    vector<Double_t> z_values;

    // Choose the pivot
    const TVector3 pivot = measuredCoordinates[0];

    // Compute phi0 and dr
    const Double_t phi0 = atan2(pivot.Y() - yC, pivot.X() - xC);

    // Compute arc lengths s
    Double_t previous_phi = phi0;
    Double_t previous_s = 0;

    for(const auto &point : measuredCoordinates)
    {
        Double_t dx = point.X() - xC;
        Double_t dy = point.Y() - yC;
        Double_t phi = atan2(dy, dx);

        // Angle unwrapping
        Double_t dphi = phi - previous_phi;
        if(dphi > TMath::Pi())
            dphi -= 2 * TMath::Pi();
        if(dphi < -TMath::Pi())
            dphi += 2 * TMath::Pi();

        Double_t s = previous_s + R * abs(dphi);

        s_values.push_back(s);
        z_values.push_back(point.Z());

        previous_phi = phi;
        previous_s = s;
    }

    // --- Fit z vs s ---
    auto *graphZvsS = new TGraphErrors(s_values.size());
    for(size_t i = 0; i < s_values.size(); ++i)
    {
        Double_t s = s_values[i];
        Double_t z = z_values[i];

        // Get the hit
        const auto &hit = measuredCoordinates[i];
        const TMatrixDSym &matCov
            = hitCov.GetMatrixCartesian(cylinders[i], atan2(hit.Y(), hit.X()));

        // Uncertainties on Z
        Double_t sigmaZ = sqrt(matCov(2, 2));

        graphZvsS->SetPoint(i, s, z);
        graphZvsS->SetPointError(i, 0, sigmaZ);
    }

    TF1 *fitZvsS = new TF1("fitZvsS", "[0] + x*[1]", -10, 10);
    fitZvsS->SetParNames("z0", "tan(lambda)");
    fitZvsS->SetParameters(0, 1);

    TFitResultPtr fitlinePtr = nullptr;
    if(Config::get().processSingle)
        fitlinePtr = graphZvsS->Fit(fitZvsS, "S");
    else
        fitlinePtr = graphZvsS->Fit(fitZvsS, "SQ");

    // Get results
    Double_t z0 = fitZvsS->GetParameter(0);
    Double_t tanLambda = fitZvsS->GetParameter(1);

    TMatrixDSym covLine = fitlinePtr->GetCovarianceMatrix();

    // Helix Fit result
    // Construct the full params vector
    fParameters = { xC, yC, R, phi0, z0, tanLambda };

    // Construct the full covariance matrix
    TMatrixDSym covHelix(5);
    for(Int_t i = 0; i < 3; ++i)
        for(Int_t j = 0; j < 3; ++j)
            covHelix(i, j) = covCircle(i, j);
    for(Int_t i = 0; i < 2; ++i)
        for(Int_t j = 0; j < 2; ++j)
            covHelix(i + 3, j + 3) = covLine(i, j);

    delete graphZvsS;
    delete fitZvsS;

    return fParameters;
}

void FITALG::AddFakeHitFromHelix(genfit::TrackCand &trackCand, Int_t hitIndex,
    Double_t sortingParameter, Double_t s, Double_t xC, Double_t yC, Double_t R, Double_t z0,
    Double_t phi0, Double_t tanLambda, TClonesArray &chetHitArray,
    vector<TVector3> &virtualCoordinates, Double_t sigmaBig)
{
    // phi(s) = phi0 + h * (s / R)
    Double_t angle_shift = s / R;
    Double_t current_phi = phi0 + muEDM::h * angle_shift;

    // Coordinates
    const Double_t x = xC + R * cos(current_phi);
    const Double_t y = yC + R * sin(current_phi);
    const Double_t z = z0 + s * tanLambda;

    TVector3 fakeHit(x, y, z);

    TMatrixDSym bigCov(3);
    bigCov.Zero();
    bigCov(0, 0) = sigmaBig * sigmaBig;
    bigCov(1, 1) = sigmaBig * sigmaBig;
    bigCov(2, 2) = sigmaBig * sigmaBig;

    virtualCoordinates.push_back(fakeHit);

    new(chetHitArray[hitIndex]) genfit::mySpacepointDetectorHit(fakeHit, bigCov);
    trackCand.addHit(0, hitIndex, -1, sortingParameter);
}

// fitTrack->addTrackRep(repHelix);

// for(auto i = 0; i < nHits; ++i)
//{
//     auto tP = fitTrack->getPoint(i);
//     auto kFI = new genfit::KalmanFitterInfo(tP, rep);
//     //kFI->setRefenceState();
//     tP->setFitterInfo(kFI);
//     cout << tP->getKalmanFitterInfo() << endl;
// }
