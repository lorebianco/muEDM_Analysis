#include "auxiliaryalgorithms.hh"

using namespace std;

void AUXALG::DrawXYView_hits(TVector3 *origin, vector<vector<Double_t>> hitsCoordinates,
    TCanvas *canvas, vector<vector<Double_t>> virtualCoordinates)
{
    canvas->cd();
    auto frame = canvas->DrawFrame(-9, -9, 9, 9);
    auto hframe = (TH1F *)gPad->GetPrimitive("hframe");
    hframe->SetLineWidth(0);
    frame->Draw();
    frame->SetTitle("X-Y View;X [cm];Y [cm]");

    Int_t nHits = hitsCoordinates.size();
    Double_t *x = new Double_t[nHits];
    Double_t *y = new Double_t[nHits];
    for(Int_t i = 0; i < nHits; i++)
    {
        x[i] = hitsCoordinates[i].at(0);
        y[i] = hitsCoordinates[i].at(1);
    }
    TGraph *gr = new TGraph(nHits, x, y);
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    Double_t x0 = origin->X();
    Double_t y0 = origin->Y();
    TGraph *ogr = new TGraph(1, &x0, &y0);
    ogr->SetMarkerStyle(20);
    ogr->SetMarkerColor(kGreen);

    auto cyls = CHeT::Config::GetCylinders();
    int nCyls = cyls.size();
    std::vector<TEllipse *> ell(nCyls);
    for(auto i = 0; i < nCyls; i++)
    {
        ell[i] = new TEllipse(0, 0, cyls[nCyls - 1 - i].inner.radius / 10.0);
        ell[i]->SetFillStyle(0);
        ell[i]->Draw("same");
    }

    gr->Draw("PC same");
    ogr->Draw("P same");

    Int_t nVirtualHits = virtualCoordinates.size();
    if(nVirtualHits > 0)
    {
        Double_t *xv = new Double_t[nVirtualHits];
        Double_t *yv = new Double_t[nVirtualHits];
        for(Int_t i = 0; i < nVirtualHits; i++)
        {
            xv[i] = virtualCoordinates[i].at(0);
            yv[i] = virtualCoordinates[i].at(1);
        }
        TGraph *grv = new TGraph(nVirtualHits, xv, yv);
        grv->SetMarkerStyle(24);
        grv->SetLineColor(kGreen);
        grv->Draw("P same");
    }

    canvas->Update();
}

void AUXALG::DrawYZView_hits(TVector3 *origin, vector<vector<Double_t>> hitsCoordinates,
    TCanvas *canvas, vector<vector<Double_t>> virtualCoordinates)
{
    canvas->cd();
    auto frame = canvas->DrawFrame(-40, -9, 40, 9);
    auto hframe = (TH1F *)gPad->GetPrimitive("hframe");
    hframe->SetLineWidth(0);
    frame->Draw();
    frame->SetTitle("Z-Y View;Z [cm]; Y [cm]");

    Int_t nHits = hitsCoordinates.size();
    Double_t *y = new Double_t[nHits];
    Double_t *z = new Double_t[nHits];
    for(Int_t i = 0; i < nHits; i++)
    {
        y[i] = hitsCoordinates[i].at(1);
        z[i] = hitsCoordinates[i].at(2);
    }
    TGraph *gr = new TGraph(nHits, z, y);
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    Double_t z0 = origin->Z();
    Double_t y0 = origin->Y();
    TGraph *ogr = new TGraph(1, &z0, &y0);
    ogr->SetMarkerStyle(20);
    ogr->SetMarkerColor(kGreen);

    auto cyls = CHeT::Config::GetCylinders();
    int nCyls = cyls.size();
    double length_cm = CHeT::Config::L_HALF / 10.0;
    std::vector<TBox *> box(nCyls);
    for(auto i = 0; i < nCyls; i++)
    {
        double r_cm = cyls[nCyls - 1 - i].inner.radius / 10.0;
        box[i] = new TBox(-length_cm, -r_cm, length_cm, r_cm);
        box[i]->SetFillStyle(0);
        box[i]->SetLineColor(kBlack);
        box[i]->SetFillStyle(0);
        box[i]->Draw("same");
    }

    gr->Draw("PC same");
    ogr->Draw("P same");

    Int_t nVirtualHits = virtualCoordinates.size();
    if(nVirtualHits > 0)
    {
        Double_t *zv = new Double_t[nVirtualHits];
        Double_t *yv = new Double_t[nVirtualHits];
        for(Int_t i = 0; i < nVirtualHits; i++)
        {
            zv[i] = virtualCoordinates[i].at(2);
            yv[i] = virtualCoordinates[i].at(1);
        }
        TGraph *grv = new TGraph(nVirtualHits, zv, yv);
        grv->SetMarkerStyle(24);
        grv->SetLineColor(kGreen);
        grv->Draw("P same");
    }

    canvas->Update();
}

void AUXALG::DrawXYZView_hits(vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas,
    vector<vector<Double_t>> virtualCoordinates)
{
    canvas->cd();

    Int_t nHits = hitsCoordinates.size();
    Double_t *x = new Double_t[nHits];
    Double_t *y = new Double_t[nHits];
    Double_t *z = new Double_t[nHits];

    for(Int_t i = 0; i < nHits; ++i)
    {
        x[i] = hitsCoordinates[i].at(0);
        y[i] = hitsCoordinates[i].at(1);
        z[i] = hitsCoordinates[i].at(2);
    }

    TGraph2D *gr = new TGraph2D(nHits, z, x, y);
    gr->SetTitle("3D View; Z [cm]; X [cm]; Y [cm]");
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    gr->Draw("P LINE");

    Int_t nVirtualHits = virtualCoordinates.size();
    if(nVirtualHits > 0)
    {
        Double_t *xv = new Double_t[nVirtualHits];
        Double_t *yv = new Double_t[nVirtualHits];
        Double_t *zv = new Double_t[nVirtualHits];

        for(Int_t i = 0; i < nVirtualHits; ++i)
        {
            xv[i] = virtualCoordinates[i].at(0);
            yv[i] = virtualCoordinates[i].at(1);
            zv[i] = virtualCoordinates[i].at(2);
        }

        TGraph2D *grv = new TGraph2D(nVirtualHits, zv, xv, yv);
        grv->SetMarkerStyle(24);
        grv->SetLineColor(kGreen);

        grv->Draw("P SAME");
    }

    canvas->Update();
}

void AUXALG::DrawXYView_arc(Double_t xC, Double_t yC, Double_t R,
    const vector<vector<Double_t>> &hitsCoordinates, TCanvas *canvas, Int_t nTurns, Float_t turnID)
{

    TVector3 first(
        hitsCoordinates.front()[0], hitsCoordinates.front()[1], hitsCoordinates.front()[2]);

    TVector3 last(hitsCoordinates.back()[0], hitsCoordinates.back()[1], hitsCoordinates.back()[2]);

    // Compute phi start and end
    Double_t phi1 = atan2(first.Y() - yC, first.X() - xC);
    Double_t phi2 = atan2(last.Y() - yC, last.X() - xC);

    // Clockwise
    Double_t dphi = phi2 - phi1;

    // Unwrapping standard
    while(dphi > TMath::Pi())
        dphi -= TMath::TwoPi();
    while(dphi < -TMath::Pi())
        dphi += TMath::TwoPi();
    // if(dphi > 0)
    // dphi -= TMath::TwoPi();

    if(muEDM::h < 0 && dphi > 0)
        dphi -= TMath::TwoPi();
    if(muEDM::h > 0 && dphi < 0)
        dphi += TMath::TwoPi();

    // Full turn
    if(nTurns > 0 && turnID >= 1)
    {
        phi1 = 0;
        dphi = muEDM::h * TMath::TwoPi();
    }

    // Use TPolyline for arc
    const Int_t nPoints = 100;
    TPolyLine *arc = new TPolyLine(nPoints);
    for(Int_t i = 0; i < nPoints; ++i)
    {
        Double_t t = (Double_t)i / (nPoints - 1);
        Double_t phi = phi1 + t * dphi;
        Double_t x = xC + R * cos(phi);
        Double_t y = yC + R * sin(phi);
        arc->SetPoint(i, x, y);
    }

    arc->SetLineColor(kRed);
    arc->SetLineWidth(2);
    canvas->cd();
    arc->Draw("L same");
    canvas->Update();
}

void AUXALG::DrawXYZView_helixFromHits(Double_t xC, Double_t yC, Double_t R, Double_t z0,
    Double_t phi0, Double_t tanLambda, vector<vector<Double_t>> &hitsCoordinates, TCanvas *canvas)
{
    // Cumulative arc length computation (Ricostruiamo s esattamente come nel fit)
    vector<Double_t> s_cumulative;
    Double_t previous_phi
        = phi0; // O l'angolo del primo hit, meglio phi0 del fit se vuoi plottare il fit

    // Per plottare il fit esattamente sui punti, ricalcoliamo s dai punti
    // Nota: assumo che phi0 passato sia quello del fit (al PCA),
    // ma qui calcolo s relativo ai punti.
    // Per coerenza visiva veloce, usiamo l'approccio differenziale sui punti:

    // Prendiamo il primo punto come riferimento per s=0 visivo relativo ai punti
    Double_t prev_phi_pt = atan2(hitsCoordinates[0][1] - yC, hitsCoordinates[0][0] - xC);
    Double_t current_s = 0;
    s_cumulative.push_back(0);

    for(size_t i = 1; i < hitsCoordinates.size(); ++i)
    {
        Double_t dx = hitsCoordinates[i][0] - xC;
        Double_t dy = hitsCoordinates[i][1] - yC;
        Double_t phi = atan2(dy, dx);

        Double_t dphi = phi - prev_phi_pt;
        while(dphi > M_PI)
            dphi -= 2 * M_PI;
        while(dphi < -M_PI)
            dphi += 2 * M_PI;

        // s cresce sempre
        current_s += R * TMath::Abs(dphi);
        s_cumulative.push_back(current_s);

        prev_phi_pt = phi;
    }

    // Range di s da plottare
    Double_t s_min = s_cumulative.front();
    Double_t s_max = s_cumulative.back();

    // Generazione punti fit
    vector<vector<Double_t>> helixPoints3D;

    // ATTENZIONE: phi0 è definito a z0.
    // Dobbiamo trovare l'offset di s per il primo punto rispetto a z0.
    // z = z0 + s * tanLambda.
    // I punti misurati hanno un loro Z.
    // Possiamo usare z per pilotare s se tanLambda è ben definito,
    // oppure usare la s geometrica calcolata sopra.
    // Usiamo la formula canonica dell'elica parametrizzata in s rispetto a (z0, phi0).

    // Dobbiamo mappare il range [s_min, s_max] dei punti relativi al primo hit
    // al sistema di coordinate del fit (dove s=0 è a z=z0).
    // Assumiamo che il fit sia buono: z_hit ~ z0 + s_hit_abs * tanL
    // Quindi s_parametro = (z_hit_start - z0)/tanL + s_relativo

    // APPROCCIO PIÙ SEMPLICE PER PLOTTARE LA LINEA ROSSA SUI PUNTI:
    // Usiamo Z come variabile indipendente se tanLambda != 0
    Double_t z_start = hitsCoordinates.front()[2];
    Double_t z_end = hitsCoordinates.back()[2];

    Int_t nSteps = 200;
    Double_t z_step = (z_end - z_start) / nSteps;

    for(int i = 0; i <= nSteps; ++i)
    {
        Double_t z = z_start + i * z_step;

        // Invertiamo: s = (z - z0) / tanLambda
        Double_t s_val = (z - z0) / tanLambda;

        // Calcoliamo phi: phi(s) = phi0 + h * (s/R)
        Double_t phi = phi0 + muEDM::h * (s_val / R);

        Double_t x = xC + R * cos(phi);
        Double_t y = yC + R * sin(phi);

        helixPoints3D.push_back({ x, y, z });
    }

    // Disegno
    const int nP = helixPoints3D.size();
    auto graph3D = new TPolyLine3D(nP);
    for(int i = 0; i < nP; ++i)
        graph3D->SetPoint(i, helixPoints3D[i][2], helixPoints3D[i][0], helixPoints3D[i][1]);

    graph3D->SetLineColor(kRed);
    graph3D->SetLineWidth(2);

    canvas->cd();
    graph3D->Draw("same");
}

void AUXALG::DrawZvsSFit(TGraph *graph, TF1 *fitFunc, TCanvas *canvas)
{
    if(!graph || !fitFunc || !canvas)
        return;

    canvas->cd();

    graph->SetTitle("Z vs arc length s; s [cm]; Z [cm]");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);
    graph->Draw("AP");

    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");

    canvas->Update();
}

tuple<vector<Double_t>, vector<Double_t>, vector<Double_t>> AUXALG::GetResults(
    genfit::Track *fitTrack, genfit::AbsTrackRep *rep, TVector3 truePos, Double_t trueMom,
    Double_t trueTheta, Double_t truePhi, Bool_t normalizedPulls)
{
    try
    {
        // Extrapolate to orbit
        const genfit::MeasuredStateOnPlane &stFirst = fitTrack->getFittedState();
        TVector3 posProj;
        TVector3 momProj;
        TMatrixDSym covProj;

        stFirst.getPosMomCov(posProj, momProj, covProj);
        genfit::MeasuredStateOnPlane stateOrbit(rep);
        rep->setPosMomCov(stateOrbit, posProj, momProj, covProj);
        rep->extrapolateToPlane(stateOrbit,
            genfit::SharedPlanePtr(new genfit::DetPlane(
                TVector3(0., 0., truePos.Z()), TVector3(1, 0, 0), TVector3(0, 1, 0))));
        stateOrbit.getPosMomCov(posProj, momProj, covProj);
        covProj = ANS::CovFromCardinalToCylindricalMom(covProj, momProj);

        // Compute angles
        Double_t x = posProj.X();
        Double_t y = posProj.Y();
        Double_t theta = momProj.Theta();
        Double_t momProjPhi = momProj.Phi();

        Double_t pz = TMath::Cos(theta);
        Double_t pr = TMath::Sin(theta) * (x * TMath::Cos(momProjPhi) + y * TMath::Sin(momProjPhi))
            / TMath::Sqrt(x * x + y * y);

        Double_t momProjTheta
            = TMath::ATan2(pz, pr); // angle in plane (e_r, z) in radiants, in (-pi, pi)

        // Pulls:
        Double_t dX, dY, dZ, dMom, dTheta, dPhi;
        if(normalizedPulls)
        {
            // Position
            dX = (posProj.X() - truePos.X()) / sqrt(covProj(0, 0));
            dY = (posProj.Y() - truePos.Y()) / sqrt(covProj(1, 1));
            dZ = (posProj.Z() - truePos.Z()) / sqrt(covProj(2, 2));
            // Momentum
            dMom = (momProj.Mag() * 1E3 - trueMom) / (sqrt(covProj(3, 3)) * 1E3);
            dTheta = (momProjTheta - trueTheta) / sqrt(covProj(4, 4));
            dPhi = TMath::ATan2(sin(momProjPhi - truePhi), cos(momProjPhi - truePhi))
                / sqrt(covProj(5, 5));
        }
        else
        {
            // Position
            dX = posProj.X() - truePos.X();
            dY = posProj.Y() - truePos.Y();
            dZ = posProj.Z() - truePos.Z();
            // Momentum
            dMom = momProj.Mag() * 1E3 - trueMom;
            dTheta = momProjTheta - trueTheta;
            dPhi = TMath::ATan2(sin(momProjPhi - truePhi), cos(momProjPhi - truePhi));
        }

        return make_tuple(vector<Double_t> { posProj.X(), posProj.Y(), posProj.Z(),
                              momProj.Mag() * 1E3, momProjTheta, momProjPhi },
            vector<Double_t> { sqrt(covProj(0, 0)), sqrt(covProj(1, 1)), sqrt(covProj(2, 2)),
                sqrt(covProj(3, 3)) * 1E3, sqrt(covProj(4, 4)), sqrt(covProj(5, 5)) },
            vector<Double_t> { dX, dY, dZ, dMom, dTheta, dPhi });
    }
    catch(genfit::Exception &e)
    {
        cerr << "Exception, next track" << endl;
        cerr << e.what();
    }

    return make_tuple(vector<Double_t>(), vector<Double_t>(), vector<Double_t>());
}

tuple<vector<Double_t>, vector<Double_t>, vector<Double_t>> AUXALG::GetResults(
    TVectorD &fittedState, TMatrixDSym &covFittedState, TVector3 truePos, Double_t trueMom,
    Double_t trueTheta, Double_t truePhi, Bool_t normalizedPulls)
{
    // Pulls:
    Double_t dX, dY, dZ, dMom, dTheta, dPhi;
    if(normalizedPulls)
    {
        // Position
        dX = (fittedState(0) - truePos.X()) / sqrt(covFittedState(0, 0));
        dY = (fittedState(1) - truePos.Y()) / sqrt(covFittedState(1, 1));
        dZ = (fittedState(2) - truePos.Z()) / sqrt(covFittedState(2, 2));
        // Momentum
        dMom = (fittedState(3) - trueMom) / (sqrt(covFittedState(3, 3)));
        dTheta = remainder(fittedState(4) - trueTheta, TMath::TwoPi()) / sqrt(covFittedState(4, 4));
        dPhi = TMath::ATan2(sin(fittedState(5) - truePhi), cos(fittedState(5) - truePhi))
            / sqrt(covFittedState(5, 5));
    }
    else
    {
        // Position
        dX = fittedState(0) - truePos.X();
        dY = fittedState(1) - truePos.Y();
        dZ = fittedState(2) - truePos.Z();
        // Momentum
        dMom = fittedState(3) - trueMom;
        dTheta = remainder(fittedState(4) - trueTheta, TMath::TwoPi());
        dPhi = TMath::ATan2(sin(fittedState(5) - truePhi), cos(fittedState(5) - truePhi));
    }

    return make_tuple(vector<Double_t> { fittedState(0), fittedState(1), fittedState(2),
                          fittedState(3), fittedState(4), fittedState(5) },
        vector<Double_t> { sqrt(covFittedState(0, 0)), sqrt(covFittedState(1, 1)),
            sqrt(covFittedState(2, 2)), sqrt(covFittedState(3, 3)), sqrt(covFittedState(4, 4)),
            sqrt(covFittedState(5, 5)) },
        vector<Double_t> { dX, dY, dZ, dMom, dTheta, dPhi });
}

vector<Double_t> AUXALG::SmearMeasurement(Int_t cylID, vector<Double_t> hitCoords)
{
    if(cylID < 0 || cylID >= 7)
        throw out_of_range("Invalid cylinder ID");

    Double_t x = hitCoords.at(0);
    Double_t y = hitCoords.at(1);
    Double_t z = hitCoords.at(2);

    auto cyls = CHeT::Config::GetCylinders();
    Double_t r_nominal = cyls[cylID].inner.radius / 10.0;
    Double_t stAngle
        = (cyls[cylID].inner.nominalStereoAngle + cyls[cylID].outer.nominalStereoAngle) / 2.0;
    Double_t fiberW = CHeT::Config::FIBER_WIDTH / 10.0;

    Double_t phi = TMath::ATan2(y, x);
    Double_t rphi = r_nominal * phi;

    // --- Radial smearing ---
    Double_t r = gRandom->Gaus(r_nominal, (2. * fiberW) / sqrt(12));

    // --- Angular and longitudinal smearing ---
    TMatrixDSym A(2);
    A(0, 0) = cos(stAngle);
    A(0, 1) = sin(stAngle);
    A(1, 0) = cos(-stAngle);
    A(1, 1) = sin(-stAngle);

    TMatrixD A_U(2, 2);
    A_U(0, 0) = A(0, 0);
    A_U(0, 1) = A(0, 1);
    A_U(1, 0) = 0.;
    A_U(1, 1) = 0.;

    TMatrixD A_L(2, 2);
    A_L(0, 0) = 0.;
    A_L(0, 1) = 0.;
    A_L(1, 0) = A(1, 0);
    A_L(1, 1) = A(1, 1);

    Double_t sigmaPitch = (ANS::smearNFibers * fiberW) / sqrt(12);
    Double_t G_1 = gRandom->Gaus(0., 1.);
    Double_t G_2 = gRandom->Gaus(0., 1.);

    TVectorD hit_1(2);
    TVectorD hit_2(2);
    hit_1(0) = rphi + sigmaPitch * G_1 * cos(stAngle);
    hit_1(1) = z + sigmaPitch * G_1 * sin(stAngle);

    hit_2(0) = rphi + sigmaPitch * G_2 * cos(-stAngle);
    hit_2(1) = z + sigmaPitch * G_2 * sin(-stAngle);

    A.Invert();
    TVectorD hitSmeared = A * A_U * hit_1 + A * A_L * hit_2;

    // --- Back to cartesian ---
    x = r * cos(hitSmeared(0) / r_nominal);
    y = r * sin(hitSmeared(0) / r_nominal);
    z = hitSmeared(1);

    return { x, y, z };
}
