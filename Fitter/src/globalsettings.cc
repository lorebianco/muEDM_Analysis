#include "globalsettings.hh"

using namespace std;


namespace muEDM
{
    Double_t h = 1.;

    namespace Fields
    {
        // Default values (fallback)
        FieldBType currentBType = FieldBType::MultiMap;
        string basePath = "/home/lorenzo/muEDM_Project/David_Files/Field24.5MeV/";
        Double_t globalBScale = 1.0;
        Bool_t useRadial = true;
        Double_t constBz = 28.9; // Default 2.9 T
    }
}



void muEDM::Fields::Configure(FieldBType type, string path, Double_t scale, Double_t bZ)
{
    currentBType = type;
    basePath = path;
    globalBScale = scale;
    constBz = bZ * scale;

    // Update h
    Double_t refB = (type == FieldBType::Constant) ? constBz : (28.9 * globalBScale);
    muEDM::h = -1.0 * TMath::Sign(1.0, 1. * refB);
    
    cout << ">>> [Configure] Updated configuration:" << endl;
    cout << "    Type: " << (int)type << " (0:Multi, 1:Single, 2:Const)" << endl;
    cout << "    Path: " << basePath << endl;
    cout << "    Scale: " << globalBScale << endl;
    cout << "    ConstBz (Bz*Scale): " << constBz << endl;
    cout << endl;
}



void muEDM::Fields::SetupBField()
{    
    // Const field
    if(currentBType == FieldBType::Constant)
    {
        cout << ">>> [SetupBField] Initializing CONSTANT field: (0, 0, " << constBz << ")" << endl;
        genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., constBz));
        return;
    }

    // Field Maps (Single or Multi)
    auto* bField = new genfit::AbsBFieldDer();
    string bFieldInfo = ">>> FieldMap Info: ";
    
    // Maps vector
    vector<MapConfig> mapsToLoad;

    if(currentBType == FieldBType::SingleMap)
    {
        mapsToLoad.push_back({ "Sol_CC_WF_newCCConfig/", globalBScale });
    }
    else
    {
        mapsToLoad =
        {
            { "WF/",    1.50005  * globalBScale },
            { "CCin/",  1.95894  * globalBScale },
            { "CCout/", 2.50350  * globalBScale },
            { "SOL/",   1.113636 * globalBScale }
        };
    }

    // Loop
    for(const auto& config : mapsToLoad)
    {
        string fullPath = basePath;
        if(fullPath.back() != '/') fullPath += "/";
        fullPath += config.filename;
        
        bField->addFieldPath(fullPath, "Spline_radial", config.scaleFactor, useRadial);

        bFieldInfo += "[" + config.filename + ", factor=" + to_string(config.scaleFactor) + "] ";
    }

    cout << ">>> [SetupBField] Reading field maps... " << endl;
    cout << bFieldInfo << endl;
    
    bField->readFieldMaps();
    
    cout << ">>> [SetupBField] Done.\n" << endl;

    genfit::FieldManager::getInstance()->init(bField);
}



TMatrixDSym ANS::CovFromCardinalToCylindricalMom(TMatrixDSym cov, TVector3 mom)
{
    // Transform a covariance matrix in x,y,z, momx, momy, momz
    // into a covariance matrix in x, y, z, mom, theta, phi
    TMatrixDSym covCyl(cov);

    TMatrixD Jac(6, 6);
    Jac.Zero();

    Double_t p = mom.Mag();
    Double_t pt = TMath::Hypot(mom.X(), mom.Y());

    if(p == 0 || pt == 0)
        return covCyl;

    // Calculate Jacobian
    Jac[0][0] = 1.;
    Jac[1][1] = 1.;
    Jac[2][2] = 1.;

    Jac[3][3] = mom.X() / p;
    Jac[3][4] = mom.Y() / p;
    Jac[3][5] = mom.Z() / p;

    Jac[4][3] = mom.X() * mom.Z() / p / p / pt;
    Jac[4][4] = mom.Y() * mom.Z() / p / p / pt;
    Jac[4][5] = - pt / p / p;

    Jac[5][3] = - mom.Y() / pt / pt;
    Jac[5][4] = mom.X() / pt / pt;
    Jac[5][5] = 0.;

    covCyl.Similarity(Jac);
    return covCyl;
}



TMatrixD ANS::ComputeHelixJacobian(
    const TVector3& pivot,
    Double_t xC, Double_t yC, Double_t R,
    Double_t z0, Double_t tanLambda,
    Double_t B)
{
    // Setup
    const Double_t k = 2.99792458; // MeV/c / (T * cm)
    
    Double_t x0 = pivot.X();
    Double_t y0 = pivot.Y();

    Double_t dx = x0 - xC;
    Double_t dy = y0 - yC;
    Double_t rho2 = dx*dx + dy*dy;
    Double_t rho = sqrt(rho2);

    Double_t cos_phi0 = dx / rho;
    Double_t sin_phi0 = dy / rho;
    Double_t phi0 = atan2(dy, dx);

    Double_t phi = - z0 / (R * tanLambda);
    Double_t theta = phi0 - phi;

    Double_t cos_theta = cos(theta);
    Double_t sin_theta = sin(theta);

    Double_t dtheta_dxC = dy / rho2;
    Double_t dtheta_dyC = -dx / rho2;
    Double_t dtheta_dR  = - z0 / (R*R * tanLambda);
    Double_t dtheta_ddz = 1.0 / (R * tanLambda);
    Double_t dtheta_dtanL = -z0 / (R * tanLambda * tanLambda);

    Double_t pT = k * B * R;

    // Jacobian: rows are (x, y, z, px, py, pz), cols are (xC, yC, R, z0, tanLambda)
    TMatrixD J(6, 5);

    // Column: xC
    J(0, 0) = 1 - R * sin_theta * dtheta_dxC;
    J(1, 0) =     - R * cos_theta * dtheta_dxC;
    J(2, 0) = 0;

    J(3, 0) = pT * cos_theta * dtheta_dxC;
    J(4, 0) = pT * sin_theta * dtheta_dxC;
    J(5, 0) = 0;

    // Column: yC
    J(0, 1) =     - R * sin_theta * dtheta_dyC;
    J(1, 1) = 1 - R * cos_theta * dtheta_dyC;
    J(2, 1) = 0;

    J(3, 1) = pT * cos_theta * dtheta_dyC;
    J(4, 1) = pT * sin_theta * dtheta_dyC;
    J(5, 1) = 0;

    // Column: R
    J(0, 2) = cos_theta - R * sin_theta * dtheta_dR;
    J(1, 2) = sin_theta + R * cos_theta * dtheta_dR;
    J(2, 2) = 0;

    J(3, 2) = k * B * (sin_theta + R * cos_theta * dtheta_dR);
    J(4, 2) = k * B * (-cos_theta + R * sin_theta * dtheta_dR);
    J(5, 2) = k * B * tanLambda; // d/dR [k B R tanλ] = k B tanλ

    // Column: dz
    J(0, 3) = R * sin_theta * dtheta_ddz;
    J(1, 3) = -R * cos_theta * dtheta_ddz;
    J(2, 3) = 0;
    
    J(3, 3) = pT * cos_theta * dtheta_ddz;
    J(4, 3) = pT * sin_theta * dtheta_ddz;
    J(5, 3) = 0; // No dz in tanLambda

    // Column: tanLambda
    J(0, 4) = R * sin_theta * dtheta_dtanL;
    J(1, 4) = -R * cos_theta * dtheta_dtanL;
    J(2, 4) = 0;
    
    J(3, 4) = pT * cos_theta * dtheta_dtanL;
    J(4, 4) = pT * sin_theta * dtheta_dtanL;
    J(5, 4) = k * B * R; // d/dtanλ [k B R tanλ]

    return J;
}