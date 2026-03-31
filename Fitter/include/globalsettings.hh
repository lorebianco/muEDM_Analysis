#ifndef GLOBALSETTINGS_HH
#define GLOBALSETTINGS_HH

#include <AbsBFieldDer.h>
#include <CHeT/CHeTGlobalSettings.hh>
#include <ConstField.h>
#include <FieldManager.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TROOT.h>
#include <TVector3.h>

namespace muEDM
{
extern Double_t h;

namespace Fields
{
enum class FieldBType
{
    MultiMap,
    SingleMap,
    Constant
};

// Config variables
extern FieldBType currentBType;
extern std::string basePath;
extern Double_t globalBScale;
extern Bool_t useRadial;

// Const field
extern Double_t constBz;

// Struct of map config
struct MapConfig
{
    std::string filename;
    Double_t scaleFactor;
};

// --- Methods ---
void Configure(FieldBType type, std::string path, Double_t scale, Double_t bZ = 0.0);

// Setup
void SetupBField();
} // namespace Fields
} // namespace muEDM

namespace CHeT
{

// Resolutions [cm] - Note that CHeT Library uses mm, so we divide by 10.
struct Resolutions
{
    Double_t sigmaPitch = (Config::FIBERS_PER_SIPM * Config::FIBER_WIDTH / 10.0) / sqrt(12);
    Double_t sigmaR = (2. * Config::FIBER_WIDTH / 10.0) / sqrt(12);

    // Fitting tricks
    Double_t scaleCov = 1.;

    // Matrix
    // (I'm using r, Rphi, z)
    inline TMatrixDSym GetMatrixCylindrical(Int_t cylID) const
    {
        TMatrixDSym C_RphiZ(3);
        const Double_t k = 0.5 * sigmaPitch * sigmaPitch;

        const auto &cyls = Config::GetCylinders();
        Double_t stereoAngle = 0.0;
        if(cylID >= 0 && cylID < cyls.size())
        {
            stereoAngle
                = (cyls[cylID].inner.nominalStereoAngle + cyls[cylID].outer.nominalStereoAngle)
                / 2.0;
        }

        C_RphiZ(0, 0) = sigmaR * sigmaR;
        C_RphiZ(0, 1) = 0.;
        C_RphiZ(0, 2) = 0.;
        C_RphiZ(1, 0) = 0.;
        C_RphiZ(1, 1) = k * pow(cos(stereoAngle), -2.);
        C_RphiZ(1, 2) = 0.;
        C_RphiZ(2, 0) = 0.;
        C_RphiZ(2, 1) = 0.;
        C_RphiZ(2, 2) = k * pow(sin(stereoAngle), -2.);

        // C_RphiZ.Print();
        return scaleCov * scaleCov * C_RphiZ;
    };

    inline TMatrixDSym GetMatrixCartesian(Int_t cylID, Double_t phi) const
    {
        TMatrixDSym C_xyz = GetMatrixCylindrical(cylID);

        // std::cout << ">>> C_RphiZ = " << std::endl;
        // C_xyz.Print();

        // Define the Jacobian matrix J
        TMatrixD J(3, 3);

        J(0, 0) = cos(phi);
        J(0, 1) = -sin(phi);
        J(0, 2) = 0;
        J(1, 0) = sin(phi);
        J(1, 1) = cos(phi);
        J(1, 2) = 0;
        J(2, 0) = 0;
        J(2, 1) = 0;
        J(2, 2) = 1;

        // Compute transformed covariance: C_xyz = J * C_RphiZ * J^T
        C_xyz.Similarity(J);

        // std::cout << ">>> C_xyz = " << std::endl;
        // C_xyz.Print();
        return C_xyz;
    };
};
}; // namespace CHeT

namespace ANS
{
// Analysis tools global values
// Detector response
constexpr Int_t smearNFibers = 4;

// Pattern recognition response
constexpr Double_t sigmaSeedPos = 0.1; // cm
constexpr Double_t sigmaSeedMomMag = 0.1; // relative
constexpr Double_t sigmaSeedMomDir = 3. / 180. * TMath::Pi(); // rad

TMatrixDSym CovFromCardinalToCylindricalMom(TMatrixDSym cov, TVector3 mom);
TMatrixD ComputeHelixJacobian(const TVector3 &pivot, Double_t xC, Double_t yC, Double_t R,
    Double_t z0, Double_t tanLambda, Double_t B = 2.2);

}; // namespace ANS

#endif // GLOBALSETTINGS_HH
