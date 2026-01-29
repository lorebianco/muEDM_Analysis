#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Activates namespace GS (GlobalSettings)
#pragma link C++ namespace GS;
#pragma link C++ class GS::LayerConfig + ;
#pragma link C++ class GS::CylinderConfig + ;
#pragma link C++ class GS::FiberProp + ;
#pragma link C++ class GS::BundlesIntersection + ;
#pragma link C++ function GS::GetCylinders;
#pragma link C++ function GS::GetFiberProp;
#pragma link C++ function GS::FindIntersections;
#pragma link C++ function GS::MapExplorer;
#pragma link C++ function GS::PrintBundleMapping;

// Activates namespace CV
#pragma link C++ namespace CV;
#pragma link C++ class CV::VisLineTrack + ;
#pragma link C++ class CV::VisPoint2D + ;
#pragma link C++ class CV::VisPoint3D + ;
#pragma link C++ function CV::Draw2D;
#pragma link C++ function CV::Draw3D;

#endif
