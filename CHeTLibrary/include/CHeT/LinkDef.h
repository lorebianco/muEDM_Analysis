#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Activates namespace CHeT
#pragma link C++ namespace CHeT;

// Activates namespace CHeT::Config (ex GS)
#pragma link C++ namespace CHeT::Config;
#pragma link C++ class CHeT::Config::LayerConfig + ;
#pragma link C++ class CHeT::Config::CylinderConfig + ;
#pragma link C++ class CHeT::Config::FiberProp + ;
#pragma link C++ class CHeT::Config::BundlesIntersection + ;
#pragma link C++ function CHeT::Config::GetOffsetExp;
#pragma link C++ function CHeT::Config::SetOffsetExp;
#pragma link C++ function CHeT::Config::GetDelta1;
#pragma link C++ function CHeT::Config::SetDelta1;
#pragma link C++ function CHeT::Config::GetDelta2;
#pragma link C++ function CHeT::Config::SetDelta2;
#pragma link C++ function CHeT::Config::SetRotation;
#pragma link C++ function CHeT::Config::GetRotation;
#pragma link C++ function CHeT::Config::SetTranslation;
#pragma link C++ function CHeT::Config::GetTranslation;
#pragma link C++ function CHeT::Config::ApplyRotation;
#pragma link C++ function CHeT::Config::ApplyInverseRotation;
#pragma link C++ function CHeT::Config::ApplyTransformation;
#pragma link C++ function CHeT::Config::ApplyInverseTransformation;
#pragma link C++ function CHeT::Config::GetBoardGlobalOffset;
#pragma link C++ function CHeT::Config::GetCylinders;
#pragma link C++ function CHeT::Config::wrap0_2pi;
#pragma link C++ function CHeT::Config::GetFiberProp;
#pragma link C++ function CHeT::Config::GetGlobalBundleId;
#pragma link C++ function CHeT::Config::FindIntersections;
#pragma link C++ function CHeT::Config::GetGlobalIdFromGeometry;
#pragma link C++ function CHeT::Config::MapExplorer;
#pragma link C++ function CHeT::Config::PrintBundleMapping;

// Activates namespace CHeT::Vis (ex CV)
#pragma link C++ namespace CHeT::Vis;
#pragma link C++ class CHeT::Vis::VisLineTrack + ;
#pragma link C++ class CHeT::Vis::VisPoint2D + ;
#pragma link C++ class CHeT::Vis::VisPoint3D + ;
#pragma link C++ function CHeT::Vis::Draw2D;
#pragma link C++ function CHeT::Vis::Draw3D;

// Activates namespace CHeT::Data (contains Reader)
#pragma link C++ namespace CHeT::Data;
#pragma link C++ class CHeT::Data::Reader + ;

#endif
