//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2bDetectorConstruction.cc
/// \brief Implementation of the B2bDetectorConstruction class
 
#include "B2bDetectorConstruction.hh"
#include "B2bDetectorMessenger.hh"
#include "B2bChamberParameterisation.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B2bDetectorConstruction::fMagFieldMessenger = 0;
 
B2bDetectorConstruction::B2bDetectorConstruction()
:G4VUserDetectorConstruction(),
 fLogicTarget(NULL), labLV(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL), 
 fCheckOverlaps(true)
{
  fMessenger = new B2bDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2bDetectorConstruction::~B2bDetectorConstruction()
{
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2bDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");
  fShieldMaterial = nistManager->FindOrBuildMaterial("G4_WATER");
  // LAB material
  G4bool isotopes = false;
  G4Element* H = nistManager->FindOrBuildElement("H", isotopes);
  G4Element* C = nistManager->FindOrBuildElement("C", isotopes);
  G4double density = 0.863*g/cm3;
  fChamberMaterial = new G4Material("LAB", density, 2);
  fChamberMaterial->AddElement(H, 30);
  fChamberMaterial->AddElement(C, 18);
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2bDetectorConstruction::DefineVolumes()
{
  G4Material* air  = G4Material::GetMaterial("G4_AIR");
  G4Material* pmtMaterial = G4Material::GetMaterial("G4_Fe");
  // Sizes of the principal geometrical components (solids)
  
  G4int NbOfChambers = 5;
  G4double chamberSpacing = 80*cm; // from chamber center to center!

  G4double chamberWidth = 20.0*cm; // width of the chambers
  G4double targetLength =  5.0*cm; // full length of Target
  
  G4double trackerLength = (NbOfChambers+1)*chamberSpacing;

  G4double labRadius = 7*m;
  G4double waterRadius = 8*m;
  G4double waterHeight = 2*waterRadius;
  G4double pbRadius1 = 8.1*m;
  G4double pbRadius2 = 8.2*m; 
  G4double pbHeight = waterHeight+20*cm;
  G4double worldLength = 2.5 * 2 * pbRadius2;


  G4double targetRadius  = 0.5*targetLength;   // Radius of Target
  targetLength = 0.5*targetLength;             // Half length of the Target  
  G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldLength/2,worldLength/2,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // pb sheild

  G4Tubs* targetS
    = new G4Tubs("target",pbRadius1,pbRadius2,pbHeight/2,0.*deg,360.*deg);
  fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
  new G4PVPlacement(0,              // no rotation
                    G4ThreeVector(), // at (x,y,z)
                    fLogicTarget,   // its logical volume
                    "Target",       // its name
                    worldLV,        // its mother volume
                    false,          // no boolean operations
                    0,              // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << 2*pbHeight/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;
  // water sheild

  G4Tubs* waterS
    = new G4Tubs("watersheild",labRadius,waterRadius,waterHeight/2,0.*deg,360.*deg);
  G4LogicalVolume* fLogicWater
    = new G4LogicalVolume(waterS, fShieldMaterial,"Watershield",0,0,0);
  
  new G4PVPlacement(0,              // no rotation
                    G4ThreeVector(), // at (x,y,z)
                    fLogicWater,   // its logical volume
                    "WaterPv",       // its name
                    worldLV,        // its mother volume
                    false,          // no boolean operations
                    0,              // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Water is " << 2*waterHeight/cm << " cm of "
         << fShieldMaterial->GetName() << G4endl;

  // Tracker
 
  G4ThreeVector positionLab = G4ThreeVector(0,0,0);

  G4Sphere* labS
    = new G4Sphere("lab",0, labRadius, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
  labLV
    = new G4LogicalVolume(labS, fChamberMaterial, "Lab",0,0,0);  
  new G4PVPlacement(0,               // no rotation
                    positionLab, // at (x,y,z)
                    labLV,       // its logical volume
                    "Lab",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  // Tracker segments
  
  // An example of Parameterised volumes
  // Dummy values for G4Tubs -- modified by parameterised volume

  G4Tubs* chamberS
    = new G4Tubs("tracker",0, 100*cm, 100*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* pmtLV 
    = new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber",0,0,0);
  
  G4double firstPosition = -trackerSize + chamberSpacing;
  G4double firstLength   = trackerLength/10;
  G4double lastLength    = trackerLength;

  G4VPVParameterisation* chamberParam =
    new B2bChamberParameterisation(
                                  NbOfChambers,   // NoChambers
                                  firstPosition,  // Z of center of first
                                  chamberSpacing, // Z spacing of centers
                                  chamberWidth,  // chamber width
                                  firstLength,    // initial length 
                                  lastLength);    // final length
                           
  // dummy value : kZAxis -- modified by parameterised volume
/*
  new G4PVParameterised("Chamber",       // their name
                        labLV,   // their logical volume
                        trackerLV,       // Mother logical volume
                        kZAxis,          // Are placed along this axis 
                        NbOfChambers,    // Number of chambers
                        chamberParam,    // The parametrisation
                        fCheckOverlaps); // checking overlaps 

  G4cout << "There are " << NbOfChambers << " chambers in the tracker region. "
         << G4endl
         << "The chambers are " << chamberWidth/cm << " cm of "
         << fChamberMaterial->GetName() << G4endl
         << "The distance between chamber is " << chamberSpacing/cm << " cm" 
         << G4endl;
*/
  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  worldLV   ->SetVisAttributes(boxVisAtt);  
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  fLogicWater ->SetVisAttributes(boxVisAtt);

  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  labLV->SetVisAttributes(chamberVisAtt);
  
  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.5*chamberWidth;
  fStepLimit = new G4UserLimits(maxStep);
  labLV->SetUserLimits(fStepLimit);
  
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2bDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector( labLV,  aTrackerSD );

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        if (labLV) labLV->SetMaterial(fChamberMaterial);
        G4cout
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2bDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
