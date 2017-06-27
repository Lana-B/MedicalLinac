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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//

#include "ML2Ph_SiDiode.hh"
#include "G4SystemOfUnits.hh"

CML2Ph_SiDiode::CML2Ph_SiDiode()
{
	// phantom size and position
	halfSize.set(150*mm,150*mm,150*mm);
	// phantom position
	centre.set(0.,0.,0.);
}

CML2Ph_SiDiode::~CML2Ph_SiDiode(void)
{
}
void CML2Ph_SiDiode::writeInfo()
{
	std::cout<<"\n\n\tcentre of the phantom: " <<centre/mm<<" [mm]"<< G4endl;
	std::cout<<"\thalf thickness of the phantom: " <<halfSize/mm<<" [mm]\n"<< G4endl;
}
bool CML2Ph_SiDiode::Construct(G4VPhysicalVolume *PWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG)
{
	PVWorld=PWorld;

	G4double A, Z;
	A = 28.09*g/mole;
	// G4Element* elSi = new G4Element("Silicon", "Si", Z=14., A);
	G4NistManager* man = G4NistManager::Instance();

	G4Material* Si = man->FindOrBuildMaterial("G4_Si");

	bool bCreated=false;
	// G4Material *WATER=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
	// G4Box *SiDiodePhantomBox = new G4Box("SiDiodePhantomBox", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4double innerRadius = 0.*mm;
	G4double outerRadius = halfSize.getX();
	G4double hz = 0.5*mm;
	G4double startAngle = 0.*deg;
	G4double spanningAngle = 360.*deg;

	G4Tubs* siliconTube
	= new G4Tubs("siTube",
	             innerRadius, 
	             outerRadius,
	             hz,
	             startAngle, 
	             spanningAngle);

    G4double  pRmin1 = 0*mm;
    G4double  pRmax1 = 1.3*mm;
    G4double  pRmin2 = 0*mm;
    G4double  pRmax2 = 1.2*mm;
    G4double  pDz = 0.49*mm;
    G4double  pSPhi = 0.*deg;
    G4double  pDPhi = 360.*deg;


	G4Cons* emptyCone 
	= new G4Cons("subCone",
                pRmin1,
                pRmax1,
                pRmin2,
                pRmax2,
                pDz,
                pSPhi,
                pDPhi);

	G4SubtractionSolid* subtraction =
    new G4SubtractionSolid("Diode", siliconTube, emptyCone);

	G4LogicalVolume *SiDiodeLV = new G4LogicalVolume(subtraction, Si, "SiDiodeLV", 0, 0, 0);
	SiDiodePV = new G4PVPlacement(0, centre, "SiDiodePV", SiDiodeLV, PVWorld, false, 0);

	// G4LogicalVolume *SiDiodePhantomLV = new G4LogicalVolume(SiDiodePhantomBox, WATER, "SiDiodePhantomLV", 0, 0, 0);
	// SiDiodePhantomPV = new G4PVPlacement(0, centre, "SiDiodePhantomPV", SiDiodePhantomLV, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol= new G4Region("SiDiodeR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*mm);
	regVol->SetProductionCuts(cuts);

	SiDiodeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(SiDiodeLV);

	// Visibility
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt->SetVisibility(true);
	// simpleAlSVisAtt->SetForceSolid(true);
	SiDiodeLV->SetVisAttributes(simpleAlSVisAtt);

	// Sensitive detector 
	sensDet=new CML2SDWithVoxels("Silicon Diode", saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG, G4ThreeVector(0.,0.,0.), halfSize, 100, 100, 100);
	G4SDManager *SDManager=G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(sensDet);
	
	// Read Out Geometry
	CML2ReadOutGeometry *ROG = new CML2ReadOutGeometry();
	ROG->setBuildData(PVWorld->GetFrameTranslation(), halfSize, 100, 100, 100);
	ROG->BuildROGeometry();
	sensDet->SetROgeometry(ROG);
	SiDiodeLV->SetSensitiveDetector(sensDet);

	bCreated=true;
	return bCreated;
}


