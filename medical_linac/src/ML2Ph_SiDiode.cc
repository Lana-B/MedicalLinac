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
	multFactor=8;
	// phantom size and position
	halfSize.set(60*1.1*multFactor*mm,60*1.1*multFactor*mm,340*mm);
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

	bool bCreated=false;
	G4double A, Z;

	G4NistManager* man = G4NistManager::Instance();  //creating material manager

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	//////////////////////////////////  Building silicon diode  /////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	G4Material* SiMat = man->FindOrBuildMaterial("G4_Si");

    pRmin1 = (1.0*multFactor)*mm;
    pRmax1 = (1.5*multFactor)*mm;
    pRmin2 = (1.2*multFactor)*mm;
    pRmax2 = (1.5*multFactor)*mm;
    pDz = (2.5*multFactor)*mm;
    pSPhi = 0.*deg;
    pDPhi = 360.*deg;

	G4Cons* SiDiodeLVSides = new G4Cons("subCone", pRmin1, pRmax1, pRmin2, pRmax2, pDz, pSPhi, pDPhi);

	innerRadius = (0.*multFactor);
	outerRadius = (1.5*multFactor);
	hz =(0.05*multFactor);
	startAngle = 0.*deg;
	spanningAngle = 360.*deg;

	G4Tubs* SiDiodeLVTop = new G4Tubs("siTube", innerRadius, outerRadius,hz, startAngle, spanningAngle);

	G4RotationMatrix *rotateMatrixEmpty = new G4RotationMatrix();
    G4ThreeVector zTrans(0, 0,-(pDz-hz));
    // G4ThreeVector zTrans(0, 0,0);

	G4UnionSolid* SiDiodeAdd = new G4UnionSolid("Diode", SiDiodeLVSides, SiDiodeLVTop, rotateMatrixEmpty, zTrans);
	G4LogicalVolume *SiDiodeLV = new G4LogicalVolume(SiDiodeAdd, SiMat, "SiDiodeLV", 0, 0, 0);
	G4RotationMatrix*  rotateMatrixDiode=new G4RotationMatrix();
	// rotateMatrixDiode->rotateX(90.0*deg);
	// rotateMatrixDiode->rotateY(-45.0*deg);////////////////////

	// SiDiodePV = new G4PVPlacement(rotateMatrixDiode, zTrans, "SiDiodePV", SiDiodeLV, PVWorld, false, 0);
	SiDiodePV = new G4PVPlacement(rotateMatrixDiode, centre, "SiDiodePV", SiDiodeLV, PVWorld, false, 0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////                   //////////////////////////////////////////
	////////////////////////////////////////   Build Phantom   //////////////////////////////////////////
	////////////////////////////////////////                   //////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	innerRadius2 = (1.7*multFactor)*mm;
	outerRadius2 = (4*multFactor)*mm;
	hz2 = (7.5*multFactor)*mm;
	startAngle2 = 0.*deg;
	spanningAngle2 = 360.*deg;

	G4Tubs* phantomSides
	= new G4Tubs("phantomSides", innerRadius2, outerRadius2, hz2, startAngle2, spanningAngle2);

	// G4Orb* phantom
	// = new G4Orb("phantom", 
	//             300*mm);
	G4double hz3 = (1.8*multFactor)*mm;

	G4Tubs* phantomTop
	= new G4Tubs("phantomTop",
	             0*mm, 
	             outerRadius2,
	             hz3,
	             startAngle2, 
	             spanningAngle2);

	A = 1.01*g/mole;
	G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);

	A = 12.011*g/mole;
	G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);  

	A = 16.00*g/mole;
	G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);
	G4Element* elSi = new G4Element("Silicon","Si" , Z= 14., A= 28.09*g/mole);

	G4double d= 1.18*g/cm3;
	G4int natoms, ncomponents;
	G4Material* PMMA = new G4Material("Polimetilmetacrilato",d,ncomponents=3);
	// G4Material *PMMA=G4NistManager::Instance()->FindOrBuildMaterial("G4_LUNG_ICRP"); // changable

	PMMA->AddElement(elC, natoms=5);
	PMMA->AddElement(elH, natoms=8);
	PMMA->AddElement(elO, natoms=2);

	G4RotationMatrix *rotateMatrixEmptyPhantom=new G4RotationMatrix();
    // G4ThreeVector zTransPhantom(0, 0, -(175)*mm);
    G4ThreeVector zTransPhantom(0, 0, -(hz2-hz3)*mm);

	 //---------rotation matrix for phantom --------

	// G4RotationMatrix*  rotateMatrixDiode=new G4RotationMatrix();
	// rotateMatrixDiode->rotateX(90.0*deg);
	// rotateMatrixDiode->rotateY(-45.0*deg);
	G4UnionSolid* phantomAdd =
    new G4UnionSolid("Diode", phantomSides, phantomTop, rotateMatrixEmptyPhantom, zTransPhantom);

	G4LogicalVolume *phantomLV = new G4LogicalVolume(phantomAdd, PMMA, "SiDiodePhantomLV", 0, 0, 0);
	phantomPV = new G4PVPlacement(rotateMatrixDiode, centre, "phantomPV", phantomLV, PVWorld, false, 0);




	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	//////////////////////////////////          PCB F4          /////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////


	//from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
	// http://www.phenix.bnl.gov/~suhanov/ncc/geant/rad-source/src/ExN03DetectorConstruction.cc

	//Epoxy (for FR4 )
	G4double density = 1.2*g/cm3;
	G4Material* Epoxy = new G4Material("Epoxy" , density, ncomponents=2);
	Epoxy->AddElement(elH, natoms=2);
	Epoxy->AddElement(elC, natoms=2);

	G4Material* SiO2 = 
	new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
	SiO2->AddElement(elSi, natoms=1);
	SiO2->AddElement(elO , natoms=2);

	//FR4 (Glass + Epoxy)
	density = 1.86*g/cm3;
	G4Material* FR4 = new G4Material("FR4"  , density, ncomponents=2);
	G4double fractionmass;
	FR4->AddMaterial(SiO2, fractionmass=0.528);
	FR4->AddMaterial(Epoxy, fractionmass=0.472);
	
	G4double PCB_xDim = (30*multFactor)*mm;
	G4double PCB_yDim = (30*multFactor)*mm;
	G4double PCB_zDim = (0.5*multFactor)*mm;

	G4Box *PCBBox = new G4Box("PCBBox", PCB_xDim, PCB_yDim, PCB_zDim);
	G4LogicalVolume *PCB_LV = new G4LogicalVolume(PCBBox, FR4, "PCB_LV", 0, 0, 0);
	G4ThreeVector PCB_TV (0,0,(hz2+PCB_zDim)*1.01);
	PCB_PV = new G4PVPlacement(0, PCB_TV, "PCB_PV", PCB_LV, PVWorld, false, 0);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	//////////////////////////////////      Gold electrodes     /////////////////////////////////////
	//////////////////////////////////   for charge collection  /////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	G4Material* AuMat = man->FindOrBuildMaterial("G4_Au");

	G4double electrode_innerR = (0.*multFactor)*mm;
	G4double electrode_outerR = pRmin1*0.8;
	G4double electrode_hz =(0.5*multFactor)*mm;
	G4double el_startAngle = 0.*deg;
	G4double el_spanningAngle = 360.*deg;

	G4Tubs* electrode_Tubs = new G4Tubs("electrode_LV", electrode_innerR, electrode_outerR, electrode_hz, el_startAngle, el_spanningAngle);

	G4LogicalVolume *electrode_LVTop = new G4LogicalVolume(electrode_Tubs, AuMat, "electrode_LVTop", 0, 0, 0);
	G4ThreeVector electrode_TVTop (0,0,(-(pDz+electrode_hz+0.01)));
	electrode_PVTop = new G4PVPlacement(0, electrode_TVTop, "electrode_PVTop", electrode_LVTop, PVWorld, false, 0);

	G4LogicalVolume *electrode_LVBot = new G4LogicalVolume(electrode_Tubs, AuMat, "electrode_LVBot", 0, 0, 0);
	G4ThreeVector electrode_TVBot (0,0,(-(pDz-hz-hz-electrode_hz-0.01)));
	electrode_PVBot = new G4PVPlacement(0, electrode_TVBot, "electrode_PVBot", electrode_LVBot, PVWorld, false, 0); 

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	//////////////////////////////////      Phantom plugged     /////////////////////////////////////
	//////////////////////////////////     with PMMA or glue    /////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////


    plugRmin1 = (0.0*multFactor)*mm;
    plugRmax1 = (1.0*multFactor*0.97)*mm;
    plugRmin2 = (0.0*multFactor)*mm;
    plugRmax2 = (1.2*multFactor*0.97)*mm;
    plugDz = ((pDz-hz-electrode_hz)*0.97);
    plugSPhi = 0.*deg;
    plugDPhi = 360.*deg;

	G4Cons* plugCon = new G4Cons("plugCone", plugRmin1, plugRmax1, plugRmin2, plugRmax2, plugDz, plugSPhi, plugDPhi);
	G4LogicalVolume *plugLV = new G4LogicalVolume(plugCon, PMMA, "plugLV", 0, 0, 0);
	G4ThreeVector plugTV (0,0,(electrode_hz+hz));
	plugPV = new G4PVPlacement(rotateMatrixDiode, plugTV, "plugPV", plugLV, PVWorld, false, 0);


	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	//////////////////////////////////     Region for cuts      /////////////////////////////////////
	//////////////////////////////////                          /////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	G4Region *regVol= new G4Region("SiDiodeR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*mm);
	regVol->SetProductionCuts(cuts);

	SiDiodeLV->SetRegion(regVol);
	phantomLV->SetRegion(regVol);
	plugLV->SetRegion(regVol);
	PCB_LV->SetRegion(regVol);
	electrode_LVTop->SetRegion(regVol);
	electrode_LVBot->SetRegion(regVol);
	regVol->AddRootLogicalVolume(SiDiodeLV);
	regVol->AddRootLogicalVolume(phantomLV);
	regVol->AddRootLogicalVolume(plugLV);
	regVol->AddRootLogicalVolume(PCB_LV);
	regVol->AddRootLogicalVolume(electrode_LVTop);
	regVol->AddRootLogicalVolume(electrode_LVBot);

	// Visibility
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceWireframe(true);
	SiDiodeLV->SetVisAttributes(simpleAlSVisAtt);

	G4VisAttributes* simpleAlSVisAttOrb= new G4VisAttributes(G4Colour::Gray());
	simpleAlSVisAttOrb->SetVisibility(true);
	simpleAlSVisAttOrb->SetForceWireframe(true);
	phantomLV->SetVisAttributes(simpleAlSVisAttOrb);

	G4VisAttributes* simpleAlSVisAttPlug= new G4VisAttributes(G4Colour::Gray());
	simpleAlSVisAttPlug->SetVisibility(true);
	simpleAlSVisAttPlug->SetForceSolid(true);
	plugLV->SetVisAttributes(simpleAlSVisAttPlug);

	G4VisAttributes* simpleAlSVisAtt_PCB= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt_PCB->SetVisibility(true);
	simpleAlSVisAtt_PCB->SetForceWireframe(true);
	PCB_LV->SetVisAttributes(simpleAlSVisAtt_PCB);

	G4VisAttributes* simpleAlSVisAtt_electrode= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt_electrode->SetVisibility(true);
	simpleAlSVisAtt_electrode->SetForceSolid(true);
	electrode_LVTop->SetVisAttributes(simpleAlSVisAtt_electrode);
	electrode_LVBot->SetVisAttributes(simpleAlSVisAtt_electrode);

	// Sensitive detector 
	sensDet=new CML2SDWithVoxels("Silicon Diode", saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG, G4ThreeVector(0.,0.,0.), halfSize, 100, 100, 100);
	// sensDet2=new CML2SDWithVoxels("phantom", saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG, G4ThreeVector(0.,0.,0.), halfSize, 100, 100, 100);
	G4SDManager *SDManager=G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(sensDet);
	// SDManager->AddNewDetector(sensDet2);
	
	// Read Out Geometry
	CML2ReadOutGeometry *ROG = new CML2ReadOutGeometry();
	ROG->setBuildData(PVWorld->GetFrameTranslation(), halfSize, 100, 100, 100);
	ROG->BuildROGeometry();
	sensDet->SetROgeometry(ROG);
	SiDiodeLV->SetSensitiveDetector(sensDet);
	// phantomLV->SetSensitiveDetector(sensDet2);

	bCreated=true;
	return bCreated;
}


