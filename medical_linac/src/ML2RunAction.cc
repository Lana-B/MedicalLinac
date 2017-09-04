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


#include "ML2RunAction.hh"

CML2RunAction::CML2RunAction(CML2Convergence *conv, G4int nB, G4bool bOV, int iSeedN)
{
    bRotationTranslationFilesNames=true;
    convergence=conv;
    nBeam=nB;
    bOnlyVisio=bOV;
    nLoop=0;

    char aSeed[10];
    sprintf(aSeed,"%d", iSeedN);
    seedName=(G4String)aSeed;
    auto analysisManager = G4AnalysisManager::Instance();
    // analysisManager->SetFileName("medLinacOutput");
    analysisManager->SetVerboseLevel(3);
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    analysisManager->CreateNtuple("medLinac", "medL");
    analysisManager->CreateNtupleDColumn("Vol");  // column Id = 0
    analysisManager->CreateNtupleDColumn("xLoc");  // column Id = 1
    analysisManager->CreateNtupleDColumn("yLoc");  // column Id = 2
    analysisManager->CreateNtupleDColumn("zLoc"); // column Id = 3
    analysisManager->CreateNtupleDColumn("ix"); // column Id = 4
    analysisManager->CreateNtupleDColumn("iy");    // column Id = 5
    analysisManager->CreateNtupleDColumn("iz");    // column Id = 6
    analysisManager->CreateNtupleDColumn("Dose"); // column Id = 7
    analysisManager->CreateNtupleDColumn("Dose2");    // column Id = 8
    analysisManager->CreateNtupleIColumn("nEvents");    // column Id = 9
    analysisManager->CreateNtupleDColumn("xPos"); // column Id = 10
    analysisManager->CreateNtupleDColumn("yPos");    // column Id = 11
    analysisManager->CreateNtupleIColumn("zPos");    // column Id = 12
    analysisManager->FinishNtuple();
}

CML2RunAction::~CML2RunAction(void)
{
    // std::cout<<"!!!\n\n !! \n\n destructor run action \n\n !! \n\n !!"<<std::endl;
    delete G4AnalysisManager::Instance();  

}
void CML2RunAction::BeginOfRunAction(const G4Run *)
{
    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();

    // std::cout<<"!!!!\n\n!! begin: "<<analysisManager<<"!!!!\n\n!!!!"<<std::endl;
    // Open an output file 
    // The default file name is set in B5RunAction::B5RunAction(),
    // it can be overwritten in a macro
    G4String fileName = "/tmp/lb8075/medLinacOutputSinglewPhantom_"+seedName;

    analysisManager->OpenFile(fileName);

    G4String fullName;
    if (bRotationTranslationFilesNames)
    {fullName=CML2AcceleratorConstruction::GetInstance()->getCurrentRotationString()+
                CML2PhantomConstruction::GetInstance()->getCurrentTranslationString();}
    else
    {fullName="";}
    CML2PhantomConstruction::GetInstance()->setNewName(fullName);

    CML2AcceleratorConstruction::GetInstance()->writeInfo();
    CML2PhantomConstruction::GetInstance()->writeInfo();

    std::cout<<"*********************************************"<<'\n';
    if (convergence->getNMaxLoops()<0 || bOnlyVisio)
    {
        std::cout << "loop n. "<<++nLoop <<'\n';
    }
    else
    {
        std::cout << "loop n. "<<++nLoop<<"/" <<convergence->getNMaxLoops() <<'\n';
    }
    if (!bOnlyVisio)
    {std::cout << "Launched "<< nBeam <<" random primary particles" << '\n';}
    std::cout<<"*********************************************"<<'\n';
    MyTime.Start();
}
void CML2RunAction::EndOfRunAction(const G4Run *)
{
    // save histograms & ntuple
    //
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    
    CML2WorldConstruction::GetInstance()->savePhantomData();
    CML2WorldConstruction::GetInstance()->savePhaseSpaceData();
    convergence->saveResults();

    MyTime.Stop();
    loopElapsedTime=MyTime.GetUserElapsed();
    std::cout << "loop elapsed time [s] : "<< loopElapsedTime << '\n';
    std::cout <<'\n';
}
