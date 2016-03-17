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
/// \file electromagnetic/TestEm6/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 83428 2014-08-21 15:46:01Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,RunAction* RuAct,EventAction* EveAct)
:G4UserSteppingAction(),fDetector(det),fRunAction(RuAct),eventAction(EveAct)
{ 
 fMuonMass = G4MuonPlus::MuonPlus()->GetPDGMass();
 G4cout << " <<<-<-->- " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 G4StepPoint* PrePoint = aStep->GetPreStepPoint();  
 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 if (process == 0) return;  
 G4String processName = process->GetProcessName();
 fRunAction->CountProcesses(processName); //count processes

 // G4cout << " Process Name " << processName << G4endl; 
// if (processName != "GammaToMuPair") return;
 
 G4double      EGamma  = PrePoint->GetTotalEnergy();
 G4ThreeVector PGamma  = PrePoint->GetMomentum();
    
 G4double      Eplus(0), Eminus(0);
 G4ThreeVector Pplus   , Pminus;
 const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
 for (size_t lp=0; lp<(*secondary).size(); lp++) {
   if ((*secondary)[lp]->GetDefinition()==G4MuonPlus::MuonPlusDefinition()) {
     Eplus  = (*secondary)[lp]->GetTotalEnergy();
     Pplus  = (*secondary)[lp]->GetMomentum();
   } else {
     Eminus = (*secondary)[lp]->GetTotalEnergy();
     Pminus = (*secondary)[lp]->GetMomentum();                 
   }
 }
               
 G4double xPlus = Eplus/EGamma, xMinus = Eminus/EGamma;
 G4double thetaPlus = PGamma.angle(Pplus), thetaMinus = PGamma.angle(Pminus);
 G4double GammaPlus = Eplus/fMuonMass;
 G4double GammaMinus= Eminus/fMuonMass;
 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 if(0.0 == thetaPlus || 0.0 == thetaMinus) {
   G4cout << "SteppingAction: "
          << "thetaPlus= " << thetaPlus << " thetaMinus= " << thetaMinus
          << " gamPlus= " << GammaPlus << " gamMinus= " <<  GammaMinus
          << "  " << thetaPlus *GammaPlus - thetaMinus*GammaMinus << G4endl;
   return;
 }
 analysisManager->FillH1(1,1./(1.+std::pow(thetaPlus*GammaPlus,2)));
 analysisManager->FillH1(2,std::log10(thetaPlus*GammaPlus));

 analysisManager->FillH1(3,std::log10(thetaMinus*GammaMinus));
 analysisManager->FillH1(4,std::log10(std::fabs(thetaPlus *GammaPlus
                                              -thetaMinus*GammaMinus)));
 
 analysisManager->FillH1(5,xPlus);
 analysisManager->FillH1(6,xMinus);
/*
 for (size_t lp=0; lp<(*secondary).size(); lp++) 
 {
    G4cout << " ******* ***** **** *** **** " << G4endl;
    G4cout << (*secondary)[lp]->GetDefinition()->GetParticleName() << " " << (*secondary)[lp]->GetTotalEnergy() << G4endl;
    G4cout << (*secondary)[lp]->GetVolume()->GetName() << G4endl;
    G4cout << (*secondary)[lp]->GetPosition() << G4endl;
    G4cout << " ******* ***** **** *** **** " << G4endl;
 }
 */

  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  G4ThreeVector pos               =      aTrack->GetPosition();
  G4ThreeVector direction = endPoint->GetMomentumDirection();
  G4ThreeVector momentum = endPoint->GetMomentum()/CLHEP::MeV;
  G4ThreeVector momDir    = aTrack->GetMomentumDirection();
  G4String particleName = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4double TotalEnergy = endPoint->GetTotalEnergy();
  G4int pdg = aTrack->GetDefinition()->GetPDGEncoding();

  std::ofstream &output = fRunAction->GetOutput();
  if(PrePoint->GetTouchableHandle()->GetVolume()->GetName() == "Target" &&
		  endPoint->GetTouchableHandle()->GetVolume()->GetName() != "Target" && pdg!=22)
  {
//output << pdg << "  " << momDir.x()/CLHEP::MeV << " " << momDir.y()/CLHEP::MeV<< " " <<  momDir.z()/CLHEP::MeV
  output << pdg << "  " << momentum.x()/CLHEP::MeV << " " << momentum.y()/CLHEP::MeV<< " " <<  momentum.z()/CLHEP::MeV
	  << " " <<  TotalEnergy/CLHEP::MeV << " "
	  << pos.y()/CLHEP::cm << " " << pos.z()/CLHEP::cm << " " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


