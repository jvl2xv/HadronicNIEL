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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class 

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "B4RunAction.hh"
#include "G4Decay.hh"
#include "G4RadioActiveDecay.hh"
#include "G4Run.hh"
#include "G4EmCalculator.hh"


#include "G4ElectronIonPair.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"

#include "G4NIELCalculator.hh"
#include "G4ICRU49NuclearStoppingModel.hh"
#include "G4LindhardPartition.hh"

#include <iostream>
#include <random>
#include <ctime>


using namespace std;

//G4Decay* theDecayProcess = new G4Decay();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction, B4RunAction* runAction, int numMaterial, int numDimension)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    fRunAction(runAction),
    fnumMaterial(numMaterial),
    fnumDimension(numDimension)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* theStep)
{


   


  //get event #
  G4int eID = 0;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
  G4int rID = run->GetRunID();
  if(evt) eID = evt->GetEventID();
  G4Track *theTrack = theStep->GetTrack(); 



  G4StepPoint* preStepPoint = theStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = theStep->GetPostStepPoint();


if ((theTrack->GetVolume()->GetName() != "World")  & (preStepPoint->GetProcessDefinedStep() != 0)) {

  if ((preStepPoint->GetPhysicalVolume()->GetName() == "Sample")  || (postStepPoint->GetPhysicalVolume()->GetName() == "Sample")) {

  
  // G4double energyDepositTotal_MeV = theStep->GetTotalEnergyDeposit()/keV;
  // NOTE all this does is sample this energy loss table modelICRU49 totally bad, can't use
  // G4double energyDepositNIEL_MeV = theStep->GetNonIonizingEnergyDeposit()/keV;

  G4String process_end_step = postStepPoint->GetProcessDefinedStep()->GetProcessName();
  G4String parent_particle_name = theTrack->GetParticleDefinition()->GetParticleName();


  fEventAction->AddTotEdepSample_keV(theStep->GetTotalEnergyDeposit()/keV);
  // fEventAction->AddNIELEdepSample_keV(theStep->GetNonIonizingEnergyDeposit()/keV);

 
  const std::vector<const G4Track*>* secondary  = theStep->GetSecondaryInCurrentStep(); 

  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    G4ParticleDefinition* particle = (*secondary)[lp]->GetDefinition(); 
    G4String secondary_name   = particle->GetParticleName();    
    G4double secondary_energy = (*secondary)[lp]->GetKineticEnergy();
    G4int secondary_Z = particle->GetAtomicNumber();
    G4int secondary_A = particle->GetAtomicMass();

    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* G4_Si =   nistManager->FindOrBuildMaterial("G4_Si");


    if (secondary_Z > 1) {


      G4String ionGenProcessName =  postStepPoint->GetProcessDefinedStep()->GetProcessName();

      if (parent_particle_name == "proton") {

      // cout << secondary_name << " " << (*secondary)[lp]->GetKineticEnergy()/keV << endl;

      G4LindhardRobinsonPartition* LindhardPartition = new G4LindhardRobinsonPartition();
      // note this just uses the dominant material atom, should define yourself taking into account the whole material (but we are only doing pure elements so it's cool)
      G4double partition = LindhardPartition->PartitionNIEL(secondary_Z, secondary_A, G4_Si, secondary_energy); 

     

     /*

       int countsIN[] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 5000000, 5000000, 5000000, 10000000, 100000000, 100000000, 100000000, 100000000, 100000000, 1000000000, 1000000000, 1000000000, 1000000000, 1000000000, 1000000000, 1000000000};
       int num_in = countsIN[fnumDimension]; 
      
      G4double PKA_E_eV = (*secondary)[lp]->GetKineticEnergy()/eV; 
      G4String fileName =  std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension) + "_PKA_E_eV.txt";
      std::ostringstream commandOS;
      commandOS << fileName;
      std::ofstream ofile;
      ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
      // in cm and MeV 
      ofile << PKA_E_eV << "\n";
      ofile.close(); 


      */
    

       // implement an recoil energy threshold manually because it causes errors to do it in the physics list it seems, note seems to work better if set to zero
       if (partition*((*secondary)[lp]->GetKineticEnergy())/eV > 10.0) {


      if (ionGenProcessName == "hadElastic") {
        fEventAction->AddNIELEdepSample_keV_HadronEl(partition*((*secondary)[lp]->GetKineticEnergy()/keV));
            

            // 2. Define the distribution for the desired range [0, 100] - because mult prob by 100 so only add every 1/100
            std::uniform_int_distribution<int> dist(0, 1000);
            std::random_device rd;
            std::mt19937 engine(rd());
            // 3. Generate the random number
            int random_num = dist(engine);
            if (random_num == 0) {
              fEventAction->AddNIELEdepSample_keV(partition*((*secondary)[lp]->GetKineticEnergy()/keV)); }
            }

      else if (ionGenProcessName == "protonInelastic") {
        fEventAction->AddNIELEdepSample_keV_HadronInel(partition*((*secondary)[lp]->GetKineticEnergy()/keV));

            // 2. Define the distribution for the desired range [0, 100] - because mult prob by 100 so only add every 1/100
            std::uniform_int_distribution<int> dist(0, 1000);
            std::random_device rd;
            std::mt19937 engine(rd());
            // 3. Generate the random number
            int random_num = dist(engine);
            if (random_num == 0) {
              fEventAction->AddNIELEdepSample_keV(partition*((*secondary)[lp]->GetKineticEnergy()/keV)); }

      }
       // used to be this: (ionGenProcessName == "ScreenedElastic") but had to add extra processes for the relativistic and I am not sure what they would be called
      else  {
        fEventAction->AddNIELEdepSample_keV_Coulomb(partition*((*secondary)[lp]->GetKineticEnergy()/keV));
        fEventAction->AddNIELEdepSample_keV(partition*((*secondary)[lp]->GetKineticEnergy()/keV)); 

    }


      
    }
    } 

  }



  //G4EmCalculator emCalc;
  //G4ParticleDefinition* particle = G4Proton::Definition();



  //cout << emCalc.ComputeNuclearDEDX(5.0*MeV, particle, G4_Si)/(MeV/cm)  << " " << emCalc.ComputeTotalDEDX(5.0*MeV, particle, G4_Si)/(MeV/cm) << endl;


     } 

   }
     

     }


  // only save the first million, that's plenty of stats
   if (eID < 100000) {
  
  // record energy when they exit the sample
  if ((theTrack->GetVolume()->GetName() != "World")  & (postStepPoint->GetProcessDefinedStep() != 0)  & (preStepPoint->GetProcessDefinedStep() != 0)) {
  if ((preStepPoint->GetPhysicalVolume()->GetName() == "Sample")  & (postStepPoint->GetPhysicalVolume()->GetName() != "Sample")) {




  G4double exit_KE_MeV = postStepPoint->GetKineticEnergy()/MeV;

  G4String fileName =  std::to_string(fnumMaterial)+ "_" +std::to_string(fnumDimension)+ "_ExitEnergyMeV.txt";
  std::ostringstream commandOS;
  commandOS << fileName;
  std::ofstream ofile;
  ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
  // in cm and MeV 
  ofile << exit_KE_MeV << "\n";
  // ofile << eID << " " << eventNIELEdepSample_keV << "\n";
  ofile.close(); 




     }    } }



}























//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......