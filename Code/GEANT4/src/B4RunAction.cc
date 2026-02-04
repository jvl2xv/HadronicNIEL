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
// $Id: B4RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class 

#include "B4RunAction.hh"
#include <fstream>


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B4PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction(int numMaterial, int numDimension)
 : G4UserRunAction(),
    fnumMaterial(numMaterial),
    fnumDimension(numDimension)
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* run)
{ 

  NIEL_Total_keV = {};
  NIEL_Coulomb_keV = {};
  NIEL_HadronInel_keV = {};
  NIEL_HadronEl_keV = {};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* run)
{ 

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

    // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B4PrimaryGeneratorAction* generatorAction
   = static_cast<const B4PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }




     int countsIN[] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 10000000, 50000000, 50000000, 50000000, 50000000, 50000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000};
     int num_in = countsIN[fnumDimension]; 

      

    if (NIEL_Total_keV.size() > 0) {
    G4String fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_Total_"+std::to_string(G4Threading::G4GetThreadId())+".txt";
    std::ostringstream commandOS;
    commandOS << fileName;
    std::ofstream ofile;
    ofile.open(G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
    for (int i = 0; i < NIEL_Total_keV.size(); ++i) {
      ofile << NIEL_Total_keV[i] << "\n"; }
    ofile.close();  }

    if (NIEL_Coulomb_keV.size() > 0) {
    G4String fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_Coulomb_"+std::to_string(G4Threading::G4GetThreadId())+".txt";
    std::ostringstream commandOS2;
    commandOS2 << fileName;
    std::ofstream ofile2;
    ofile2.open(G4String(commandOS2.str()), ios::out | ios::app);     // ascii file   
    for (int i = 0; i < NIEL_Coulomb_keV.size(); ++i) {
      ofile2 << NIEL_Coulomb_keV[i] << "\n"; }
    ofile2.close(); } 

    if (NIEL_HadronEl_keV.size() > 0) {
    G4String fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_HadronEl_"+std::to_string(G4Threading::G4GetThreadId())+".txt";
    std::ostringstream commandOS3;
    commandOS3 << fileName;
    std::ofstream ofile3;
    ofile3.open(G4String(commandOS3.str()), ios::out | ios::app);     // ascii file   
    for (int i = 0; i < NIEL_HadronEl_keV.size(); ++i) {
      ofile3 << NIEL_HadronEl_keV[i] << "\n"; }
    ofile3.close();  } 

    if (NIEL_HadronInel_keV.size() > 0) {
    G4String fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_HadronInel_"+std::to_string(G4Threading::G4GetThreadId())+".txt";
    std::ostringstream commandOS4;
    commandOS4 << fileName;
    std::ofstream ofile4;
    ofile4.open(G4String(commandOS4.str()), ios::out | ios::app);     // ascii file   
    for (int i = 0; i < NIEL_HadronInel_keV.size(); ++i) {
      ofile4 << NIEL_HadronInel_keV[i] << "\n";  }
    ofile4.close(); } 


}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
