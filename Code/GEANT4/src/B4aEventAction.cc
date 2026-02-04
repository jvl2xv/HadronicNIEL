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
// $Id: B4aEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class 

#include "B4aEventAction.hh"
#include "B4RunAction.hh"


#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction(B4RunAction* runActionIn, int numMaterial, int numDimension)
 : G4UserEventAction(),
   runAction(runActionIn),
   fnumMaterial(numMaterial),
   fnumDimension(numDimension)

{ 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{ 
//eventTotEdepMask = 0.; 
 eventTotEdepSample_keV = 0.0;
 eventNIELEdepSample_keV = 0.0;

 eventNIELEdepSample_keV_Coulomb = 0.0;
 eventNIELEdepSample_keV_HadronEl = 0.0;
 eventNIELEdepSample_keV_HadronInel = 0.0;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
	//get event #
    G4int eID = 0;
    G4int runNum = 0;
	  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
    const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
    if(evt) eID = evt->GetEventID();
    if(run) runNum = run->GetRunID();


    // don't save the zeros because it takes forever to read them all in python, just pad later if you need to

  if (eventNIELEdepSample_keV > 0.0) {
  runAction->NIEL_Total_keV.push_back(eventNIELEdepSample_keV); }
  if (eventNIELEdepSample_keV_Coulomb > 0.0) {
  runAction->NIEL_Coulomb_keV.push_back(eventNIELEdepSample_keV_Coulomb); }
  if (eventNIELEdepSample_keV_HadronInel > 0.0) {
  runAction->NIEL_HadronInel_keV.push_back(eventNIELEdepSample_keV_HadronInel); }
  if (eventNIELEdepSample_keV_HadronEl > 0.0) { 
  runAction->NIEL_HadronEl_keV.push_back(eventNIELEdepSample_keV_HadronEl); }




  /*

  int countsIN[] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 5000000, 5000000, 5000000, 10000000, 10000000, 10000000, 10000000, 10000000, 10000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000};
  int num_in = countsIN[fnumDimension]; 

  G4String fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV.txt";
  std::ostringstream commandOS;
  commandOS << fileName;
  std::ofstream ofile;
  ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
  // in cm and MeV 
  ofile << eventTotEdepSample_keV << " " << eventNIELEdepSample_keV << "\n";
  // ofile << eID << " " << eventNIELEdepSample_keV << "\n";
  ofile.close(); 

  commandOS.str("");
  fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_Coulomb.txt";
  commandOS << fileName;
  ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
  // in cm and MeV 
  ofile << eventTotEdepSample_keV << " " << eventNIELEdepSample_keV_Coulomb << "\n";
  // ofile << eID << " " << eventNIELEdepSample_keV << "\n";
  ofile.close(); 

  commandOS.str("");
  fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_HadronEl.txt";
  commandOS << fileName;
  ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
  // in cm and MeV 
  ofile << eventTotEdepSample_keV << " " << eventNIELEdepSample_keV_HadronEl << "\n";
  // ofile << eID << " " << eventNIELEdepSample_keV << "\n";
  ofile.close(); 

  commandOS.str("");
  fileName = std::to_string(fnumMaterial) + "_" + std::to_string(num_in) + "_" + std::to_string(fnumDimension)+ "_TotEdepkeV_NIELEdepkeV_HadronInel.txt";
  commandOS << fileName;
  ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
  // in cm and MeV 
  ofile << eventTotEdepSample_keV << " " << eventNIELEdepSample_keV_HadronInel << "\n";
  // ofile << eID << " " << eventNIELEdepSample_keV << "\n";
  ofile.close(); 



*/


   
  }
  
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
