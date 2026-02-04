#include "TrackingAction.hh"
#include "B4aEventAction.hh"
#include "B4aSteppingAction.hh"
#include "B4RunAction.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
//#include "AIDA/AIDA.h"



#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "B4RunAction.hh"
#include "G4Decay.hh"
#include "G4RadioActiveDecay.hh"


#include "G4ElectronIonPair.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>
#include "G4LindhardPartition.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"


using namespace std;



TrackingAction::TrackingAction(B4aEventAction* EvAct, B4RunAction* run)
:evAction(EvAct), Run(run)
{ }


TrackingAction::~TrackingAction()
{ }

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{



	const G4ParticleDefinition * particleDef = aTrack->GetParticleDefinition();
	G4String particleName = particleDef->GetParticleName();
    
   // if it isn't a proton kill it, we record the energy of the secondaries properly anyway
  if (particleName == "e-") {   
      const_cast<G4Track*>(aTrack)->SetTrackStatus(fStopAndKill);
    }

 

          
        }




void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

}



