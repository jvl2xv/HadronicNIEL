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
/// \file electromagnetic/TestEm7/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "PhysListEmStandardNR.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "PhysListEmStandardSSM.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"

#include "G4DecayPhysics.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"
#include "G4HadronPhysicsFTFQGSP_BERT.hh"

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4UnitsTable.hh"

#include "G4ProcessManager.hh"
#include "G4Decay.hh"

#include "StepMax.hh"

#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4StepLimiterPhysics.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList(),
  fStepMaxProcess(nullptr)
{
  fHelIsRegisted  = false;
  fBicIsRegisted  = false;
  fBiciIsRegisted = false;

  // protected member of the base class
  verboseLevel = 1;

  fMessenger = new PhysicsListMessenger(this);

  // EM physics
  fEmName = G4String("emstandard_opt0");  
  fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);

  // G4VUserPhysicsList::SetProductionCut(1e-9 * mm, G4Proton::Proton());

  // allow setting the max step limit defined for the logical solid in detector construction
  RegisterPhysics(new G4StepLimiterPhysics());


  // Deacy physics and all particles
  fDecPhysicsList = new G4DecayPhysics(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fEmPhysicsList;
  delete fDecPhysicsList;
  for(size_t i=0; i<fHadronPhys.size(); i++) {delete fHadronPhys[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  // transportation
  //
  AddTransportation();
  
  // electromagnetic physics list
  //
  fEmPhysicsList->ConstructProcess();

  // decay physics list
  //
  fDecPhysicsList->ConstructProcess();
  
  // hadronic physics lists
  for(size_t i=0; i<fHadronPhys.size(); i++) {
    fHadronPhys[i]->ConstructProcess();
  }
  
  // step limitation (as a full process)
  //  
  AddStepMax();  

  // Bias hadron elastic and inelastic processes for protons so I don't have to simulate as many particles

  G4HadronicProcess* muInelasticProcess = nullptr;
  G4ParticleDefinition* particle_proton = G4Proton::Definition();
  G4ProcessManager* pmanager = particle_proton->GetProcessManager();

  G4ProcessVector* pvec = pmanager->GetProcessList();
  for (std::size_t i = 0; i < pvec->size(); i++) {
    if ((*pvec)[i]->GetProcessName() == "hadElastic") {
      muInelasticProcess = static_cast<G4HadronicProcess*>((*pvec)[i]);
      muInelasticProcess->BiasCrossSectionByFactor(1000);
    }
     if ((*pvec)[i]->GetProcessName() == "protonInelastic") {
      muInelasticProcess = static_cast<G4HadronicProcess*>((*pvec)[i]);
      muInelasticProcess->BiasCrossSectionByFactor(1000);
    }

  }




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "local") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandard(name);

  } else if (name == "emstandard_opt0") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);

  } else if (name == "emstandard_opt1") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1(verboseLevel);

  } else if (name == "emstandard_opt2") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2(verboseLevel);
    
  } else if (name == "emstandard_opt3") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3(verboseLevel);
    
  } else if (name == "emstandard_opt4") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);

  } else if (name == "ionGasModels") {

    AddPhysicsList("emstandard_opt0");
    fEmName = name;
    AddIonGasModels();

  } else if (name == "standardNR") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandardNR(name);

  } else if (name == "emlivermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);

  } else if (name == "empenelope") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);

  } else if (name == "emlowenergy") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics(verboseLevel);

  } else if (name == "emstandardSSM") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandardSSM("standardSSM");

  } else if (name == "emstandardWVI") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsWVI(verboseLevel);

  } else if (name == "emstandardGS") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsGS(verboseLevel);

  } else if (name == "elastic" && !fHelIsRegisted) {
    // G4HadronElasticProcess with added dataset G4BGGNucleonElasticXS (Barashenkov-Glauber-Gribov parameterization) and registered G4ChipsElasticModel mult by XSFactorNucleonElastic 
    // Barashenkov parameterisation is used below 91 GeV and Glauber-Gribov
    // note Chips/Diffuse this determines the angle elastic scatter - Chips does better for protons because parameterized from experimental data https://repository.cern/records/vqrgw-wmb90/preview/CERN_report_AYeltokov.pdf
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
    fHelIsRegisted = true;
  

  } else if (name == "DElastic" && !fHelIsRegisted) {
    // G4HadronElasticProcess with added dataset G4BGGNucleonElasticXS (Barashenkov-Glauber-Gribov parameterization) and registered G4DiffuseElastic mult by XSFactorNucleonElastic 
    fHadronPhys.push_back( new G4HadronDElasticPhysics(verboseLevel));
    fHelIsRegisted = true;

  } else if (name == "HElastic" && !fHelIsRegisted) {
    // G4HadronElasticProcess with added dataset G4BGGNucleonElasticXS (Barashenkov-Glauber-Gribov parameterization) and registered G4DiffuseElastic mult by XSFactorNucleonElastic but use Chips for hydrogen (doesn't matter)
    fHadronPhys.push_back( new G4HadronHElasticPhysics(verboseLevel));
    fHelIsRegisted = true;
  

   // https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/reference_PL/FTFP_BERT.html
  } else if (name == "binary" && !fBicIsRegisted) {
    fHadronPhys.push_back(new G4HadronInelasticQBBC(verboseLevel));
    //fHadronPhys.push_back(new G4HadronPhysicsFTFQGSP_BERT(verboseLevel)); // this didn't help
    fBicIsRegisted = true;
    
  } else if (name == "binary_ion" && !fBiciIsRegisted) {
    fHadronPhys.push_back(new G4IonPhysics(verboseLevel));
    fBiciIsRegisted = true;
    

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  fStepMaxProcess = new StepMax();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fStepMaxProcess->IsApplicable(*particle) && pmanager)
      {
        pmanager ->AddDiscreteProcess(fStepMaxProcess);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddIonGasModels()
{
  G4EmConfigurator* em_config = 
    G4LossTableManager::Instance()->EmConfigurator();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    G4String partname = particle->GetParticleName();
    if(partname == "alpha" || partname == "He3" || partname == "GenericIon") {
      G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
      G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
      G4double eth = 2.*MeV*particle->GetPDGMass()/proton_mass_c2;
      em_config->SetExtraEmModel(partname,"ionIoni",mod1,"",0.0,eth,
                                 new G4IonFluctuations());
      em_config->SetExtraEmModel(partname,"ionIoni",mod2,"",eth,100*TeV,
                                 new G4UniversalFluctuation());

    }
  }
}


void PhysicsList::SetCuts()
{ 
 // fixe lower limit for cut
 G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

 // call base class method to set cuts which default value can be
 // modified via /run/setCut/* commands
 G4VUserPhysicsList::SetCuts();

 DumpCutValuesTable();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

