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

// REQUIRED GEANT HEADER FILES

#include "B4DetectorConstruction.hh"
#include "B4aActionInitialization.hh"
#include "TrackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4ScoringManager.hh"
#include "G4PhysListFactory.hh"
#include "G4Box.hh"
//#include "FTFP_BIC.hh"




#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4ScreenedNuclearRecoil.hh"

#include "G4DecayPhysics.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"

#include "G4StepLimiterPhysics.hh"

#include "QGSP_BERT.hh"

#include "PhysicsList.hh"





// REQUIRED C++ HEADER FILES FOR IO AND MATH



#include <math.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <stdlib.h>
//#include <time.h>
#include <cmath>
#include <vector>

#include <ctime>

#include <iostream>
#include <random>

using namespace std;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " myMesh [-m macro ] [-c numMaterial] [-d numDimension] [-p EMphysics] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments

  //if ( argc > 7 ) {
  //  PrintUsage();
  //  return 1;
  //}
  
  G4String macro;
  G4String session;
  G4String numMaterial;
  G4String numDimension;
  G4String EMphysics; 

#ifdef G4MULTITHREADED
  G4int nThreads = 5;
#endif


 for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-c" ) numMaterial = argv[i+1];
    else if ( G4String(argv[i]) == "-d" ) numDimension = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) EMphysics = argv[i+1];
 #ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]); }
 #endif
    else {
      PrintUsage();
      return 1;
    }
   }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager (multithread)
  
  #ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    if ( nThreads > 0 ) { 
      runManager->SetNumberOfThreads(nThreads);
    }  
  #else
    G4RunManager * runManager = new G4RunManager;
  #endif



  // Initialize the required clases 
  
  // UI create a mesh if you want to use this
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
  // Register the physics list
  G4int verbose = 4;
  G4PhysListFactory factory;

  // initialize building a string (for c++ UI stuff)
  G4String command = "";
  std::ostringstream commandOS;
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("Shielding");

  //G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("FTFP_BIC_EMY");
  //physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics());


  // PHYSICS 1
  // G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("Shielding");

  // PHYSICS 2
  // G4PhysicsListHelper* helper = G4PhysicsListHelper::GetPhysicsListHelper(); 


  




  // PHYSICS 3
  // G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QGSP_INCLXX_HP");
  //physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics() );
  




  //physicsList->ReplacePhysics(new PhysListEmStandardISS);

  //G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QGSP_BIC_EMY");
  //physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics() );


  //physicsList->RegisterPhysics( new PhysListEmStandardNR );
  //physicsList->RegisterPhysics( new PhysListEmStandardISS );

  
  //G4VModularPhysicsList*  physicsList  = factory.GetReferencePhysList("QBBC");
  
  //physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
  //physicsList->RegisterPhysics(new PhysListEmStandardSS());
  //physicsList->RegisterPhysics(new PhysListEmStandardNR());


  //PhysicsList* physicsList = new PhysicsList();
  //physicsList->SetCuts();
  //physicsList->AddPhysicsList("standardNR");

    //QGSP_BERT *physicsList = new QGSP_BERT;
    //LBE *physicsList = new LBE;
    //physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics );
    //physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
    //physicsList->RegisterPhysics(new PhysListEmStandardSS());
    //physicsList->RegisterPhysics(new PhysListEmStandardNR());

    //G4VPhysicsConstructor*   physicsList = new PhysListEmStandardNR();



  PhysicsList* physicsList = new PhysicsList();
  // Coulomb
  physicsList->AddPhysicsList(EMphysics);
  cout << "applying physics: " << EMphysics << endl;
  // Hadronic elastic - 
  // G4HadronElasticProcess with added dataset G4BGGNucleonElasticXS (Barashenkov-Glauber-Gribov parameterization) and registered G4ChipsElasticModel mult by XSFactorNucleonElastic 
  // Barashenkov parameterisation is used below 91 GeV and Glauber-Gribov
   // note Chips/Diffuse this determines the angle elastic scatter - Chips does better for protons because parameterized from experimental data
  // G4ChipsElasticModel from 0 to 100 TeV. This model uses the Kossov parameterised cross sections.
  physicsList->AddPhysicsList("elastic");
  cout << "applying physics: hadron elastic" << endl;
  // hadron inelastic - adds for protons FTFP-BERT-BIC - Bertini intranuclear cascade and precompound model to de-excite the remnant nucleus after the high energy interaction, which then calls the Fermi breakup, neutron and light ion evaporation and photon evaporation models as needed. 
  // Inelastic nucleus-nucleus scattering for all incident A is handled by the Binary Light Ion Cascade (BIC) 
  // https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/reference_PL/FTFP_BERT.html 
  physicsList->AddPhysicsList("binary");
  cout << "applying physics: hadron binary" << endl;
  // binary cascade and FTFP to handle the pre-compound, 
  physicsList->AddPhysicsList("binary_ion");
  cout << "applying physics: hadron binary ion" << endl;
  runManager->SetUserInitialization(physicsList);
  //physicsList->DumpList();



  // Register the pointer to the visualization manager
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();



  // Process macro or start UI session
  // To run use this instead of a macro, use the -m specifier and say that the macro name is julieRun
  if ( macro.size() ) {
    if (macro == "julieRun") {

      // initialize a random seed
      srand(time(NULL));

          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_int_distribution<> dis(0, 1000000000);


          double energies_MeV[] = {0.0005, 0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 200.0, 300.0, 500.0, 700.0, 1000.0};
          double thicknesses_mm[] = {1e-5, 1e-5, 1e-5, 1e-5,1e-5, 1e-5, 1e-5, 1e-5, 5e-5, 5e-5, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 5e-3, 5e-3, 5e-3 , 5e-3 , 5e-3, 5e-3}; 
          int counts_in[] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 10000000, 50000000, 50000000, 50000000, 50000000, 50000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000, 100000000};


          G4int rand1 = dis(gen);
          G4int rand2 = dis(gen);

          cout << "My random numbers were " << rand1 << " " << rand2 << endl;
          commandOS << "/random/setSeeds " << rand1 << " " << rand2;
          UImanager->ApplyCommand(G4String(commandOS.str()));
          commandOS.str("");

          B4DetectorConstruction* detConstruction1 = new B4DetectorConstruction(std::stoi(numMaterial), std::stoi(numDimension));
          // detConstruction1->SetParameters(this_energy_num);
          detConstruction1->Construct();
          runManager->SetUserInitialization(detConstruction1);
          //runManager->GeometryHasBeenModified();

          // seems like only run 1 onward work so run a teeny one first
          int count_in[] = {2, counts_in[std::stoi(numDimension)]};
          int num_count_in = sizeof(count_in)/sizeof(count_in[0]);           
          for (int this_count_num = 0; this_count_num < num_count_in; ++this_count_num){


          B4aActionInitialization* actionInitialization1 = new B4aActionInitialization(detConstruction1, std::stoi(numMaterial), std::stoi(numDimension));
          runManager->SetUserInitialization(actionInitialization1);
          UImanager->ApplyCommand("/run/reinitializeGeometry");
          // Force all physics tables recalculated again - maybe this will help the seg fault from the geometry rebuild
          UImanager->ApplyCommand("/run/physicsModified");

          // initialization commands for geant (just applying what would be in the macro file)
          UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1 ");
          UImanager->ApplyCommand("/vis/disable");
          UImanager->ApplyCommand("/run/verbose 0");
          UImanager->ApplyCommand("/event/verbose 0");
          UImanager->ApplyCommand("/process/verbose 0");
          UImanager->ApplyCommand("/vis/verbose 0");
          UImanager->ApplyCommand("/tracking/verbose 0");


        command = "/run/setCut 0.1 nm";
        UImanager->ApplyCommand(command);
        cout << command << endl;

        command = "/cuts/setLowEdge 5.0 eV";
        UImanager->ApplyCommand(command);
        cout << command << endl;   



          // proton
          command = "/gps/particle proton";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          command = "/gps/pos/centre 0.0 0.0 0.25 cm";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          command = "/gps/direction 0 0 -1";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          // energy range of their paper = 0.001 MeV to 1000 MeV
          command =  "/gps/energy "+std::to_string(energies_MeV[std::stoi(numDimension)])+" MeV";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          command = "/run/particle/dumpCutValues";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          UImanager->ApplyCommand("/run/initialize");

          // run
          commandOS << "/run/beamOn " +std::to_string(count_in[this_count_num]);
          cout << G4String(commandOS.str()) << endl;
          UImanager->ApplyCommand(G4String(commandOS.str()));
          commandOS.str("");


       }
    
  }
    else {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro); }
  }
  else  {  

      G4String command = "";
      std::ostringstream commandOS;


        // if I want to see the geometry, run without the command line argument

        B4DetectorConstruction* detConstruction = new B4DetectorConstruction(0,0);
        detConstruction->DefineVolumes();
        // register the detector
        runManager->SetUserInitialization(detConstruction);
        // set the action initialization
        B4aActionInitialization* actionInitialization = new B4aActionInitialization(detConstruction, 0, 0);
        // register the action initialization
        runManager->SetUserInitialization(actionInitialization);

        G4double posX = 0.0;
          G4double posY = 0.0;
          G4double posZ = 0.0;

          command = "/gps/particle alpha";
          UImanager->ApplyCommand(command);

          command = "/gps/pos/type Surface";
          UImanager->ApplyCommand(command);

          command = "/gps/pos/shape Cylinder";
          UImanager->ApplyCommand(command);

          commandOS << "/gps/pos/centre " << posX << " " << posY << " " << posZ << " " << " cm";
          UImanager->ApplyCommand(G4String(commandOS.str()));
          commandOS.str("");

          command = "/gps/pos/radius 2.0 cm";
          UImanager->ApplyCommand(command);

          command = "/gps/pos/halfz 0.5 mm";
          UImanager->ApplyCommand(command);

          command = "/gps/ang/type iso";
          UImanager->ApplyCommand(command);

          command = "/gps/energy 5.486 MeV";
          UImanager->ApplyCommand(command);

          commandOS << "/run/beamOn 10000000";
          UImanager->ApplyCommand(G4String(commandOS.str()));
          commandOS.str("");

          command = "/run/setCutForAGivenParticle proton 0.001 mm";
          UImanager->ApplyCommand(command);
          cout << command << endl;

          command = "/run/initialize";
          UImanager->ApplyCommand(command);
          cout << command << endl;



    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
