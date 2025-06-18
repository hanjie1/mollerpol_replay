#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"


void replay(){

  gSystem->Load("libMollerPol");

  MollerPolApparatus *mol = new MollerPolApparatus("M", "Moller polarimeter apparatus");

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  gHaApps->Add(mol);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;
  
  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
//  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  THaRun* run = new THaRun( "data/mollerpol_test_199.evio.0" );
  run->SetLastEvent(-1);

  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(4);
  analyzer->SetOdefFile("replay.odef" );
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( "Rootfiles/mollerpol_test_199.root" );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("summary_199.log"); // optional
  
  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
