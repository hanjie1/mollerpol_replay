/* Replay script for GEM analysis data - Optimized for Parallel Processing */
#if defined(__CLING__)
R__ADD_INCLUDE_PATH(/work/halla/moller12gev/asar/Moller_Polarimeter/PolGEM/mollerpol_analyzer)
R__LOAD_LIBRARY(/work/halla/moller12gev/asar/Moller_Polarimeter/PolGEM/mollerpol_analyzer/build/libMollerPol.so)
#endif

#include "gem/MollerPolGEMSpectrometer.h"
#include "gem/MollerPolGEMSpectrometerTracker.h"


#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "TString.h"
#include "TDatime.h"

#include "THaDetector.h"
// Note: Function signature matches the shell script call
void replay_moller_gem_multi(UInt_t runnum=1859, UInt_t segment=0, UInt_t firstevent=0, Int_t nevents=-1){

// --- 1. Setup apparatus ---
// A tracking apparatus must derive from THaSpectrometer
auto* moller =
    new MollerPolGEMSpectrometer(
        "mollerPol",
        "Moller Polarimeter GEM spectrometer"
    );

auto* polgem =
    new MollerPolGEMSpectrometerTracker(
        "polgem",
        "8-module PolGEM setup",
        moller
    );

if (moller->AddDetector(polgem) != 0) {
    printf("ERROR: Failed to add PolGEM tracker\n");
    delete moller;
    return;
}

gHaApps->Add(moller);

printf("Number of detectors in mollerPol: %d\n",
       moller->GetNumDets());
moller->Print("DET");
        // --- 2. Setup Analyzer ---
    THaAnalyzer* analyzer = new THaAnalyzer;
    THaEvent* event = new THaEvent;

    // --- 3. Handle Input File (Single Segment) ---
    TString data_prefix = gSystem->Getenv("DATA_DIR");
    TString codafilename;
    codafilename.Form("%s/moller_ssp_%d.evio.%d", data_prefix.Data(), runnum, segment);

    if( gSystem->AccessPathName(codafilename.Data()) ){
        printf("ERROR: CODA file not found: %s\n", codafilename.Data());
        return;
    }

    THaRun* run = new THaRun(codafilename.Data());
    run->SetDate(TDatime());
    run->SetNumber(runnum);
    run->SetFirstEvent(firstevent);
    if(nevents > 0) run->SetLastEvent(nevents);
    run->SetDataRequired(0);

    // --- 4. Handle Output File ---
    TString out_prefix = gSystem->Getenv("OUT_DIR");
    TString outfilename;
    // Unique name so parallel jobs don't clash
    outfilename.Form("%s/moller_polgem_run%d_seg%d.root", out_prefix.Data(), runnum, segment);

    analyzer->SetVerbosity(2);
    analyzer->SetMarkInterval(1000);
    analyzer->EnableBenchmarks();

    analyzer->SetEvent(event);
    analyzer->SetOutFile(outfilename.Data());

    // Optional summary file per segment
    TString summaryName;
    summaryName.Form("summary_%d_%d.log", runnum, segment);
    analyzer->SetSummaryFile(summaryName.Data());

    // --- 5. ODEF Setup ---
    TString replay_prefix = gSystem->Getenv("MOLLER_REPLAY");
    TString odef_filename = "DB/replay_mollerPol_polgem.odef";
    TString odef_path = replay_prefix + "/" + odef_filename;

    analyzer->SetOdefFile(odef_path.Data());

    // --- 6. Process ---
    printf("Processing Run %d, Segment %d...\n", runnum, segment);
    analyzer->Process(run);

    // Cleanup to prevent memory leaks in long loops
    delete analyzer;
    delete event;
    delete run;
}
