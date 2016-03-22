/**
 * @brief Main executable to run over input mini flat ntuples,
 * read events in an object-oriented way and
 * run the event analysis class that derives from Analysis.
 * The main objective is to be minimal and to allow users to run this part
 * of the code in their laptops.
 *
 * @author Danilo Enoque Ferreira de Lima <dferreir@mail.cern.ch>
 */

#include "TChain.h"
#include <iostream>

#include "STop_Ntuple/CollectionTree.h"
#include "STop_Ntuple/ChainHelper.h"

#include <fstream>

#include "TROOT.h"
#include "TInterpreter.h"

#include <sstream>

void help()
{
  cout << "  Options:"                          << endl;
  cout << "  -n number of events to process"    << endl;
  cout << "     defaults: -1 (all events)"      << endl;

  cout << "  -k number of events to skip"       << endl;
  cout << "     defaults: 0"                    << endl;

  cout << "  -i input (file, list, or dir)"     << endl;
  cout << "     defaults: ''"                   << endl;

  cout << "  -o output (file.root)"             << endl;
  cout << "     defaults: 'output.root'"        << endl;

  cout << "  -v verbose ( 0 or 1 )"             << endl;
  cout << "     default: 0 (false)  "           << endl;

  cout << "  -e event display file (file.eps)"  << endl;
  cout << "     default: 'event_display.eps'"   << endl;

  cout << "  -isMultiJet (0 or 1)"              << endl;
  cout << "     default: '0 (false)'"           << endl;

  cout << "  -h print this help"                << endl;
}

int main(int argc, char **argv) {

  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  string sample;
  string input;
  string output = "output.root";
  string event_display = "event_display.eps";
  int verbose = 0;
  bool verb = false;
  int isMultiJetMC = 0;
  cout << "SusyDiag0L" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if      (strcmp(argv[i], "-n") == 0) nEvt = atoi(argv[++i]);
    else if (strcmp(argv[i], "-k") == 0) nSkip = atoi(argv[++i]);
    else if (strcmp(argv[i], "-i") == 0) input = argv[++i];
    else if (strcmp(argv[i], "-o") == 0) output = argv[++i];
    else if (strcmp(argv[i], "-e") == 0) event_display = argv[++i];
    else if (strcmp(argv[i], "-v") == 0) verbose = atoi(argv[++i]);
    else if (strcmp(argv[i], "-isMultiJet") == 0) isMultiJetMC = atoi(argv[++i]);
    else {
      help();
      return 0;
    }
  }

  if ( verbose == 1 ) verb = true;

  TChain* chain = new TChain("HFntupleNONE");
  ChainHelper::addInput(chain, input, dbg>0);
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  //  ROOT::Cintex::Cintex::Enable();

  gROOT->ProcessLine("#include <vector>");
  //gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  //  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  //  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");

  // input files

  // get output files
  CollectionTree ana(chain);
  ana.SetEvtDisplayName(event_display);
  ana.SetOutputName(output);
  ana.BookHistos();
  ana.SetVerbose(verb);
  ana.Loop();
  ana.WriteHistos();

  
  return 0;
}

