#define CollectionTree_cxx
#include "STop_Ntuple/CollectionTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CollectionTree::Loop()
{
  if (verbose) std::cout << "looping " << std::endl;
        
//   In a ROOT session, you can do:
//      Root > .L CollectionTree.C
//      Root > CollectionTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
 
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if ( jentry%(nentries/20) == 0 ) std::cout << "Processed " << jentry << " out of " << nentries << std::endl;

     if(verbose) std::cout << "entry " << jentry << std::endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      ProcessEvent();
   }
}

void CollectionTree::ProcessEvent() {

  //--------------------------------------------------//
  //    Build up Vectors of TLorentzVectors for Jets  //
  //    and other basic objects                       //
  //--------------------------------------------------//

  build_objects();
 
  if (verbose) std::cout << "event met " << MET_TLV.Pt() << std::endl;

  //--------------------------------------------//
  //  Event cleaning, pre-selection and trigger
  //--------------------------------------------//

  if(!passEventCleaning()) return;
  
  if( MET_TLV.Pt() < 200.) return;
  
  //---------------------------------------//
  //          cut out multijets
  //---------------------------------------//

  if(!PassFiducial()) return;

  //---------------------------------------//
  
  int nmuon = 0;
  int nel = 0;

  mu_index.resize(0);
  mu_index.resize(1);

  for ( uint i=0; i<nMu; i++ ) {
    if ( Mu_TLV.at(i).Pt() > 10.0 ) {
      nmuon++;
      mu_index.push_back(i);
    }
  }
  for ( uint i=0; i<nEl; i++ ) {
    if ( El_TLV.at(i).Pt() > 10.0 ) {
      nel++;
      el_index.push_back(i);
    }
  }

  //--------------------------------------//

  bjets_index.resize(0);
  for ( uint i=0; i<nJets; i++ ) {
    if ( Jet_isB.at(i) ) {
      bjets_index.push_back(i);
    }
  }
  int nbjets = bjets_index.size();

  fiducialjets.resize(0);
  for ( uint i=0; i<nJets; i++ ) {
    if ( fabs(Jet_TLV.at(i).Eta()) < 2.8 ) {
      fiducialjets.push_back(i);
    }
  }

  //--------------------------------------//
  //      see if jets pass pt cut         //
  //--------------------------------------//

  highptjet.resize(0);
  vector<float> jet_ptcut;
  jet_ptcut.resize(0);
  jet_ptcut.push_back(60.);
  jet_ptcut.push_back(60.);
  jet_ptcut.push_back(20.);
  jet_ptcut.push_back(20.);
  jet_ptcut.push_back(20.);
  jet_ptcut.push_back(20.);

  for ( uint i=0; i<fiducialjets.size(); i++ ) {
    int ijet = fiducialjets.at(i);
    if ( highptjet.size() < jet_ptcut.size() ) {
      if ( Jet_TLV.at(ijet).Pt() > jet_ptcut.at(highptjet.size()) ) {
        highptjet.push_back(ijet);
      }
    }
    else if ( Jet_TLV.at(ijet).Pt() > 20. ) {
      highptjet.push_back(ijet);
    }
  }

  int n_highpt_jets = highptjet.size();

  //-----------------------------------------//
  //            all hadronic case            //
  //-----------------------------------------//

  if ( nEl_baseline+nMu_baseline == 0 ) { 

    if ( nbjets >= 2 && n_highpt_jets >= 6 ) {
      run2b_AllHad_ThrustAna();
    }

  }
  return;
}

void CollectionTree::run2b_AllHad_ThrustAna() {

  //  double weight = 1.0;
  double weight = pileupweight*AnalysisWeight*XSecWeight;
  isMC = true;
  
  //----------------------------------//
  //    Calculate Thrust Direction    //
  //----------------------------------//

  FindThrustDir();

  //----------------------------------//
  //   SortB jets by P along Thrust   //
  //----------------------------------//

  if ( verbose ) std::cout << "sort b by p thrust " << std::endl;
  SortB_by_PThrust();

  //----------------------------------//
  //          Find ISR Jets           //
  //----------------------------------//

  if ( verbose ) std::cout << "find isr by p thrust " << std::endl;
  FindISR_by_PThrust();

  //----------------------------------//
  //   Look for Hadronic Top with     //
  //  large negative p along thrust   //
  //----------------------------------//

  if ( verbose ) std::cout << "ttbar had reco " << std::endl;
  Found_Neg_PThrust_Had_Top = TTbar_Had_Top_Reconstruction();
  if ( verbose ) std::cout << "ttbar had reco " << std::endl;

  if ( Found_Neg_PThrust_Had_Top ) {
    TTbar_Had_Top_vec = Jet_TLV.at(bjets_index_sorted.back())+whad_vec;
  }

  if ( verbose ) std::cout << "plot 0 " << std::endl;
  plot2b_AllHad_ThrustAna(weight, 0);

  if ( !isMC ) return;

  //---------------------------------------------//
  //    find what event decay mode was in truth  //
  //---------------------------------------------//

  bool top1_ishad = false;
  if ( abs(pdgId1_1) < 6 ) {
    top1_ishad = true;
  }
  bool top1_istau = false;
  if ( abs(pdgId1_1) == 15 ) {
    top1_istau = true;
  }

  bool top2_ishad = false;
  if ( abs(pdgId2_1) < 6 ) {
    top2_ishad = true;
  }
  bool top2_istau = false;
  if ( abs(pdgId2_1) == 15 ) {
    top2_istau = true;
  }

  //-------------------------------------------------------//
  //  plot plots corresponding to decay mode for reference //
  //-------------------------------------------------------//
  
  if ( top1_ishad && top2_ishad ) {
    plot2b_AllHad_ThrustAna(weight, 8);
  }
  else {
    plot2b_AllHad_ThrustAna(weight, 9);
  }

  //-------------------------------------------------------------------------//
  //  draw crude 2D event displays for both truth and reconstructed objects  //
  //-------------------------------------------------------------------------//

  if ( ndisplay < 200 ) DrawEvtDisplay(weight);
  ndisplay++;

  if (verbose) std::cout << "number of ttjets " << n_TTjets << std::endl;

  //---------------------------------------------------//
  //  Begin actual cutflow on reconstructed variables  //
  //---------------------------------------------------//

  // at least 4 jets on TT side
  if ( n_TTjets <= 3 ) return;
  plot2b_AllHad_ThrustAna(weight, 1);

  // at least 1 bjet on TT side
  if ( n_b_TTjets <= 0 ) return;
  plot2b_AllHad_ThrustAna(weight, 2);

  // ISR pt cut
  if ( ISRjets_vec.Pt() < 400. ) return;
  plot2b_AllHad_ThrustAna(weight, 3);

  // ISR, met deltaphi cut
  if ( fabs(ISRjets_vec.DeltaPhi(MET_TLV)) < 3.0 ) return;
  plot2b_AllHad_ThrustAna(weight, 4);

  // no top candidate reconstructed with large negative p-along-thrust
  if ( bjets_index.size() <= 3     && Found_Neg_PThrust_Had_Top && 
       TTbar_Had_Top_vec.M() < 220 && 
       CalculatePAlongThrust( &TTbar_Had_Top_vec) < -200. ) return;
  plot2b_AllHad_ThrustAna(weight, 5);

  // 300 GeV optimization  based on TT system / ISR system ratios                           
  if ( TTjets_vec.Pt()/ISRjets_vec.Pt() > 0.4 && TTjets_vec.Pt()/ISRjets_vec.Pt() < 0.65) {
    plot2b_AllHad_ThrustAna(weight, 6);
  }
  // 350 GeV optimization                                     

  if ( TTjets_vec.Pt()/ISRjets_vec.Pt() > 0.35 && TTjets_vec.Pt()/ISRjets_vec.Pt() < 0.55){
    // plot2b_AllHad_ThrustAna(weight, 5);
  }
  // 500 GeV optimization                                        
  if ( TTjets_vec.Pt()/ISRjets_vec.Pt() < 0.4 ){
    // plot2b_AllHad_ThrustAna(weight, 6);
  }

  return;

}

double CollectionTree::CalcHT(double jetPtCut)
{
  double ht = 0;
  for (uint i = 0; i < Jet_TLV.size(); i++) {
    double pt = Jet_TLV.at(i).Pt();
    if (pt > jetPtCut) ht += pt;
  }
  return ht;
}

double CollectionTree::CalcThrust(){

  if (verbose) std::cout<<"calculating thrust"<<std::endl;

  double totpt = 0.0;
  for ( uint i=0; i < fiducialjets.size(); i++ ) {
    int ijet = fiducialjets.at(i);
    totpt += Jet_TLV.at(ijet).Pt();
  }

  totpt += MET_TLV.Pt();

  double Thrust = thrustFunction(thrust_vec.X(),
                                 thrust_vec.Y())/totpt;

  return Thrust;

}

double CollectionTree::thrustFunction(double x, double y) {
  double tot = 0.0;

  for ( uint i=0; i < fiducialjets.size(); i++ ) {
    int ijet = fiducialjets.at(i);
    tot += fabs(Jet_TLV.at(ijet).Px()*x + Jet_TLV.at(ijet).Py()*y);
  }

  tot += fabs(MET_TLV.Px()*x + MET_TLV.Py()*y);
  return tot/sqrt(x*x+y*y);
}

void CollectionTree::FindThrustDir(){

  //----------------------------------------------------//                                                  
  //   Calculate the thrust vector in the px,py plane
  //   The thrust vector is the direction that maximizes the momentum along it
  //   find 2D-vector "n" such that Sum ( 3momentum dot n ), summed over all jets and met is maximized

  //  inputs: std::vector<TLorentzVector> Jet_TLV    // list of jet TLV, size njets                          
  //          std::vector<int> fiducialjets;         // list of jets with |eta|<2.8, size nfiducialjets      
  //          TLorentzVector MET_TLV;                // only the 2D vector of the MET, pz and e are 0

  //  outputs: TLorentzVector thrust_vec             // only the 2D vector of the thrust, positive is in the same direction as met
  //           TLorentzVector thrust_perp_vec        // 2D vector of the direction perpendicular to thrust, positive side is defined by 
  //                                                    positive cross product with thrust_vec

  if (verbose) std::cout <<"finding thrust"<<std::endl;

  //----------------------------------//
  //   Write sum as a TF1 formula     //
  //----------------------------------//

  TString Formula("abs(cos(x)*[0]+sin(x)*[1])+");
  TString tmp;
  for( uint i=1; i < fiducialjets.size(); i++ ) {
    tmp = TString(Form("abs(cos(x)*[%u]+sin(x)*[%u])+",2*i,2*i+1));
    Formula = Formula+tmp;
  }

  uint sizejets = fiducialjets.size();
  tmp = TString(Form("abs(cos(x)*[%u]+sin(x)*[%u])",
                     2*sizejets,
                     2*sizejets+1));

  Formula = Formula+tmp;

  TF1 fcn("ThrustFunction",Formula.Data(),0.0,TMath::Pi());

  for ( uint i=0; i < fiducialjets.size(); i++ ) {
    fcn.FixParameter(2*i,   Jet_TLV.at(fiducialjets.at(i)).Px());
    fcn.FixParameter(2*i+1, Jet_TLV.at(fiducialjets.at(i)).Py());
  }

  fcn.FixParameter(2*fiducialjets.size(),   MET_TLV.Px());
  fcn.FixParameter(2*fiducialjets.size()+1, MET_TLV.Py());

  //---------------------------------//
  //      Find Maximum of TF1        //
  //---------------------------------//

  double bestX, bestY, bestPhi;
  fcn.Eval(0.5);
  bestPhi = fcn.GetMaximumX();
  bestX = cos(bestPhi);
  bestY = sin(bestPhi);

  //--------------------------------------------------------------//
  //  define positive thrust as positive dot product with MET     //
  //--------------------------------------------------------------//

  if ( (MET_TLV.Px()*bestX + MET_TLV.Py()*bestY) < 0 ) {
    bestX = -bestX;
    bestY = -bestY;
  }

  thrust_vec.SetXYZT(bestX,bestY,0,0);

  //--------------------------------------------------------------//   
  //            Find Direction Perpendicular to Thrust            //   
  //--------------------------------------------------------------//   

  double bestPerpX, bestPerpY;
  bestPerpX = 1.0 - bestX*bestX;
  bestPerpY = 0.0 - bestX*bestY;
  double totPerp = sqrt(bestPerpX*bestPerpX+bestPerpY*bestPerpY);
  bestPerpX = bestPerpX/totPerp;
  bestPerpY = bestPerpY/totPerp;

  //------------------------------------------------------------------------//   
  //  define positive perp thrust as positive cross product with thrust     //   
  //------------------------------------------------------------------------//   

  if ( (bestX*bestPerpY - bestPerpX*bestY) < 0 ) {
    bestPerpX = - bestPerpX;
    bestPerpY = - bestPerpY;
  }

  thrust_perp_vec.SetXYZT(bestPerpX,bestPerpY,0,0);

  if (verbose) std::cout <<"founding thrust"<<std::endl;

  return;
}

float CollectionTree::CalculatePAlongThrust(double px, double py) {
  return thrust_vec.Px()*px + thrust_vec.Py()*py;
}

float CollectionTree::CalculatePAlongThrust( TLorentzVector *TLV ) {
  return CalculatePAlongThrust(TLV->Px(),TLV->Py());
}

float CollectionTree::CalculatePtAlongThrust(double px, double py) {
  return thrust_perp_vec.Px()*px + thrust_perp_vec.Py()*py;
}

float CollectionTree::CalculatePtAlongThrust( TLorentzVector *TLV ) {
  return CalculatePtAlongThrust(TLV->Px(),TLV->Py());
}

void CollectionTree::SortB_by_PThrust() {
  
  //----------------------------------------------------//
  //  Sorts the b-tagged jets by p-along-thrust
  //  if there are >=3 bjets the 2 with the highest p-along-thrust
  //  are probably the b's coming from the TT system
  //  b's with the large negative p-along thrust is probably ISR
  //
  //  inputs: std::vector<TLorentzVector> Jet_TLV    // list of jet TLV, size njets                    
  //          std::vector<bool>           Jet_isB    // is jet a bjets,  size njets                    
  //          std::vector<int> bjets_index;          // list of bjet_indexes 
  //          std::vector<int> fiducialjets;         // list of jets with |eta|<2.8, size nfiducialjets
  //                                                                                                   
  //  outputs: std::vector<int> bjets_index_sorted   // list of bjet_indexes sorted by p along thrust
  //           std::vector<int> bjets_pthrust_sorted // list of bjet p along thrust sorted by p along thrust

  bjets_index_sorted.resize(0);
  bjets_pthrust_sorted.resize(0);

  for ( uint i=0; i < bjets_index.size(); i++ ) {
    if ( verbose ) std::cout << "looping sort b" << std::endl;

    //-------------------------------------------------------------------//
    //              Calculate bjet p along thrust                        //
    //-------------------------------------------------------------------//

    int    bjet_i     = bjets_index.at(i);
    double b_pthrust  = CalculatePAlongThrust( Jet_TLV.at(bjet_i).Px(), Jet_TLV.at(bjet_i).Py() );

    if ( verbose ) std::cout <<"looping sort b " << b_pthrust << std::endl;

    //-------------------------------------------------------------------//
    //  insert b jet into sorted vector according to its p along thrust  //
    //-------------------------------------------------------------------//

    std::vector<float>::iterator it;
    it = bjets_pthrust_sorted.begin();
    std::vector<int>::iterator it2;
    it2 = bjets_index_sorted.begin();

    if ( bjets_pthrust_sorted.size() != 0 ) {
      for ( ; it != bjets_pthrust_sorted.end(); it++ ) {
	if ( *it <= b_pthrust ) {
	  it2 = bjets_index_sorted.insert(  it2,bjet_i);
	  it  = bjets_pthrust_sorted.insert(it, b_pthrust);
	  break;
	}
	else if ( it+1 == bjets_pthrust_sorted.end() ) {
	  bjets_index_sorted.push_back(  bjet_i );
	  bjets_pthrust_sorted.push_back(b_pthrust);
	  break;
	}
	it2++;
      }
    }
    else {
      bjets_index_sorted.push_back(  bjet_i );
      bjets_pthrust_sorted.push_back(b_pthrust);
    }
  }

  // output list to see if the order was correct if on debug mode
  
  if ( verbose ) {
    for ( uint i=0; i < bjets_pthrust_sorted.size(); i++ ) {
      std::cout << bjets_pthrust_sorted.at(i) << " " ;
    }
    for ( uint i=0; i < bjets_index_sorted.size(); i++ ) {
      std::cout << bjets_index_sorted.at(i) << " " ;
    }
  }

  return;
  
}

void CollectionTree::FindISR_by_PThrust() {
  
  //----------------------------------------------------// 
  //  Decide which jets come from the ISR system and which 
  //  comes from the TT system.   Jets with large negative p along thrust
  //  and jets which form large mass systems with b-jets are considered ISR
  //
  //  inputs: std::vector<TLorentzVector> Jet_TLV    // list of jet TLV, size njets
  //          std::vector<bool>           Jet_isB    // is jet a bjets,  size njets
  //          std::vector<int> bjets_index_sorted;   // list of bjet_indexes sorted by p along thrust, size nbjets
  //          std::vector<int> fiducialjets;         // list of jets with |eta|<2.8, size nfiducialjets           
  //                                                                                                              
  //  outputs: int n_ISRjets, n_TTjets               // number of ISR jets and number of TT jets found
  //           int n_b_ISRjets, n_b_ISRjets          // number of bjets considered ISR and number of bjets considered TT 
  //           TLorentzVector ISRjets_vec            // total ISR TLV (all ISR jets combined)
  //           TLorentzVector TTjets_vec             // total TT system TLV ( all TT jets combined)
  //           std::vector<int> potentialISR_index   // index of every ISR jet. feeds into vector<int> fiducialjets

  potentialISR_index.resize(0);

  // get b-jets with the most positive p-along-thrust (positive side is same side as MET)
  int b_pthrust_max_i = bjets_index_sorted.at(0);
  int b_pthrust_2nd_max_i = bjets_index_sorted.at(1);

  ISRjets_vec.SetPxPyPzE(0,0,0,0);
  TTjets_vec.SetPxPyPzE(0,0,0,0);

  n_TTjets=0;
  n_b_TTjets=0;
  n_ISRjets=0;
  n_b_ISRjets=0;

  //-----------------------------------------------------------//
  //      find ISR jets from list of fiducial jets             //
  //-----------------------------------------------------------//

  for ( uint i=0; i < fiducialjets.size(); i++ ) {

    int ijet = fiducialjets.at(i);

    bool jet_isb   = false;
    bool jet_isISR = false;

    jet_isb = Jet_isB.at(ijet);
    
    //-----------------------------------------------------------//
    //  if jet is too negative in p along thrust then it is ISR  //
    //  if jet+bjet's mass is too large then it is an ISR        //
    //-----------------------------------------------------------//

    double pthrust  = CalculatePAlongThrust( Jet_TLV.at(ijet).Px(), Jet_TLV.at(ijet).Py() );

    double B1mass   = (Jet_TLV.at(ijet)+Jet_TLV.at(b_pthrust_max_i)).M();
    double B2mass   = (Jet_TLV.at(ijet)+Jet_TLV.at(b_pthrust_2nd_max_i)).M();

    double maxBmass = TMath::Max(B1mass,B2mass);
    double minBmass = TMath::Min(B1mass,B2mass);

    //    if ( pthrust < -100.0 || maxBmass > 320. || minBmass > 170.) { //alternative definition, not sure which is best
    if ( pthrust < -100.0 || minBmass > 170.) {
      jet_isISR = true;
    }

    if ( jet_isISR ) {
      potentialISR_index.push_back(i);

      n_ISRjets++;
      ISRjets_vec += Jet_TLV.at(ijet);
      if ( jet_isb ) n_b_ISRjets++;
    }
    else             {

      TTjets_vec +=  Jet_TLV.at(ijet);
      n_TTjets++;
      if ( jet_isb ) n_b_TTjets++;
    }

  }

  return;
}

bool CollectionTree::TTbar_Had_Top_Reconstruction() {

  //----------------------------------------------------//
  //  Look to reconstruct a fully hadronic top in a cone around
  //  the b with the most negative p-along-thrust.  TTbar bkg events
  //  can have a full, fairly boosted hadronic top with large negative 
  //  p-along-thrust.  Stop events on the other hand cannot as the ISR
  //  push both tops to the positive hemisphere
  //  i.e. if we find a top here and the top p-along-thrust is very negative
  //  then the event is probably ttbar bkg and not stop.  
  //  This reconstruction method does not work for stop because stop tops are not
  //  boosted enough.  For stop, the results will probably be some random soft jets get
  //  clustered together and the top peak is not reconstructed
  //
  //  inputs: std::vector<TLorentzVector> Jet_TLV    // list of jet TLV, size njets
  //          std::vector<bool>           Jet_isB    // is jet a bjets,  size njets
  //          std::vector<int> bjets_index;          // list of bjet indexes, size nbjets
  //          std::vector<int> bjets_index_sorted;   // list of bjet_indexes sorted by p along thrust, size nbjets
  //          std::vector<int> fiducialjets;         // list of jets with |eta|<2.8, size nfiducialjets
  //          
  //  outputs: found_W                               // whether we found hadronic W candidate within cone
  //           TLorentzVector TTbar_Had_Top_vec      // TLV of best hadronic top candidate
  //           TLorentzVector whad_vec               // TLV of best hadronic W candidate

  TTbar_Had_Top_vec = TLorentzVector(0,0,0,0);
  
  //-------------------------------------------------//
  //           Look for W near top B                 //
  //-------------------------------------------------//

  if ( verbose ) std::cout << "look for W near B " << std::endl;

  whad_vec = TLorentzVector(0,0,0,0);

  double Top_W_deltaR_cut = 1.5;
  bool found_W = false;
  
  double w_polemass = 80.39;

  double minWDeltaM=1000000.;
  int minWDeltaM_i = -1;
  int minWDeltaM_j = -1;
  
  double minWDeltaR=1000000.;
  int minWDeltaR_i = -1;
  int minWDeltaR_j = -1;

  if ( verbose ) std::cout << bjets_index.size() << " " << bjets_index_sorted.size() << std::endl;
  if ( verbose ) std::cout << fiducialjets.size() << std::endl;

  //-----------------------------------------------------//
  //   Get the two most negative p-along-thrust b jets   //
  //-----------------------------------------------------//

  int b_min_pthrust_i     = bjets_index_sorted.at(bjets_index_sorted.size()-1);
  int b_2nd_min_pthrust_i = bjets_index_sorted.at(bjets_index_sorted.size()-2);

  if ( verbose ) std::cout << b_min_pthrust_i << std::endl;
  if ( verbose ) std::cout << b_2nd_min_pthrust_i << std::endl;

  int count = 0;

  double bjet_mass = Jet_TLV.at(b_min_pthrust_i).M();

  if ( verbose ) std::cout << "look for W near B " << std::endl;
  if ( verbose ) std::cout << bjets_index_sorted.size() << std::endl;
  if ( verbose ) std::cout << bjet_mass << std::endl;

  //-------------------------------------//
  //            2 btag version           //
  //-------------------------------------// 

  if ( bjets_index.size() == 2 ) {


    while ( !found_W && count <= 1 ) {
      
      for ( uint i=0; i<fiducialjets.size(); i++ ) {
	
	int ijet = fiducialjets.at(i);

	if ( verbose ) std::cout << ijet << std::endl;

	if ( Jet_isB.at(ijet) ) continue;
	
	//----------------------------------------------------------------------------------------//
	// if b jet is massive only add one another jet as bjet maybe 2 jets from a boosted top   //
	//----------------------------------------------------------------------------------------//    

	if ( bjet_mass > 60.0 ) {

	  double B_ij_deltaR = Jet_TLV.at(b_min_pthrust_i).DeltaR(Jet_TLV.at(ijet));
	  
	  if ( B_ij_deltaR < Top_W_deltaR_cut && minWDeltaR > B_ij_deltaR ) {
	    found_W = true;

	    minWDeltaR   = B_ij_deltaR;
	    minWDeltaR_i = ijet;
	  }
	}
	else {
	  if ( i+1 == fiducialjets.size() ) continue;
	  for ( uint j=i+1; j<fiducialjets.size(); j++ ) {
	    
	    int jjet = fiducialjets.at(j);
	    if ( verbose ) std::cout << jjet << std::endl;
	    if ( Jet_isB.at(jjet) ) continue;
	    
	    double B_ij_deltaR = Jet_TLV.at(b_min_pthrust_i).DeltaR(Jet_TLV.at(ijet)+Jet_TLV.at(jjet));

	    //----------------------------------------------------//
	    //   Find all jet pairs that fall within deltaR cone  //
	    //----------------------------------------------------//
	    
	    // Note, one or both individual jets may be outside the cone.  
	    //       Only the combined 4vector of the jet pair must be within the cone
	    
	    if ( B_ij_deltaR < Top_W_deltaR_cut ) {
	      found_W = true;
	      
	      double ij_mass   = (Jet_TLV.at(ijet)+Jet_TLV.at(jjet)).M();
	      
	      //  pick out jet pair that has the mass that is closest to the w mass

	      if ( minWDeltaM > fabs( ij_mass - w_polemass ) ) {
		minWDeltaM   = fabs( ij_mass - w_polemass );
		minWDeltaM_i = ijet;
		minWDeltaM_j = jjet;
	      }
	    }
	  }
	}
      }

      //--------------------------------------------//
      // if cannot find top, search in a wider cone //
      //--------------------------------------------//

      if ( !found_W ) {
	Top_W_deltaR_cut+=0.5;
	count++;
      }
    }

    if ( found_W ) {
      if ( verbose ) std::cout << "found W" << std::endl;

      if ( bjet_mass > 60.0 ) {
	whad_vec = Jet_TLV.at(minWDeltaR_i);
      }
      else {
	whad_vec = Jet_TLV.at(minWDeltaM_i)+Jet_TLV.at(minWDeltaM_j);
      }

      if ( verbose ) std::cout << "returning top reco" << std::endl;
      return true;
    }
    else { return false; }
  }

  //-------------------------------------// 
  //           3 btag version            //
  //-------------------------------------// 

  if ( bjets_index.size() == 3 ) {

    while ( !found_W && count <= 1 ) {

      if ( bjet_mass > 60.0 ) {
	double B_ij_deltaR = Jet_TLV.at(b_min_pthrust_i).DeltaR(Jet_TLV.at(b_2nd_min_pthrust_i));

	if ( B_ij_deltaR < Top_W_deltaR_cut && minWDeltaR > B_ij_deltaR ) {
	  found_W = true;
	  minWDeltaR   = B_ij_deltaR;
	  minWDeltaR_i = b_2nd_min_pthrust_i;
	}
      }
      else {
	for ( uint i=0; i<fiducialjets.size(); i++ ) {
	  
	  int ijet = fiducialjets.at(i);
	  if ( ijet == bjets_index_sorted.at(0) ) continue;
	  
	  double B_ij_deltaR = (Jet_TLV.at(b_min_pthrust_i)+Jet_TLV.at(b_2nd_min_pthrust_i)).DeltaR(Jet_TLV.at(ijet));

	  if ( B_ij_deltaR < Top_W_deltaR_cut ) {
	    found_W = true;
	    
	    if ( minWDeltaR > B_ij_deltaR ) {
	      minWDeltaM   = B_ij_deltaR;
	      minWDeltaM_i = b_2nd_min_pthrust_i;
	      minWDeltaM_j = ijet;
	    }
	  }
	}
      }
      if ( !found_W ) {
        Top_W_deltaR_cut+=0.5;
        count++;
      }
    }
    if ( found_W ) {
      if ( bjet_mass > 60.0 ) {
	whad_vec = Jet_TLV.at(minWDeltaR_i);
      }
      else {
	whad_vec = Jet_TLV.at(minWDeltaM_i)+Jet_TLV.at(minWDeltaM_j);
      }
      return true;
    }
    else {
      return false;
    }
  }

  // if more than 3 bjets just accept the event as stop
  return false;
}

void CollectionTree::build_objects() {

  //--------------------------//
  //        outputs
  //--------------------------//

  Jet_TLV.resize(0); // signal jets
  Jet_isB.resize(0); // whether jet is a b or not

  El_TLV.resize(0); // signal electrons
  Mu_TLV.resize(0); // signal muons

  MET_TLV.SetPxPyPzE(0,0,0,0);

  //-------------------------//

  //  double btag_cut = -0.0436; // 70% working point
  //  double btag_cut = -0.4434; // 77% working point
  double btag_cut = -0.7887; // 85% working point

  // MET is in GeV
  // jet, muon, el are in MeV

  for ( int i=0; i<nJets; i++ ) {
    Jet_TLV.push_back(TLorentzVector(jet_px->at(i)/1000.,jet_py->at(i)/1000.,
				     jet_pz->at(i)/1000.,jet_e->at(i)/1000.));
    if ( jet_MV2c20->at(i) > btag_cut ) {
      Jet_isB.push_back(true);
    }
    else {
      Jet_isB.push_back(false);
    }
  }

  for ( int i=0; i<nEl; i++ ) {
    if ( verbose ) std::cout << "electron px " << el_px->at(i) << std::endl;
    El_TLV.push_back(TLorentzVector(el_px->at(i)/1000.,el_py->at(i)/1000.,
				    el_pz->at(i)/1000.,el_e->at(i)/1000.));
    
  }
  for ( int i=0; i<nMu;i++ ) {
    if ( verbose ) std::cout << "muon px " << mu_px->at(i) << std::endl;
    Mu_TLV.push_back(TLorentzVector(mu_px->at(i)/1000.,mu_py->at(i)/1000.,
				    mu_pz->at(i)/1000.,mu_e->at(i)/1000.));

  }

  MET_TLV.SetPxPyPzE(MET_px,MET_py,0,0);

  return;
}

bool CollectionTree::passEventCleaning() {

  /*
    444                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(2, "GRL");
    445                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(3, "LAr and Tile error");
    446                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(4, "Trigger");
    447                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(5, "Primary vertex exists");
    448                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(6, "Clean jet");
    449                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(7, "Cosmic");
    450                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(9, "Lepton veto");
    451                m_cutflow[0][systI]->GetXaxis()->SetBinLabel(8, "Clean Muon");
  */

  if( !passGRL )               return false;
  n_pass_grl++;

  if ( !Cleaning ) return false;

  //  if( !LArError )            return false;
  //  n_pass_LAr++;

  //  if( !TileTrip )           return false;
  //  n_pass_tileErr++;

  if ( !HLT_xe70_tc_lcw ) return false;

  if ( numVtx <= 0 ) return false;

  if ( !JetCleaning ) return false;

  if ( nMu_cosmic > 0 ) return false;

  if ( !MuonCleaning ) return false;

  return true;

}

bool CollectionTree::PassFiducial() {

  if ( nJets < 3 ) return false;
  
  //--------------------------------------------//
  //            Cut on Delta Phi Met            //
  //--------------------------------------------//

  double deltaPhi;

  if ( eT_miss_track < 30 ) return false;

  deltaPhi = fabs(dPhi_met_trackmet);
  if ( deltaPhi > TMath::Pi()/2.0 ) return false;

  deltaPhi = fabs(Jet_TLV.at(0).DeltaPhi(MET_TLV));
  if ( deltaPhi < 0.50 ) return false;

  deltaPhi = fabs(Jet_TLV.at(1).DeltaPhi(MET_TLV));
  if ( deltaPhi < 0.40 ) return false;

  deltaPhi = fabs(Jet_TLV.at(2).DeltaPhi(MET_TLV));
  //  if ( deltaPhi < 0.40 ) return;               

  return true;

  //---------------------------------------------//

}

void CollectionTree::plot2b_AllHad_ThrustAna(double weight, int icut) {

  if (verbose) std::cout << "plot2b all had thrust" <<std::endl;
  
  double AllThrust_T = CalcThrust();
  Ana_AllHad_2b_Thrust_Thrust[icut].Fill(AllThrust_T,weight);
  
  double deltaPhi_Thrust_Met = fabs(thrust_vec.DeltaPhi(MET_TLV));
  Ana_AllHad_2b_Thrust_Thrust_Met_dPhi[icut].Fill(deltaPhi_Thrust_Met,weight);

  for ( uint i=0; i < fiducialjets.size(); i++ ) {
    
    bool jet_isISR = false;
    for ( uint j=0; j < potentialISR_index.size(); j++ ) {
      if ( i == potentialISR_index.at(j) ) {
	jet_isISR = true;
      }
    }
    
    int ijet = fiducialjets.at(i);
    
    double JetB1Mass = (Jet_TLV.at(ijet)+Jet_TLV.at(bjets_index_sorted.at(0))).M();
    double JetB2Mass = (Jet_TLV.at(ijet)+Jet_TLV.at(bjets_index_sorted.at(1))).M();

    double MaxJetBMass = TMath::Max(JetB1Mass,JetB2Mass);
    double MinJetBMass = TMath::Min(JetB1Mass,JetB2Mass);
    
    double pthrust  = CalculatePAlongThrust( Jet_TLV.at(ijet).Px(), Jet_TLV.at(ijet).Py());
    double ptthrust = CalculatePtAlongThrust(Jet_TLV.at(ijet).Px(), Jet_TLV.at(ijet).Py());
    
    Ana_AllHad_2b_Thrust_MaxJetBMass[icut].Fill(MaxJetBMass, weight );
    Ana_AllHad_2b_Thrust_MinJetBMass[icut].Fill(MinJetBMass, weight );
    
    Ana_AllHad_2b_Thrust_Jet_PThrust[icut].Fill(pthrust, weight );
    Ana_AllHad_2b_Thrust_Jet_PtThrust[icut].Fill(ptthrust,weight );
    Ana_AllHad_2b_Thrust_Jet_PThrust2D[icut].Fill( ptthrust, pthrust, weight );
    
    if ( jet_isISR ) {
      Ana_AllHad_2b_Thrust_ISRJet_PThrust[icut].Fill(    pthrust, weight );
      Ana_AllHad_2b_Thrust_ISRJet_PtThrust[icut].Fill(ptthrust,weight );
      Ana_AllHad_2b_Thrust_ISRJet_PThrust2D[icut].Fill(  ptthrust, pthrust, weight );
    }
    
    for( uint j=0; j < bjets_index.size(); j++ ) {
      if ( fiducialjets.at(i) == bjets_index.at(j) ) {
	Ana_AllHad_2b_Thrust_bJet_PThrust[icut].Fill(pthrust, weight);
	Ana_AllHad_2b_Thrust_bJet_PtThrust[icut].Fill(ptthrust, weight);
	Ana_AllHad_2b_Thrust_bJet_PThrust2D[icut].Fill( ptthrust, pthrust, weight );
      }
    }
    
  }

  //---------------------------------------------------//
  
  Ana_AllHad_2b_Thrust_nISRJets[icut].Fill( potentialISR_index.size(),weight );
  Ana_AllHad_2b_Thrust_n_TTjets[icut].Fill(  n_TTjets,weight );
  Ana_AllHad_2b_Thrust_n_bTTjets[icut].Fill( n_b_TTjets,weight );
  
  if ( potentialISR_index.size() > 0 ) {
    
    double MET_dirRatio_PTISR = ( MET_TLV.Px()*ISRjets_vec.Px() + MET_TLV.Py()*ISRjets_vec.Py() )/(ISRjets_vec.Pt()*ISRjets_vec.Pt());
    
    double ISRMetDeltaPhi = fabs( MET_TLV.DeltaPhi(  ISRjets_vec ) );
    double ISRTTDeltaPhi  = fabs( TTjets_vec.DeltaPhi( ISRjets_vec ) );
    double TTMetDeltaPhi  = fabs( TTjets_vec.DeltaPhi( MET_TLV ) );
    
    double ISR_pthrust = CalculatePAlongThrust(ISRjets_vec.Px() , ISRjets_vec.Py());
    double MET_pthrust = CalculatePAlongThrust( MET_TLV.Px(), MET_TLV.Py() );
    double TT_pthrust  = CalculatePAlongThrust(TTjets_vec.Px(), TTjets_vec.Py());
    
    double TT_overISR  = TTjets_vec.Pt() / ISRjets_vec.Pt();
    double MET_overISR = MET_TLV.Pt() / ISRjets_vec.Pt();
    double MET_overTT  = MET_TLV.Pt() / TTjets_vec.Pt();
    
    double ISRMass = ISRjets_vec.M();
    double TTMass  = TTjets_vec.M();
    
    Ana_AllHad_2b_Thrust_METISR_PtDirRatio[icut].Fill( MET_dirRatio_PTISR, weight);
    
    Ana_AllHad_2b_Thrust_ISRMass[icut].Fill( ISRMass, weight);
    Ana_AllHad_2b_Thrust_TTMass[icut].Fill(  TTMass,  weight);
    
    Ana_AllHad_2b_Thrust_ISRMet_DPhi[icut].Fill( ISRMetDeltaPhi, weight);
    
    Ana_AllHad_2b_Thrust_ISR_Pt[icut].Fill(ISRjets_vec.Pt(), weight);
    Ana_AllHad_2b_Thrust_ISR_PThrust[icut].Fill(ISR_pthrust, weight);
    Ana_AllHad_2b_Thrust_METoverISR_pt[icut].Fill(MET_overISR, weight);
    Ana_AllHad_2b_Thrust_METoverISR_PThrust[icut].Fill(MET_pthrust/ISR_pthrust, weight);
    
    Ana_AllHad_2b_Thrust_ISRTT_DPhi[icut].Fill( ISRTTDeltaPhi, weight);
    
    Ana_AllHad_2b_Thrust_TT_Pt[icut].Fill(TTjets_vec.Pt(), weight);
    Ana_AllHad_2b_Thrust_TT_PThrust[icut].Fill(TT_pthrust, weight);
    Ana_AllHad_2b_Thrust_TToverISR_pt[icut].Fill(     TT_overISR, weight);
    Ana_AllHad_2b_Thrust_TToverISR_PThrust[icut].Fill(TT_pthrust/ISR_pthrust, weight);
    
    Ana_AllHad_2b_Thrust_TTMET_DPhi[icut].Fill( TTMetDeltaPhi, weight);
    
    Ana_AllHad_2b_Thrust_METoverTT_pt[icut].Fill(MET_overTT, weight);
    Ana_AllHad_2b_Thrust_METoverTT_PThrust[icut].Fill(MET_pthrust/TT_pthrust, weight);
    
    Ana_AllHad_2b_Thrust_ISR_MET_pt[icut].Fill(ISRjets_vec.Pt(), MET_TLV.Pt(), weight);
    Ana_AllHad_2b_Thrust_MET_TT_pt[icut].Fill( MET_TLV.Pt(),        TTjets_vec.Pt(), weight);
    Ana_AllHad_2b_Thrust_ISR_TT_pt[icut].Fill( ISRjets_vec.Pt(), TTjets_vec.Pt(), weight);
    
    Ana_AllHad_2b_Thrust_ISR_MET_pthrust[icut].Fill(ISR_pthrust, MET_pthrust, weight);
    Ana_AllHad_2b_Thrust_MET_TT_pthrust[icut].Fill(  MET_pthrust, TT_pthrust, weight);
    Ana_AllHad_2b_Thrust_ISR_TT_pthrust[icut].Fill(  ISR_pthrust, TT_pthrust, weight);
    
    if ( potentialISR_index.size() > 1 ) {
      
      Ana_AllHad_2b_Thrust_multISR_ISR_Pt[icut].Fill(ISRjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_multISR_ISR_PThrust[icut].Fill(ISR_pthrust, weight);
      Ana_AllHad_2b_Thrust_multISR_TT_Pt[icut].Fill(TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_multISR_TT_PThrust[icut].Fill(TT_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_multISR_METoverISR_pt[icut].Fill(     MET_overISR, weight);
      Ana_AllHad_2b_Thrust_multISR_METoverISR_PThrust[icut].Fill(MET_pthrust/ISR_pthrust, weight);
      Ana_AllHad_2b_Thrust_multISR_TToverISR_pt[icut].Fill(      TT_overISR, weight);
      Ana_AllHad_2b_Thrust_multISR_TToverISR_PThrust[icut].Fill( TT_pthrust/ISR_pthrust, weight);
      Ana_AllHad_2b_Thrust_multISR_METoverTT_pt[icut].Fill(      MET_overTT, weight);
      Ana_AllHad_2b_Thrust_multISR_METoverTT_PThrust[icut].Fill( MET_pthrust/TT_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_multISR_ISRMet_DPhi[icut].Fill( ISRMetDeltaPhi, weight);
      Ana_AllHad_2b_Thrust_multISR_ISRTT_DPhi[icut].Fill(  ISRTTDeltaPhi, weight);
      Ana_AllHad_2b_Thrust_multISR_TTMET_DPhi[icut].Fill( TTMetDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_multISR_ISR_MET_pt[icut].Fill(ISRjets_vec.Pt(), MET_TLV.Pt(), weight);
      Ana_AllHad_2b_Thrust_multISR_MET_TT_pt[icut].Fill( MET_TLV.Pt(),        TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_multISR_ISR_TT_pt[icut].Fill( ISRjets_vec.Pt(), TTjets_vec.Pt(), weight);
      
      Ana_AllHad_2b_Thrust_multISR_ISR_MET_pthrust[icut].Fill(ISR_pthrust, MET_pthrust, weight);
      Ana_AllHad_2b_Thrust_multISR_MET_TT_pthrust[icut].Fill(  MET_pthrust, TT_pthrust, weight);
      Ana_AllHad_2b_Thrust_multISR_ISR_TT_pthrust[icut].Fill(  ISR_pthrust, TT_pthrust, weight);
    }
  }

  //---------------------------------------------------//
    
  double HT = CalcHT( 20. );
  double CT = TMath::Min((double) (MET_TLV.Pt()), (double) HT - 100.);
  
  Ana_AllHad_2b_Thrust_MET[icut].Fill(MET_TLV.Pt(),weight);
  Ana_AllHad_2b_Thrust_HT[icut].Fill(HT,weight);
  Ana_AllHad_2b_Thrust_CT[icut].Fill(CT,weight);
  Ana_AllHad_2b_Thrust_METsig[icut].Fill(MET_TLV.Pt()/sqrt(HT), weight);
  Ana_AllHad_2b_Thrust_MET_over_ISR[icut].Fill(MET_TLV.Pt()/ISRjets_vec.Pt(), weight);
  
  Ana_AllHad_2b_Thrust_nJets[icut].Fill(fiducialjets.size(), weight);
  Ana_AllHad_2b_Thrust_nHighPtJets[icut].Fill(highptjet.size(), weight);
  
  if (verbose) std::cout << "plotting jets" << std::endl;
  
  //---------------------------------------------------//
  
  if( fiducialjets.size() > 0 ) Ana_AllHad_2b_Thrust_Jet1Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(0)).Pt(), weight);
  if( fiducialjets.size() > 0 ) Ana_AllHad_2b_Thrust_Jet1Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(0)).Eta(), weight);
  if( fiducialjets.size() > 1 ) Ana_AllHad_2b_Thrust_Jet2Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(1)).Pt(), weight);
  if( fiducialjets.size() > 1 ) Ana_AllHad_2b_Thrust_Jet2Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(1)).Eta(), weight);
  if( fiducialjets.size() > 2 ) Ana_AllHad_2b_Thrust_Jet3Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(2)).Pt(), weight);
  if( fiducialjets.size() > 2 ) Ana_AllHad_2b_Thrust_Jet3Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(2)).Eta(), weight);
  if( fiducialjets.size() > 3 ) Ana_AllHad_2b_Thrust_Jet4Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(3)).Pt(), weight);
  if( fiducialjets.size() > 3 ) Ana_AllHad_2b_Thrust_Jet4Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(3)).Eta(), weight);
  if( fiducialjets.size() > 4 ) Ana_AllHad_2b_Thrust_Jet5Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(4)).Pt(), weight);
  if( fiducialjets.size() > 4 ) Ana_AllHad_2b_Thrust_Jet5Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(4)).Eta(), weight);
  if( fiducialjets.size() > 5 ) Ana_AllHad_2b_Thrust_Jet6Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(5)).Pt(), weight);
  if( fiducialjets.size() > 5 ) Ana_AllHad_2b_Thrust_Jet6Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(5)).Eta(), weight);
  if( fiducialjets.size() > 6 ) Ana_AllHad_2b_Thrust_Jet7Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(6)).Pt(), weight);
  if( fiducialjets.size() > 6 ) Ana_AllHad_2b_Thrust_Jet7Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(6)).Eta(), weight);
  if( fiducialjets.size() > 7 ) Ana_AllHad_2b_Thrust_Jet8Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(7)).Pt(), weight);
  if( fiducialjets.size() > 7 ) Ana_AllHad_2b_Thrust_Jet8Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(7)).Eta(), weight);
  if( fiducialjets.size() > 8 ) Ana_AllHad_2b_Thrust_Jet9Pt[icut].Fill(Jet_TLV.at(fiducialjets.at(8)).Pt(), weight);
  if( fiducialjets.size() > 8 ) Ana_AllHad_2b_Thrust_Jet9Eta[icut].Fill(Jet_TLV.at(fiducialjets.at(8)).Eta(), weight);
  
  //----------------------------------------------------//

  Ana_AllHad_2b_Thrust_bJet1Pt[icut].Fill(Jet_TLV.at(bjets_index.at(0)).Pt(), weight);
  Ana_AllHad_2b_Thrust_bJet1Eta[icut].Fill(Jet_TLV.at(bjets_index.at(0)).Eta(), weight);
  Ana_AllHad_2b_Thrust_bJet1_ijet[icut].Fill(bjets_index.at(0), weight);

  Ana_AllHad_2b_Thrust_bJet2Pt[icut].Fill(Jet_TLV.at(bjets_index.at(1)).Pt(), weight);
  Ana_AllHad_2b_Thrust_bJet2Eta[icut].Fill(Jet_TLV.at(bjets_index.at(1)).Eta(), weight);
  Ana_AllHad_2b_Thrust_bJet2_ijet[icut].Fill(bjets_index.at(1), weight);

  if (bjets_index.size()>2) {
    Ana_AllHad_2b_Thrust_bJet3Pt[icut].Fill(Jet_TLV.at(bjets_index.at(2)).Pt(), weight);
    Ana_AllHad_2b_Thrust_bJet3Eta[icut].Fill(Jet_TLV.at(bjets_index.at(2)).Eta(), weight);
    Ana_AllHad_2b_Thrust_bJet3_ijet[icut].Fill(bjets_index.at(2), weight);
  }

  //-------------------------------------------------//
  //         Look for W near the min pthrust b       //
  //-------------------------------------------------//

  if (verbose) std::cout << "plotting W jets" << std::endl;

  double top_hadronic_pthrust = CalculatePAlongThrust(TTbar_Had_Top_vec.Px(), TTbar_Had_Top_vec.Py());
  double top_delta_pthrust   = top_hadronic_pthrust + CalculatePAlongThrust(MET_TLV.Px(), MET_TLV.Py());

  //----------------------------------------------------//
  //       Look for a full W in the negative hemi       //
  //----------------------------------------------------//

  double top_met_deltaphi    = fabs(TTbar_Had_Top_vec.DeltaPhi(MET_TLV));

  if ( Found_Neg_PThrust_Had_Top ) {

    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_FoundTop[icut].Fill( 1.0, weight );
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_Mass[icut].Fill(               TTbar_Had_Top_vec.M(), weight);
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust[icut].Fill(            top_hadronic_pthrust, weight);
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_met_deltaphi[icut].Fill(   top_met_deltaphi, weight );
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_DeltaR[icut].Fill(             Jet_TLV.at(bjets_index_sorted.back()).DeltaR(whad_vec), 
										weight);
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust_mass[icut].Fill(       top_hadronic_pthrust, TTbar_Had_Top_vec.M(), weight);
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_delta_pthrust_mass[icut].Fill( top_delta_pthrust,    TTbar_Had_Top_vec.M(), weight);

  }
  else {
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_FoundTop[icut].Fill( 0.0, weight );
  }

  //--------------------------------------------------//
  //            Plot Truth information                //
  //--------------------------------------------------//

  bool top1_ishad = false;
  if ( abs(pdgId1_1) < 6 ) {
    top1_ishad = true;
  }
  bool top1_istau = false;
  if ( abs(pdgId1_1) == 15 ) {
    top1_istau = true;
  }

  bool top2_ishad = false;
  if ( abs(pdgId2_1) < 6 ) {
    top2_ishad = true;
  }
  bool top2_istau = false;
  if ( abs(pdgId2_1) == 15 ) {
    top2_istau = true;
  }

  if ( top1_ishad && top2_ishad ) {  // all had
    Ana_AllHad_2b_Thrust_TruthTopID[icut].Fill( 0.0, weight );
  }
  else if ( (top1_ishad || top2_ishad) && ( top1_istau || top2_istau ) ) { // single tau
    Ana_AllHad_2b_Thrust_TruthTopID[icut].Fill( 1.0, weight );
  }
  else if ( (top1_ishad || top2_ishad) && ( !top1_istau && !top2_istau ) ) { // single lepton
    Ana_AllHad_2b_Thrust_TruthTopID[icut].Fill( 2.0, weight );
  }
  else { // dilepton/ditau equivalent to !top1_ishad && !top2_ishad
    Ana_AllHad_2b_Thrust_TruthTopID[icut].Fill( 3.0, weight );
  }

  if ( verbose ) {
    if ( top1_istau || top2_istau ) {
      std::cout << "top->tau event" << std::endl;
    }
    else if ( top1_ishad && top2_ishad ) {
      std::cout << "all hadronic evt" << std::endl;
    }
    else {
      std::cout << "leptonic evt" << std::endl;
    }
  }

  //----------------------------------------------------------//

  double truth_met_px = stop1_xi_px+stop2_xi_px;
  double truth_met_py = stop1_xi_py+stop2_xi_py;

  if ( abs(pdgId1_1) == 12 || abs(pdgId1_1) == 14 || abs(pdgId1_1) == 16 ) {
    truth_met_px += stop1_w1_px;
    truth_met_py += stop1_w1_py;
  }
  if ( abs(pdgId1_2) == 12 || abs(pdgId1_2) == 14 || abs(pdgId1_2) == 16 ) {
    truth_met_px += stop1_w2_px;
    truth_met_py += stop1_w2_py;
  }
  if ( abs(pdgId2_1) == 12 || abs(pdgId2_1) == 14 || abs(pdgId2_1) == 16 ) {
    truth_met_px += stop2_w1_px;
    truth_met_py += stop2_w1_py;
  }
  if ( abs(pdgId2_2) == 12 || abs(pdgId2_2) == 14 || abs(pdgId2_2) == 16 ) {
    truth_met_px += stop2_w2_px;
    truth_met_py += stop2_w2_py;
  }

  double truth_met_pt = sqrt( truth_met_px*truth_met_px +
                              truth_met_py*truth_met_py );

  truth_tot_ISR_TLV.SetPxPyPzE(ISR1_px+ISR2_px, ISR1_py+ISR2_py,
			       ISR1_pz+ISR2_pz, ISR1_e +ISR2_e);

  //----------------------------------------------------//

  if ( top1_ishad && top2_ishad ) {

    if ( potentialISR_index.size() > 0 ) {

      double MET_dirRatio_PTISR = ( MET_TLV.Px()*ISRjets_vec.Px() + MET_TLV.Py()*ISRjets_vec.Py() )/(ISRjets_vec.Pt()*ISRjets_vec.Pt());
      
      double ISRMetDeltaPhi = fabs( MET_TLV.DeltaPhi(  ISRjets_vec ) );
      double ISRTTDeltaPhi  = fabs( TTjets_vec.DeltaPhi( ISRjets_vec ) );
      double TTMetDeltaPhi  = fabs( TTjets_vec.DeltaPhi( MET_TLV ) );
      
      double ISR_pthrust = CalculatePAlongThrust(ISRjets_vec.Px() , ISRjets_vec.Py());
      double MET_pthrust = CalculatePAlongThrust( MET_TLV.Px(), MET_TLV.Py() );
      double TT_pthrust  = CalculatePAlongThrust(TTjets_vec.Px(), TTjets_vec.Py());
      
      double TT_overISR  = TTjets_vec.Pt() / ISRjets_vec.Pt();
      double MET_overISR = MET_TLV.Pt() / ISRjets_vec.Pt();
      double MET_overTT  = MET_TLV.Pt() / TTjets_vec.Pt();
      
      double ISRMass = ISRjets_vec.M();
      double TTMass  = TTjets_vec.M();
      
      Ana_AllHad_2b_Thrust_TruthAllHad_METISR_PtDirRatio[icut].Fill( MET_dirRatio_PTISR, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_ISRMass[icut].Fill( ISRMass, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_TTMass[icut].Fill(  TTMass,  weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_ISRMet_DPhi[icut].Fill( ISRMetDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_Pt[icut].Fill(ISRjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_PThrust[icut].Fill(ISR_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_pt[icut].Fill(MET_overISR, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_PThrust[icut].Fill(MET_pthrust/ISR_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_ISRTT_DPhi[icut].Fill( ISRTTDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_TT_Pt[icut].Fill(TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_TT_PThrust[icut].Fill(TT_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_pt[icut].Fill(     TT_overISR, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_PThrust[icut].Fill(TT_pthrust/ISR_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_TTMET_DPhi[icut].Fill( TTMetDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_pt[icut].Fill(MET_overTT, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_PThrust[icut].Fill(MET_pthrust/TT_pthrust, weight);

      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pt[icut].Fill(ISRjets_vec.Pt(), MET_TLV.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pt[icut].Fill( MET_TLV.Pt(),        TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pt[icut].Fill( ISRjets_vec.Pt(), TTjets_vec.Pt(), weight);

      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pthrust[icut].Fill(ISR_pthrust, MET_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pthrust[icut].Fill(  MET_pthrust, TT_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pthrust[icut].Fill(  ISR_pthrust, TT_pthrust, weight);

      Ana_AllHad_2b_Thrust_TruthAllHad_TruthISR_RecoISR_pt[icut].Fill(truth_tot_ISR_TLV.Pt(), ISRjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthAllHad_TruthMET_RecoMET_pt[icut].Fill(truth_met_pt, MET_TLV.Pt(), weight );
    }
  }
  else {
    if ( potentialISR_index.size() > 0 ) {
      
      double MET_dirRatio_PTISR = ( MET_TLV.Px()*ISRjets_vec.Px() + MET_TLV.Py()*ISRjets_vec.Py() )/(ISRjets_vec.Pt()*ISRjets_vec.Pt());
      
      double ISRMetDeltaPhi = fabs( MET_TLV.DeltaPhi(  ISRjets_vec ) );
      double ISRTTDeltaPhi  = fabs( TTjets_vec.DeltaPhi( ISRjets_vec ) );
      double TTMetDeltaPhi  = fabs( TTjets_vec.DeltaPhi( MET_TLV ) );
      
      double ISR_pthrust = CalculatePAlongThrust(ISRjets_vec.Px() , ISRjets_vec.Py());
      double MET_pthrust = CalculatePAlongThrust( MET_TLV.Px(), MET_TLV.Py() );
      double TT_pthrust  = CalculatePAlongThrust(TTjets_vec.Px(), TTjets_vec.Py());
      
      double TT_overISR  = TTjets_vec.Pt() / ISRjets_vec.Pt();
      double MET_overISR = MET_TLV.Pt() / ISRjets_vec.Pt();
      double MET_overTT  = MET_TLV.Pt() / TTjets_vec.Pt();
      
      double ISRMass = ISRjets_vec.M();
      double TTMass  = TTjets_vec.M();
      
      Ana_AllHad_2b_Thrust_TruthNotHad_METISR_PtDirRatio[icut].Fill( MET_dirRatio_PTISR, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISRMass[icut].Fill( ISRMass, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_TTMass[icut].Fill(  TTMass,  weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISRMet_DPhi[icut].Fill( ISRMetDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_Pt[icut].Fill(ISRjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_PThrust[icut].Fill(ISR_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_pt[icut].Fill(MET_overISR, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_PThrust[icut].Fill(MET_pthrust/ISR_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISRTT_DPhi[icut].Fill( ISRTTDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_TT_Pt[icut].Fill(TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_TT_PThrust[icut].Fill(TT_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_pt[icut].Fill(     TT_overISR, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_PThrust[icut].Fill(TT_pthrust/ISR_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_TTMET_DPhi[icut].Fill( TTMetDeltaPhi, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_pt[icut].Fill(MET_overTT, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_PThrust[icut].Fill(MET_pthrust/TT_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pt[icut].Fill(ISRjets_vec.Pt(), MET_TLV.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pt[icut].Fill( MET_TLV.Pt(),        TTjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pt[icut].Fill( ISRjets_vec.Pt(), TTjets_vec.Pt(), weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pthrust[icut].Fill(ISR_pthrust, MET_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pthrust[icut].Fill(  MET_pthrust, TT_pthrust, weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pthrust[icut].Fill(  ISR_pthrust, TT_pthrust, weight);
      
      Ana_AllHad_2b_Thrust_TruthNotHad_TruthISR_RecoISR_pt[icut].Fill(truth_tot_ISR_TLV.Pt(), ISRjets_vec.Pt(), weight);
      Ana_AllHad_2b_Thrust_TruthNotHad_TruthMET_RecoMET_pt[icut].Fill(truth_met_pt, MET_TLV.Pt(), weight );
    }
  }

  return;
}

void CollectionTree::SetVerbose(bool v) {
  verbose = v;
  return;
}

void CollectionTree::SetOutputName(std::string name) {
  output_name = name;
  return;
}

void CollectionTree::SetEvtDisplayName(std::string name)
{
  evt_display_name = name;
  return;
}

void CollectionTree::DrawEvtDisplay(double weight){

  TCanvas *c2 = new TCanvas("c2");
  c2->Range(0,0,1,1);

  TArrow *metArrow = new TArrow();
  TArrow *arrow = new TArrow();

  TLatex p;
  p.SetTextSize(0.05);
  p.SetTextFont(42);

  metArrow->SetLineColor(kBlack);

  double maximum = 2.2*TMath::Max((double) MET_TLV.Pt(), (double) Jet_TLV.at(0).Pt());

  double length  = MET_TLV.Pt()/maximum;
  metArrow->DrawArrow(0.5, 0.5, 0.5-length*cos(MET_TLV.Phi()), 0.5+length*sin(MET_TLV.Phi()),0.01,"|>");

  double place = 0.6;

  for ( uint i=0; i< fiducialjets.size(); i++ ) {
    
    int ijet = fiducialjets.at(i);

    length = Jet_TLV.at(ijet).Pt()/maximum;

    bool isB = false;
    for ( uint j=0; j < bjets_index.size(); j++ ) {
      if ( fiducialjets.at(i) == bjets_index.at(j) ) {
        isB = true;
      }
    }
    
    bool isISR = false;
    for ( uint j=0; j < potentialISR_index.size(); j++ ) {
      if ( i == potentialISR_index.at(j) ) {
        isISR = true;
      }
    }


    if ( isB ) {
      arrow->SetLineColor(kBlue);
      if ( isISR ) {
        arrow->SetLineColor(kViolet);
      }
    }
    else if ( isISR ) {
      arrow->SetLineColor(kRed);
    }
    else {
      arrow->SetLineColor(kOrange);
    }


    if ( i < 3 ) {
      arrow->SetFillColor(kGreen);
    }
    else {
      arrow->SetFillColor(kYellow);
    }

    arrow->DrawArrow(0.5, 0.5, 0.5-length*cos(Jet_TLV.at(ijet).Phi()), 0.5+length*sin(Jet_TLV.at(ijet).Phi()),0.01,"|>");

    if ( isB ) {
      p.DrawLatex(0.1,place,Form("#color[4]{bjet pt = %3.2f}",Jet_TLV.at(ijet).Pt()));
      place -= 0.1;
    }
  }

  p.DrawLatex(0.1,0.9,"Reconstructed");
  p.DrawLatex(0.1,0.8,Form("#color[1]{MET     = %3.2f}",MET_TLV.Pt()));
  p.DrawLatex(0.1,0.7,Form("#color[2]{ISR pt  = %3.2f}",ISRjets_vec.Pt()) );

  c2->Update();  ps->NewPage();
  c2->Clear();

  //----------------------------------------------------//
  //                  draw truth info                   //
  //----------------------------------------------------//

  if ( isMC ) {
    
    arrow->SetFillColor(kWhite);
    arrow->SetLineColor(kBlue);
    arrow->DrawArrow(0.5, 0.5, 0.5-stop1_b_px/maximum, 0.5+stop1_b_py/maximum,0.01,"|>");
    arrow->DrawArrow(0.5, 0.5, 0.5-stop2_b_px/maximum, 0.5+stop2_b_py/maximum,0.01,"|>");

    if ( abs(pdgId1_1) < 6 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kYellow);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w1_px/maximum, 0.5+stop1_w1_py/maximum,0.01,"|>");
    }
    else if ( abs(pdgId1_1) == 11 || abs(pdgId1_1) == 13 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kGreen);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w1_px/maximum, 0.5+stop1_w1_py/maximum,0.01,"|>");
    }
    else if ( abs(pdgId1_1) == 15 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kMagenta);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w1_px/maximum, 0.5+stop1_w1_py/maximum,0.01,"|>");
    }

    if ( abs(pdgId1_2) < 6 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kYellow);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w2_px/maximum, 0.5+stop1_w2_py/maximum,0.01,"|>");
    }   
    else if ( abs(pdgId1_2) == 11 || abs(pdgId1_2) == 13 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kGreen);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w2_px/maximum, 0.5+stop1_w2_py/maximum,0.01,"|>");
    }
    else if ( abs(pdgId1_2) == 15 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kMagenta);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop1_w2_px/maximum, 0.5+stop1_w2_py/maximum,0.01,"|>");
    }

    if ( abs(pdgId2_1) < 6 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kYellow);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w1_px/maximum, 0.5+stop2_w1_py/maximum,0.01,"|>");
    }   
    else if ( abs(pdgId2_1) == 11 || abs(pdgId2_1) == 13 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kGreen);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w1_px/maximum, 0.5+stop2_w1_py/maximum,0.01,"|>");
    }
    else if ( abs(pdgId2_1) == 15 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kMagenta);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w1_px/maximum, 0.5+stop2_w1_py/maximum,0.01,"|>");
    }
      
    if ( abs(pdgId2_2) < 6 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kYellow);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w2_px/maximum, 0.5+stop2_w2_py/maximum,0.01,"|>");
    }   
    else if ( abs(pdgId2_2) == 11 || abs(pdgId2_2) == 13 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kGreen);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w2_px/maximum, 0.5+stop2_w2_py/maximum,0.01,"|>");
    }
    else if ( abs(pdgId2_2) == 15 ) {
      arrow->SetFillColor(kWhite);
      arrow->SetLineColor(kMagenta);
      arrow->DrawArrow(0.5, 0.5, 0.5-stop2_w2_px/maximum, 0.5+stop2_w2_py/maximum,0.01,"|>");
    }

    arrow->SetFillColor(kWhite);
    arrow->SetLineColor(kRed);
    arrow->DrawArrow(0.5, 0.5, 0.5-ISR1_px/maximum, 0.5+ISR1_py/maximum,0.01,"|>");
    arrow->DrawArrow(0.5, 0.5, 0.5-ISR2_px/maximum, 0.5+ISR2_py/maximum,0.01,"|>");

    double truth_met_px = stop1_xi_px+stop2_xi_px;
    double truth_met_py = stop1_xi_py+stop2_xi_py;

    if ( abs(pdgId1_1) == 12 || abs(pdgId1_1) == 14 || abs(pdgId1_1) == 16 ) {
      truth_met_px += stop1_w1_px;
      truth_met_py += stop1_w1_py;
    }
    if ( abs(pdgId1_2) == 12 || abs(pdgId1_2) == 14 || abs(pdgId1_2) == 16 ) {
      truth_met_px += stop1_w2_px;
      truth_met_py += stop1_w2_py;
    }
    if ( abs(pdgId2_1) == 12 || abs(pdgId2_1) == 14 || abs(pdgId2_1) == 16 ) {
      truth_met_px += stop2_w1_px;
      truth_met_py += stop2_w1_py;
    }
    if ( abs(pdgId2_2) == 12 || abs(pdgId2_2) == 14 || abs(pdgId2_2) == 16 ) {
      truth_met_px += stop2_w2_px;
      truth_met_py += stop2_w2_py;
    }
    
    arrow->SetFillColor(kWhite);
    arrow->SetLineColor(kBlack);
    arrow->DrawArrow(0.5, 0.5, 0.5-truth_met_px/maximum, 0.5+truth_met_py/maximum,0.01,"|>");

    p.DrawLatex(0.1,0.9,"Truth");
    p.DrawLatex(0.1,0.8,Form("#color[1]{MET     = %3.2f}",sqrt(truth_met_px*truth_met_px+
							       truth_met_py*truth_met_py)) );
    p.DrawLatex(0.1,0.7,Form("#color[2]{ISR pt  = %3.2f}",sqrt((ISR1_px+ISR2_px)*
							       (ISR1_px+ISR2_px)+
							       (ISR1_py+ISR2_py)*
							       (ISR1_py+ISR2_py))) );
    //  p.DrawLatex(0.1,0.6,Form("#color[4]{bjet pt = %3.2f}",jets.at(bjet1)->Pt()));                                                                               
    //  p.DrawLatex(0.1,0.5,Form("#color[4]{bjet pt = %3.2f}",jets.at(bjet2)->Pt()));                                                                               
    
    c2->Update();  ps->NewPage();

  }

  delete c2;
  delete arrow;
  delete metArrow;
}

void CollectionTree::BookHistos()
{

  if( verbose ) std::cout << "booking histos" << std::endl;

  ps= new TPostScript(evt_display_name.c_str(), 112);

  //--------------------------------------------------------//

  All_nlep       = new TH1F("All_nlep", "All n leptons", 10, -0.5, 9.5);
  All_nHighPtLep = new TH1F("All_nHighPtLep", "All n high pt leptons", 10, -0.5, 9.5);

  All_lep1Pt  = new TH1F("All_lep1Pt",  "All lead lepton Pt",  200, 0.0, 200. );
  All_lep1Eta = new TH1F("All_lep1Eta", "All lead lepton Eta", 200, -5.0, 5.0 );
  All_lep1Iso = new TH1F("All_lep1Iso", "All lead lepton Iso", 200, 0.0, 1.0 );

  All_lep2Pt  = new TH1F("All_lep2Pt",  "All sublead lepton Pt",  200, 0.0, 100. );
  All_lep2Eta = new TH1F("All_lep2Eta", "All sublead lepton Eta", 200, -5.0, 5.0 );
  All_lep2Iso = new TH1F("All_lep2Iso", "All sublead lepton Iso", 200, 0.0, 1.0 );

  All_MET    = new TH1F("All_MET", "All MET", 200, 0.0, 800.);
  All_METPhi = new TH1F("All_METPhi", "All METPhi", 200, -TMath::Pi(), TMath::Pi());

  All_HT  = new TH1F("All_HT",  "All HT",  200, 0.0, 1300.);
  All_CT  = new TH1F("All_CT",  "All CT",  200, 0.0, 500.);

  All_nJets = new TH1F("All_nJets", "All n jets", 21, -0.5, 20.5);

  All_Jet1Pt  = new TH1F("All_Jet1Pt", "All Leading Jet Pt", 200, 0.0, 1200.);
  All_Jet1Eta = new TH1F("All_Jet1Eta", "All Leading Jet Eta", 200, -5.0, 5.0 );

  All_Jet2Pt  = new TH1F("All_Jet2Pt", "All subLeading Jet Pt", 200, 0.0, 400.);
  All_Jet2Eta = new TH1F("All_Jet2Eta", "All subLeading Jet Eta", 200,-5.0, 5.0 );

  All_DiHiPtJet_DeltaPhi = new TH1F("All_DiJet_DeltaPhi", "All DiJet DeltaPhi", 200, 0.0, TMath::Pi());

  All_Mt = new TH1F("All_Mt", "All Mt", 200, 0.0, 100.);

  //--------------------------------------------------------//


  AfterTrigger_nlep       = new TH1F("AfterTrigger_nlep", "AfterTrigger n leptons", 10, -0.5, 9.5);
  AfterTrigger_nHighPtLep = new TH1F("AfterTrigger_nHighPtLep", "AfterTrigger n high pt leptons", 10, -0.5, 9.5);

  AfterTrigger_lep1Pt  = new TH1F("AfterTrigger_lep1Pt",  "AfterTrigger lead lepton Pt",  200, 0.0, 200. );
  AfterTrigger_lep1Eta = new TH1F("AfterTrigger_lep1Eta", "AfterTrigger lead lepton Eta", 200, -5.0, 5.0 );
  AfterTrigger_lep1Iso = new TH1F("AfterTrigger_lep1Iso", "AfterTrigger lead lepton Iso", 200, 0.0, 1.0 );

  AfterTrigger_lep2Pt  = new TH1F("AfterTrigger_lep2Pt",  "AfterTrigger sublead lepton Pt",  200, 0.0, 100. );
  AfterTrigger_lep2Eta = new TH1F("AfterTrigger_lep2Eta", "AfterTrigger sublead lepton Eta", 200, -5.0, 5.0 );
  AfterTrigger_lep2Iso = new TH1F("AfterTrigger_lep2Iso", "AfterTrigger sublead lepton Iso", 200, 0.0, 1.0 );

  AfterTrigger_MET    = new TH1F("AfterTrigger_MET", "AfterTrigger MET", 200, 0.0, 800.);
  AfterTrigger_METPhi = new TH1F("AfterTrigger_METPhi", "AfterTrigger METPhi", 200, -TMath::Pi(), TMath::Pi());

  AfterTrigger_HT  = new TH1F("AfterTrigger_HT",  "AfterTrigger HT",  200, 0.0, 1400.);
  AfterTrigger_CT  = new TH1F("AfterTrigger_CT",  "AfterTrigger CT",  200, 0.0, 1400.);

  AfterTrigger_nJets = new TH1F("AfterTrigger_nJets", "AfterTrigger n jets", 21, -0.5, 20.5);

  AfterTrigger_Jet1Pt  = new TH1F("AfterTrigger_Jet1Pt", "AfterTrigger Leading Jet Pt", 200, 0.0, 1300.);
  AfterTrigger_Jet1Eta = new TH1F("AfterTrigger_Jet1Eta", "AfterTrigger Leading Jet Eta", 200, -5.0, 5.0 );

  AfterTrigger_Jet2Pt  = new TH1F("AfterTrigger_Jet2Pt", "AfterTrigger subLeading Jet Pt", 200, 0.0, 400.);
  AfterTrigger_Jet2Eta = new TH1F("AfterTrigger_Jet2Eta", "AfterTrigger subLeading Jet Eta", 200,-5.0, 5.0 );

  AfterTrigger_DiHiPtJet_DeltaPhi = new TH1F("AfterTrigger_DiJet_DeltaPhi", "AfterTrigger DiJet DeltaPhi", 200, 0.0, TMath::Pi());

  AfterTrigger_Mt = new TH1F("AfterTrigger_Mt", "AfterTrigger Mt", 200, 0.0, 1000.);

  //---------------------------------------------------------//

  PreSel_nlep       = new TH1F("PreSel_nlep", "PreSel n leptons", 10, -0.5, 9.5);
  PreSel_nHighPtLep = new TH1F("PreSel_nHighPtLep", "PreSel n high pt leptons", 10, -0.5, 9.5);

  PreSel_Mu1Pt  = new TH1F("PreSel_Mu1Pt",  "PreSel lead Muon Pt",  200, 0.0, 100. );
  PreSel_Mu1Eta = new TH1F("PreSel_Mu1Eta", "PreSel lead Muon Eta", 200, -5.0, 5.0 );
  PreSel_Mu1Iso = new TH1F("PreSel_Mu1Iso", "PreSel lead Muon Iso", 200, 0.0, 1.0 );

  PreSel_Mu2Pt  = new TH1F("PreSel_Mu2Pt",  "PreSel sublead muon Pt",  200, 0.0, 100. );
  PreSel_Mu2Eta = new TH1F("PreSel_Mu2Eta", "PreSel sublead muon Eta", 200, -5.0, 5.0 );
  PreSel_Mu2Iso = new TH1F("PreSel_Mu2Iso", "PreSel sublead muon Iso", 200, 0.0, 1.0 );

  PreSel_Ele1Pt  = new TH1F("PreSel_Ele1Pt",  "PreSel lead Ele Pt",  200, 0.0, 100. );
  PreSel_Ele1Eta = new TH1F("PreSel_Ele1Eta", "PreSel lead Ele Eta", 200, -5.0, 5.0 );
  PreSel_Ele1Iso = new TH1F("PreSel_Ele1Iso", "PreSel lead Ele Iso", 200, 0.0, 1.0 );

  PreSel_Ele2Pt  = new TH1F("PreSel_Ele2Pt",  "PreSel sublead Ele Pt",  200, 0.0, 100. );
  PreSel_Ele2Eta = new TH1F("PreSel_Ele2Eta", "PreSel sublead Ele Eta", 200, -5.0, 5.0 );
  PreSel_Ele2Iso = new TH1F("PreSel_Ele2Iso", "PreSel sublead Ele Iso", 200, 0.0, 1.0 );

  PreSel_MET    = new TH1F("PreSel_MET", "PreSel MET", 200, 0.0, 1000.);
  PreSel_METPhi = new TH1F("PreSel_METPhi", "PreSel METPhi", 200, -TMath::Pi(), TMath::Pi());

  PreSel_HT  = new TH1F("PreSel_HT",  "PreSel HT",  200, 0.0, 1400.);
  PreSel_CT  = new TH1F("PreSel_CT",  "PreSel CT",  200, 0.0, 1400.);

  PreSel_nJets = new TH1F("PreSel_nJets", "PreSel n jets", 21, -0.5, 20.5);

  PreSel_JetISR_Pt  = new TH1F("PreSel_JetISR_Pt", "PreSel ISR Jet Pt", 200, 0.0, 1300.);
  PreSel_JetISR_Eta = new TH1F("PreSel_JetISR_Eta", "PreSel ISR Jet Eta", 200, -5.0, 5.0 );

  PreSel_Jet1Pt  = new TH1F("PreSel_Jet1Pt", "PreSel Leading Jet Pt", 200, 0.0, 1300.);
  PreSel_Jet1Eta = new TH1F("PreSel_Jet1Eta", "PreSel Leading Jet Eta", 200, -5.0, 5.0 );

  PreSel_Jet2Pt  = new TH1F("PreSel_Jet2Pt", "PreSel subLeading Jet Pt", 200, 0.0, 400.);
  PreSel_Jet2Eta = new TH1F("PreSel_Jet2Eta", "PreSel subLeading Jet Eta", 200,-5.0, 5.0 );

  PreSel_DiHiPtJet_DeltaPhi = new TH1F("PreSel_DiJet_DeltaPhi", "PreSel DiJet DeltaPhi", 200, 0.0, TMath::Pi());

  PreSel_Mt = new TH1F("PreSel_Mt", "PreSel Mt", 200, 0.0, 1000.);

  //-----------------------------------------------------------------------------// 

  TString name;

  for ( int i=0; i<10; i++ ) {

    stringstream ssi;

    ssi << i;
    string istr = ssi.str();

    name = "Ana_AllHad_2b_Thrust_Thrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Thrust.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1.0) );
    name = "Ana_AllHad_2b_Thrust_Thrust_Met_dPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Thrust_Met_dPhi.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, TMath::Pi()) );

    name = "Ana_AllHad_2b_Thrust_Jet_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -700.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_Jet_PtThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet_PtThrust.push_back( TH1F(name.Data(),name.Data(), 400, -500.0, 500.0) );
    name = "Ana_AllHad_2b_Thrust_Jet_PThrust2D_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet_PThrust2D.push_back( TH2F(name.Data(),name.Data(), 40, -500.,500., 40, -700.0, 600.0) );
    
    name = "Ana_AllHad_2b_Thrust_ISRJet_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRJet_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -700.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_ISRJet_PtThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRJet_PtThrust.push_back( TH1F(name.Data(),name.Data(), 400, -500.0, 500.0) );
    name = "Ana_AllHad_2b_Thrust_ISRJet_PThrust2D_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRJet_PThrust2D.push_back( TH2F(name.Data(),name.Data(), 40, -500.,500., 40, -700.0, 600.0) );
    
    name = "Ana_AllHad_2b_Thrust_bJet_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -700.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_bJet_PtThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet_PtThrust.push_back( TH1F(name.Data(),name.Data(), 400, -500.0, 500.0) );
    name = "Ana_AllHad_2b_Thrust_bJet_PThrust2D_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet_PThrust2D.push_back( TH2F(name.Data(),name.Data(), 40, -500.,500., 40, -700.0, 600.0) );

    //----------------------------------------------------------//

    name = "Ana_AllHad_2b_Thrust_METISR_PtDirRatio_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METISR_PtDirRatio.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_MaxJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MaxJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_MinJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MinJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    
    name = "Ana_AllHad_2b_Thrust_n_TTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_n_TTjets.push_back( TH1F(name.Data(),name.Data(), 8, -0.5, 7.5) );
    name = "Ana_AllHad_2b_Thrust_n_bTTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_n_bTTjets.push_back( TH1F(name.Data(),name.Data(), 4, -0.5, 3.5) );
    
    name = "Ana_AllHad_2b_Thrust_nISRJets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_nISRJets.push_back( TH1F(name.Data(),name.Data(), 6, -0.5, 5.5) );
    
    name = "Ana_AllHad_2b_Thrust_ISRMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 800.) );
    name = "Ana_AllHad_2b_Thrust_TTMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TTMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.) );

    name = "Ana_AllHad_2b_Thrust_ISR_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISR_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_ISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 800, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_ISRMet_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRMet_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_METoverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METoverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_METoverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METoverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );
    
    
    name = "Ana_AllHad_2b_Thrust_TT_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TT_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TToverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TToverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TToverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TToverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_METoverTT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METoverTT_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_METoverTT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METoverTT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_ISRTT_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISRTT_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_TTMET_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TTMET_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    
    name = "Ana_AllHad_2b_Thrust_ISR_MET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISR_MET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_MET_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MET_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
name = "Ana_AllHad_2b_Thrust_ISR_TT_pt_"+istr+"cut";
Ana_AllHad_2b_Thrust_ISR_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
name = "Ana_AllHad_2b_Thrust_ISR_MET_pthrust_"+istr+"cut";
Ana_AllHad_2b_Thrust_ISR_MET_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, 0.0, 600.0) );
name = "Ana_AllHad_2b_Thrust_MET_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MET_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, 0., 600., 100, -900.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_ISR_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_ISR_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, -900.0, 900.0) );

    //-------------------------------------------------------//

    name = "Ana_AllHad_2b_Thrust_multISR_ISR_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_ISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_ISRMet_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISRMet_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_multISR_METoverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_METoverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_METoverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_METoverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_multISR_TT_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_TT_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_TT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_TT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_TToverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_TToverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 6.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_TToverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_TToverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 6.0) );

    name = "Ana_AllHad_2b_Thrust_multISR_METoverTT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_METoverTT_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 6.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_METoverTT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_METoverTT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 6.0) );

    name = "Ana_AllHad_2b_Thrust_multISR_ISRTT_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISRTT_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_multISR_TTMET_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_TTMET_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );

    name = "Ana_AllHad_2b_Thrust_multISR_ISR_MET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_MET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_MET_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_MET_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_ISR_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_ISR_MET_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_MET_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, 0.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_MET_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_MET_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, 0., 600., 100, -900.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_multISR_ISR_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_multISR_ISR_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, -900.0, 900.0) );

    //----------------------------------------------------------// 

    name = "Ana_AllHad_2b_Thrust_Met_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MET.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1200.) );
    name = "Ana_AllHad_2b_Thrust_HT_"+istr+"cut";
    Ana_AllHad_2b_Thrust_HT.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1500.) );
    name = "Ana_AllHad_2b_Thrust_CT_"+istr+"cut";
    Ana_AllHad_2b_Thrust_CT.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1500.) );

    name = "Ana_AllHad_2b_Thrust_nJets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_nJets.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );
    name = "Ana_AllHad_2b_Thrust_nHighPtJet_"+istr+"cut";
    Ana_AllHad_2b_Thrust_nHighPtJets.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );

    name = "Ana_AllHad_2b_Thrust_NJets_NegHemi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_NJets_NegHemi.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );
    name = "Ana_AllHad_2b_Thrust_NJets_PosHemi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_NJets_PosHemi.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );

    name = "Ana_AllHad_2b_Thrust_Jet1Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet1Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1200.) );
    name = "Ana_AllHad_2b_Thrust_Jet1Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet1Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet2Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet2Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 400. ) );
    name = "Ana_AllHad_2b_Thrust_Jet2Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet2Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet3Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet3Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 300. ) );
    name = "Ana_AllHad_2b_Thrust_Jet3Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet3Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet4Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet4Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 300. ) );
    name = "Ana_AllHad_2b_Thrust_Jet4Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet4Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet5Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet5Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 200. ) );
    name = "Ana_AllHad_2b_Thrust_Jet5Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet5Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet6Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet6Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 200. ) );
    name = "Ana_AllHad_2b_Thrust_Jet6Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet6Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet7Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet7Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 200. ) );
    name = "Ana_AllHad_2b_Thrust_Jet7Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet7Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet8Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet8Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 200. ) );
    name = "Ana_AllHad_2b_Thrust_Jet8Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet8Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_Jet9Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet9Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 200. ) );
    name = "Ana_AllHad_2b_Thrust_Jet9Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_Jet9Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_minB_Mt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_minB_Mt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 500 ) );
    name = "Ana_AllHad_2b_Thrust_METsig_"+istr+"cut";
    Ana_AllHad_2b_Thrust_METsig.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 20 ) );
    name = "Ana_AllHad_2b_Thrust_MET_over_ISR_"+istr+"cut";
    Ana_AllHad_2b_Thrust_MET_over_ISR.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 10. ) );
    name = "Ana_AllHad_2b_Thrust_bJet1Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet1Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 600. ) );
    name = "Ana_AllHad_2b_Thrust_bJet1Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet1Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_bJet1_ijet_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet1_ijet.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5 ) );
    name = "Ana_AllHad_2b_Thrust_bJet2Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet2Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 400. ) );
    name = "Ana_AllHad_2b_Thrust_bJet2Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet2Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_bJet2_ijet_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet2_ijet.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );
    name = "Ana_AllHad_2b_Thrust_bJet3Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet3Pt.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 300. ) );
    name = "Ana_AllHad_2b_Thrust_bJet3Eta_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet3Eta.push_back( TH1F(name.Data(),name.Data(), 200, -5.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_bJet3_ijet_"+istr+"cut";
    Ana_AllHad_2b_Thrust_bJet3_ijet.push_back( TH1F(name.Data(),name.Data(), 21, -0.5, 20.5) );

    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_FoundTop_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_FoundTop.push_back( TH1F(name.Data(),name.Data(), 2, -0.5, 1.5 ) );

    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_Mass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_Mass.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 1000.0 ) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust.push_back( TH1F(name.Data(),name.Data(), 500, -700.0, 600.0 ) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_met_deltaphi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_met_deltaphi.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, TMath::Pi() ) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_thrust_deltaphi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_thrust_deltaphi.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, TMath::Pi() ) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_DeltaR_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_DeltaR.push_back( TH1F(name.Data(),name.Data(), 200, 0.0, 5.0 ) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust_mass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust_mass.push_back( TH2F(name.Data(),name.Data(), 40, -700.0, 600.0, 40, 0.0, 1000.) );
    name = "Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_delta_pthrust_mass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_delta_pthrust_mass.push_back( TH2F(name.Data(),name.Data(), 40, -700.0, 600.0, 40, 0.0, 1000.) );

    //----------------------------------------------------------------//

    name = "Ana_AllHad_2b_Thrust_TruthTopID_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthTopID.push_back( TH1F(name.Data(),name.Data(),4,-0.5,3.5) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_METISR_PtDirRatio_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_METISR_PtDirRatio.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_MaxJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_MaxJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_MinJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_MinJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_n_TTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_n_TTjets.push_back( TH1F(name.Data(),name.Data(), 8, -0.5, 7.5) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_n_bTTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_n_bTTjets.push_back( TH1F(name.Data(),name.Data(), 4, -0.5, 3.5) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_nISRJets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_nISRJets.push_back( TH1F(name.Data(),name.Data(), 6, -0.5, 5.5) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISRMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISRMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 800.) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TTMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TTMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 800, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISRMet_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISRMet_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TT_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TT_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_PThrust_"+istr+"cut";
Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISRTT_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISRTT_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TTMET_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TTMET_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, 0.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, 0., 600., 100, -900.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, -900.0, 900.0) );

    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TruthISR_RecoISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TruthISR_RecoISR_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthAllHad_TruthMET_RecoMET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthAllHad_TruthMET_RecoMET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );

    //-------------------------------------------------------//

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_METISR_PtDirRatio_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_METISR_PtDirRatio.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_MaxJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_MaxJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_MinJetBMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_MinJetBMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_n_TTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_n_TTjets.push_back( TH1F(name.Data(),name.Data(), 8, -0.5, 7.5) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_n_bTTjets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_n_bTTjets.push_back( TH1F(name.Data(),name.Data(), 4, -0.5, 3.5) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_nISRJets_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_nISRJets.push_back( TH1F(name.Data(),name.Data(), 6, -0.5, 5.5) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISRMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISRMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 800.) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TTMass_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TTMass.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 800, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISRMet_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISRMet_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TT_Pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TT_Pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -1000.0, 1000.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_pt.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, 5.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_PThrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_PThrust.push_back( TH1F(name.Data(),name.Data(), 400, -5.0, 5.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISRTT_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISRTT_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TTMET_DPhi_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TTMET_DPhi.push_back( TH1F(name.Data(),name.Data(), 400, 0.0, TMath::Pi()) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, 0.0, 600.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, 0., 600., 100, -900.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pthrust_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pthrust.push_back( TH2F(name.Data(),name.Data(), 100, -900., 900., 100, -900.0, 900.0) );

    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TruthISR_RecoISR_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TruthISR_RecoISR_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
    name = "Ana_AllHad_2b_Thrust_TruthNotHad_TruthMET_RecoMET_pt_"+istr+"cut";
    Ana_AllHad_2b_Thrust_TruthNotHad_TruthMET_RecoMET_pt.push_back( TH2F(name.Data(),name.Data(), 100, 0., 900., 100, 0.0, 900.0) );
  }

  if (verbose) std::cout << "booked histos" << std::endl;

  return;

}

void CollectionTree::WriteHistos()
{

  if ( verbose ) std::cout << "writing" << std::endl;

  output = new TFile(output_name.c_str(), "RECREATE");

  All_nlep->Write();
  All_nHighPtLep->Write();
  All_lep1Pt->Write();
  All_lep1Eta->Write();
  All_lep1Iso->Write();
  All_lep2Pt->Write();
  All_lep2Eta->Write();
  All_lep2Iso->Write();
  All_MET->Write();
  All_METPhi->Write();
  All_HT->Write();
  All_CT->Write();
  All_nJets->Write();
  All_Jet1Pt->Write();
  All_Jet1Eta->Write();
  All_Jet2Pt->Write();
  All_Jet2Eta->Write();
  All_DiHiPtJet_DeltaPhi->Write();
  All_Mt->Write();

  AfterTrigger_nlep->Write();
  AfterTrigger_nHighPtLep->Write();
  AfterTrigger_lep1Pt->Write();
  AfterTrigger_lep1Eta->Write();
  AfterTrigger_lep1Iso->Write();
  AfterTrigger_lep2Pt->Write();
  AfterTrigger_lep2Eta->Write();
  AfterTrigger_lep2Iso->Write();
  AfterTrigger_MET->Write();
  AfterTrigger_METPhi->Write();
  AfterTrigger_HT->Write();
  AfterTrigger_CT->Write();
  AfterTrigger_nJets->Write();
  AfterTrigger_Jet1Pt->Write();
  AfterTrigger_Jet1Eta->Write();
  AfterTrigger_Jet2Pt->Write();
  AfterTrigger_Jet2Eta->Write();
  AfterTrigger_DiHiPtJet_DeltaPhi->Write();
  AfterTrigger_Mt->Write();

  PreSel_nlep->Write();
  PreSel_nHighPtLep->Write();
  PreSel_Mu1Pt->Write();
  PreSel_Mu1Eta->Write();
  PreSel_Mu1Iso->Write();
  PreSel_Mu2Pt->Write();
  PreSel_Mu2Eta->Write();
  PreSel_Mu2Iso->Write();
  PreSel_Ele1Pt->Write();
  PreSel_Ele1Eta->Write();
  PreSel_Ele1Iso->Write();
  PreSel_Ele2Pt->Write();
  PreSel_Ele2Eta->Write();
  PreSel_Ele2Iso->Write();
  PreSel_MET->Write();
  PreSel_METPhi->Write();
  PreSel_HT->Write();
  PreSel_CT->Write();
  PreSel_nJets->Write();
  PreSel_JetISR_Pt->Write();
  PreSel_JetISR_Eta->Write();
  PreSel_Jet1Pt->Write();
  PreSel_Jet1Eta->Write();
  PreSel_Jet2Pt->Write();
  PreSel_Jet2Eta->Write();
  PreSel_DiHiPtJet_DeltaPhi->Write();
  PreSel_Mt->Write();

  for ( int i=0; i<Ana_AllHad_2b_Thrust_Thrust.size(); i++ ) {

    Ana_AllHad_2b_Thrust_Thrust[i].Write();
    Ana_AllHad_2b_Thrust_Thrust_Met_dPhi[i].Write();

    Ana_AllHad_2b_Thrust_Jet_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_Jet_PtThrust[i].Write();
    Ana_AllHad_2b_Thrust_Jet_PThrust2D[i].Write();
    Ana_AllHad_2b_Thrust_ISRJet_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_ISRJet_PtThrust[i].Write();
    Ana_AllHad_2b_Thrust_ISRJet_PThrust2D[i].Write();
    Ana_AllHad_2b_Thrust_bJet_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_bJet_PtThrust[i].Write();
    Ana_AllHad_2b_Thrust_bJet_PThrust2D[i].Write();

    Ana_AllHad_2b_Thrust_METISR_PtDirRatio[i].Write();
    Ana_AllHad_2b_Thrust_n_TTjets[i].Write();
    Ana_AllHad_2b_Thrust_n_bTTjets[i].Write();

    Ana_AllHad_2b_Thrust_MaxJetBMass[i].Write();
    Ana_AllHad_2b_Thrust_MinJetBMass[i].Write();

    Ana_AllHad_2b_Thrust_nISRJets[i].Write();

    Ana_AllHad_2b_Thrust_ISRMass[i].Write();
    Ana_AllHad_2b_Thrust_TTMass[i].Write();

    Ana_AllHad_2b_Thrust_ISR_Pt[i].Write();
    Ana_AllHad_2b_Thrust_ISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_ISRMet_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_METoverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_METoverISR_PThrust[i].Write();

    Ana_AllHad_2b_Thrust_TT_Pt[i].Write();
    Ana_AllHad_2b_Thrust_TT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TToverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TToverISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_METoverTT_pt[i].Write();
    Ana_AllHad_2b_Thrust_METoverTT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_ISRTT_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_TTMET_DPhi[i].Write();

    Ana_AllHad_2b_Thrust_ISR_MET_pt[i].Write();
    Ana_AllHad_2b_Thrust_MET_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_ISR_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_ISR_MET_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_MET_TT_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_ISR_TT_pthrust[i].Write();

    Ana_AllHad_2b_Thrust_multISR_ISR_Pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISRMet_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_multISR_METoverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_METoverISR_PThrust[i].Write();

    Ana_AllHad_2b_Thrust_multISR_TT_Pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_TT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_TToverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_TToverISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_METoverTT_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_METoverTT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISRTT_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_multISR_TTMET_DPhi[i].Write();

    Ana_AllHad_2b_Thrust_multISR_ISR_MET_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_MET_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISR_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISR_MET_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_MET_TT_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_multISR_ISR_TT_pthrust[i].Write();

    Ana_AllHad_2b_Thrust_NJets_NegHemi[i].Write();
    Ana_AllHad_2b_Thrust_NJets_PosHemi[i].Write();

    Ana_AllHad_2b_Thrust_HT[i].Write();
    Ana_AllHad_2b_Thrust_CT[i].Write();
    Ana_AllHad_2b_Thrust_nJets[i].Write();
    Ana_AllHad_2b_Thrust_nHighPtJets[i].Write();
    Ana_AllHad_2b_Thrust_Jet1Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet1Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet2Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet2Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet3Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet3Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet4Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet4Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet5Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet5Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet6Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet6Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet7Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet7Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet8Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet8Eta[i].Write();
    Ana_AllHad_2b_Thrust_Jet9Pt[i].Write();
    Ana_AllHad_2b_Thrust_Jet9Eta[i].Write();
    Ana_AllHad_2b_Thrust_minB_Mt[i].Write();
    Ana_AllHad_2b_Thrust_METsig[i].Write();
    Ana_AllHad_2b_Thrust_MET_over_ISR[i].Write();
    Ana_AllHad_2b_Thrust_bJet1Pt[i].Write();
    Ana_AllHad_2b_Thrust_bJet1Eta[i].Write();
    Ana_AllHad_2b_Thrust_bJet1_ijet[i].Write();
    Ana_AllHad_2b_Thrust_bJet2Pt[i].Write();
    Ana_AllHad_2b_Thrust_bJet2Eta[i].Write();
    Ana_AllHad_2b_Thrust_bJet2_ijet[i].Write();
    Ana_AllHad_2b_Thrust_bJet3Pt[i].Write();
    Ana_AllHad_2b_Thrust_bJet3Eta[i].Write();
    Ana_AllHad_2b_Thrust_bJet3_ijet[i].Write();

    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_FoundTop[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_Mass[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_met_deltaphi[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_top_thrust_deltaphi[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_DeltaR[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_pthrust_mass[i].Write();
    Ana_AllHad_2b_Thrust_T_MaxPThrust_WMinDeltaR_delta_pthrust_mass[i].Write();

    Ana_AllHad_2b_Thrust_TruthTopID[i].GetXaxis()->SetBinLabel(1, "Truth All Hadronic");
    Ana_AllHad_2b_Thrust_TruthTopID[i].GetXaxis()->SetBinLabel(2, "Truth Single Tau");
    Ana_AllHad_2b_Thrust_TruthTopID[i].GetXaxis()->SetBinLabel(3, "Truth Single lepton");
    Ana_AllHad_2b_Thrust_TruthTopID[i].GetXaxis()->SetBinLabel(4, "Truth Di-leptonic/Di tau");

    Ana_AllHad_2b_Thrust_TruthTopID[i].Write();

    //--------------------------------------------------//

    Ana_AllHad_2b_Thrust_TruthAllHad_METISR_PtDirRatio[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_n_TTjets[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_n_bTTjets[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_MaxJetBMass[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_MinJetBMass[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_nISRJets[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_ISRMass[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TTMass[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_Pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISRMet_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverISR_PThrust[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_TT_Pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TToverISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_METoverTT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISRTT_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TTMET_DPhi[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_MET_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_MET_TT_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_ISR_TT_pthrust[i].Write();

    Ana_AllHad_2b_Thrust_TruthAllHad_TruthISR_RecoISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthAllHad_TruthMET_RecoMET_pt[i].Write();

    //--------------------------------------------------//         

    Ana_AllHad_2b_Thrust_TruthNotHad_METISR_PtDirRatio[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_n_TTjets[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_n_bTTjets[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_MaxJetBMass[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_MinJetBMass[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_nISRJets[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_ISRMass[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TTMass[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_Pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISRMet_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverISR_PThrust[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_TT_Pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TToverISR_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_METoverTT_PThrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISRTT_DPhi[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TTMET_DPhi[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_MET_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_MET_TT_pthrust[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_ISR_TT_pthrust[i].Write();

    Ana_AllHad_2b_Thrust_TruthNotHad_TruthISR_RecoISR_pt[i].Write();
    Ana_AllHad_2b_Thrust_TruthNotHad_TruthMET_RecoMET_pt[i].Write();
  }

  output->Close();

  return;

}
