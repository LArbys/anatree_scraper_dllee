#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

std::vector<std::string> parseargs( const int nargs, const char** argv ) {
  std::vector<std::string> filelist;
  
  if (nargs==1) {
    std::cout << "arguments" << std::endl;
    std::cout << " -f [single root file]" << std::endl;
    std::cout << " -d [directory with root files]" << std::endl;
    std::cout << " -l [new line separated file paths]" << std::endl;
    std::cout << "[single root file]" << std::endl;
    return filelist;
  }
  
  for ( int iarg=1; iarg<nargs; iarg++ ) {
    //std::cout << "arg: " << argv[iarg] << std::endl;
    
    if (argv[iarg][0]=='-' ) {
      if ( argv[iarg][1]=='f' ) {
	filelist.push_back( argv[++iarg] );
      }
      else if ( argv[iarg][1]=='d' ) {
	std::string dirwildcard = std::string(argv[++iarg]) + "/*.root";
	filelist.push_back(  dirwildcard );
      }
      else if ( argv[iarg][1]=='l' ) {
	std::ifstream inputlist( argv[++iarg] );
	char buffer[1028];
	std::string sbuf;
	std::string last = "__start__";
	while (inputlist.good()) {
	  inputlist >> buffer;
	  sbuf = buffer;
	  if ( sbuf!="" && sbuf!=last ) {
	    filelist.push_back( sbuf );
	    last = sbuf;
	  }
	}
	inputlist.close();
      }
      else {
	std::stringstream ss;
	ss << "Unrecognized flag/argument: " << argv[iarg] << std::endl;
	std::runtime_error( ss.str() );
      }
    }
    else {
      filelist.push_back( argv[++iarg] );
    }
    
  }// end of argument loop
  std::cout << "number of items in filelist: " << filelist.size() << std::endl;
  
  return filelist;
}

int main( const int nargs, const char** argv ) {
  
  std::cout << "===========================================" << std::endl;
  std::cout << "  Scrape MicroBooNE Anatree " << std::endl;
  std::cout << "===========================================" << std::endl;

  std::vector< std::string > inputlist = parseargs( nargs, argv );
  if ( inputlist.size()==0 )
    return 0;

  TChain* tree = new TChain("analysistree/anatree");
  TChain* pottree = new TChain("analysistree/pottree");

  for ( auto& input : inputlist ) {
    std::cout << "adding " << input << std::endl;
    tree->Add( input.c_str() );
    pottree->Add( input.c_str() );
  }

  // Variables

  // MC TRUTH
  int   mcevts_truth;
  float enu_truth[100];
  int   mode_truth[100];
  int   ccnc_truth[100];
  int   nuPDG_truth[100];
  float nuvtxx_truth[100];
  float nuvtxy_truth[100];
  float nuvtxz_truth[100];
  int   ptype_flux[100];  
  float evtwgt_weight;
  float Q2_truth;
  float W_truth;
  float lep_mom_truth;
  std::vector< std::vector<double> >* pevtwgt_weight = NULL;
  tree->SetBranchAddress("mcevts_truth", &mcevts_truth);
  tree->SetBranchAddress("enu_truth", enu_truth);
  tree->SetBranchAddress("mode_truth", mode_truth);
  tree->SetBranchAddress("ccnc_truth", ccnc_truth);
  tree->SetBranchAddress("nuPDG_truth", nuPDG_truth);
  tree->SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
  tree->SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
  tree->SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
  tree->SetBranchAddress("ptype_flux", ptype_flux);
  tree->SetBranchAddress("evtwgt_weight", &pevtwgt_weight );
  tree->SetBranchAddress("Q2_truth", &Q2_truth );
  tree->SetBranchAddress("W_truth", &W_truth );
  tree->SetBranchAddress("lep_mom_truth", &lep_mom_truth);

  // MC TRACK
  int no_mctracks;
  int mctrk_pdg[1000];
  int mctrk_TrackId[1000];
  int mctrk_Process[1000];
  float mctrk_startX[1000];
  float mctrk_startY[1000];
  float mctrk_startZ[1000];  
  float mctrk_endX[1000];
  float mctrk_endY[1000];
  float mctrk_endZ[1000];
  int mctrk_Ancestorpdg[1000];
  int mctrk_AncestorTrkId[1000];
  int mctrk_AncestorProcess[1000];
  float mctrk_p_drifted[1000];
  tree->SetBranchAddress("no_mctracks", &no_mctracks);
  tree->SetBranchAddress("mctrk_pdg", mctrk_pdg);
  tree->SetBranchAddress("mctrk_TrackId", mctrk_TrackId);
  tree->SetBranchAddress("mctrk_Process", mctrk_Process);
  tree->SetBranchAddress("mctrk_startX", mctrk_startX);  
  tree->SetBranchAddress("mctrk_startY", mctrk_startY);  
  tree->SetBranchAddress("mctrk_startZ", mctrk_startZ);
  tree->SetBranchAddress("mctrk_endX", mctrk_endX);  
  tree->SetBranchAddress("mctrk_endY", mctrk_endY);
  tree->SetBranchAddress("mctrk_endZ", mctrk_endZ);
  tree->SetBranchAddress("mctrk_endZ", mctrk_endZ);  
  tree->SetBranchAddress("mctrk_Ancestorpdg", mctrk_Ancestorpdg);
  tree->SetBranchAddress("mctrk_AncesotorTrkId", mctrk_AncestorTrkId);
  tree->SetBranchAddress("mctrk_Process", mctrk_Process);
  tree->SetBranchAddress("mctrk_p_drifted", &mctrk_p_drifted );
    
  // MC SHOWER
  int no_mcshowers;
  int mcshwr_pdg[1000];
  int mcshwr_TrackId[1000];
  int mcshwr_Process[1000];
  int mcshwr_startX[1000];
  int mcshwr_startY[1000];
  int mcshwr_startZ[1000];  
  int mcshwr_endX[1000];
  int mcshwr_endY[1000];
  int mcshwr_endZ[1000];
  int mcshwr_Ancestorpdg[1000];
  int mcshwr_AncestorTrkId[1000];
  int mcshwr_AncestorProcess[1000];
  tree->SetBranchAddress("no_mcshowers", &no_mcshowers);
  tree->SetBranchAddress("mcshwr_pdg", mcshwr_pdg);
  tree->SetBranchAddress("mcshwr_TrackId", mcshwr_TrackId);
  tree->SetBranchAddress("mcshwr_Process", mcshwr_Process);
  tree->SetBranchAddress("mcshwr_startX", mcshwr_startX);  
  tree->SetBranchAddress("mcshwr_startY", mcshwr_startY);  
  tree->SetBranchAddress("mcshwr_startZ", mcshwr_startZ);
  tree->SetBranchAddress("mcshwr_endX", mcshwr_endX);  
  tree->SetBranchAddress("mcshwr_endY", mcshwr_endY);
  tree->SetBranchAddress("mcshwr_endZ", mcshwr_endZ);
  tree->SetBranchAddress("mcshwr_Ancestorpdg", mcshwr_Ancestorpdg);
  tree->SetBranchAddress("mcshwr_AncesotorTrkId", mcshwr_AncestorTrkId);
  tree->SetBranchAddress("mcshwr_Process", mcshwr_Process);  
  
  double pot;
  pottree->SetBranchAddress("pot",&pot);

  TFile* outfile = new TFile("output_scrapedtree.root", "RECREATE");
  TTree* scraped = new TTree("scrapedana", "Scraped Anatree");
  float enugev;
  int mode;
  int ccnc;
  int nufluxpdg;
  int nufluxparentpdg;
  int nuxsecpdg;  
  float fluxweight;
  float nuvtx[3];
  float nuvtxsce[3];
  float q2truth;
  float wtruth;
  float lmom;
  float protonmaxke;
  int nprotons60mev;
  int nshowers;
  int npi0;
  float closestshowerdist;
  float lepdwall;
  scraped->Branch("enugev",&enugev,"enugev/F"); // enu_truth
  scraped->Branch("mode",&mode,"mode/I");       // mode_truth
  scraped->Branch("ccnc",&ccnc,"ccnc/I");       // ccnc_truth
  scraped->Branch("nufluxpdg",&nufluxpdg,"nufluxpdg/I"); // nuPDG_truth
  scraped->Branch("nufluxparentpdg",&nufluxparentpdg,"nufluxparentpdg/I"); // 
  scraped->Branch("nuxsecpdg",&nuxsecpdg,"nuxsecpdg/I"); 
  scraped->Branch("fluxweight",&fluxweight,"fluxweight/F"); 
  scraped->Branch("nuvtx", nuvtx, "nuvtx[3]/F");     // nuvtxx_truth
  scraped->Branch("nuvtxsce", nuvtxsce, "nuvtxsce[3]/F");     // nuvtxx_truth  
  scraped->Branch("q2truth", &q2truth, "q2truth/F");  // Q2_truth
  scraped->Branch("wtruth", &wtruth, "wtruth/F");     // W_truth
  scraped->Branch("lmom", &lmom, "lmom/F" );          // lep_mom_truth
  scraped->Branch("protonmaxke", &protonmaxke, "protonmaxke/F" ); // derived from mctrk
  scraped->Branch("nprotons60mev", &nprotons60mev, "nprotons60mev/I" );
  scraped->Branch("nshowers", &nshowers, "nshowers/I" );
  scraped->Branch("npi0", &npi0, "npi0/I" );
  scraped->Branch("closestshowerdist", &closestshowerdist, "closestshowerdist/F" );
  scraped->Branch("lepdwall", &lepdwall, "lepdwall/F" );
  

  TTree* scrapedpot = new TTree("pot", "Scraped POT");
  scrapedpot->Branch("pot",&pot,"pot/D");

  ULong_t ientry = 0;
  ULong_t bytes = tree->GetEntry(ientry);
  while (bytes!=0) {

    if ( ientry%100==0 )
      std::cout << "Anatree entry " << ientry << std::endl;
    
    for (int inu=0; inu<mcevts_truth; inu++) {

      enugev = enu_truth[inu];
      mode   = mode_truth[inu];
      ccnc   = ccnc_truth[inu];
      nufluxpdg = nuPDG_truth[inu];
      nufluxparentpdg = ptype_flux[inu];
      nuvtx[0] = nuvtxx_truth[inu];
      nuvtx[1] = nuvtxy_truth[inu];
      nuvtx[2] = nuvtxz_truth[inu];
      nuxsecpdg = nufluxpdg;
      fluxweight = 1.0;
      q2truth = Q2_truth;
      wtruth  = W_truth;
      lmom = lep_mom_truth;
      for ( auto& evtwgt_types : (*pevtwgt_weight) ) {
	for ( auto& evtwgt : evtwgt_types ) {
	  if ( !std::isnan(evtwgt) && !std::isinf(evtwgt) )
	    fluxweight *= evtwgt;
	}
      }

      // final state information
      // get lepton end point
      // get proton ke and count
      protonmaxke = 0;
      nprotons60mev = 0;
      for (int itrk=0; itrk<no_mctracks; itrk++) {
	float dist=0;
	dist += (nuvtx[0]-mctrk_startX[itrk])*(nuvtx[0]-mctrk_startX[itrk]);
	dist += (nuvtx[1]-mctrk_startY[itrk])*(nuvtx[1]-mctrk_startY[itrk]);
	dist += (nuvtx[2]-mctrk_startZ[itrk])*(nuvtx[2]-mctrk_startZ[itrk]);
	dist = sqrt(dist);

	if ( dist<0.1 )
	
	if ( mctrk_pdg[itrk]==2212 ) {
	  
	  float ke = ( sqrt( mctrk_p_drifted[itrk]*mctrk_p_drifted[itrk] + 938*938 )-938 );
	  std::cout << " primary proton: ke=" << ke << " process=" << mctrk_Process[itrk] << " ancestor=" << mctrk_Ancestorpdg[itrk] << std::endl;
	  if ( ke>60 )
	    nprotons60mev++;
	  if ( ke>protonmaxke ) {
	    protonmaxke = ke;
	  }
	}
      }
      std::cout << "number of >60 mev protons: " << nprotons60mev << ", maxke=" << protonmaxke << std::endl;
      
      scraped->Fill();
    }

    ientry++;
    bytes = tree->GetEntry(ientry);
  }

  // POT TREE LOOP
  ientry = 0;
  bytes = pottree->GetEntry(ientry);
  while ( bytes!=0 ) {
    scrapedpot->Fill();
    ientry++;
    bytes = pottree->GetEntry(ientry);    
  }

  outfile->Write();
  
  return 0;
}