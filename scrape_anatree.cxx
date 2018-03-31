#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "dwall.h"

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
  float lepke;
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
  int mctrk_origin[1000];
  int mctrk_TrackId[1000];
  int mctrk_MotherTrackId[1000];
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
  tree->SetBranchAddress("mctrk_origin", mctrk_origin);
  tree->SetBranchAddress("mctrk_TrackId", mctrk_TrackId);
  tree->SetBranchAddress("mctrk_MotherTrackId", mctrk_MotherTrackId);
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
  int mcshwr_origin[1000];
  int mcshwr_TrackId[1000];
  int mcshwr_MotherTrkId[1000];
  int mcshwr_Process[1000];
  float mcshwr_startX[1000];
  float mcshwr_startY[1000];
  float mcshwr_startZ[1000];  
  float mcshwr_endX[1000];
  float mcshwr_endY[1000];
  float mcshwr_endZ[1000];
  int mcshwr_Ancestorpdg[1000];
  int mcshwr_AncestorTrkId[1000];
  int mcshwr_AncestorProcess[1000];
  int mcshwr_isEngDeposited[1000];
  float mcshwr_CombEngX[1000];
  float mcshwr_CombEngY[1000];
  float mcshwr_CombEngZ[1000];  
  tree->SetBranchAddress("no_mcshowers", &no_mcshowers);
  tree->SetBranchAddress("mcshwr_pdg", mcshwr_pdg);
  tree->SetBranchAddress("mcshwr_origin", mcshwr_origin);
  tree->SetBranchAddress("mcshwr_TrackId", mcshwr_TrackId);
  tree->SetBranchAddress("mcshwr_MotherTrkId", mcshwr_MotherTrkId);
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
  tree->SetBranchAddress("mcshwr_isEngDeposited", mcshwr_isEngDeposited);
  tree->SetBranchAddress("mcshwr_CombEngX", mcshwr_CombEngX );
  tree->SetBranchAddress("mcshwr_CombEngY", mcshwr_CombEngY );
  tree->SetBranchAddress("mcshwr_CombEngZ", mcshwr_CombEngZ );  
  
  double pot;
  pottree->SetBranchAddress("pot",&pot);

  TFile* outfile = new TFile("output_scrapedtree.root", "RECREATE");
  TTree* scraped = new TTree("scrapedana", "Scraped Anatree");
  // mode/flavors
  int mode;
  int ccnc;
  int nufluxpdg;
  int nufluxparentpdg;
  int nuxsecpdg;  
  float fluxweight;
  float nuvtx[3];
  float nuvtxsce[3];
  // kinematics
  float enugev;
  float q2truth;
  float wtruth;
  float lmom;
  // primary counters
  int nmuon1;
  int nelectron1;
  int npi01;
  int nchargedpi1;
  int nchargedpi35mev1;
  int nprotons1;
  int nprotons60mev1;
  int nshowers; // all (means gamma)
  int nobsgamma;
  // secondary counters
  int nmuon2;
  int nelectron2;
  int npi02;
  int nchargedpi2;  
  int nchargedpi35mev2;  
  int nprotons2;
  int nprotons60mev2;
  // truth cuts
  float protonmaxke;           // proton threshold
  float closestshowerdist;     // all showers
  float closestelectron1dist;  // primary electron
  float closestelectron2dist;  // secondary electron
  float closestpi0shower1dist; // primary pi0 showers
  int nobspi0gamma;
  int nmissedpi0gamma;  
  float lepdwall;
  float vtxdwall;

  // mode/flavors
  scraped->Branch("mode",&mode,"mode/I");       // mode_truth
  scraped->Branch("ccnc",&ccnc,"ccnc/I");       // ccnc_truth
  scraped->Branch("nufluxpdg",&nufluxpdg,"nufluxpdg/I"); // nuPDG_truth
  scraped->Branch("nufluxparentpdg",&nufluxparentpdg,"nufluxparentpdg/I"); // 
  scraped->Branch("nuxsecpdg",&nuxsecpdg,"nuxsecpdg/I"); 
  scraped->Branch("fluxweight",&fluxweight,"fluxweight/F"); 
  scraped->Branch("nuvtx", nuvtx, "nuvtx[3]/F");     // nuvtxx_truth
  scraped->Branch("nuvtxsce", nuvtxsce, "nuvtxsce[3]/F");     // nuvtxx_truth  
  // kinematics
  scraped->Branch("enugev",&enugev,"enugev/F"); // enu_truth
  scraped->Branch("q2truth", &q2truth, "q2truth/F");  // Q2_truth
  scraped->Branch("wtruth", &wtruth, "wtruth/F");     // W_truth
  scraped->Branch("lmom", &lmom, "lmom/F" );          // lep_mom_truth
  scraped->Branch("lepke", &lepke, "lepke/F");
  // primary counters
  scraped->Branch("nmuon1", &nmuon1, "nmuon1/I");
  scraped->Branch("nelectron1", &nelectron1, "nelectron1/I");
  scraped->Branch("npi01", &npi01, "npi01/I");
  scraped->Branch("nchargedpi1", &nchargedpi1, "nchargedpi1/I");
  scraped->Branch("nchargedpi35mev1", &nchargedpi35mev1, "nchargedpi35mev1/I");
  scraped->Branch("nprotons1", &nprotons1, "nprotons1/I");
  scraped->Branch("nprotons60mev1", &nprotons60mev1, "nprotons60mev1/I" ); // derived from mctrk
  scraped->Branch("nshowers", &nshowers, "nshowers/I" ); // derived from mcshower
  scraped->Branch("nobsgamma", &nobsgamma, "nobsgamma/I" ); // derived from mcshower
  // secondary counters
  scraped->Branch("nelectron2", &nelectron2, "nelectron2/I");
  // truth cuts
  scraped->Branch("protonmaxke", &protonmaxke, "protonmaxke/F" ); // derived from mctrk
  scraped->Branch("closestshowerdist", &closestshowerdist, "closestshowerdist/F" );  
  scraped->Branch("closestelectron1dist", &closestelectron1dist, "closestelectron1dist/F" );  
  scraped->Branch("closestelectron2dist", &closestelectron2dist, "closestelectron2dist/F" );  
  scraped->Branch("lepdwall", &lepdwall, "lepdwall/F" );
  scraped->Branch("vtxdwall", &vtxdwall, "vtxdwall/F" );
  

  TTree* scrapedpot = new TTree("pot", "Scraped POT");
  scrapedpot->Branch("pot",&pot,"pot/D");

  ULong_t ientry = 0;
  ULong_t bytes = tree->GetEntry(ientry);
  while (bytes!=0) {

    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Anatree entry " << ientry << std::endl;


    // primary particles
    bool haslepton = false;
    int primleptonid = -1;
    nmuon1 = 0;
    nelectron1 = 0;
    closestelectron1dist = 1.0e6;
    for (int inu=0; inu<1; inu++) {
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
      lepke = -1.0;
      if ( ccnc==0 )
	lmom = lep_mom_truth*1000.0;
      else
	lmom = 0;

      if ( abs(nufluxpdg)==12 ) { 
	lepke = sqrt(lmom*lmom + 0.911*0.911)-0.911;
	if ( ccnc==0 ) {
	  nelectron1++;
	  closestelectron1dist = 0.0;
	}
	
      }
      else if ( abs(nufluxpdg)==14 ) {
	lepke = sqrt(lmom*lmom + 105.0*105.0)-105.0;
	if ( ccnc==0 )
	  nmuon1++;
      }

      for ( auto& evtwgt_types : (*pevtwgt_weight) ) {
	for ( auto& evtwgt : evtwgt_types ) {
	  if ( !std::isnan(evtwgt) && !std::isinf(evtwgt) )
	    fluxweight *= evtwgt;
	}
      }
      int boundary_type = 0;
      std::vector<float> nuvtx_v(3,0);
      for (int i=0; i<3; i++)
	nuvtx_v[i] = nuvtx[i];
      vtxdwall = anascraper::dwall( nuvtx_v, boundary_type );
    }
    std::cout << "  Mode: " << mode << " CCNC=" << ccnc << std::endl;

    // final state information
    // get lepton end point
    // get proton ke and count

    // primary counters
    npi01 = 0;
    nchargedpi1 = 0;
    nchargedpi35mev1 = 0;
    nprotons1 = 0;
    nprotons60mev1 = 0;
    nshowers = 0;
    nobsgamma = 0;
    // secondary counters
    nelectron2 = 0;
    npi02 = 0;
    nchargedpi2 = 0;
    nchargedpi35mev2 = 0;
    nprotons2 = 0;
    nprotons60mev2 = 0;
    // truth cuts
    protonmaxke = 0.0;
    closestshowerdist = 1.0e6;
    closestelectron2dist = 1.0e6;
    closestpi0shower1dist = 1.0e6;    
    lepdwall = 1000;

    // loop over tracks
    for (int itrk=0; itrk<no_mctracks; itrk++) {

      if ( mctrk_origin[itrk]!=1 )
	continue;

      float dist=0;
      dist += (nuvtx[0]-mctrk_startX[itrk])*(nuvtx[0]-mctrk_startX[itrk]);
      dist += (nuvtx[1]-mctrk_startY[itrk])*(nuvtx[1]-mctrk_startY[itrk]);
      dist += (nuvtx[2]-mctrk_startZ[itrk])*(nuvtx[2]-mctrk_startZ[itrk]);
      dist = sqrt(dist);
      
      if ( mctrk_pdg[itrk]==2212 )
	std::cout << "  nu-origin proton:"
		  << " id=" << mctrk_TrackId[itrk]
		  << " mother=" << mctrk_MotherTrackId[itrk]
		  << " process=" << mctrk_Process[itrk] 
		  << " ancestor=" << mctrk_Ancestorpdg[itrk] 
		  << " dist2vtx=" << dist << std::endl;


      // look at tracks coming from the vertex
      if ( dist>1.0e-3 )
	continue;
	
      if ( mctrk_pdg[itrk]==2212 ) {
	
	float ke = ( sqrt( mctrk_p_drifted[itrk]*mctrk_p_drifted[itrk] + 938*938 )-938 );
	std::cout << "  primary proton:"
		  << " ke=" << ke 
		  << " id=" << mctrk_TrackId[itrk]
		  << " mother=" << mctrk_MotherTrackId[itrk]
		  << " process=" << mctrk_Process[itrk] 
		  << " ancestor=" << mctrk_Ancestorpdg[itrk] 
		  << " dist2vtx=" << dist << std::endl;
	if ( ke>60 )
	  nprotons60mev1++;
	if ( ke>protonmaxke ) {
	  protonmaxke = ke;
	}
	nprotons1++;
      }
      else if ( mctrk_pdg[itrk]==13 || mctrk_pdg[itrk]==-13 ) {
	// primary muon
	std::vector<float> endpos(3);
	endpos[0] = mctrk_endX[itrk];
	endpos[1] = mctrk_endY[itrk];
	endpos[2] = mctrk_endZ[itrk];
	int boundary_type = 0;
	lepdwall = anascraper::dwall( endpos, boundary_type );
      }
      else if ( mctrk_pdg[itrk]==111 ) {
	npi01++;

	// shower info
	std::cout << "  Chase down pi0 gammas: pi0 id=" << mctrk_TrackId[itrk] << std::endl;
	for ( int ishw=0; ishw<no_mcshowers; ishw++) {
	  if ( mcshwr_AncestorTrkId[ishw]==mctrk_TrackId[itrk] ) {

	    float shwdist = -1;
	    if ( mcshwr_isEngDeposited[ishw]>0 ) {
	      nobspi0gamma++;

	      float shwdir[3];
	      shwdir[0] = mcshwr_CombEngX[ishw]-nuvtx[0];
	      shwdir[1] = mcshwr_CombEngY[ishw]-nuvtx[1];
	      shwdir[2] = mcshwr_CombEngZ[ishw]-nuvtx[2];
	      shwdist = sqrt( shwdir[0]*shwdir[0] + shwdir[1]*shwdir[1] + shwdir[2]*shwdir[2] );
	      if ( (closestpi0shower1dist<0 || shwdist<closestpi0shower1dist) ) {
		closestpi0shower1dist = shwdist;
	      }
	    }// if energy deposited
	    
	    std::cout << "  pi0-mcshower: " << mcshwr_TrackId[ishw]
		      << " pdg=" << mcshwr_pdg[ishw]
		      << " ancestor=" << mcshwr_AncestorTrkId[ishw]
		      << " shwdist=" << shwdist
		      << " is_edep=" << mcshwr_isEngDeposited[ishw]
		      << " start=(" << mcshwr_startX[ishw] << "," << mcshwr_startY[ishw] << "," << mcshwr_startZ[ishw] << ")"
		      << " eng=(" << mcshwr_CombEngX[ishw] << "," << mcshwr_CombEngY[ishw] << "," << mcshwr_CombEngZ[ishw] << ")"
		      << std::endl;
	    if ( shwdist<0 )
	      nmissedpi0gamma++;
	  }//end of if gamma from pi0
	}//end of shwr loop
	//if ( nobspi0gamma>0 )
	//	  std::cin.get();
      }
      else if ( mctrk_pdg[itrk]==211 || mctrk_pdg[itrk]==-211 ) {
	// charged pions
	nchargedpi1++;
	float ke_pi = ( sqrt( mctrk_p_drifted[itrk]*mctrk_p_drifted[itrk] + 139.5*139.5 )-139.5 );
	if ( ke_pi>35.0 )
	  nchargedpi35mev1++;
      }
      else {
	std::cout << "  other pdg: " << mctrk_pdg[itrk] << std::endl;
      }
    }//end of mc track
    std::cout << "  number of primary >60 mev protons: " << nprotons60mev1 << ", maxke=" << protonmaxke << " nprotons=" << nprotons1 << std::endl;
    std::cout << "  number of primary >35 mev chargedpi: " << nchargedpi35mev1 << " nchargedpi=" << nchargedpi1 << std::endl;

    
    // shower info
    for ( int ishw=0; ishw<no_mcshowers; ishw++) {

      if (mcshwr_origin[ishw]!=1)
	continue;

      std::cout << "  shower at nu-origin=" << mcshwr_origin[ishw] << " pdg=" << mcshwr_pdg[ishw] 
		<< " id=" << mcshwr_TrackId[ishw]
		<< " motherid=" << mcshwr_MotherTrkId[ishw] << std::endl;

      float shwdist = -1;
      if ( mcshwr_isEngDeposited[ishw]>0 ) {
	if ( mcshwr_pdg[ishw]==22 )
	  nobsgamma++;
	if ( abs(mcshwr_pdg[ishw])==11 ) {
	  if ( mcshwr_TrackId[ishw]==mcshwr_MotherTrkId[ishw] ) {
	    //nobselectron1++;
	  }
	  else {
	    nelectron2++;
	  }
	}

	float shwdir[3];
	shwdir[0] = mcshwr_CombEngX[ishw]-nuvtx[0];
	shwdir[1] = mcshwr_CombEngY[ishw]-nuvtx[1];
	shwdir[2] = mcshwr_CombEngZ[ishw]-nuvtx[2];
	shwdist = sqrt( shwdir[0]*shwdir[0] + shwdir[1]*shwdir[1] + shwdir[2]*shwdir[2] );
	if ( mcshwr_pdg[ishw]==22 && (closestshowerdist<0 || shwdist<closestshowerdist) ) {
	  closestshowerdist = shwdist;
	}

	if ( abs(mcshwr_pdg[ishw])==11 ) {
	  if (mcshwr_TrackId[ishw]==mcshwr_MotherTrkId[ishw] ) {
	    if ( closestelectron1dist<0 || shwdist<closestelectron1dist ) {
	      closestelectron1dist = shwdist;
	    }
	  }
	  else {
	    if ( closestelectron2dist<0 || shwdist<closestelectron2dist ) {
	      closestelectron2dist = shwdist;
	    }
	  }
	}

	// if ( shwdist<0.6 )
	//   nobsgamma2pix++;
      }
    }
    std::cout << "  Primary counters" << std::endl;
    std::cout << "    muon=" << nmuon1 << std::endl;
    std::cout << "    electron=" << nelectron1 << " closestdist=" << closestelectron1dist <<  std::endl;
    std::cout << "    nshowers=" << nshowers   << " closestdist=" << closestshowerdist << std::endl;
    std::cout << "    nobsgamma=" << nobsgamma  << std::endl;
    std::cout << "    npi01=" << npi01  << std::endl;
    std::cout << "    nchargedpi1=" << nchargedpi1  << std::endl;
    std::cout << "    nchargedpi35mev1=" << nchargedpi35mev1  << std::endl;
    std::cout << "    nprotons1=" << nprotons1  << std::endl;
    std::cout << "    nprotons60mev1=" << nprotons60mev1  << std::endl;
    std::cout << "  Secondary counters" << std::endl;
    std::cout << "    electron=" << nelectron2 << " closestdist=" << closestelectron2dist << std::endl;

      
    scraped->Fill();

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
