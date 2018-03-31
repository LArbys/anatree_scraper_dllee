import os,sys
import ROOT as rt
from math import sqrt

h = rt.TH1D("tmp","",1000,0,10.0)

cclist = ["cc","nc"]
ccdef = {"cc":"(ccnc==0)",
         "nc":"(ccnc==1)"}

cutlist = ["qe","res","dis","coh","mec","other"]
cutdef = {"qe":"(mode==0)",
          "res":"(mode==1)",
          "dis":"(mode==2)",
          "coh":"(mode==3 || mode==4)",
          "mec":"(mode==10)",
          "other":"(mode>=5 && mode!=10)"}
cutcolor = { "ccqe":rt.kRed,
             "ccmec":rt.kRed-6,
             "ccres":rt.kMagenta+3,
             "ccdis":rt.kOrange+4,
             "cccoh":rt.kMagenta-5,
             "ccother":rt.kOrange,
             "ncqe":rt.kBlue+2,
             "ncmec":rt.kBlue-7,
             "ncres":rt.kCyan+2,
             "ncdis":rt.kCyan-5,
             "nccoh":rt.kGreen+3,
             "ncother":rt.kGreen-8 }

enurange_mev = [100,3000.0]
nbins = 29
ana_enurange = [200.0,800.0]

def calculate_pot( tpot ):
    ientry = 0
    bites = tpot.GetEntry(ientry)
    totpot = 0.
    while bites>0:
        totpot += tpot.pot
        ientry += 1
        bites = tpot.GetEntry(ientry)
    return totpot

def get_entries( cutcmd, tana, pot, potmc ):
    potratio = pot/potmc
    entry_cut = "fluxweight*(%s)"%(cutcmd)
    h.Reset()
    nevts = float(tana.Draw("enugev>>tmp",entry_cut))
    if nevts>0:
        fracerr = sqrt(nevts)/nevts
    else:
        fracerr = 0.0
    entries = h.Integral()*potratio
    print "Cut: ",cutcmd
    print "  Total: ",entries," +/- ",entries*fracerr
    for cc in cclist:
        for mode in cutlist:
            cut = "fluxweight*(%s && %s && %s)"%(cutcmd,cutdef[mode],ccdef[cc])
            #cut = "(%s && %s && %s)"%(cutcmd,cutdef[mode],ccdef[cc])
            nevts = float(tana.Draw("enugev>>tmp",cut))
            if nevts>0:
                fracerr = sqrt(nevts)/nevts
            else:
                fracerr = 0.0
            print "  [",cc,",",mode,"] ",h.Integral()*potratio," +/- ",h.Integral()*potratio*fracerr
    return entries


if __name__ == "__main__":

    #fname = "/home/twongjirad/working/data/larbys/mcc8v4/anatree/mcc8v4_cocktail_scraped_anatree.root"
    fname = "cocktail_mcc8.4_anatree_scraped_merged.root"
    #fname = "intrinsicnue_mcc8.4_anatree_scraped_merged.root"
    #fname = "intrinsicnue_mcc8.4_anatree_scraped_merged_primcounters.root"
    
    tana = rt.TChain("scrapedana")
    tana.Add( fname )

    tpot = rt.TChain("pot")
    tpot.Add( fname )

    #pot = 5.0e19
    pot = 6.6e20    
    potmc = calculate_pot( tpot )

    print "Target POT: ",pot
    print "MC POT: ",potmc
    print "POT ratio: ",pot/potmc

    # -------------------------------------------------------------------------
    # ALL NEUTRINOS

    #get_entries("vtxdwall>17.0 && abs(nufluxpdg)==14",tana,pot,potmc)
    print
    
    #get_entries("vtxdwall>17.0 && abs(nufluxpdg)==12",tana,pot,potmc)
    print
    
    #get_entries("vtxdwall>10.0",tana,pot,potmc)
    print

    print "CC-NuMu background: 1mu1p"
    cut = "vtxdwall>10.0 && abs(nufluxpdg)==14"
    cut += " && (nprotons60mev1+nmuon1+nchargedpi35mev1)==1 && (nprotons60mev1==1 || (nmuon1==1 && lepke>35.0))" # 1 track
    cut += " && ccnc==1"
    #cut += " && enugev>0.2 && enugev<0.3"
    #cut += " && enugev>0.4 && enugev<0.6"
    #cut += " && enugev>0.3 && enugev<0.475"
    #cut += " && enugev>0.475 && enugev<0.600"
    #cut += " && enugev>0.600"
    cut += " && ( nelectron1+nobsgamma+nelectron2 )==1"
    cut += " && ( (nelectron1==1 && lepke>35.0) || (nobsgamma==1 && closestshowerdist<3.0) || (nelectron2==1 && closestelectron2dist<3.0 && lepke/2.2<100.0) )"
    get_entries(cut,tana,pot,potmc)    
    print

    print "CC-nue: 1e1p"
    cut = "vtxdwall>10.0 && abs(nufluxpdg)==12"
    cut += " && (nprotons60mev1+nmuon1+nchargedpi35mev1)==1 && (nprotons60mev1==1 || (nmuon1==1 && lepke>35.0))" # 1 track
    cut += " && ccnc==0"
    #cut += " && enugev>0.2 && enugev<0.3"
    #cut += " && enugev>0.4 && enugev<0.6"
    #cut += " && enugev>0.3 && enugev<0.475"
    #cut += " && enugev>0.475 && enugev<0.600"
    cut += " && enugev>0.600"
    cut += " && ( nelectron1+nobsgamma+nelectron2 )==1"
    cut += " && ( (nelectron1==1 && lepke>35.0) || (nobsgamma==1 && closestshowerdist<3.0) || (nelectron2==1 && closestelectron2dist<3.0) )"
    #get_entries(cut,tana,pot,potmc)
    print
    
    print "NC background"
    cut = "vtxdwall>10.0 && nchargedpi==0 && nprotons60mev==1 && ccnc==1"
    #cut += " && enugev>0.2 && enugev<0.3"
    #cut += " && enugev>0.3 && enugev<0.475"
    #cut += " && enugev>0.475 && enugev<0.600"
    #cut += " && enugev>0.600"
    cut += " && (closestshowerdist<3.0 || closestelectrondist<3.0) && (nobsgamma+nobselectron>=1)"
    #get_entries(cut,tana,pot,potmc)
    print
    

    
    raw_input()
        
