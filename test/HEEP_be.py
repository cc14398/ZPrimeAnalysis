
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from array import array
import argparse
import sys
import ROOT
import json
import re
import os
import itertools
import math
import time

import EgammaUser.EgammaDAS2020.CoreTools as CoreTools
import EgammaUser.EgammaDAS2020.MathTools as MathTools

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection

def get_bincontent_2D(hist,x,y):
    x_bin = hist.GetXaxis().FindBin(x)
    y_bin = hist.GetYaxis().FindBin(y)
    if x_bin<0 or y_bin<0: 
        return None
    else:
        return hist.GetBinContent(x_bin,y_bin)

def trigger_match(ele,trig_objs,min_pt=32,max_dr=0.1):
    max_dr2 = max_dr*max_dr
    for trig_obj in trig_objs:
        if (trig_obj.id ==11 and #i.e. is trigger object an electron?
            trig_obj.pt > min_pt and 
            MathTools.cal_delta_r2(ele.eta+ele.deltaEtaSC,ele.phi,trig_obj.eta,trig_obj.phi) < max_dr2 and
            (trig_obj.filterBits & 0x2 ) !=0 ): # need () due to C++ operator precedence, checks if passed any WPTight trigger
            return True
    return False


def is_good_lumi(grl_dict,run,lumi):
    if str(run) in grl_dict:
        lumi_ranges = grl_dict[str(run)]
        for lumi_range in lumi_ranges:
            if lumi>=lumi_range[0] and lumi<=lumi_range[1]:
                return True
    return False


if __name__ == "__main__":


    #now read the cmd line to find our input filename
    #local file "file:<filename>
    #remote file "root:// 
    #note unlike CMSSW proper its not smart enough to resolve to fall over to xrootd automatically
    #so it has to be specified
    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('in_filenames',nargs="+",help='input filenames')
    parser.add_argument('--prefix','-p',default='',help='file prefix')
    parser.add_argument('--evtweight','-w',default=1.0,type=float,help='weight for histograms, 1.0 for data, ')
    parser.add_argument('--maxevents','-n',type=int,default=-1,help='max events <0, process all')
    parser.add_argument('--idsf_file','-s',type=str,default=None,help='id sf file')
    parser.add_argument('--recosf_file','-r',type=str,default=None,help='reco sf file')
    parser.add_argument('--out','-o',default="output.root",help='output filename')
    parser.add_argument('--grl','-g',default=None,help='good lumi json')
    args = parser.parse_args()

    max_events = args.maxevents if args.maxevents>0 else ROOT.TTree.kMaxEntries

    in_filenames_with_prefix = CoreTools.get_filenames(args.in_filenames,args.prefix)

    Events = ROOT.TChain("Events","chain");
    for x in in_filenames_with_prefix:
        print("adding {}".format(x))
        Events.Add(x)
    
    #ID scale factor histogram
    if args.idsf_file:
        idsf_file  = ROOT.TFile(args.idsf_file,"READ") 
        idsf_hist = idsf_file.EGamma_SF2D 
    else:
        idsf_hist = None

    #reco scale factor histogram
    if args.recosf_file:
        recosf_file  = ROOT.TFile(args.recosf_file,"READ") 
        recosf_hist = recosf_file.EGamma_SF2D 
    else:
        recosf_hist = None


    #now open the output root file
    out_file = ROOT.TFile(args.out,"RECREATE")
    #create a histogram
    ROOT.massAllHist = ROOT.TH1D("massAllHist",";M_{ee} (GeV);#entries;",40,70,110)
    ROOT.massTrigHist = ROOT.TH1D("massTrigHist",";M_{ee} (GeV);#entries;",40,70,110)
    ROOT.massIDSelHist = ROOT.TH1D("massIDSelHist",";M_{ee} (GeV);#entries;",40,70,110)
    ROOT.massPassingProbesHist = ROOT.TH1D("massPassingProbesHist",";M_{ee} (GeV);#entries;",40,70,110)
    ROOT.massFailingProbesHist = ROOT.TH1D("massFailingProbesHist",";M_{ee} (GeV);#entries;",40,70,110)
    #ROOT.massIDSFCorrHist = ROOT.TH1D("massIDSFCorrHist",";M_{ee} (GeV);#entries;",40,70,110)
    #ROOT.massESSCorrHist =  ROOT.TH1D("massESSCorrHist",";M_{ee} (GeV);#entries;",40,70,110)

    ROOT.pTAllPairsHist = ROOT.TH1D("pTAllPairsHist",";p_{T} (GeV);#entries;",50,0,100)
    ROOT.pTIDSelPairsHist = ROOT.TH1D("pTIDSelPairsHist",";p_{T} (GeV);#entries;",50,0,100)
    ROOT.pTAllProbesHist = ROOT.TH1D("pTAllProbesHist",";p_{T} (GeV);#entries;",50,0,100)
    ROOT.pTPassingProbesHist = ROOT.TH1D("pTPassingProbesHist",";p_{T} (GeV);#entries;",50,0,100)

    ROOT.etaAllPairsHist = ROOT.TH1D("EtaAllPairsHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.etaIDSelPairsHist = ROOT.TH1D("EtaIDSelPairsHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.etaAllProbesHist = ROOT.TH1D("EtaAllProbesHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.etaPassingProbesHist = ROOT.TH1D("EtaPassingProbesHist",";eta;#entries;",40,-2.6,2.6)

    ROOT.phiAllPairsHist = ROOT.TH1D("PhiAllPairsHist",";phi;#entries;",25,-3.2,3.2)
    ROOT.phiIDSelPairsHist = ROOT.TH1D("PhiIDSelPairsHist",";phi;#entries;",25,-3.2,3.2)
    ROOT.phiAllProbesHist = ROOT.TH1D("PhiAllProbesHist",";phi;#entries;",25,-3.2,3.2)
    ROOT.phiPassingProbesHist = ROOT.TH1D("PhiPassingProbesHist",";phi;#entries;",25,-3.2,3.2)
    
    ROOT.pileupAllPairsHist = ROOT.TH1D("PileupAllPairsHist",";eta;#entries;",25,-3.2,3.2)
    ROOT.pileupIDSelPairsHist = ROOT.TH1D("PileupIDSelPairsHist",";eta;#entries;",25,-3.2,3.2)
    ROOT.pileupAllProbesHist = ROOT.TH1D("PileupAllProbesHist",";eta;#entries;",25,-3.2,3.2)
    ROOT.pileupPassingProbesHist = ROOT.TH1D("PileupPassingProbesHist",";eta;#entries;",25,-3.2,3.2)

    ROOT.EtaAllPairsNoEtaReqTagHist = ROOT.TH1D("EtaAllPairsNoEtaReqTagHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaAllPairsTagHist = ROOT.TH1D("EtaAllPairsTagHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaTrigTagHist = ROOT.TH1D("EtaTrigTagHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaIDSelTagHist = ROOT.TH1D("EtaIDSelTagHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaPassingProbesTagHist = ROOT.TH1D("EtaPassingProbesTagHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaFailingProbesTagHist = ROOT.TH1D("EtaFailingProbesTagHist",";eta;#entries;",40,-2.6,2.6)

    ROOT.EtaAllPairsNoEtaReqProbeHist = ROOT.TH1D("EtaAllPairsNoEtaReqProbeHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaAllPairsNoEtaReq2DHist = ROOT.TH2D("EtaAllPairsNoEtaReq2DHist",";probe_eta;tag_eta;",40,-2.6,2.6,40,-2.6,2.6)
    ROOT.EtaAllPairsProbeHist = ROOT.TH1D("EtaAllPairsProbeHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaTrigProbeHist = ROOT.TH1D("EtaTrigProbeHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaIDSelProbeHist = ROOT.TH1D("EtaIDSelProbeHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaPassingProbesProbeHist = ROOT.TH1D("EtaPassingProbesProbeHist",";eta;#entries;",40,-2.6,2.6)
    ROOT.EtaFailingProbesProbeHist = ROOT.TH1D("EtaFailingProbesProbeHist",";eta;#entries;",40,-2.6,2.6)
    
    trig_name = "HLT_Ele32_WPTight_Gsf"
        
    trig_psetid = None
    trig_indx = None

    Events.Draw(">>elist","Sum$(Electron_pt>25 && abs(Electron_eta+Electron_deltaEtaSC)<1.4442)>=2","entrylist goff",max_events)
    elist = ROOT.gDirectory.Get('elist')
    elist.SetDirectory(0) #removing it from the file not to write it out
    Events.SetEntryList(elist)

    branchsel = BranchSelection("EgammaUser/EgammaDAS2020/data/HEEP_branches.txt")
    branchsel.selectBranches(Events)

    grl_dict = None
    if args.grl:
        with open(args.grl,'r') as f:
            grl_dict = json.load(f)

    
    nr_events = elist.GetN()
    for event_nr in range(0,nr_events):

        if args.maxevents >=0 and event_nr>=args.maxevents:
            break
        if event_nr%10000==0:
            print("processing event {} / {} {}".format(event_nr,nr_events,time.ctime()))
        
        entry_nr = Events.GetEntryNumber(event_nr)
        Events.GetEntry(entry_nr)
        event = Events

        good_lumi = is_good_lumi(grl_dict,event.run,event.luminosityBlock) if grl_dict else True
                
        eles = Collection(event,"Electron")
        trig_objs = Collection(Events,"TrigObj")

        trig_pass = getattr(event,trig_name)
        pileup = getattr(event,"PV_npvsGood")
        
        for tag,probe in itertools.permutations(eles,2):
            ROOT.EtaAllPairsNoEtaReqTagHist.Fill(tag.eta+tag.deltaEtaSC,args.evtweight)
            ROOT.EtaAllPairsNoEtaReqProbeHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)
            ROOT.EtaAllPairsNoEtaReq2DHist.Fill(tag.eta+tag.deltaEtaSC,probe.eta+probe.deltaEtaSC)
            #remember tag must be a barrel electron
            if (abs(tag.eta+tag.deltaEtaSC)<1.4442 and abs(probe.eta+probe.deltaEtaSC)>1.566):
                # and abs(probe.eta+probe.deltaEtaSC)<2.5 
                #note we have a problem with p4() due to how its declared in the c++
                #easier to just use polarP4 
                mass = math.sqrt(2*tag.pt*probe.pt*(math.cosh(tag.eta-probe.eta) - math.cos(tag.phi-probe.phi)))
                ROOT.massAllHist.Fill(mass,args.evtweight)
                ROOT.etaAllPairsHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)
                ROOT.phiAllPairsHist.Fill(probe.phi,args.evtweight)
                ROOT.pTAllPairsHist.Fill(probe.pt,args.evtweight)
                ROOT.pileupAllPairsHist.Fill(pileup,args.evtweight)

                ROOT.EtaAllPairsTagHist.Fill(tag.eta+tag.deltaEtaSC,args.evtweight)
                ROOT.EtaAllPairsProbeHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)

                # check tag passes the HEEP ID V7.0-2018Prompt
                if not (tag.cutBased_HEEP == 1):
                    continue
                ROOT.massIDSelHist.Fill(mass,args.evtweight) 
                ROOT.etaIDSelPairsHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)
                ROOT.phiIDSelPairsHist.Fill(probe.phi,args.evtweight)
                ROOT.pTIDSelPairsHist.Fill(probe.pt,args.evtweight)
                #ROOT.pileupIDSelHist(pileup,args.evtweight)

                ROOT.EtaIDSelTagHist.Fill(tag.eta+tag.deltaEtaSC,args.evtweight)
                ROOT.EtaIDSelProbeHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)

                tag_trig_match = trigger_match(tag,trig_objs)
                mass_pass_pair = True if (mass > 70 and mass < 110) else False
                
                # check tag matched to HLT_Ele32_WPTight trigger and invariant mass of pair is between 70 and 110 GeV/c^2
                if not (tag_trig_match and mass_pass_pair):
                    continue
                ROOT.massTrigHist.Fill(mass,args.evtweight) 

                ROOT.EtaTrigTagHist.Fill(tag.eta+tag.deltaEtaSC,args.evtweight)
                ROOT.EtaTrigProbeHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)

                if (probe.pt > 35):
                    ROOT.etaAllProbesHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)
                    ROOT.phiAllProbesHist.Fill(probe.phi,args.evtweight)
                    #ROOT.pileupAllProbesHist(pileup,args.evtweight)
                ROOT.pTAllProbesHist.Fill(probe.pt,args.evtweight)
                

                if not (probe.cutBased_HEEP == 1):
                    ROOT.massFailingProbesHist.Fill(mass,args.evtweight)
                    continue
                ROOT.etaPassingProbesHist.Fill(probe.eta+probe.deltaEtaSC,args.evtweight)
                if (probe.pt > 35):
                    ROOT.pTPassingProbesHist.Fill(probe.pt,args.evtweight)
                    ROOT.phiPassingProbesHist.Fill(probe.phi,args.evtweight)
                    #ROOT.pileupPassingProbesHist(pileup,args.evtweight)
                ROOT.massPassingProbesHist.Fill(mass,args.evtweight)
                
                """
                id_sf_weight = 1.
                reco_sf_weight = 1.
                #now do ID scale factors
                if idsf_hist:
                    tag_idsf = get_bincontent_2D(idsf_hist,tag.eta+tag.deltaEtaSC,tag.pt)
                    probe_idsf = get_bincontent_2D(idsf_hist,probe.eta+probe.deltaEtaSC,probe.pt)
                    id_sf_weight = tag_idsf*probe_idsf
                if recosf_hist:
                    tag_recosf = get_bincontent_2D(recosf_hist,tag.eta+tag.deltaEtaSC,tag.pt)
                    probe_recosf = get_bincontent_2D(recosf_hist,probe.eta+probe.deltaEtaSC,probe.pt)
                    reco_sf_weight = tag_recosf*probe_recosf

                ROOT.massIDSFCorrHist.Fill(mass,args.evtweight*id_sf_weight*reco_sf_weight)

                #now energy corrections are already applied
                ROOT.massESSCorrHist.Fill(mass,args.evtweight*id_sf_weight)
                """

        #make efficiency plots - only works for a single thread, use HEEPplotter for multithreading
        ROOT.pTHEEPEfficiencyHist = ROOT.pTPassingProbesHist.Clone("pTHEEPEfficiencyHist")
        ROOT.pTHEEPEfficiencyHist.Sumw2()
        ROOT.pTHEEPEfficiencyHist.Divide(ROOT.pTPassingProbesHist,ROOT.pTAllProbesHist,1,1,"")
        
        ROOT.etaHEEPEfficiencyHist = ROOT.etaPassingProbesHist.Clone("etaHEEPEfficiencyHist")
        ROOT.etaHEEPEfficiencyHist.Sumw2()
        ROOT.etaHEEPEfficiencyHist.Divide(ROOT.etaPassingProbesHist,ROOT.etaAllProbesHist,1,1,"")

        ROOT.phiHEEPEfficiencyHist = ROOT.phiPassingProbesHist.Clone("phiHEEPEfficiencyHist")
        ROOT.phiHEEPEfficiencyHist.Sumw2()
        ROOT.phiHEEPEfficiencyHist.Divide(ROOT.phiPassingProbesHist,ROOT.phiAllProbesHist,1,1,"")
        

    out_file.Write()
