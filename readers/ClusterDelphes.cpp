#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <stdlib.h>
#include <functional>
#include <time.h>
#include <math.h>
#include <numeric>

#include "TROOT.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

static int NMAX = 150;

//---------------------------------------------------------------------------


struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float charge = 1;
  float hardfrac = 1;  
  float vtxid = -1;
  float pdgid = 0;
};

template <typename T>
void
fill(vector<float> &vattr, vector<PFCand> &particles, T fn_attr)
{
  vattr.clear();
  for (auto& p : particles)
    vattr.push_back(fn_attr(p));
}


//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  srand(time(NULL));

  if(argc < 3) {
    cout << " Usage: " << "PapuDelphes" << " input_file"
         << " output_file" << endl;
    cout << " input_file - input file in ROOT format," << endl;
    cout << " output_file - output file in ROOT format" << endl;
    return 1;
  }

  // figure out how to read the file here 
  //

  TFile* ifile = TFile::Open(argv[1], "READ");
  TTree* itree = (TTree*)ifile->Get("Delphes;1");

  auto* fout = TFile::Open(argv[2], "RECREATE");
  auto* tout = new TTree("events", "events");

  unsigned int nevt = itree->GetEntries();
  TBranch* partbranch = (TBranch*)itree->GetBranch("Particle");
  TBranch* pfbranch = (TBranch*)itree->GetBranch("ParticleFlowCandidate");
  TBranch* genjetbranch = (TBranch*)itree->GetBranch("GenJet");
  TBranch* pumbranch = (TBranch*)itree->GetBranch("PileUpMix");
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;
  vector<TLorentzVector> pum_particles;
  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac, vvtxid;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpuppi.reserve(NMAX); vpdgid.reserve(NMAX); 
  vhardfrac.reserve(NMAX); vvtxid.reserve(NMAX);
  
  float pumjetpt=-99., pumjeteta=-99., pumjetphi=-99.,pumjetm=-99.;
  float genjetpt=-99., genjeteta=-99., genjetphi=-99., genjetm=-99.;
  float puppijetpt=-99., puppijeteta=-99., puppijetphi=-99., puppijetm=-99.;
  float truthjetpt=-99., truthjeteta=-99., truthjetphi=-99., truthjetm=-99.;
  float pfjetpt=-99., pfjeteta=-99., pfjetphi=-99., pfjetm=-99.;

  // jet branches
  TBranch* b_genjetpt = tout->Branch("genjetpt",&genjetpt, "genjetpt/F");
  TBranch* b_genjeteta = tout->Branch("genjeteta",&genjeteta, "genjeteta/F");
  TBranch* b_genjetphi = tout->Branch("genjetphi",&genjetphi, "genjetphi/F");
  TBranch* b_genjetm = tout->Branch("genjetm",&genjetm, "genjetm/F");
  TBranch* b_puppijetpt = tout->Branch("puppijetpt",&puppijetpt, "puppijetpt/F");
  TBranch* b_puppijeteta = tout->Branch("puppijeteta",&puppijeteta, "puppijeteta/F");
  TBranch* b_puppijetphi = tout->Branch("puppijetphi",&puppijetphi, "puppijetphi/F");
  TBranch* b_puppijetm = tout->Branch("puppijetm",&puppijetm, "puppijetm/F");
  TBranch* b_truthjetpt = tout->Branch("truthjetpt",&truthjetpt, "truthjetpt/F");
  TBranch* b_truthjeteta = tout->Branch("truthjeteta",&truthjeteta, "truthjeteta/F");
  TBranch* b_truthjetphi = tout->Branch("truthjetphi",&truthjetphi, "truthjetphi/F");
  TBranch* b_truthjetm = tout->Branch("truthjetm",&truthjetm, "truthjetm/F");
  TBranch* b_pfjetpt = tout->Branch("pfjetpt",&pfjetpt, "pfjetpt/F");
  TBranch* b_pfjeteta = tout->Branch("pfjeteta",&pfjeteta, "pfjeteta/F");
  TBranch* b_pfjetphi = tout->Branch("pfjetphi",&pfjetphi, "pfjetphi/F");
  TBranch* b_pfjetm = tout->Branch("pfjetm",&pfjetm, "pfjetm/F");
  // PF cand branches
  tout->Branch("pt", &vpt);
  tout->Branch("eta", &veta);
  tout->Branch("phi", &vphi);
  tout->Branch("e", &ve);
  tout->Branch("puppi", &vpuppi);
  tout->Branch("pdgid", &vpdgid);
  tout->Branch("hardfrac", &vhardfrac);
  tout->Branch("vtxid", &vvtxid);


  ExRootProgressBar progressBar(nevt);
  
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);

  // Soft drop -- only needed for large radius jets
  double radius = 0.4;
  double sdZcut = 0.1;
  double sdBeta = 0.;
  
  fastjet::contrib::SoftDrop softDrop = fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);
  int counter_yeet = 0;
  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);
    pfjetpt = -99;
    puppijetpt = -99;
    truthjetpt = -99;
    if (k%100==0)
      std::cout << k << " / " << nevt << std::endl;

    pum_particles.clear();
    unsigned int npumjets = pumbranch->GetEntries();
    npumjets = itree->GetLeaf("PileUpMix_size")->GetValue(0);
    for (unsigned int j=0; j<npumjets; j++){
      TLorentzVector pumjet;
      pumjet.SetPtEtaPhiM(itree->GetLeaf("PileUpMix.PT")->GetValue(j),itree->GetLeaf("PileUpMix.Eta")->GetValue(j),itree->GetLeaf("PileUpMix.Phi")->GetValue(j),itree->GetLeaf("PileUpMix.Mass")->GetValue(j));
      pumjetpt = pumjet.Pt();
      pumjeteta = pumjet.Eta();
      pumjetphi = pumjet.Phi();
      pumjetm = pumjet.M();
      if ((abs(itree->GetLeaf("PileUpMix.PID")->GetValue(j))==11||abs(itree->GetLeaf("PileUpMix.PID")->GetValue(j))==11) && pumjetpt>10){
	pum_particles.push_back(pumjet);
      }
      else{
	continue;
      }
    }

    unsigned int ngenjets = genjetbranch->GetEntries();
    ngenjets = itree->GetLeaf("GenJet_size")->GetValue(0);
    TLorentzVector genjet;
    for (unsigned int j=0; j<ngenjets; j++){
      int c = 0;
      genjet.SetPtEtaPhiM(itree->GetLeaf("GenJet.PT")->GetValue(j),itree->GetLeaf("GenJet.Eta")->GetValue(j),itree->GetLeaf("GenJet.Phi")->GetValue(j),itree->GetLeaf("GenJet.Mass")->GetValue(j));
      genjetpt = genjet.Pt();
      genjeteta = genjet.Eta();
      genjetphi = genjet.Phi();
      genjetm = genjet.M();
      for (unsigned int jj = 0; jj<pum_particles.size(); jj++){
	if(pum_particles[jj].DeltaR(genjet)<.4){
	  c = 1;
	  break;
	}
	else{
	  continue;
	}
      }
      if (c == 1){
	continue;
      }
      else{
	break;// only get leading gen jet
      }
    }

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);
    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.puppi = itree->GetLeaf("ParticleFlowCandidate.PuppiW")->GetValue(j);
      tmppf.hardfrac = itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j);
      tmppf.pdgid = itree->GetLeaf("ParticleFlowCandidate.PID")->GetValue(j);
      if (itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j)!=0){
        if (itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j)==1)
          tmppf.vtxid = 0;
        else
          tmppf.vtxid = 1;
      }
      else
        tmppf.vtxid = -1;
      if((abs(tmppf.pdgid)==11||abs(tmppf.pdgid)==13) && tmppf.pt > 10){
	continue;
      }
      else{
	input_particles.push_back(tmppf);
      }
    }
    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);

    vector<fastjet::PseudoJet> finalStates_puppi;
    vector<fastjet::PseudoJet> finalStates_truth;
    vector<fastjet::PseudoJet> finalStates_pf;
    int pfid = 0;
    for(auto &p : input_particles){
      if (p.vtxid==1){
	pfid+=1;
	continue; // Charged Hadron Subtraction
      }
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);

      // PF
      fastjet::PseudoJet curjet_pf(tmp.Px(), tmp.Py(), tmp.Pz(), tmp.E());
      curjet_pf.set_user_index(pfid);    
      finalStates_pf.emplace_back(curjet_pf);
      
      // Puppi
      if (p.puppi>0.){
	fastjet::PseudoJet curjet_puppi(p.puppi*tmp.Px(), p.puppi*tmp.Py(), p.puppi*tmp.Pz(), p.puppi*tmp.E());
	curjet_puppi.set_user_index(pfid);
	finalStates_puppi.emplace_back(curjet_puppi);
      }

      // Truth
      if (p.hardfrac>0.){
	fastjet::PseudoJet curjet_truth(p.hardfrac*tmp.Px(), p.hardfrac*tmp.Py(), p.hardfrac*tmp.Pz(), p.hardfrac*tmp.E());
	curjet_truth.set_user_index(pfid);
	finalStates_truth.emplace_back(curjet_truth);
      }
        pfid += 1;
    }


    fastjet::ClusterSequence seq_puppi(sorted_by_pt(finalStates_puppi), *jetDef);
    fastjet::ClusterSequence seq_truth(sorted_by_pt(finalStates_truth), *jetDef);
    fastjet::ClusterSequence seq_pf(sorted_by_pt(finalStates_pf), *jetDef);

    float minpt = 15;

    vector<fastjet::PseudoJet> allJets_puppi(sorted_by_pt(seq_puppi.inclusive_jets(minpt)));
    vector<fastjet::PseudoJet> allJets_truth(sorted_by_pt(seq_truth.inclusive_jets(minpt)));
    vector<fastjet::PseudoJet> allJets_pf(sorted_by_pt(seq_pf.inclusive_jets(minpt)));
    for (auto& puppiJet : allJets_puppi) {
      if (puppiJet.perp() < minpt)
	break;

      //fastjet::PseudoJet sdJet = (softDrop)(puppiJet);

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(puppiJet.perp(),puppiJet.eta(),puppiJet.phi(),puppiJet.m());
      if (tmp.DeltaR(genjet)<0.4){
	puppijetpt = tmp.Pt();
	puppijeteta = tmp.Eta();
	puppijetphi = tmp.Phi();
	puppijetm = tmp.M();
	break;
      }
    }

    for (auto& truthJet : allJets_truth) {
      if (truthJet.perp() < minpt)
	break;

      //fastjet::PseudoJet sdJet = (softDrop)(truthJet);

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(truthJet.perp(),truthJet.eta(),truthJet.phi(),truthJet.m());
      if (tmp.DeltaR(genjet)<0.4){
	truthjetpt = tmp.Pt();
	truthjeteta = tmp.Eta();
	truthjetphi = tmp.Phi();
	truthjetm = tmp.M();
	break;
      }
    }
    //    cout<<k<<endl;
    counter_yeet = 0;
    // cout<<k<<endl;
    for (auto& pfJet : allJets_pf) {
      if (pfJet.perp() < minpt)
	break;

      //fastjet::PseudoJet sdJet = (softDrop)(puppiJet);

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(pfJet.perp(),pfJet.eta(),pfJet.phi(),pfJet.m());
      if (tmp.DeltaR(genjet)<0.4){
	pfjetpt = tmp.Pt();
	pfjeteta = tmp.Eta();
	pfjetphi = tmp.Phi();
	pfjetm = tmp.M();
	break;
      }
      counter_yeet++;
    }
    if(pfjetpt < 0 || puppijetpt < 0 || truthjetpt < 0)
      continue;
    //cout<<counter_yeet<<endl;
    //cout<<allJets_pf.size()<<endl;
    // cout<<k<<endl;
    output_particles.clear();
    // cout<<k<<endl;
    std::vector< PseudoJet > constituents = allJets_pf[counter_yeet].constituents();
    // cout<<k<<endl;
    TLorentzVector jet(0,0,0,0);
    //cout<<k<<endl;
    for(unsigned int j = 0; j<constituents.size(); j++){
      // TLorentzVector tmp;
      //tmp.SetPtEtaPhiM(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
      //jet += tmp;
      output_particles.push_back( input_particles[constituents[j].user_index()]);
    }

    output_particles.resize(NMAX);
    fill(vpt, output_particles, [](PFCand& p) { return p.pt; });
    fill(veta, output_particles, [](PFCand& p) { return p.eta; });
    fill(vphi, output_particles, [](PFCand& p) { return p.phi; });
    fill(ve, output_particles, [](PFCand& p) { return p.e; });
    fill(vpuppi, output_particles, [](PFCand& p) { return p.puppi; });
    fill(vpdgid, output_particles, [](PFCand& p) { return p.pdgid; });
    fill(vhardfrac, output_particles, [](PFCand& p) { return p.hardfrac; });
    fill(vvtxid, output_particles, [](PFCand& p) { return p.vtxid; }); 
    tout->Fill();

    progressBar.Update(k, k);
  }

  fout->Write();
  fout->Close();

}
