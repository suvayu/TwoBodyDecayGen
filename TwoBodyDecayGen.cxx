/**
 * @file   TwoBodyDecayGen.cxx
 * @author Suvayu Ali <Suvayu.Ali@cern.ch>
 * @date   Thu Oct 11 13:45:20 2012
 *
 * @brief  Implementation of TwoBodyDecayGen
 *
 *
 */

#include <iostream>

#include <TRandom3.h>

#include "TwoBodyDecayGen.hxx"


double TwoBodyDecayGen::generate(TLorentzVector &momp,
				 std::vector<TLorentzVector> &particle_lvs)
{
  // setup decay and generate
  _generator.SetDecay(momp, NDAUS, _daumasses);
  double evt_wt = _generator.Generate();

  // retrieve decays
  for (unsigned j = 0; j < NDAUS; ++j) {
    particle_lvs.push_back(*(_generator.GetDecay(j)));
  }

  for (unsigned j = 0; j < NDAUS; ++j) {
    if (_daus[j]) {
      // FIXME: Assumption: particle_lvs is of size 2
      evt_wt += _daus[j]->generate(particle_lvs[j+1], particle_lvs);
      evt_wt /= 2.0;
    } // FIXME: the handling of weights is probably wrong
  }

  return evt_wt;
}


TTree* TwoBodyDecayGen::get_event_tree(unsigned nevents, TH1 *hmomp)
{
  std::vector<TLorentzVector> particle_lvs;
  double evt_wt(1.0);

  TTree *decaytree =
    new TTree("TwoBodyDecayGen_decaytree", "Vector of decay product "
	      "TLorentzVectors");
  decaytree->Branch("particle_lvs", &particle_lvs);
  decaytree->Branch("evt_wt", &evt_wt, "evt_wt/D");

  // Reset gRandom to TRandom3
  gRandom = new TRandom3();
  std::cout << "Generating " << nevents << " events." << std::endl;

  TLorentzVector momp(0.0, 0.0, 4.0, _mommass);
  for (unsigned i = 0; i < nevents; ++i) {
    particle_lvs.clear();

    // generate event and fill tree
    momp.SetXYZM( 0.0, 0.0, hmomp->GetRandom(), 5.367);
    particle_lvs.push_back(momp);
    evt_wt = this->generate(momp, particle_lvs);
    decaytree->Fill();
  }

  return decaytree;
}
