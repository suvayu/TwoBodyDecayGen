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

#include <boost/foreach.hpp>

#include <TRandom3.h>

#include "TwoBodyDecayGen.hxx"


#define DEBUG(COUNT, MSG)                                  \
  std::cout << "SA-DEBUG: [" << COUNT << "] (" << __func__ << ") " \
  << MSG << std::endl; \
  COUNT++;


TwoBodyDecayGen::TwoBodyDecayGen(double mommass,
				 double dau1mass,
				 double dau2mass,
				 TwoBodyDecayGen *dau1,
				 TwoBodyDecayGen *dau2) :
  _generator(TGenPhaseSpace()), _mommass(mommass)
{
  _daumasses[0] = dau1mass;
  _daumasses[1] = dau2mass;

  std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);
  daus[0] = dau1;
  daus[1] = dau2;

  DauNode priNode(daus, 1.0);
  _dauchannels.push_back(priNode);
}


TwoBodyDecayGen::TwoBodyDecayGen(double mommass, double *daumasses,
				 TwoBodyDecayGen *dau1,
				 TwoBodyDecayGen *dau2) :
  _generator(TGenPhaseSpace()), _mommass(mommass)//, _daumasses(daumasses)
{
  _daumasses[0] = daumasses[0];
  _daumasses[1] = daumasses[1];

  std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);
  daus[0] = dau1;
  daus[1] = dau2;

  DauNode priChannel(daus, 1.0);
  _dauchannels.push_back(priChannel);
}


TwoBodyDecayGen::TwoBodyDecayGen(double *masses, unsigned nparts) :
  _generator(TGenPhaseSpace()), _mommass(masses[0])
{
  unsigned msgcount(0);

  if (nparts > 7) {
    std::cout << "Greater than two levels of decay is not tested. "
      " Expect the unexpected!" << std::endl;
  }
  unsigned nodes = (nparts - 1) / 2;
  for (unsigned i = 0; i < nodes; ++i) {
    if (masses[i] < 0.0) continue;
    // FIXME: gibberish is entered for last particle
    double daumasses[NDAUS] = {masses[2*i + 1], masses[2*i + 2]};
    std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);

    DEBUG(msgcount, "mom: " << masses[i] << " dau[" << 2*i+1 << ","
	  << 2*i+2 <<  "]: (" << daumasses[0] << "," << daumasses[1]
	  << ")");

    if (0 == i) {
      _daumasses[0] = daumasses[0];
      _daumasses[1] = daumasses[1];
    } else {
      daus[i-1] = new TwoBodyDecayGen(masses[i], daumasses[0], daumasses[1]);
    }

    DauNode priChannel(daus, 1.0);
    _dauchannels.push_back(priChannel);
  }
}


double TwoBodyDecayGen::generate(TLorentzVector &momp,
				 std::vector<TLorentzVector> &particle_lvs,
				 unsigned ich)
{
  // setup decay and generate
  if (not _generator.SetDecay(momp, NDAUS, _daumasses)) {
    return -1.0;
  }
  double evt_wt = _generator.Generate();

  // retrieve decays
  for (unsigned j = 0; j < NDAUS; ++j) {
    particle_lvs.push_back(*(_generator.GetDecay(j)));
  }

  for (unsigned j = 0; j < NDAUS; ++j) {
    if (_dauchannels[ich].first[j]) {
      // FIXME: Ignoring branching fraction for now
      evt_wt += _dauchannels[ich].first[j]->generate(particle_lvs[j+1],
						     particle_lvs);
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
  unsigned eff_nevents(0);

  for (unsigned j = 0; j < _dauchannels.size(); ++j) {
    eff_nevents = _dauchannels[j].second * nevents;

    for (unsigned i = 0; i < eff_nevents; ++i) {
      particle_lvs.clear();

      // generate event and fill tree
      momp.SetXYZM( 0.0, 0.0, hmomp->GetRandom(), 5.367);
      particle_lvs.push_back(momp);
      evt_wt = this->generate(momp, particle_lvs, j);
      if (evt_wt < 0) {
	std::cout << "Decay not permitted by kinematics, skipping!"
		  << std::endl;
	continue;
      }
      decaytree->Fill();
    }
  }

  return decaytree;
}


void TwoBodyDecayGen::Print(unsigned indent) {
  std::string prefix(indent * 2, ' ');
  std::cout << prefix << "mommass: " << _mommass << ", daumass: ("
	    << _daumasses[0] << "," << _daumasses[1] << ")"
	    << std::endl;

  BOOST_FOREACH(DauNode node, _dauchannels) {
    for (unsigned j = 0; j < NDAUS; ++j) {
      if (node.first[j]) {
	std::cout << "Dau " << j << ":" << std::endl;
	node.first[j]->Print(indent + 1);
      }
    }
  }
  return;
}
