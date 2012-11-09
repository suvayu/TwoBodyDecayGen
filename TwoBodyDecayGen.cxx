/**
 * @file   TwoBodyDecayGen.cxx
 * @author Suvayu Ali <Suvayu.Ali@cernNOSPAM.ch>
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


/**
 * \def DEBUG(COUNT, MSG)
 * Debug statement with a counter
 */

#define DEBUG(COUNT, MSG)                                  \
  std::cout << "DEBUG: [" << COUNT << "] (" << __func__ << ") " \
  << MSG << std::endl; \
  COUNT++;


/**
 * \def WARNING(COUNT, MSG)
 * Warning with a counter
 */

#define WARNING(COUNT, MSG)                                  \
  std::cout << "WARNING: [" << COUNT << "] (" << __func__ << ") " \
  << MSG << std::endl; \
  COUNT++;


/**
 * \def ERROR(COUNT, MSG)
 * Error message with a counter
 */

#define ERROR(COUNT, MSG)                                  \
  std::cout << "ERROR: [" << COUNT << "] (" << __func__ << ") " \
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
  // c++11 only, compile with -std=c++11 or -std=gnu++11
  // _daumasses{dau1, dau2} {}
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
  this->add_decay_channel(masses, nparts);
}


bool TwoBodyDecayGen::add_decay_channel(double *masses, unsigned nparts,
					double brfr)
{
  unsigned msgcount(0);

  if (masses[0] - _mommass > 1E-4) {
    ERROR(msgcount, "Mass of the mothers do not match!"
	  " Skipping new decay channel.");
    return false;
  }
  if (nparts > 7) {
    WARNING(msgcount, "Greater than two levels of decay is not tested."
	    " Expect the unexpected!");
  }

  unsigned nodes = (nparts - 1) / 2;
  for (unsigned i = 0; i < nodes; ++i) {
    if (masses[i] < 0.0) continue;
    double daumasses[NDAUS] = {masses[2*i + 1], masses[2*i + 2]};
    std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);

    DEBUG(msgcount, "mom: " << masses[i] << " dau[" << 2*i+1 << ","
	  << 2*i+2 <<  "]: (" << daumasses[0] << "," << daumasses[1] << ")");

    if (0 == i) {
      _daumasses[0] = daumasses[0];
      _daumasses[1] = daumasses[1];
    } else {
      daus[i-1] = new TwoBodyDecayGen(masses[i], daumasses[0], daumasses[1]);
    }

    if (_dauchannels.empty()) { // Ignore provided B.F. when there are
      brfr = 1.0;		// no primary channels
    } else {
      _dauchannels[0].second -= brfr; // Correct primary channel B.F.
    }
    DauNode channel(daus, brfr);
    _dauchannels.push_back(channel);
  }

  return true;
}


double TwoBodyDecayGen::generate(TLorentzVector &momp,
				 std::vector<TLorentzVector> &particle_lvs,
				 int ich)
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

  // propagate generate to daughters
  for (unsigned j = 0; j < NDAUS; ++j) {
    if (ich < 0) {		// For daughter of daughters
      for (unsigned i = 0; i < _dauchannels.size(); ++i) {
	if (_dauchannels[i].first[j]) {
	  evt_wt += _dauchannels[i].first[j]->generate(particle_lvs[j+1],
						       particle_lvs, -1);
	  evt_wt /= 2.0;
	}
      }
    } else {
      if (_dauchannels[ich].first[j]) {
	evt_wt += _dauchannels[ich].first[j]->generate(particle_lvs[j+1],
						       particle_lvs, -1);
	evt_wt /= 2.0;
      }
    }
  } // FIXME: the handling of weights is probably wrong

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


void TwoBodyDecayGen::print(unsigned indent) {
  std::string prefix(indent * 2, ' ');
  std::cout << prefix << "mommass: " << _mommass << ", daumass: ("
	    << _daumasses[0] << "," << _daumasses[1] << ")"
	    << std::endl;

  BOOST_FOREACH(DauNode node, _dauchannels) {
    for (unsigned j = 0; j < NDAUS; ++j) {
      if (node.first[j]) {
	std::cout << "Dau " << j << ":" << std::endl;
	node.first[j]->print(indent + 1);
      }
    }
  }
  return;
}
