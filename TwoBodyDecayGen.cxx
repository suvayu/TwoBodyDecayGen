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
#include <iomanip>

#include <boost/foreach.hpp>

#include <TRandom3.h>

#include "TwoBodyDecayGen.hxx"


/**
 * \def DEBUG(MSG)
 * Debug statement with a counter
 */

#define DEBUG(MSG)  \
  std::cout << "DEBUG: [" << std::setw(4) << std::setfill('0') << _count \
  << "] (" << __func__ << ") " << MSG << std::endl; \
  _count++;


/**
 * \def WARNING(MSG)
 * Warning with a counter
 */

#define WARNING(MSG)                                  \
  std::cout << "WARNING: [" << std::setw(4) << std::setfill('0') << _count \
  << "] (" << __func__ << ") " << MSG << std::endl; \
  _count++;


/**
 * \def ERROR(MSG)
 * Error message with a counter
 */

#define ERROR(MSG)                                  \
  std::cout << "ERROR: [" << std::setw(4) << std::setfill('0') << _count \
  << "] (" << __func__ << ") "	<< MSG << std::endl; \
  _count++;


void TwoBodyDecayGen::printQ(std::string prefix, std::deque<chBFpair> queue)
{
  DEBUG(prefix << "Q size: " << queue.size());
  BOOST_FOREACH(chBFpair pair, queue) {
    DEBUG(prefix << "(" << pair.first << ", " << pair.second << ")");
  }
  return;
}


void TwoBodyDecayGen::printQ(std::string prefix, std::vector<std::deque<chBFpair> > QVec)
{
  DEBUG(prefix << "V size: " << QVec.size());
  BOOST_FOREACH(std::deque<chBFpair> queue, QVec) {
    printQ(prefix, queue);
  }
  return;
}


unsigned long long TwoBodyDecayGen::_count(0);


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

  DauNode priChannel(daus, 1.0);
  _dauchannels.push_back(priChannel);
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
  _daumasses[0] = masses[1];
  _daumasses[1] = masses[2];

  if (nparts > 3) {
    this->add_decay_channel(masses, nparts, 1.0);
  } else {
    std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);
    DauNode priChannel(daus, 1.0);
    _dauchannels.push_back(priChannel);
  }
}


bool TwoBodyDecayGen::add_decay_channel(double *masses, unsigned nparts,
					double brfr)
{
  if ((std::fabs(masses[0] - _mommass) > 1E-4) or
      (std::fabs(masses[1] - _daumasses[0]) > 1E-4) or
      (std::fabs(masses[2] - _daumasses[1]) > 1E-4)) {
    ERROR("Mass of the mothers do not match!"
	  " Skipping new decay channel.");
    return false;
  }

  if (brfr - 1.0 > 0.0) {
    ERROR("Branching fraction cannot be > 1.0,"
	  " skipping new decay channel.");
    return false;
  }

  if (nparts > 7) {
    WARNING("Greater than two levels of decay is not tested."
	    " Expect the unexpected!");
  }

  std::vector<double> dau1tree, dau2tree;

  // FIXME: spurious ghost daughters for stable daughter
  const unsigned nodes = (nparts - 1) / 2;
  for (unsigned i = 1; i < nodes; ++i) {
    if (i % 2) {
      if (i == 1) dau1tree.push_back(masses[i]);
      dau1tree.push_back(masses[2*i + 1]);
      dau1tree.push_back(masses[2*i + 2]);
    } else {
      if (i == 2) dau2tree.push_back(masses[i]);
      dau2tree.push_back(masses[2*i + 1]);
      dau2tree.push_back(masses[2*i + 2]);
    }
  }

  std::vector<TwoBodyDecayGen*> daus(NDAUS, NULL);
  if (not dau1tree.empty()) {
    daus[0] = new TwoBodyDecayGen(&dau1tree[0], dau1tree.size());
  }
  if (not dau2tree.empty()) {
    daus[1] = new TwoBodyDecayGen(&dau2tree[0], dau2tree.size());
  }

  if (not _dauchannels.empty()) {
    _dauchannels[0].second -= brfr; // Correct primary channel B.F.
  }
  DauNode channel(daus, brfr);
  _dauchannels.push_back(channel);
  return true;
}


TwoBodyDecayGen* TwoBodyDecayGen::get_daughter(unsigned chid, unsigned dauid)
{
  if (dauid > 1) {
    ERROR("Only " << NDAUS << " daughters.  dauid cannot be greater than 1.");
    return NULL;
  }
  return _dauchannels[chid].first[dauid];
}


double TwoBodyDecayGen::get_brfr(unsigned chid)
{
  return _dauchannels[chid].second;
}


void TwoBodyDecayGen::find_leaf_nodes(std::vector<std::deque<chBFpair> > &brfrVec,
				      std::deque<chBFpair> &brfrQ)
{
  DEBUG("Passed Qptr: " << &brfrQ);
  printQ("Passed ", brfrQ);
  DEBUG("_dauchannels.size(): " << _dauchannels.size());
  std::deque<chBFpair> brfrQcopy(brfrQ);

  // loop over channels
  for (unsigned chid = 0; chid < _dauchannels.size(); ++chid) {
    // FIXME: copy not working, copying gibberish
    // need copy when both daughters are _not_ leaves
    brfrQ.push_back(std::make_pair(chid, this->get_brfr(chid))); // channel BF
    DEBUG("Ch id: " << chid);

    unsigned leafcounter(0);
    // loop over daughters for each channel
    for (unsigned dauid = 0; dauid < NDAUS; ++dauid) {
      TwoBodyDecayGen *dau = this->get_daughter(chid, dauid);

      // recursive calls
      if (dau) { // not a leaf node, propagate call
	dau->find_leaf_nodes(brfrVec, brfrQ);
      } else { // leaf branch
	++leafcounter;
	DEBUG("leafcounter: " << leafcounter);
      }
      if (leafcounter == 2) { // leaf decay node
	brfrVec.push_back(brfrQ);
	// FIXME: probably this does not work
	brfrQ = brfrQcopy;
	DEBUG("Leaf node found!");
      }
    } // end daughter loop
  } // end channel loop

  if (_dauchannels.empty()) { // leaf decay node
    brfrVec.push_back(brfrQ);
    DEBUG("Leaf node found!");
  }

  printQ("End ", brfrQ);

  DEBUG("Vector size: " << brfrVec.size());
  return;
}


double TwoBodyDecayGen::generate(TLorentzVector &momp,
				 std::vector<TLorentzVector> &particle_lvs,
				 std::deque<chBFpair> chQ)
{
  unsigned ich(chQ.front().first);
  chQ.pop_front();
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
    DEBUG("_dauchannels size: " << _dauchannels.size());
    if (_dauchannels.empty()) continue;
    if (_dauchannels[ich].first.empty()) {
      DEBUG("Empty vector!");
      continue;
    }
    DEBUG("pointer: " << _dauchannels[ich].first[j]);
    if (chQ.empty()) {
      DEBUG("Channel Q empty!");
      continue;
    } else {
      DEBUG("Channel Q size: " << chQ.size());
    }
    if (_dauchannels[ich].first[j]) {
      _dauchannels[ich].first[j]->print(4);
      evt_wt += _dauchannels[ich].first[j]->generate(particle_lvs[j+1],
						     particle_lvs, chQ);
      evt_wt /= 2.0;
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

  std::vector<std::deque<chBFpair> > brfrVec;
  std::deque<chBFpair> brfrQ;
  // brfrQ.push_back(std::make_pair(1, 0.81));
  DEBUG("Initial Qptr: " << &brfrQ);
  this->find_leaf_nodes(brfrVec, brfrQ);
  this->printQ("Final ", brfrVec);

  BOOST_FOREACH(std::deque<chBFpair> chQ, brfrVec) {
    double eff_brfr(1.0);
    unsigned eff_nevents(0);

    DEBUG("chQ size: " << chQ.size());
    BOOST_FOREACH(chBFpair ch, chQ) {
      eff_brfr *= ch.second;
      DEBUG("Ch id: " << ch.first << ", BF: " << ch.second);
    }

    eff_nevents = eff_brfr * nevents;
    DEBUG("Effective BF: " << eff_brfr << ", effective events: " << eff_nevents);
    for (unsigned i = 0; i < eff_nevents; ++i) {
      particle_lvs.clear();
      // generate event and fill tree
      momp.SetXYZM( 0.0, 0.0, hmomp->GetRandom(), _mommass);
      particle_lvs.push_back(momp);
      if (chQ.empty()) {
	DEBUG("Empty channel queue!");
	continue;
      }
      evt_wt = this->generate(momp, particle_lvs, chQ);
      if (evt_wt < 0) {
	std::cout << "Decay not permitted by kinematics, skipping!"
		  << std::endl;
	continue;
      }
      decaytree->Fill();
    } // end of loop over events per leaf branch/decay node
  }   // end of loop over leaves

  return decaytree;
}


void TwoBodyDecayGen::print(unsigned indent) {
  std::string prefix(indent * 2, ' ');
  std::cout << prefix << "mommass: " << _mommass << ", daumass: ("
	    << _daumasses[0] << "," << _daumasses[1] << ") with "
	    << _dauchannels.size() << " daughter channel(s)." << std::endl;

  BOOST_FOREACH(DauNode node, _dauchannels) {
    if (not node.first.empty()) {
      std::cout << prefix << "Channel BF: " << node.second << std::endl;
    }
    for (unsigned j = 0; j < NDAUS; ++j) {
      if (node.first[j]) {
	std::cout << prefix << "Node " << j << ":" << std::endl;
	node.first[j]->print(indent + 1);
      }
    }
  }
  return;
}
