/**
 * @file   TwoBodyDecayGen.hxx
 * @author Suvayu Ali <Suvayu.Ali@cern.ch>
 * @date   Wed Oct 10 22:39:13 2012
 *
 * @brief  This class implements an interface to generate 2-body decays
 *
 *
 */

#ifndef TWOBODYDECAYGEN_HXX
#define TWOBODYDECAYGEN_HXX

#include <string>
#include <vector>

#include <TTree.h>
#include <TGenPhaseSpace.h>

#define NDAUS 2


class TwoBodyDecayGen {
public:

  /**
   *
   *
   */

  TwoBodyDecayGen(double *daumasses, TwoBodyDecayGen *dau1=NULL,
		  TwoBodyDecayGen *dau2=NULL) :
    _generator(TGenPhaseSpace()), _daumasses(daumasses)
  { _daus[0] = dau1; _daus[1] = dau2; }
    // c++11 only, compile with -std=c++11 or -std=gnu++11
    // _daus{dau1, dau2} {}
  ~TwoBodyDecayGen() {}

  double generate(TLorentzVector &momp,
		  std::vector<TLorentzVector*> &particle_lvs);
  TTree* get_event_tree(unsigned nevents, TTree *tmomp=NULL);

private:

  /**
   *
   *
   */

  TGenPhaseSpace _generator;
  double *_daumasses;
  TwoBodyDecayGen *_daus[NDAUS];
};

#endif	// TWOBODYDECAYGEN_HXX
