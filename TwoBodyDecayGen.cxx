/**
 * @file   TwoBodyDecayGen.cxx
 * @author Suvayu Ali <Suvayu.Ali@cern.ch>
 * @date   Thu Oct 11 13:45:20 2012
 * 
 * @brief  Implementation of TwoBodyDecayGen
 * 
 * 
 */

#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "TwoBodyDecayGen.hxx"


bool TwoBodyDecayGen::generate(unsigned nevents)
{
  TClonesArray dau_lvs(TLorentzVector::Class());
  double evt_wt(1.0);
  _decaytree->Branch("dau_lvs", &dau_lvs);
  _decaytree->Branch("evt_wt", &evt_wt, "evt_wt/D");

  // set generator
  TLorentzVector initial_lv(0.0, 0.0, 4.0, 5.367);
  _generator.SetDecay(initial_lv, NDAUS, _daumasses);

  for (unsigned i = 0; i < nevents; ++i) {
    dau_lvs.Clear();
    evt_wt = _generator.Generate();
    for (unsigned j = 0; j < NDAUS; ++j) {
      TLorentzVector *temp = _generator.GetDecay(j);
      new(dau_lvs[j]) TLorentzVector(*temp);
    }
    _decaytree->Fill();
  }

  return true;
}
