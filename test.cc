#include <iostream>
#include <cstdlib>
#include <cassert>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TPad.h>


int main()
{
  std::string fname = "eventtree.root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("TwoBodyDecayGen_decaytree"));

  TH1D Bsmomp("Bsmomp", "", 100, 0.0, 300.0);
  intree->Draw("particle_lvs[0].P()>>Bsmomp");
  gPad->Print("Bs_mom_sav.png");

  TH1D dau1momp("dau1momp", "", 100, 0.0, 300.0);
  intree->Draw("particle_lvs[1].P()>>dau1momp");
  gPad->Print("dau1_mom_sav.png");

  TH1D dau2momp("dau2momp", "", 100, 0.0, 300.0);
  intree->Draw("particle_lvs[2].P()>>dau2momp");
  gPad->Print("dau2_mom_sav.png");

  return 0;
}
