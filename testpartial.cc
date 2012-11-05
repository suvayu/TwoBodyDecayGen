#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TLegend.h>


// some constants
static const double BSMASS(5366.3), DSMASS(1968.49), KMASS(493.677),
  BDMASS(5279.53), DMASS(1869.62), PIMASS(139.57018), DSSTMASS(2112.34),
  KSTMASS(891.66), LBMASS(5620.2), LCMASS(2286.46), PMASS(938.27203);
  /*, RHOMASS(775.49), DSTMASS(2010.25), DSTMASS2(2460.1);*/


int kfactorp(TTree *intree, TTree *outtree);


int main(int argc, char* argv[])
{
  // program arguments
  std::string mode;
  if (argc == 2) {
    mode = argv[1];
  } else {
    std::cout << "Not enough arguments!" << std::endl;
    return -1;
  }

  // read ntuple from file
  std::string fname = "smalltree-" + mode + ".root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("ftree"));

  TLatex *label = new TLatex();
  label->SetTextSize(0.04);
  label->SetNDC(true);

  // B momentum distribution from MC
  TH1D Bsmomp("Bsmomp", "", 100, 0.0, 300.0);
  TH1D dau1momp("dau1momp", "", 100, 0.0, 300.0);
  TH1D dau2momp("dau2momp", "", 100, 0.0, 300.0);

  intree->Draw("1E-3*BsMom.P()>>Bsmomp");
  label->DrawLatex(0.6, 0.5, "B_{s} momentum (MC)");
  gPad->Update();
  gPad->Print("Bs_mom_MC.png");


  // read generated 4-vectors
  fname = "eventtree-" + mode + ".root";
  TFile outfile(fname.c_str(), "read");
  TTree *outtree = dynamic_cast<TTree*>(outfile.Get("TwoBodyDecayGen_decaytree"));

  kfactorp(intree, outtree);

  return 0;
}


int kfactorp(TTree *intree, TTree *outtree)
{
  std::vector<TLorentzVector> *particle_lvs = NULL;
  TLorentzVector BsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector DsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector hMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector Bs_rec(0.0, 0.0, 0.0, 0.0);

  outtree->SetBranchAddress("particle_lvs", &particle_lvs);

  long oentries(outtree->GetEntries());

  TH1D kfactor("kfactor", "", 100, 0., 1.15);

  double kfactorp(0.0), kfactorm(0.0), kfactorpm(0.0);

  for (unsigned i = 0; i < oentries; ++i) {
    outtree->GetEntry(i);
    if (particle_lvs->size() < 3) continue;
    (*particle_lvs)[2].SetVectM((*particle_lvs)[2].Vect(), KMASS*1E-3);
    Bs_rec = (*particle_lvs)[3] + (*particle_lvs)[2];
    
    kfactorp = Bs_rec.P() / (*particle_lvs)[0].P();
    kfactorm = (*particle_lvs)[0].M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactor.Fill(kfactorpm);
  }

  kfactor.SetLineColor(kAzure);
  kfactor.Draw("hist");

  TLatex *label = new TLatex();
  label->DrawLatex(0.6, 0.5, "k-factor (gen)");

  gPad->Update();
  gPad->Print("kfactor_DsstPi_gen.png");

  TFile hdump("hdump-DsstPi.root", "update");
  hdump.WriteTObject(&kfactor, NULL, "WriteDelete");
  hdump.Close();

  return 0;
}
