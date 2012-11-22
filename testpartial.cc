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


int kfactorp(TTree *outtree, std::string mode);


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

  fname = mode + "_Bs_mom_MC.png";
  intree->Draw("1E-3*BsMom.P()>>Bsmomp");
  label->DrawLatex(0.6, 0.5, "B_{s} momentum (MC)");
  gPad->Update();
  gPad->Print(fname.c_str());


  // read generated 4-vectors
  fname = "eventtree-" + mode + ".root";
  TFile outfile(fname.c_str(), "read");
  TTree *outtree = dynamic_cast<TTree*>(outfile.Get("TwoBodyDecayGen_decaytree"));

  kfactorp(outtree, mode);

  return 0;
}


int kfactorp(TTree *outtree, std::string mode)
{
  std::string fname;
  std::vector<TLorentzVector> *particle_lvs = NULL;
  TLorentzVector BsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector DsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector hMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector Bs_rec(0.0, 0.0, 0.0, 0.0);

  outtree->SetBranchAddress("particle_lvs", &particle_lvs);

  long oentries(outtree->GetEntries());

  TH1D hkfactorm("hkfactorm", "", 100, 0., 1.15);
  TH1D hkfactorp("hkfactorp", "", 100, 0., 1.15);
  TH1D hkfactorpm("hkfactorpm", "", 100, 0., 1.15);

  double kfactorp(0.0), kfactorm(0.0), kfactorpm(0.0);

  for (unsigned i = 0; i < oentries; ++i) {
    outtree->GetEntry(i);
    if (particle_lvs->size() < 3) continue;
    (*particle_lvs)[2].SetVectM((*particle_lvs)[2].Vect(), KMASS*1E-3);
    Bs_rec = (*particle_lvs)[3] + (*particle_lvs)[2];
    
    kfactorp = Bs_rec.P() / (*particle_lvs)[0].P();
    kfactorm = (*particle_lvs)[0].M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    hkfactorm.Fill(kfactorm);
    hkfactorp.Fill(kfactorp);
    hkfactorpm.Fill(kfactorpm);
  }

  hkfactorm.SetLineColor(kAzure);
  hkfactorp.SetLineColor(kAzure);
  hkfactorpm.SetLineColor(kAzure);

  TLatex *label = new TLatex();
  label->SetTextSize(0.04);
  label->SetNDC(true);

  hkfactorm.Draw("hist");
  label->DrawLatex(0.4, 0.5, "k-factor (m)");
  fname = mode + "_kfactorm_gen.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  hkfactorp.Draw("hist");
  label->DrawLatex(0.4, 0.5, "k-factor (p)");
  fname = mode + "_kfactorp_gen.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  hkfactorpm.Draw("hist");
  label->DrawLatex(0.4, 0.5, "k-factor (m/p)");
  fname = mode + "_kfactorpm_gen.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  fname = "hdump-" + mode + ".root";
  TFile hdump(fname.c_str(), "update");
  hdump.WriteTObject(&hkfactorm, NULL, "WriteDelete");
  hdump.WriteTObject(&hkfactorp, NULL, "WriteDelete");
  hdump.WriteTObject(&hkfactorpm, NULL, "WriteDelete");
  hdump.Close();

  return 0;
}
