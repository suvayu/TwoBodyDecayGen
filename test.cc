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

  intree->Draw("1E-3*DsMom.P()>>dau1momp");
  label->DrawLatex(0.6, 0.5, "dau1 momentum (MC)");
  gPad->Update();
  gPad->Print("dau1_mom_MC.png");

  intree->Draw("1E-3*hMom.P()>>dau2momp");
  label->DrawLatex(0.6, 0.5, "dau2 momentum (MC)");
  gPad->Update();
  gPad->Print("dau2_mom_MC.png");


  // read generated 4-vectors
  fname = "eventtree-" + mode + ".root";
  TFile outfile(fname.c_str(), "read");
  TTree *outtree = dynamic_cast<TTree*>(outfile.Get("TwoBodyDecayGen_decaytree"));

  // generated B momentum distribution
  outtree->Draw("particle_lvs[0].P()>>Bsmomp");
  label->DrawLatex(0.6, 0.5, "B_{s} momentum (gen)");
  gPad->Update();
  gPad->Print("Bs_mom_gen.png");

  outtree->Draw("particle_lvs[1].P()>>dau1momp");
  label->DrawLatex(0.6, 0.5, "dau1 momentum (gen)");
  gPad->Update();
  gPad->Print("dau1_mom_gen.png");

  outtree->Draw("particle_lvs[2].P()>>dau2momp");
  label->DrawLatex(0.6, 0.5, "dau2 momentum (gen)");
  gPad->Update();
  gPad->Print("dau2_mom_gen.png");

  kfactorp(intree, outtree);

  return 0;
}


int kfactorp(TTree *intree, TTree *outtree)
{
  TLorentzVector *tru_BsMom = NULL;
  TLorentzVector *tru_hMom = NULL;
  TLorentzVector *tru_DsMom = NULL;
  TLorentzVector Bs_rec(0.0, 0.0, 0.0, 0.0);

  intree->SetBranchAddress("tru_BsMom", &tru_BsMom);
  intree->SetBranchAddress("tru_hMom" , &tru_hMom);
  intree->SetBranchAddress("tru_DsMom", &tru_DsMom);

  long ientries(intree->GetEntries());

  TH1D kfactorMC("kfactorMC", "", 1000, 0.9, 1.1);

  // TLatex *label = new TLatex();
  // label->SetTextSize(0.04);
  // label->SetNDC(true);

  double kfactorp(0.0), kfactorm(0.0), kfactorpm(0.0);

  for (unsigned i = 0; i < ientries; ++i) {
    intree->GetEntry(i);

    if (tru_BsMom->P() > 3E5) continue;

    tru_hMom->SetVectM(tru_hMom->Vect(), KMASS);
    Bs_rec = *tru_hMom + *tru_DsMom;

    kfactorp = Bs_rec.P() / tru_BsMom->P();
    kfactorm = tru_BsMom->M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactorMC.Fill(kfactorm);
  }

  std::vector<TLorentzVector> *particle_lvs = NULL;
  TLorentzVector BsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector DsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector hMom(0.0, 0.0, 0.0, 0.0);

  outtree->SetBranchAddress("particle_lvs", &particle_lvs);

  long oentries(outtree->GetEntries());

  // kfactor.Reset("ICESM");
  TH1D kfactortru("kfactortru", "", 1000, 0.9, 1.1);
  TH1D kfactortruf("kfactortruf", "", 1000, 0.9, 1.1);

  for (unsigned i = 0; i < oentries; ++i) {
    outtree->GetEntry(i);
    if (particle_lvs->size() < 3) continue;
    (*particle_lvs)[2].SetVectM((*particle_lvs)[2].Vect(), KMASS*1E-3);
    Bs_rec = (*particle_lvs)[1] + (*particle_lvs)[2];

    kfactorp = Bs_rec.P() / (*particle_lvs)[0].P();
    kfactorm = (*particle_lvs)[0].M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactortru.Fill(kfactorm);

    // 10s of keV
    volatile int iBs[3] = {int((*particle_lvs)[0].X()*1E5), int((*particle_lvs)[0].Y()*1E5),
			   int((*particle_lvs)[0].Z()*1E5)};
    volatile int iDs[3] = {int((*particle_lvs)[1].X()*1E5), int((*particle_lvs)[1].Y()*1E5),
			   int((*particle_lvs)[1].Z()*1E5)};
    volatile int ih[3] = {int((*particle_lvs)[2].X()*1E5), int((*particle_lvs)[2].Y()*1E5),
			  int((*particle_lvs)[2].Z()*1E5)};

    float mBs((*particle_lvs)[0].M()*1E3),
      mDs((*particle_lvs)[1].M()*1E3),
      mh((*particle_lvs)[2].M()*1E3);

    // MeV
    double dBs[4] = {iBs[0]*1E-2, iBs[1]*1E-2, iBs[2]*1E-2, 0.};
    double dDs[4] = {iDs[0]*1E-2, iDs[1]*1E-2, iDs[2]*1E-2, 0.};
    double dh[4] = {ih[0]*1E-2, ih[1]*1E-2, ih[2]*1E-2, 0.};
    dBs[3] = std::sqrt(double(mBs) * double(mBs) + dBs[0] *dBs[0] +
		       dBs[1] * dBs[1] + dBs[2] * dBs[2]);
    dDs[3] = std::sqrt(double(mDs) * double(mDs) + dDs[0] *dDs[0] +
		       dDs[1] * dDs[1] + dDs[2] * dDs[2]);
    dh[3] = std::sqrt(double(mh) * double(mh) + dh[0] *dh[0] +
		       dh[1] * dh[1] + dh[2] * dh[2]);

    float fBs[4] = {float(dBs[0]), float(dBs[1]), float(dBs[2]), float(dBs[3])};
    float fDs[4] = {float(dDs[0]), float(dDs[1]), float(dDs[2]), float(dDs[3])};
    float fh[4] = {float(dh[0]), float(dh[1]), float(dh[2]), float(dh[3])};

    // GeV
    double ddBs[4] = {fBs[0]*1E-3, fBs[1]*1E-3, fBs[2]*1E-3, fBs[3]*1E-3};
    double ddDs[4] = {fDs[0]*1E-3, fDs[1]*1E-3, fDs[2]*1E-3, fDs[3]*1E-3};
    double ddh[4] =  {fh[0]*1E-3, fh[1]*1E-3, fh[2]*1E-3, fh[3]*1E-3};

    BsMom.SetXYZT(ddBs[0], ddBs[1], ddBs[2], ddBs[3]);
    DsMom.SetXYZT(ddDs[0], ddDs[1], ddDs[2], ddDs[3]);
    hMom.SetXYZT(ddh[0], ddh[1], ddh[2], ddh[3]);

    Bs_rec = DsMom + hMom;

    kfactorp = Bs_rec.P() / BsMom.P();
    kfactorm = BsMom.M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactortruf.Fill(kfactorm);
  }

  // kfactor.Draw("hist");
  // label->DrawLatex(0.6, 0.5, "k-factor (MC)");
  // gPad->Update();
  // gPad->Print("kfactor_MC.png");
  // kfactor.Print("all");

  // kfactor.Draw("hist");
  // label->DrawLatex(0.6, 0.5, "k-factor (gen)");
  // gPad->Update();
  // gPad->Print("kfactor_gen.png");

  kfactorMC.SetLineColor(kAzure);
  kfactortru.SetLineColor(kRed);
  kfactortruf.SetLineColor(kGreen);

  kfactorMC.Draw("hist ");
  kfactortru.Draw("hist same");

  TLegend legend(0.6, 0.6, 0.85, 0.4);
  legend.AddEntry(&kfactorMC, "Monte Carlo", "l");
  legend.AddEntry(&kfactortru, "Generated", "l");
  legend.SetFillColor(0);
  legend.SetLineColor(0);
  legend.Draw();
  gPad->Update();
  gPad->Print("kfactor_both.png");

  TFile hdump("hdump.root", "update");
  hdump.WriteTObject(&kfactorMC, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortru, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortruf, NULL, "WriteDelete");
  hdump.Close();

  return 0;
}
