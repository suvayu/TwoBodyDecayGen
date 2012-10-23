#include <iostream>
#include <cstdlib>
#include <cassert>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TPad.h>

#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooPlot.h>

#include "TwoBodyDecayGen.hxx"


void usage(char * prog)
{
  std::cout << "Usage: $ " << prog << " <nevents> <mode>"
    " # args are case sensitive" << std::endl;
}


int main(int argc, char* argv[])
{
  // program arguments
  if (argc > 3) {
    std::cout << "Too many arguments!" << std::endl;
    usage(argv[0]);
    return -1;
  }

  int nevents(100);
  std::string mode;
  if (argc == 3) {
    nevents = atol(argv[1]);
    mode = argv[2];
  } else {
    std::cout << "Not enough arguments!" << std::endl;
    usage(argv[0]);
    return -1;
  }

  // read ntuple from file
  std::string fname = "smalltree-" + mode + ".root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("ftree"));

  // get B momentum distribution
  TH1D Bsmomp("Bsmomp", "", 100, 0.0, 300.0);
  intree->Draw("1E-3*BsMom.P()>>Bsmomp");
  gPad->Print("Bs_mom.png");

  TH1D dau1momp("dau1momp", "", 100, 0.0, 300.0);
  intree->Draw("1E-3*DsMom.P()>>dau1momp");
  gPad->Print("dau1_mom.png");

  TH1D dau2momp("dau2momp", "", 100, 0.0, 300.0);
  intree->Draw("1E-3*hMom.P()>>dau2momp");
  gPad->Print("dau2_mom.png");

  // make dataset from histogram
  RooRealVar momentum("momentum", "momentum", 0.0, 300.0);
  RooDataHist datahist("datahist", "datahist", RooArgList(momentum), &Bsmomp);
  RooHistPdf histpdf("histpdf", "histpdf", RooArgSet(momentum), datahist);
  RooDataSet *dataset = histpdf.generate(RooArgSet(momentum),
					 RooFit::NumEvents(nevents),
					 RooFit::AutoBinned(false));

  // xcheck: plot generated B momenta
  RooPlot *momframe = momentum.frame();
  dataset->plotOn(momframe);
  momframe->Draw();
  gPad->Print("Bs_mom_gen.png");

  // ROOT file dump
  TFile *file = new TFile("eventtree.root", "recreate");
  file->cd();

  // fill tree from dataset
  double mom(0.0);
  TTree* momtree = new TTree("momtree", "");
  momtree->Branch("momentum", &mom);
  for (int i = 0; i < nevents; ++i) {
    const RooArgSet *entry = dataset->get(i);
    mom = entry->getRealValue("momentum");
    momtree->Fill();
  }

  // generator config
  double masses[NDAUS] = {1.968, 0.494};
  TwoBodyDecayGen generator(masses);

  // generate, print summary and dump to ROOT file
  TTree* eventtree = generator.get_event_tree(nevents, momtree, mom);
  eventtree->Print("all");
  file->WriteTObject(momtree);
  file->WriteTObject(eventtree);
  file->Close();

  return 0;
}
