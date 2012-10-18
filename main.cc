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

#include "TwoBodyDecayGen.hxx"


int main(int argc, char* argv[])
{
  if (argc > 3) {
    std::cout << "Too many arguments!" << std::endl;
    return -1;
  }

  int nevents(100);
  std::string mode;
  if (argc == 3) {
    nevents = atol(argv[1]);
    mode = argv[2];
  } else {
    std::cout << "Not enough arguments!" << std::endl;
    return -1;
  }

  std::string fname = "smalltree-" + mode + ".root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("ftree"));
  
  TH1D hmomp("hmomp", "", 100, 0.0, 300.0);
  intree->Draw("1E-3*BsMom.P()>>hmomp");
  gPad->Print("Bs_mom.png");

  RooRealVar momentum("momentum", "momentum", 0.0, 300.0);
  RooDataHist datahist("datahist", "datahist", RooArgList(momentum), &hmomp);
  RooHistPdf histpdf("histpdf", "histpdf", RooArgSet(momentum), datahist);
  RooDataSet *dataset = histpdf.generate(RooArgSet(momentum),
					 RooFit::NumEvents(nevents),
					 RooFit::AutoBinned(false));

  double mom(0.0);
  TTree* momtree = new TTree("momtree", "");
  momtree->Branch("momentum", &mom);

  std::cout << dataset->numEntries() << " events generated." << std::endl;
  for (int i = 0; i < nevents; ++i) {
    const RooArgSet *entry = dataset->get(i);
    mom = entry->getRealValue("momentum");
    momtree->Fill();
  }

  double masses[NDAUS] = {1.968, 0.494};
  TwoBodyDecayGen generator(masses);

  TFile file("eventtree.root", "recreate");
  file.cd();

  TTree* eventtree = generator.get_event_tree(nevents, momtree);
  std::cout << eventtree->GetEntries() << " events generated" << std::endl;
  // eventtree->Scan("evt_wt");
  eventtree->Print("all");

  file.WriteTObject(eventtree);
  file.Close();

  return 0;
}
