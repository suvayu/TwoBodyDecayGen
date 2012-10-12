#include <iostream>

#include "TwoBodyDecayGen.hxx"

int main()
{
  double masses[NDAUS] = {1.968, 0.494};
  TwoBodyDecayGen generator(masses);
  generator.generate(100);
  TTree* eventtree = generator.get_tree();
  std::cout << eventtree->GetEntries() << " events generated" << std::endl;
  eventtree->Scan("evt_wt");
  eventtree->Print("all");
  return 0;
}
