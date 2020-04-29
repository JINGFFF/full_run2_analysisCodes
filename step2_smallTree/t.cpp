#include "smallTree.h"
#include "smallTree.C"
using namespace std;

int main(){

TFile *file1 =new TFile("treePKU_137.root");
TDirectory * dir1 = (TDirectory*)file1->Get("treeDumper");
TTree *tree1 = (TTree*) dir1->Get("PKUCandidates");
TString outname = "out1.root";
smallTree m(tree1,outname);
//m.m_dataset = "out1.root";
m.Loop();
m.endJob();
return 0;
}
