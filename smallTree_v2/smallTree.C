#define smallTree_cxx
#include "smallTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void smallTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L smallTree.C
//      root> smallTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   int num = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t npp = fChain->GetEntries("theWeight>0.");
   Long64_t nmm = fChain->GetEntries("theWeight<0.");
   cout<< "numberofnp:" << npp << "  numberofnm:" <<nmm << endl;

   TFile * input1 = new TFile ("puweight.root");
   TH1D* h = (TH1D*)input1->Get("h2");

   TFile * input2 = new TFile ("photonmediumID.root");
   TH2D* h2 = (TH2D*)input2->Get("EGamma_SF2D");

   TFile * input3 = new TFile ("muonTRIGGERbf.root");
   TDirectory * dir3 = (TDirectory*)input3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins");
   TH2D* h3 = (TH2D*)dir3->Get("abseta_pt_ratio");

   TFile * input4 = new TFile ("muonTRIGGERgh.root");
   TDirectory * dir4 = (TDirectory*)input4->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins");
   TH2D* h4 = (TH2D*)dir4->Get("abseta_pt_ratio");

   TFile * input5 = new TFile ("muonIDbf.root");
   TDirectory * dir5 = (TDirectory*)input5->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
   TH2D* h5 = (TH2D*)dir5->Get("abseta_pt_ratio");

   TFile * input6 = new TFile ("muonIDgh.root");
   TDirectory * dir6 = (TDirectory*)input6->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
   TH2D* h6 = (TH2D*)dir6->Get("abseta_pt_ratio");

   TFile * input7 = new TFile ("muonISObf.root");
   TDirectory * dir7 = (TDirectory*)input7->Get("TightISO_TightID_pt_eta");
   TH2D* h7 = (TH2D*)dir7->Get("abseta_pt_ratio");

   TFile * input8 = new TFile ("muonISOgh.root");
   TDirectory * dir8 = (TDirectory*)input8->Get("TightISO_TightID_pt_eta");
   TH2D* h8 = (TH2D*)dir8->Get("abseta_pt_ratio");

   TFile * input9 = new TFile ("elereco.root");
   TH2D* h9 = (TH2D*)input9->Get("EGamma_SF2D");

   TFile * input10 = new TFile ("eletightID.root");
   TH2D* h10 = (TH2D*)input10->Get("EGamma_SF2D");

   TFile * input11 = new TFile ("eletrigger.root");
   TH2D* h11 = (TH2D*)input11->Get("Ele27_WPTight_Gsf");

   Long64_t nbytes = 0, nb = 0;
//-----------------Loop
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0||jentry>10000) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      scalef = 1;
      use_f(type);
      //------cross section factor
      if(m_dataset.Contains("aqgc.root")){ scalef=1000.*1.864/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WAJJ.root")){ scalef=1000.*0.776/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WA.root")){ scalef=1000.*178.6/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WJets-2.root")){ scalef=1000.*61526.7/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WJets-1.root")){ scalef=1000.*61526.7/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("ZJets.root")){ scalef=1000.*5765.4/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("ZA.root")){ scalef=1000.*47.46/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("TTA.root")){ scalef=1000.*3.697/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("TTJets.root")){ scalef=1000.*831.76/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("STs.root")){ scalef=1000.*3.36/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("STtbarw.root")){ scalef=1000.*35.85/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("STtw.root")){ scalef=1000.*35.85/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("STt.root")){ scalef=1000.*136.02/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("STtbar.root")){ scalef=1000.*80.95/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WW.root")){ scalef=1000.*118.7/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WZ.root")){ scalef=1000.*47.13/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("ZZ.root")){ scalef=1000.*16.523/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("DY.root")){ scalef=1000.*5765.4/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("DY1.root")){ scalef=1000.*1016/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("DY2.root")){ scalef=1000.*331.4/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("DY3.root")){ scalef=1000.*96.36/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("DY4.root")){ scalef=1000.*51.4/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("W3Jets.root")){ scalef=1000.*942.3/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("QCD.root")){ scalef=1000.*162060000/float(npp-nmm)*fabs(theWeight)/theWeight; }
      if(m_dataset.Contains("WA-mlm.root")){ scalef=1000.*405.271/float(npp-nmm)*fabs(theWeight)/theWeight; }
      //------cross section factor

      //------photon electron veto scalefactor
      if(PHOTONET>10 && PHOTONET<200)
      {
         if(fabs(PHOTONETA)>0 && fabs(PHOTONETA)<1.4442){scalef=scalef*0.9978;}
         if(fabs(PHOTONETA)>1.566 && fabs(PHOTONETA)<2.5){scalef=scalef*0.9931;}
      }//------photon electron veto scalefactor


      //------add factor for MC
      if (!m_dataset.Contains("SMu") && !m_dataset.Contains("SEle") ) {
         pileupWeight=h->GetBinContent(h->GetXaxis()->FindBin(npT));
         scale_sys_up=scalef;
         scale_sys_low=scalef;

         TLorentzVector p4, v4;
         p4.SetPtEtaPhiM(PHOTONET, PHOTONETA, PHOTONPHI,0);
         v4.SetPtEtaPhiM(ptVlepJEC, yVlepJEC, phiVlepJEC, massVlepJEC);
         WGmass=(p4+v4).M();

         //------photon medium ID scalefactor
         if(PHOTONET<500)
         {
            scalef=scalef*h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET));
            scale_sys_up=scale_sys_up*(h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET))+h2->GetBinErrorUp(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)));
            scale_sys_low=scale_sys_low*(h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET))-h2->GetBinErrorLow(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)));
         }

         if(PHOTONET>=500)
         {
            scalef=scalef*h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)-1);
            scale_sys_up=scale_sys_up*(h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)-1)+h2->GetBinErrorUp(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)-1));
            scale_sys_low=scale_sys_low*(h2->GetBinContent(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)-1)-h2->GetBinErrorLow(h2->GetXaxis()->FindBin(PHOTONSCETA),h2->GetYaxis()->FindBin(PHOTONET)-1));
         } //------photon medium ID scalefactor

         double bf=0.5497;   double gh=0.4503;
         //------add factor for muon
         if(lep==13)
         {
            //------muon tracking scalefactor
            if(fabs(etalep1)<0.20000004)
            {scalef=scalef*0.9969965; scale_sys_low=scale_sys_low*(0.9969965-7.293973e-05); scale_sys_up=scale_sys_up*(0.9969965+7.271734e-05);}
            if(fabs(etalep1)>0.20000004 && fabs(etalep1)<0.40000002)
            {scalef=scalef*0.9977118; scale_sys_low=scale_sys_low*(0.9977118-8.230405e-05); scale_sys_up=scale_sys_up*(0.9977118+8.204939e-05);}
            if(fabs(etalep1)>0.40000002 && fabs(etalep1)<0.59999999)
            {scalef=scalef*0.9980776; scale_sys_low=scale_sys_low*(0.9980776-7.097223e-05); scale_sys_up=scale_sys_up*(0.9980776+7.158818e-05);}
            if(fabs(etalep1)>0.59999999 && fabs(etalep1)<0.80000004)
            {scalef=scalef*0.9978039; scale_sys_low=scale_sys_low*(0.9978039-7.447371e-05); scale_sys_up=scale_sys_up*(0.9978039+7.435091e-05);}
            if(fabs(etalep1)>0.80000004 && fabs(etalep1)<1.00000003)
            {scalef=scalef*0.9979708; scale_sys_low=scale_sys_low*(0.9979708-0.0001073152); scale_sys_up=scale_sys_up*(0.9979708+0.0001061586);}
            if(fabs(etalep1)>1.00000003 && fabs(etalep1)<1.20000008)
            {scalef=scalef*0.9971477; scale_sys_low=scale_sys_low*(0.9971477-0.0001852945); scale_sys_up=scale_sys_up*(0.9971477+0.0001864667);}
            if(fabs(etalep1)>1.20000008 && fabs(etalep1)<1.40000008)
            {scalef=scalef*0.9962274; scale_sys_low=scale_sys_low*(0.9962274-0.000180911); scale_sys_up=scale_sys_up*(0.9962274+0.0001816305);}
            if(fabs(etalep1)>1.40000008 && fabs(etalep1)<1.60000003)
            {scalef=scalef*0.9954786; scale_sys_low=scale_sys_low*(0.9954786-0.0001704205); scale_sys_up=scale_sys_up*(0.9954786+0.000170625);}
            if(fabs(etalep1)>1.60000003 && fabs(etalep1)<1.80000004)
            {scalef=scalef*0.9957808; scale_sys_low=scale_sys_low*(0.9957808-0.0001752651); scale_sys_up=scale_sys_up*(0.9957808+0.0001767324);}
            if(fabs(etalep1)>1.80000004 && fabs(etalep1)<2.00000004)
            {scalef=scalef*0.9938919; scale_sys_low=scale_sys_low*(0.9938919-0.0002350023); scale_sys_up=scale_sys_up*(0.9938919+0.0002347843);}
            if(fabs(etalep1)>2.00000004 && fabs(etalep1)<2.20000001)
            {scalef=scalef*0.9929427; scale_sys_low=scale_sys_low*(0.9929427-0.0003260826); scale_sys_up=scale_sys_up*(0.9929427+0.000330644);}
            if(fabs(etalep1)>2.20000001 && fabs(etalep1)<2.4)
            {scalef=scalef*0.9873133; scale_sys_low=scale_sys_low*(0.9873133-0.0008527153); scale_sys_up=scale_sys_up*(0.9873133+0.0008593464);}
            //------muon tracking scalefactor

            //------muon trigger scalefactor
            if(ptlep1<800)
            {
               scalef=scalef*(bf*h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1))+gh*h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)));
               scale_sys_up=scale_sys_up*(bf*(h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1))+h3->GetBinErrorUp(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)))+gh*(h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1))+h4->GetBinErrorUp(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1))));
               scale_sys_low=scale_sys_low*(bf*(h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1))-h3->GetBinErrorLow(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)))+gh*(h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1))-h4->GetBinErrorLow(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1))));
            }

            if(ptlep1>=800)
            {
               scalef=scalef*(bf*h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)-1)+gh*h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)-1));
               scale_sys_up=scale_sys_up*(bf*(h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)-1)+h3->GetBinErrorUp(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)-1))+gh*(h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)-1)+h4->GetBinErrorUp(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)-1)));
               scale_sys_low=scale_sys_low*(bf*(h3->GetBinContent(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)-1)-h3->GetBinErrorLow(h3->GetXaxis()->FindBin(fabs(etalep1)),h3->GetYaxis()->FindBin(ptlep1)-1))+gh*(h4->GetBinContent(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)-1)-h4->GetBinErrorLow(h4->GetXaxis()->FindBin(fabs(etalep1)),h4->GetYaxis()->FindBin(ptlep1)-1)));
            }
            //------muon trigger scalefactor

            //------muon tight ID  and tight ISO scalefactor
            if(ptlep1<120)
            {
               scalef=scalef*(bf*h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1))+gh*h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)))*(bf*h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1))+gh*h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)));

               scale_sys_up=scale_sys_up*(bf*(h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1))+h5->GetBinErrorUp(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)))+gh*(h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1))+h6->GetBinErrorUp(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1))))*(bf*(h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1))+h7->GetBinErrorUp(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)))+gh*(h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1))+h8->GetBinErrorUp(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1))));

               scale_sys_low=scale_sys_low*(bf*(h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1))-h5->GetBinErrorLow(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)))+gh*(h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1))-h6->GetBinErrorLow(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1))))*(bf*(h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1))-h7->GetBinErrorLow(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)))+gh*(h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1))-h8->GetBinErrorLow(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1))));
            }

            if(ptlep1>=120)
            {
               scalef=scalef*(bf*h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)-1)+gh*h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)-1))*(bf*h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)-1)+gh*h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)-1));

               scale_sys_up=scale_sys_up*(bf*(h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)-1)+h5->GetBinErrorUp(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)-1))+gh*(h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)-1)+h6->GetBinErrorUp(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)-1)))*(bf*(h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)-1)+h7->GetBinErrorUp(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)-1))+gh*(h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)-1)+h8->GetBinErrorUp(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)-1)));

              scale_sys_low=scale_sys_low*(bf*(h5->GetBinContent(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)-1)-h5->GetBinErrorLow(h5->GetXaxis()->FindBin(fabs(etalep1)),h5->GetYaxis()->FindBin(ptlep1)-1))+gh*(h6->GetBinContent(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)-1)-h6->GetBinErrorLow(h6->GetXaxis()->FindBin(fabs(etalep1)),h6->GetYaxis()->FindBin(ptlep1)-1)))*(bf*(h7->GetBinContent(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)-1)-h7->GetBinErrorLow(h7->GetXaxis()->FindBin(fabs(etalep1)),h7->GetYaxis()->FindBin(ptlep1)-1))+gh*(h8->GetBinContent(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)-1)-h8->GetBinErrorLow(h8->GetXaxis()->FindBin(fabs(etalep1)),h8->GetYaxis()->FindBin(ptlep1)-1)));
            }
            //------muon tight ID  and tight ISO scalefactor
         }//------add factor for muon

         //------add factor for electron
         if(lep==11)
         {
            //------electron reco scalefactor
            if(ptlep1<500)
            {
               scalef=scalef*h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1));
               scale_sys_up=scale_sys_up*(h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1))+h9->GetBinErrorUp(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)));
               scale_sys_low=scale_sys_low*(h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1))-h9->GetBinErrorLow(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)));
            }

            if(ptlep1>=500)
            {
               scalef=scalef*h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)-1);
               scale_sys_up=scale_sys_up*(h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)-1)+h9->GetBinErrorUp(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)-1));
               scale_sys_low=scale_sys_low*(h9->GetBinContent(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)-1)-h9->GetBinErrorLow(h9->GetXaxis()->FindBin(etalep1),h9->GetYaxis()->FindBin(ptlep1)-1));
            }
            //------electron reco scalefactor

            //------electron tight ID scalefactor
            if(ptlep1<500)
            {
               scalef=scalef*h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1));
               scale_sys_up=scale_sys_up*(h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1))+h10->GetBinErrorUp(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)));
               scale_sys_low=scale_sys_low*(h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1))-h10->GetBinErrorLow(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)));
            }

            if(ptlep1>=500)
            {
               scalef=scalef*h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)-1);
               scale_sys_up=scale_sys_up*(h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)-1)+h10->GetBinErrorUp(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)-1));
               scale_sys_low=scale_sys_low*(h10->GetBinContent(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)-1)-h10->GetBinErrorLow(h10->GetXaxis()->FindBin(etalep1),h10->GetYaxis()->FindBin(ptlep1)-1));
            }
            //------electron tight ID scalefactor

            //------electron trigger scalefactor
            if(etalep1<200)
            {
               scalef=scalef*h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1),h11->GetYaxis()->FindBin(etalep1));
               scale_sys_up=scale_sys_up*(h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1),h11->GetYaxis()->FindBin(etalep1))+h11->GetBinErrorUp(h11->GetXaxis()->FindBin(ptlep1),h11->GetYaxis()->FindBin(etalep1)));
               scale_sys_low=scale_sys_low*(h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1),h11->GetYaxis()->FindBin(etalep1))-h11->GetBinErrorLow(h11->GetXaxis()->FindBin(ptlep1),h11->GetYaxis()->FindBin(etalep1)));
            }

            if(etalep1>=200)
            {
               scalef=scalef*h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1)-1,h11->GetYaxis()->FindBin(etalep1));
               scale_sys_up=scale_sys_up*(h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1)-1,h11->GetYaxis()->FindBin(etalep1))+h11->GetBinErrorUp(h11->GetXaxis()->FindBin(ptlep1)-1,h11->GetYaxis()->FindBin(etalep1)));
               scale_sys_low=scale_sys_low*(h11->GetBinContent(h11->GetXaxis()->FindBin(ptlep1)-1,h11->GetYaxis()->FindBin(etalep1))-h11->GetBinErrorLow(h11->GetXaxis()->FindBin(ptlep1)-1,h11->GetYaxis()->FindBin(etalep1)));
            }
            //------electron trigger scalefactor

         }//------add factor for electron

         //------add factor for jets b, c
    if(fabs(JET1ETA)<2.4 && fabs(JET2ETA)<2.4)
    {
        if(fabs(JET1PF)==5 && fabs(JET2PF)==5)  //partonflavour b=5  c=4 
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*central_b_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*up_b_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*down_b_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-beff(JET1PT)*central_b_scale(JET1PT))/(1-beff(JET1PT))*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*(1-beff(JET1PT)*up_b_scale(JET1PT))/(1-beff(JET1PT))*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*(1-beff(JET1PT)*down_b_scale(JET1PT))/(1-beff(JET1PT))*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
        }
        if(fabs(JET1PF)==5 && fabs(JET2PF)==4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*central_c_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*up_c_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*down_c_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-beff(JET1PT)*central_b_scale(JET1PT))/(1-beff(JET1PT))*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*(1-beff(JET1PT)*up_b_scale(JET1PT))/(1-beff(JET1PT))*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*(1-beff(JET1PT)*down_b_scale(JET1PT))/(1-beff(JET1PT))*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
        }

        if(fabs(JET1PF)==5 && fabs(JET2PF)!=5 && fabs(JET2PF)!=4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*central_l_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*up_l_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*down_l_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT)*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT)*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT)*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-beff(JET1PT)*central_b_scale(JET1PT))/(1-beff(JET1PT))*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*(1-beff(JET1PT)*up_b_scale(JET1PT))/(1-beff(JET1PT))*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*(1-beff(JET1PT)*down_b_scale(JET1PT))/(1-beff(JET1PT))*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
        }
        if(fabs(JET1PF)==4 && fabs(JET2PF)==5)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*central_b_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*up_b_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*down_b_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-ceff(JET1PT)*central_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*(1-ceff(JET1PT)*up_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*(1-ceff(JET1PT)*down_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
        }

        if(fabs(JET1PF)==4 && fabs(JET2PF)==4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*central_c_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*up_c_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*down_c_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-ceff(JET1PT)*central_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*(1-ceff(JET1PT)*up_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*(1-ceff(JET1PT)*down_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
        }
        if(fabs(JET1PF)==4 && fabs(JET2PF)!=5 && fabs(JET2PF)!=4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*central_l_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*up_l_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*down_l_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT)*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT)*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT)*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-ceff(JET1PT)*central_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*(1-ceff(JET1PT)*up_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*(1-ceff(JET1PT)*down_c_scale(JET1PT))/(1-ceff(JET1PT))*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
        }

        if(fabs(JET1PF)!=5 && fabs(JET1PF)!=4 && fabs(JET2PF)==5)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*central_b_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*up_b_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*down_b_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-leff(JET1PT)*central_l_scale(JET1PT))/(1-leff(JET1PT))*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*(1-leff(JET1PT)*up_l_scale(JET1PT))/(1-leff(JET1PT))*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*(1-leff(JET1PT)*down_l_scale(JET1PT))/(1-leff(JET1PT))*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
        }
        if(fabs(JET1PF)!=5 && fabs(JET1PF)!=4 && fabs(JET2PF)==4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*central_c_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*up_c_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*down_c_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-leff(JET1PT)*central_l_scale(JET1PT))/(1-leff(JET1PT))*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*(1-leff(JET1PT)*up_l_scale(JET1PT))/(1-leff(JET1PT))*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*(1-leff(JET1PT)*down_l_scale(JET1PT))/(1-leff(JET1PT))*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
        }

        if(fabs(JET1PF)!=5 && fabs(JET1PF)!=4 && fabs(JET2PF)!=5 && fabs(JET2PF)!=4)
        {
            if(JET1ICSV>0.8484 && JET2ICSV>0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*central_l_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*up_l_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*down_l_scale(JET2PT);
            }
            if(JET1ICSV>0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT)*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT)*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT)*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
            if(JET1ICSV<0.8484 && JET2ICSV<0.8484)
            {
            scalef=scalef*(1-leff(JET1PT)*central_l_scale(JET1PT))/(1-leff(JET1PT))*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*(1-leff(JET1PT)*up_l_scale(JET1PT))/(1-leff(JET1PT))*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*(1-leff(JET1PT)*down_l_scale(JET1PT))/(1-leff(JET1PT))*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
        }
    }
    if(fabs(JET1ETA)<2.4 && fabs(JET2ETA)>2.4)
    {
        if(fabs(JET1PF)==5)
        {
            if(JET1ICSV>0.8484)
            {
            scalef=scalef*central_b_scale(JET1PT);
            scale_sys_up=scale_sys_up*up_b_scale(JET1PT);
            scale_sys_low=scale_sys_low*down_b_scale(JET1PT);
            }
            if(JET1ICSV<0.8484)
            {
            scalef=scalef*(1-beff(JET1PT)*central_b_scale(JET1PT))/(1-beff(JET1PT));
            scale_sys_up=scale_sys_up*(1-beff(JET1PT)*up_b_scale(JET1PT))/(1-beff(JET1PT));
            scale_sys_low=scale_sys_low*(1-beff(JET1PT)*down_b_scale(JET1PT))/(1-beff(JET1PT));
            }
        }
        if(fabs(JET1PF)==4)
        {
            if(JET1ICSV>0.8484)
            {
            scalef=scalef*central_c_scale(JET1PT);
            scale_sys_up=scale_sys_up*up_c_scale(JET1PT);
            scale_sys_low=scale_sys_low*down_c_scale(JET1PT);
            }
            if(JET1ICSV<0.8484)
            {
            scalef=scalef*(1-ceff(JET1PT)*central_c_scale(JET1PT))/(1-ceff(JET1PT));
            scale_sys_up=scale_sys_up*(1-ceff(JET1PT)*up_c_scale(JET1PT))/(1-ceff(JET1PT));
            scale_sys_low=scale_sys_low*(1-ceff(JET1PT)*down_c_scale(JET1PT))/(1-ceff(JET1PT));
            }
        }
        if(fabs(JET1PF)!=5 && fabs(JET1PF)!=4)
        {
            if(JET1ICSV>0.8484)
            {
            scalef=scalef*central_l_scale(JET1PT);
            scale_sys_up=scale_sys_up*up_l_scale(JET1PT);
            scale_sys_low=scale_sys_low*down_l_scale(JET1PT);
            }
            if(JET1ICSV<0.8484)
            {
            scalef=scalef*(1-leff(JET1PT)*central_l_scale(JET1PT))/(1-leff(JET1PT));
            scale_sys_up=scale_sys_up*(1-leff(JET1PT)*up_l_scale(JET1PT))/(1-leff(JET1PT));
            scale_sys_low=scale_sys_low*(1-leff(JET1PT)*down_l_scale(JET1PT))/(1-leff(JET1PT));
            }
        }
    }

    if(fabs(JET1ETA)>2.4 && fabs(JET2ETA)<2.4)
    {
        if(fabs(JET2PF)==5)
        {
            if(JET2ICSV>0.8484)
            {
            scalef=scalef*central_b_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_b_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_b_scale(JET2PT);
            }
            if(JET2ICSV<0.8484)
            {
            scalef=scalef*(1-beff(JET2PT)*central_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_up=scale_sys_up*(1-beff(JET2PT)*up_b_scale(JET2PT))/(1-beff(JET2PT));
            scale_sys_low=scale_sys_low*(1-beff(JET2PT)*down_b_scale(JET2PT))/(1-beff(JET2PT));
            }
        }
        if(fabs(JET2PF)==4)
        {
            if(JET2ICSV>0.8484)
            {
            scalef=scalef*central_c_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_c_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_c_scale(JET2PT);
            }
            if(JET2ICSV<0.8484)
            {
            scalef=scalef*(1-ceff(JET2PT)*central_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_up=scale_sys_up*(1-ceff(JET2PT)*up_c_scale(JET2PT))/(1-ceff(JET2PT));
            scale_sys_low=scale_sys_low*(1-ceff(JET2PT)*down_c_scale(JET2PT))/(1-ceff(JET2PT));
            }
        }
        if(fabs(JET2PF)!=5 && fabs(JET1PF)!=4)
        {
            if(JET2ICSV>0.8484)
            {
            scalef=scalef*central_l_scale(JET2PT);
            scale_sys_up=scale_sys_up*up_l_scale(JET2PT);
            scale_sys_low=scale_sys_low*down_l_scale(JET2PT);
            }
            if(JET2ICSV<0.8484)
            {
            scalef=scalef*(1-leff(JET2PT)*central_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_up=scale_sys_up*(1-leff(JET2PT)*up_l_scale(JET2PT))/(1-leff(JET2PT));
            scale_sys_low=scale_sys_low*(1-leff(JET2PT)*down_l_scale(JET2PT))/(1-leff(JET2PT));
            }
        }
    }

         //------add factor for jets b, c


      }//------add factor for MC



      //--------fill tree
      if(iphoton>=0 || iphoton_f>=0){
          num++;
          //scalef = 2;
          cout<<"iphoton: "<<iphoton<<"  iphoton_f: "<<iphoton_f<<endl;
          ExTree->Fill();
      }//--------fill tree



      // if (Cut(ientry) < 0) continue;
   }//----------------Loop
}
