#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;
struct Evt{
    Double_t nEnergy;
    Double_t pEnergy;
    Double_t pAniEnergy;
    Double_t proton_x;
    Double_t proton_y;
    Double_t proton_z;
};
void bgCut(TString backgroundfilename="../Bg-build/background_bak.root"){
    TFile* file = new TFile(backgroundfilename);
    TTree* tr =(TTree*) file->Get("Energy");
    Evt evt;
    Int_t N=tr->GetEntries();
    TBranch* b = tr->GetBranch("nEnergy");
    TBranch* b2 = tr->GetBranch("pEnergy");
    TBranch* b3 = tr->GetBranch("pAniEnergy");
    TBranch* bx = tr->GetBranch("proton_x");
    TBranch* by = tr->GetBranch("proton_z");
    TBranch* bz = tr->GetBranch("proton_y");
    //tr->SetBranchAddress("row_wise_branch",&evt);
    b->SetAddress(&evt.nEnergy);
    b2->SetAddress(&evt.pEnergy);
    b3->SetAddress(&evt.pAniEnergy);
    bx->SetAddress(&evt.proton_x);
    by->SetAddress(&evt.proton_y);
    bz->SetAddress(&evt.proton_z);
    TH1D* nE = new TH1D("n", "n initial Energy",100,0,10);
    TH1D* pE = new TH1D("p", "p initial Energy",100,0,10);
    TH1D* pAniE = new TH1D("pAni", "p deposit Energy",100,0,80);
    TH1D* R = new TH1D("R", "radius", 100,0,12);
    cout<<N<<endl;
    Int_t bkEnNum=0, bkRaNum=0, bkCut=0;
    Double_t r;
    for(int i=0;i<N;i++){
        tr->GetEntry(i);
        if(evt.pAniEnergy<0.6 && evt.pAniEnergy>0.1){
            bkEnNum += 1;
            r=sqrt(evt.proton_x*evt.proton_x+evt.proton_y*evt.proton_y+evt.proton_z*evt.proton_z)/1000;
            if(r<11){
                bkCut += 1;
                pAniE->Fill(evt.pAniEnergy*100);
                R->Fill(r);
            }
        }
    }
    cout<<"Energy Cut"<<Double_t(bkEnNum)/N<<"Radius cut"<<Double_t(bkCut)/bkEnNum<<"total cut"<<Double_t(bkCut)/N<<endl;

    //return;
    TCanvas* c1 = new TCanvas("c1", "mu decay",1000, 400);
    bool logstyle = false;
    if(logstyle)
        gStyle->SetOptLogy();
    c1->Divide(2,1);
    c1->cd(1);
    pAniE->SetXTitle("Energy(MeV)");
    pAniE->SetTitle("p Annilate spectra");
    pAniE->Draw();
    c1->cd(2);
    R->SetXTitle("/m");
    R->SetTitle("Radius p generate");
    R->Draw();
    if(logstyle)
        c1->SaveAs(TString::Format("background%s_2.png","log"));
    else
    {
        c1->SaveAs(TString::Format("background%s_2.png","linear"));
    }
    
}
