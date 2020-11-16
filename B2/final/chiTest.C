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
void chiTest(TString backgroundfilename="../Bg-build/background_bak.root"){
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
    TH1D* pAniE = new TH1D("n", "n initial Energy",100,10,60);
    Double_t r;
    for(int i=0;i<N;i++){
        tr->GetEntry(i);
        if(evt.pAniEnergy<0.6 && evt.pAniEnergy>0.1){
            r=sqrt(evt.proton_x*evt.proton_x+evt.proton_y*evt.proton_y+evt.proton_z*evt.proton_z)/1000;
            if(r<11){
                pAniE->Fill(evt.pAniEnergy*100);
            }
        }
    }
    TF1* fit = new TF1("f1","[0]*exp(x/1000)",10,60);
    TH1D* h2= new TH1D("fit", "fit the mc data", 100, 10,60);
    TH1D* chi = new TH1D("chi", "chi of Mc data", 100, 0, 100);
    for(Int_t i=0;i<1000;i++){
        h2->FillRandom(pAniE, 1000);
        //fit->SetParameters(100);
        h2->Fit("f1");
    }
}
void bayes(Double_t beta=0.05, Double_t mub=37, Double_t obs=37, Double_t mus=0){
    Double_t bcacu = 0;
    Double_t left =0, right =100;
    Double_t step=1;
    Int_t times =0;
    mus = 50;
    while(abs(beta-bcacu)>0.000001 && times<100){
         bcacu = (1-ROOT::Math::chisquared_cdf(2*(mub+mus),2*(obs+1)))/(1-ROOT::Math::chisquared_cdf(2*(mub),2*(obs+1)));

        if(beta > bcacu){
            right = mus;
        }
        else{
            left = mus;
        }
        mus = (right+left)/2;
        times++;
        //cout<<bcacu<<" "<<beta-bcacu<<endl;
    }
    cout<<"bc:"<<bcacu<<"mus"<<mus<<"times"<<times<<endl;
    // 0.051 13,87
}
void clsCacu(Double_t beta=0.05, Double_t mub=37, Double_t obs=37, Double_t mus=0){
    Double_t bcacu = 0;
    Double_t left =0, right =100;
    Double_t step=1;
    Int_t times =0;
    while(abs(beta-bcacu)>0.000001 && times<100){
         bcacu = (1-ROOT::Math::chisquared_cdf(2*(mub+mus),2*(obs)))/(1-ROOT::Math::chisquared_cdf(2*(mub),2*(obs)));

        if(beta > bcacu){
            right = mus;
        }
        else{
            left = mus;
        }
        mus = (right+left)/2;
        times++;
        //cout<<bcacu<<" "<<beta-bcacu<<endl;
    }
    cout<<"bc:"<<bcacu<<"mus"<<mus<<"times"<<times<<endl;
    // 0.051 13,87
}
void cls(Double_t beta=0.05, Double_t mub=37, Double_t mus=0){

}