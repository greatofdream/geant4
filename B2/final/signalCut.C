#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;
struct Evt{
    Double_t nu_eEnergy;
    Double_t eEnergy;
    Double_t nu_muEnergy;
};
void signalCut(TString filename="../B2b-build/signal.root"){
    TFile* file = new TFile(filename);
    TTree* tree =(TTree*) file->Get("Energy");
    Evt evt;
    Int_t N=tree->GetEntries();
    TBranch* b = tree->GetBranch("nu_eEnergy");
    TBranch* b2 = tree->GetBranch("eEnergy");
    TBranch* b3 = tree->GetBranch("nu_muEnergy");
    //tree->SetBranchAddress("row_wise_branch",&evt);
    b->SetAddress(&evt.nu_eEnergy);
    b2->SetAddress(&evt.eEnergy);
    b3->SetAddress(&evt.nu_muEnergy);
    TH1D* nu_eE = new TH1D("anti_nu", "anti nu_e spectrum",100,0,70);
    TH1D* nu_muE = new TH1D("nu", "nu_mu spectrum",100,0,70);
    TH1D* eE = new TH1D("e", "e spectrum",100,0,70);
    TH1D* total = new TH1D("total", "total energy", 100,0,300);
    cout<<N<<endl;
    Int_t ncut=0;
    for(int i=0;i<N;i++){
        tree->GetEntry(i);
        if(evt.nu_eEnergy>10){
            ncut += 1;
            nu_eE->Fill(evt.nu_eEnergy);
        }
        //nu_muE->Fill(evt.nu_muEnergy);
        //eE->Fill(evt.eEnergy);
        //total->Fill(evt.nu_eEnergy+evt.nu_muEnergy+evt.eEnergy);
        //cout<<evt.nu_muEnergy<<" "<<evt.nu_eEnergy<<" "<<evt.eEnergy<<endl;
    }
    cout<<"total cut"<<Double_t(ncut)/N<<endl;
    //return;
    TCanvas* c1 = new TCanvas("c1", "mu decay",1000, 800);
    bool logstyle = false;
    if(logstyle)
        gStyle->SetOptLogy();
    ///c1->Divide(2,2);
    //c1->cd(1);
    nu_eE->SetXTitle("Energy(MeV)");
    nu_eE->SetTitle("nu_e spectra");
    nu_eE->Draw();
    /*c1->cd(2);
    nu_muE->SetXTitle("Energy(MeV)");
    nu_muE->SetTitle("nu_mu spectra");
    nu_muE->Draw();
    c1->cd(3);
    eE->SetXTitle("Energy(MeV)");
    eE->SetTitle("e- spectra");
    eE->Draw();
    c1->cd(4);
    total->SetXTitle("Energy(MeV)");
    total->SetTitle("total spectra");
    total->Draw();
    */
    if(logstyle)
        c1->SaveAs(TString::Format("signal%s.png","log"));
    else
    {
        c1->SaveAs(TString::Format("signal%s.png","linear"));
    }
    
}