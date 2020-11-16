#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "iostream"
#include "TLegend.h"
#include "TLine.h"
using namespace std;

Double_t Q(Double_t nobs, Double_t ns,Double_t nb)
{

    Double_t a=TMath::Log((ns+nb)/nb);
    return -2.*(-ns+nobs*a);
}

Double_t CLs(Double_t Nobs,Double_t ns,Double_t nb)
{
    ///toy
    Int_t N=10000;
    const Int_t Nbins=50;
    const Double_t Start=-20;
    const Double_t End=-Start;
    TH1D QNsb("QNSB","Q Distribution;Q;count",Nbins,Start,End);
    TH1D QNb("QNB","F QNB",Nbins,Start,End);
    for(Int_t i=0;i<N;i++)
    {
        Double_t nobs=Double_t(gRandom->Poisson(ns+nb));
        //cout<<Q(nobs,ns,nb)<<endl;
        QNsb.Fill(Q(nobs,ns,nb));
    }
    for(Int_t i=0;i<N;i++)
    {
        Double_t nobs=Double_t(gRandom->Poisson(nb));
        QNb.Fill(Q(nobs,ns,nb));
    }

    Double_t Qobs=Q(Nobs,ns,nb);
    /*////////////////                         把*号去掉就可以画图了
    TCanvas mycanvas("mycanvas","my");
    gStyle->SetOptStat(0);
    QNsb.SetFillStyle(3003);
    QNsb.SetFillColor(kRed);

    QNsb.Draw();
    QNb.SetFillStyle(3008);
    QNb.SetFillColor(kBlue);
    QNb.Draw("same");
    /////////////////
    TLine Line(Qobs,0.5,Qobs,400);
    Line.SetLineColor(kRed);
    Line.SetLineWidth(2);
    Line.Draw("same");
    TLegend * leg1=new TLegend(0.7,0.7,0.9,0.9);
    leg1->AddEntry(&QNsb,"Q Ns+Nb");
    leg1->AddEntry(&QNb, "Q Nb");
    leg1->AddEntry(&Line,"Qobs");
    leg1->Draw("same");


    mycanvas.Print("Qdist_.png");
    ///////////////////*/

    Double_t Psb;
    Double_t Pb;
    Double_t norm=1;

    Int_t nbinQ=Int_t((Qobs-Start)/((End-Start)/Nbins));
    QNsb.Scale(norm/QNsb.Integral());
    QNb.Scale(norm/QNb.Integral());

    Psb=QNsb.Integral(nbinQ,Nbins);
    Pb=QNb.Integral(0,nbinQ);
    cout<<ns<<" "<<Qobs<<" "<<nbinQ<<" "<<Psb<<" "<<Pb<<" "<<Psb/(1-Pb)<<endl;
    //delete leg1;
    return Psb/(1-Pb);

}
Double_t CLsplot(Double_t Nobs,Double_t ns,Double_t nb)
{
    ///toy
    Int_t N=10000;
    const Int_t Nbins=50;
    const Double_t Start=-20;
    const Double_t End=-Start;
    TH1D QNsb("QNSB","Q Distribution;Q;count",Nbins,Start,End);
    TH1D QNb("QNB","F QNB",Nbins,Start,End);
    for(Int_t i=0;i<N;i++)
    {
        Double_t nobs=Double_t(gRandom->Poisson(ns+nb));
        //cout<<Q(nobs,ns,nb)<<endl;
        QNsb.Fill(Q(nobs,ns,nb));
    }
    for(Int_t i=0;i<N;i++)
    {
        Double_t nobs=Double_t(gRandom->Poisson(nb));
        QNb.Fill(Q(nobs,ns,nb));
    }

    Double_t Qobs=Q(Nobs,ns,nb);
    /////////////////                         把*号去掉就可以画图了
    TCanvas mycanvas("mycanvas","my");
    gStyle->SetOptStat(0);
    QNsb.SetFillStyle(3003);
    QNsb.SetFillColor(kRed);

    QNsb.Draw();
    QNb.SetFillStyle(3008);
    QNb.SetFillColor(kBlue);
    QNb.Draw("same");
    /////////////////
    TLine Line(Qobs,0.5,Qobs,400);
    Line.SetLineColor(kRed);
    Line.SetLineWidth(2);
    Line.Draw("same");
    TLegend * leg1=new TLegend(0.7,0.7,0.9,0.9);
    leg1->AddEntry(&QNsb,"Q Ns+Nb");
    leg1->AddEntry(&QNb, "Q Nb");
    leg1->AddEntry(&Line,"Qobs");
    leg1->Draw("same");


    mycanvas.Print(TString::Format("Qdist_%.2f_0.05.png",ns));
    ///////////////////*/

    Double_t Psb;
    Double_t Pb;
    Double_t norm=1;

    Int_t nbinQ=Int_t((Qobs-Start)/((End-Start)/Nbins));
    QNsb.Scale(norm/QNsb.Integral());
    QNb.Scale(norm/QNb.Integral());

    Psb=QNsb.Integral(nbinQ,Nbins);
    Pb=QNb.Integral(0,nbinQ);
    cout<<ns<<" "<<Qobs<<" "<<nbinQ<<" "<<Psb<<" "<<Pb<<" "<<Psb/(1-Pb)<<endl;
    //delete leg1;
    return Psb/(1-Pb);

}

void CLsMethod(Double_t beta=0.05)
{
    Double_t CLS=0;
    Double_t Ns=30;
    Double_t left =0, right=50;
    Int_t times =0;
    while(abs(CLS-beta)>0.00001 &&times<100){
        times++;
        CLS=CLs(365.,Ns,365.);
        if(CLS<beta){
            right = Ns;
        }else{
            left = Ns;
        }
        Ns=(right+left)/2;
    }
    CLS=CLsplot(365.,Ns,365.);

    cout<<Ns<<endl;
    cout<<CLS<<endl;
    cout<<times<<endl;

}
/*
void CLsMethod(){
    Double_t CLS=0;
    Double_t Ns=30;
    do{
        Ns+=0.1;
    CLS=CLs(365.,Ns,365.);

    }while(CLS>0.1);
    CLS=CLs(365.,Ns,365.);

    cout<<Ns<<endl;
    cout<<CLS<<endl;
}*/