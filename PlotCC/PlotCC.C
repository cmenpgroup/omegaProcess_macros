// PlotCC.C
//
// macro to plot eg2a histograms for CC
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

const Int_t NSECTORS = 6;
const Int_t NPHOTONS = 2;
const Int_t NDETPART = 5;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};
char *Plabel[3] = {"#pi^{+}","#pi^{-}","#pi^{+}#pi^{-}"};

char hname[50];
char htitle[500];
char title[500];
char cname[50];
char ctitle[500];
char xtitle[100];
char ytitle[100];
char OutCan[100];
char OutText[100];

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

//
// Plot_CCnphe - plot histogram for the number of photo-electrons by particle
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_CCnphe(char *fAna="Ana.root", char *target="C", int iPart=0)
{
    Int_t i;
    TH2D *h2D;
    TH1D *h1D;
    
    if(iPart>=0 && iPart<NDETPART){
        // Canvas to plot histogram
        TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
        c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        c1->SetBorderSize(5);
        gStyle->SetOptStat(1111);
        c1->SetFillStyle(4000);
    
        // data files contain the trees
        printf("Analyzing file %s\n",fAna);
        TFile *fm = new TFile(fAna,"READ");
        TDirectory *tmp = fm->GetDirectory("Detectors");
    
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
    
        sprintf(hname,"CCnphe");
        h2D = (TH2D*)tmp->Get(hname);

        h1D = (TH1D*)h2D->ProjectionX("h1D",iPart+1,iPart+1,"");
        sprintf(htitle,"%s Runs",target);
        h1D->SetTitle(htitle);
        h1D->GetXaxis()->CenterTitle();
        h1D->GetYaxis()->CenterTitle();
        h1D->GetYaxis()->SetTitleOffset(yoff);
        h1D->SetLineWidth(2);
        h1D->Draw();
    
        sprintf(OutCan,"Plot_CCnphe_%s_%i.gif",target,iPart);
        c1->Print(OutCan);
        sprintf(OutCan,"Plot_CCnphe_%s_%i.eps",target,iPart);
        c1->Print(OutCan);
    }
    else{
        cout<<"Wrong particle number: "<<iPart<<endl;
    }
}


