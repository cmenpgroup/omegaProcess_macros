// PlotTopology.C
//
// macro to plot eg2a histograms for particle topology
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};
char *Plabel[4] = {"e^{-}","#pi^{-}","#pi^{+}","Photons"};

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
// PlotTopo_Grid - plot histogram for the number of detected particles
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void PlotTopo_Grid(char *fAna="Ana.root", char *target="C")
{
    Int_t i;
    TH2D *h2D;
    TH1D *h1D[4];
    
    char strname[50];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(2,2);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
    
    sprintf(hname,"NumDetPart");
    h2D = (TH2D*)tmp->Get(hname);

    for(int i=0; i<4; i++){
        c1->cd(i+1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        sprintf(strname,"h1D_%i",i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        sprintf(htitle,"%s Runs, %s",target,Plabel[i]);
        h1D[i]->SetTitle(htitle);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->Draw();
    }
    
    sprintf(OutCan,"PlotTopo_Grid_%s.gif",target);
    c1->Print(OutCan);
    sprintf(OutCan,"PlotTopo_Grid_%s.eps",target);
    c1->Print(OutCan);
}


