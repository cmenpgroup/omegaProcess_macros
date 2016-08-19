// ProjectHists_dmsOmega.C
//
// macro to project 1D eg2a histograms from 2D histograms
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.75;

// 
// ProjectHists_dmsOmega - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  fDir = directory name in ROOT file
//                  hname = histogram name
//                  title = histogram title
//                  tgtOut = target name
//                  chanLo = lower bin
//                  chanHi = upper bin
//
void ProjectHists_dmsOmega(char *fAna,  char *fDir, char *hname, char *title, char *tgtOut, Int_t chanLo=0, Int_t chanHi=0, char* projAxis="x")
{
	char OutCan[100];
    char strname[100];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
	
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(fDir);
	
	c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH2D *h2D = (TH2D*)tmp->Get(hname);
    TH1D *h1D;
    
    switch(projAxis){
        case "x":
        case "X":
            sprintf(strname,"%s_px_%i_%i",hname,chanLo,chanHi);
            h1D = (TH1D*)h2D->ProjectionX(strname,chanLo,chanHi,"");
            break;
        case "y":
        case "Y":
            sprintf(strname,"%s_py_%i_%i",hname,chanLo,chanHi);
            h1D = (TH1D*)h2D->ProjectionY(strname,chanLo,chanHi,"");
            break;
        default:
            cout << "Incorrect projection axis " << projAxis << endl;
            exit(0);
            break;
    }
	h1D->SetTitle(title);
	h1D->GetXaxis()->CenterTitle();
	h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
	h1D->GetYaxis()->SetTitleOffset(yoff);
    h1D->SetLineWidth(2);
    h1D->Draw();

	sprintf(OutCan,"dmsOmega_%s_%s.gif",tgtOut,strname);
	c1->Print(OutCan);
	sprintf(OutCan,"dmsOmega_%s_%s.eps",tgtOut,strname);
	c1->Print(OutCan);
}

