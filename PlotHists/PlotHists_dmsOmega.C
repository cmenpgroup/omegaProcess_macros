// PlotHists_dmsOmega.C
//
// macro to plot eg2a histograms
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
// PlotHists_dmsOmega - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  hname = histogram name
//                  target = target name
//
void PlotHists_dmsOmega(char *fAna, char *fDir, char *hname, char *target, char *tgtOut, int iDim = 1)
{
	char OutCan[100];
    
	// data files contain the trees

    printf("Analyzing file %s:/%s\n",fAna,fDir);
    
    
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(fDir);
    
	TH1F *hist = (TH1F*)tmp->Get(hname);
	hist->SetTitle(target);
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleOffset(yoff);

    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    switch(iDim){
        case 1: hist->SetLineWidth(2); hist->Draw(); break;
        case 2:
//            gPad->SetLogz();
            hist->Draw("colz");
            break;
        default:
            cout << "Incorrect dimension for histogram." << endl;
            exit(0);
            break;
    }

	sprintf(OutCan,"dmsOmega_%s_%s.gif",tgtOut,hname);
	c1->Print(OutCan);
	sprintf(OutCan,"dmsOmega_%s_%s.eps",tgtOut,hname);
	c1->Print(OutCan);
}

