// 
// PlotMixedEvent - make plots of mixed event spectra
//
// M. H. Wood, Canisius College
//--------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();

Int_t lcol[10] = {1,2,4,6,7,8,9,12,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

char *legName[5] = {"Photon 1","Photon 2","#pi^{+}","#pi^{-}","#pi^{0}"};
char *RunName[4] = {"C12","Fe56","Sn","Pb208"};
char *TgtName[3] = {"NoTarget","LD2","Nuc"};
char *HistName[6] = {"IMOmega","IMOmega_OpAng_Electron_Cut","IMOmega_MassPi0Cut","IMOmega_ZVertCut","IMOmega_QsqCut","IMOmega_AllCuts"};
                     
void PlotMixedEvent(char *rootFile, Int_t iRun=0, Int_t iTgt=0, Int_t iCut=0, Int_t chanLo=1, Int_t chanHi=1)
{
	Int_t i;

    TH2D *hME2D;
    TH1D *hME1D;

    char title[100];
    char hnameME2D[50];
    char hnameME1D[50];
    
	char plotFilePrefix[100];
    
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
    
    sprintf(hnameME2D,"%s_ME_%s",HistName[iCut],TgtName[iTgt]);
	hME2D = (TH2D*)fd->Get(hnameME2D); // get the histogram from the ROOT file
    
    sprintf(hnameME1D,"%s-%i_%i",hnameME2D,chanLo,chanHi);
    hME1D = (TH1D*)hME2D->ProjectionX(hnameME1D,chanLo,chanHi,"");
    
    sprintf(plotFilePrefix,"PlotME_%s_%s",RunName[iRun],hnameME1D);
   
    sprintf(title,"#omega meson, Mixed Event Background");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    hME1D->GetYaxis()->SetTitleOffset(yoff);
    hME1D->Draw();
    
    char OutCan[100];
    sprintf(OutCan,"%s.gif",plotFilePrefix);
    can1->Print(OutCan);
    sprintf(OutCan,"%s.eps",plotFilePrefix);
    can1->Print(OutCan);
    
//	fd->Close();  // close ROOT file
}

void CompMixedEvent(char *rootFile, Int_t iRun=0, Int_t iTgt=0, Int_t iCut=0)
{
    Int_t i;
    
    TH2D *hME2D;
    TH1D *hME1D[5];
    
    char title[100];
    char hnameME2D[50];
    char hnameME1D[50];
    
    char plotFilePrefix[100];
    char legLabel[50];
    
    // data files contain the trees
    printf("Analyzing file %s\n",rootFile);
    TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
    
    sprintf(hnameME2D,"%s_ME_%s",HistName[iCut],TgtName[iTgt]);
    hME2D = (TH2D*)fd->Get(hnameME2D); // get the histogram from the ROOT file
 
    sprintf(plotFilePrefix,"CompME_%s_%s",RunName[iRun],hnameME2D);

    sprintf(title,"#omega meson, Mixed Event Background");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.55,0.5,0.95,0.875);
    
    for(i=0; i<5; i++){
        sprintf(hnameME1D,"%s-%i_%i",hnameME2D,i+1,i+1);
        hME1D[i] = (TH1D*)hME2D->ProjectionX(hnameME1D,i+1,i+1,"");
        hME1D[i]->GetYaxis()->SetTitleOffset(yoff);
        hME1D[i]->SetLineWidth(2);
        hME1D[i]->SetLineColor(lcol[i]);
        hME1D[i]->Draw(fSame[i]);
        
        hME1D[i]->GetYaxis()->SetTitle("Counts");
        hME1D[i]->GetYaxis()->CenterTitle();
        hME1D[i]->GetYaxis()->SetTitleOffset(yoff);
        
        sprintf(legLabel,"%s",legName[i]);
        leg->AddEntry(hME1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Mixing Particle");
    leg->Draw();
    
    char OutCan[100];
    sprintf(OutCan,"%s.gif",plotFilePrefix);
    can1->Print(OutCan);
    sprintf(OutCan,"%s.eps",plotFilePrefix);
    can1->Print(OutCan);
    
    //	fd->Close();  // close ROOT file
}




