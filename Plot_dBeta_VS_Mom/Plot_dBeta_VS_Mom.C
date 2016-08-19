// Plot_dBeta_VS_Mom.C
//
// macro to plot eg2a histograms for beta difference vs momentum
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

const Int_t NPART = 5;
const Int_t MAX_HIST = 2;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};
char *Plabel[NPART] = {"Electron","Pi-","Pi+","Photon1","Photon2"};
char *HistName[MAX_HIST] = {"dBeta_VS_Momentum","dBeta_VS_Momentum_EPC"};
char *BeforeAfter[MAX_HIST] = {"Before","After"};

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
// Plot_dBeta_VS_Mom - plot histogram of beta difference vs momentum
//
//                  fAna = output from eg2a DMS
//                  histIndex = index for histogram name
//                  target = target name
//
//
void Plot_dBeta_VS_Mom(char *fAna="Ana.root", Int_t histIndex, char *target)
{
    Int_t i;
    
    TH2D *hist[NPART];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(3,2);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
    
    for(i=0; i<NPART; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        c1->cd(i+1);
        sprintf(hname,"%s_%s",HistName[histIndex],Plabel[i]);
        cout<<hname<<endl;
        hist[i] = (TH2D*)tmp->Get(hname);
        sprintf(htitle,"%s, %s",Plabel[i],target);
        hist[i]->SetTitle(htitle);
        hist[i]->GetXaxis()->CenterTitle();
        hist[i]->GetYaxis()->CenterTitle();
        hist[i]->GetYaxis()->SetTitleOffset(yoff);
        hist[i]->Draw("colz");

    }
    
    sprintf(OutCan,"Plot_%s_%s.gif",HistName[histIndex],target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_%s_%s.eps",HistName[histIndex],target);
    c1->Print(OutCan);
}

// BeforeAfter_dBeta_VS_Mom_Particle - plot before and after the particle ID cuts
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
//
void BeforeAfter_dBeta_VS_Mom_Particle(char *fAna="Ana.root", char *target)
{
    Int_t i, j;
    
    TH2D *hist[NPART][MAX_HIST];
    TCanvas *can[NPART];

    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
    
    for(j=0; j<NPART; j++){
        // Canvas to plot histogram
        sprintf(cname,"can%i",j);
        sprintf(ctitle,"%s Canvas",Plabel[j]);
        can[j] = new TCanvas(cname,ctitle,50*j,50*j,800,400);
        can[j]->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        can[j]->SetBorderSize(5);
        gStyle->SetOptStat(1111);
        can[j]->SetFillStyle(4000);
        can[j]->Divide(2,1);
    
        for(i=0; i<MAX_HIST; i++){
            gPad->SetLeftMargin(Lmar);
            gPad->SetRightMargin(Rmar);
            gPad->SetFillColor(0);
        
            can[j]->cd(i+1);
            sprintf(hname,"%s_%s",HistName[i],Plabel[j]);
        
            hist[j][i] = (TH2D*)tmp->Get(hname);
            sprintf(htitle,"%s, %s, %s ID cuts",Plabel[j],target,BeforeAfter[i]);
            hist[j][i]->SetTitle(htitle);
            hist[j][i]->GetXaxis()->CenterTitle();
            hist[j][i]->GetYaxis()->CenterTitle();
            hist[j][i]->GetYaxis()->SetTitleOffset(yoff);
            hist[j][i]->Draw("colz");
        }
        sprintf(OutCan,"BeforeAfter_dBvsP_%s_%s.gif",Plabel[j],target);
        can[j]->Print(OutCan);
        sprintf(OutCan,"BeforeAfter_dBvsP_%s_%s.eps",Plabel[j],target);
        can[j]->Print(OutCan);
    }
}
