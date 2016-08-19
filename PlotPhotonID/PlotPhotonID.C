// PlotPhotonID.C
//
// macro to plot eg2a histograms for PhotonID
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
// PlotPhotonID - plot histogram with labels for photons
//
//                  fAna = output from eg2a DMS
//                  HistName = histogram name
//                  nDim = histogram dimension
//                  target = target name
//
//
void PlotPhotonID(char *fAna="Ana.root", char *HistName, Int_t nDim, char *target)
{
    Int_t i;
    char drawType[50];
    
    switch (nDim){
        case 1: TH1D *hist[NPHOTONS]; break;
        case 2: TH2D *hist[NPHOTONS]; break;
        default: cout<<"Wrong histogram dimension: "<<nDim<<endl; exit(0); break;
    }
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,400);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(2,1);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("PhotonID");
    
    for(i=0; i<NPHOTONS; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        c1->cd(i+1);
        sprintf(hname,"%s%i",HistName,i+1);

        switch (nDim){
            case 1: hist[i] = (TH1D*)tmp->Get(hname); sprintf(drawType,""); break;
            case 2: hist[i] = (TH2D*)tmp->Get(hname); sprintf(drawType,"colz"); break;
            default: cout<<"Wrong histogram dimension: "<<nDim<<endl; exit(0); break;
        }

        sprintf(htitle,"Photon %i, %s",i+1,target);
        hist[i]->SetTitle(htitle);
        hist[i]->GetXaxis()->CenterTitle();
        hist[i]->GetYaxis()->CenterTitle();
        hist[i]->GetYaxis()->SetTitleOffset(yoff);
        hist[i]->Draw(drawType);

    }
    
    sprintf(OutCan,"Plot_%s_%s.gif",HistName,target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_%s_%s.eps",HistName,target);
    c1->Print(OutCan);
}

//
// OverlayPhotonIDcuts - overlay histograms with before and after cuts
//
//                  fAna = output from eg2a DMS
//                  HistName = histogram name
//                  target = target name
//
//
void OverlayPhotonIDcuts(char *fAna="Ana.root", char *HistName, char *target)
{
    Int_t i;
    
    TH1D *hist[NPHOTONS];
    TH2D *histCuts[NPHOTONS];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,400);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(2,1);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("PhotonID");
    
    for(i=0; i<NPHOTONS; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        c1->cd(i+1);
        sprintf(hname,"%s%i",HistName,i+1);
        hist[i] = (TH1D*)tmp->Get(hname);
        
        sprintf(htitle,"Photon %i, %s",i+1,target);
        hist[i]->SetTitle(htitle);
        hist[i]->GetXaxis()->CenterTitle();
        hist[i]->GetYaxis()->CenterTitle();
        hist[i]->GetYaxis()->SetTitleOffset(yoff);
        hist[i]->SetLineWidth(2);
        hist[i]->Draw();
        
        sprintf(hname,"%s%i_cut",HistName,i+1);
        histCuts[i] = (TH2D*)tmp->Get(hname);
        histCuts[i]->SetLineColor(2);
        histCuts[i]->SetLineWidth(2);
        histCuts[i]->SetFillColor(2);
        histCuts[i]->SetFillStyle(3002);
        histCuts[i]->Draw("same");
    }
    
    sprintf(OutCan,"OL_PhotonIDcuts_%s_%s.gif",HistName,target);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_PhotonIDcuts_%s_%s.eps",HistName,target);
    c1->Print(OutCan);
}

//
// OverlayEC_xy_Photons_1Sector - overlay histogram of EC x vs y for one sector
//
//                  fAna = output from eg2a DMS
//                  iSector = sector number
//                  target = target name
//
void OverlayEC_xy_Photons_1Sector(char *fAna="Ana.root", int iSector,char *target)
{
    Int_t i;
    TH2D *hist1[NPHOTONS];
    TH2D *hist2[NPHOTONS];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,400);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(2,1);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("PhotonID");
    
    for(i=0; i<NPHOTONS; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        c1->cd(i+1);
        sprintf(hname,"EC_XvsY_local_FidCut%i_Photon%i",iSector,i+1);
        hist1[i] = (TH2D*)tmp->Get(hname);
        sprintf(htitle,"Sector %i, %s",iSector,target);
        hist1[i]->SetTitle(htitle);
        hist1[i]->GetXaxis()->CenterTitle();
        hist1[i]->GetYaxis()->CenterTitle();
        hist1[i]->GetYaxis()->SetTitleOffset(yoff);
        hist1[i]->Draw("colz");
    
        sprintf(hname,"EC_XvsY_local_AntiFidCut%i_Photon%i",iSector,i+1);
        hist2[i] = (TH2D*)tmp->Get(hname);
        hist2[i]->SetMarkerColor(1);
        hist2[i]->Draw("same");
    }
    
    sprintf(OutCan,"OL_EC_XvsY_Photons_FidCut_VS_AntiFidCut_S%i_%s.gif",iSector,target);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_EC_XvsY_Photons_FidCut_VS_AntiFidCut_S%i_%s.eps",iSector,target);
    c1->Print(OutCan);
}

