// PlotOmega_dmsOmega.C
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

Int_t lcol[10] = {1,2,4,6,7,8,9,13,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

const Int_t MAX_HIST = 3;
char *HistName[MAX_HIST] = {"IMOmega_","IMOmega_woCut_","IMOmega_antiCut_"};
char *legHeader[MAX_HIST] = {"Cuts: ","All Cuts Except:","Anti-Cuts:"};

const Int_t MAX_RUN = 4;
char *RunName[MAX_RUN] = {"C12","Fe56","Sn","Pb208"};

const Int_t MAX_TGT = 3;
char *TgtName[MAX_TGT] = {"NoTarget","LD2","Nuc"};

const Int_t MAX_CUT = 17;
char *CutName[MAX_CUT] = {"None","All #omega cuts","M(#pi^{0})","Q^{2}","W","V_{z} Matching","Part. Topology","#theta_{e-,#gamma}","M(#pi^{+}#pi^{-}","EC 2nd Moment for #gamma's,  Region 1","EC 2nd Moment for #gamma's, Region 2","EC 2nd Moment for #gamma's, Region 3","EC 3rd Moment for #gamma's, Region 1","EC 3rd Moment for #gamma's, Region 2","EC 3rd Moment for #gamma's, Region 3","Photon TOF M^{2}","Dalitz 1"};

// 
// PlotOmega_CutIndex - plot omega inv. mass for a specific cut selection
//                  
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//                  chanLo = lower bin
//                  chanHi = upper bin
//
void PlotOmega_CutIndex(char *fAna, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chanLo = 0, Int_t chanHi=0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char strname[100];
    char legLabel[50];
    
    Int_t iColor = 0;
    
    TH1D *h1D[MAX_CUT];
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chanLo);
    Check_CutIndex(chanHi);
    Check_CutLoHi(chanLo,chanHi);
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
	
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(TgtName[tgtIndex]);
    
	c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",HistName[histIndex],TgtName[tgtIndex]);
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D->ProjectionX(strname,2,2,"");
    
    sprintf(title,"Target: %s",TgtName[tgtIndex]);
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();
    
    sprintf(legLabel,"%s",CutName[1]);
    leg->AddEntry(h1D_noCut,legLabel,"l");
    
    iColor++;
    
    sprintf(hname,"%s%s",HistName[histIndex],TgtName[tgtIndex]);
    
    for(i=chanLo; i<chanHi+1; i++){
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(lcol[iColor]);
        h1D[i]->Draw("same");
        
        sprintf(legLabel,"%s",CutName[i]);
        leg->AddEntry(h1D[i],legLabel,"l");
        
        iColor++;
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(legHeader[histIndex]);
    leg->Draw();
    
	sprintf(OutCan,"PlotOmega_%s_%i_%i.gif",hname,chanLo,chanHi);
	c1->Print(OutCan);
	sprintf(OutCan,"PlotOmega_%s_%i_%i.eps",hname,chanLo,chanHi);
	c1->Print(OutCan);
}

//
// PlotOmega_AllGrid - plot omega inv. mass for all cut selections
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//
void PlotOmega_AllGrid(char *fAna, Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char strname[100];
    
    TH1D *h1D[MAX_CUT];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    int nrow = ceil(sqrt(MAX_CUT));
    c1->Divide(nrow,nrow);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(TgtName[tgtIndex]);
    
    sprintf(hname,"%s%s",HistName[histIndex],TgtName[tgtIndex]);
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    for(i=0; i<MAX_CUT; i++){
        c1->cd(i+1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        
        sprintf(title,"Target: %s, Cut: %s",TgtName[tgtIndex],CutName[i]);
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->Draw();
    }
    
    sprintf(OutCan,"PlotOmega_AllGrid_%s.gif",hname);
    c1->Print(OutCan);
    sprintf(OutCan,"PlotOmega_AllGrid_%s.eps",hname);
    c1->Print(OutCan);
}

//
// PlotOmega_CutVSAntiCut - plot omega inv. mass for a specific cut selection
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//                  chanLo = lower bin
//                  chanHi = upper bin
//
void PlotOmega_CutVsAntiCut(char *fAna, Int_t tgtIndex = 0, Int_t chan = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char strname[100];
    char legLabel[50];
    
    Int_t iColor = 0;
    
    TH1D *h1D_Cut[MAX_CUT];
    TH1D *h1D_AntiCut[MAX_CUT];
    
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chan);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(TgtName[tgtIndex]);
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",HistName[0],TgtName[tgtIndex]);
    TH2D *h2D_Cut = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D_Cut->ProjectionX(strname,2,2,"");
    
    sprintf(title,"Target: %s",TgtName[tgtIndex]);
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();
    
    sprintf(legLabel,"%s",CutName[1]);
    leg->AddEntry(h1D_noCut,legLabel,"l");
    
    sprintf(strname,"%s_%i",hname,chan);
    h1D_Cut[i] = (TH1D*)h2D_Cut->ProjectionX(strname,chan+1,chan+1,"");
    h1D_Cut[i]->SetLineWidth(2);
    h1D_Cut[i]->SetLineColor(2);
    h1D_Cut[i]->Draw("same");
        
    sprintf(legLabel,"%s",CutName[chan]);
    leg->AddEntry(h1D_Cut[i],legLabel,"l");

    sprintf(hname,"%s%s",HistName[2],TgtName[tgtIndex]);
    TH2D *h2D_AntiCut = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,chan);
    h1D_AntiCut[i] = (TH1D*)h2D_AntiCut->ProjectionX(strname,chan+1,chan+1,"");
    h1D_AntiCut[i]->SetLineWidth(2);
    h1D_AntiCut[i]->SetLineColor(4);
    h1D_AntiCut[i]->Draw("same");
    
    sprintf(legLabel,"%s (anti)",CutName[chan]);
    leg->AddEntry(h1D_AntiCut[i],legLabel,"l");
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(legHeader[0]);
    leg->Draw();
    
    sprintf(OutCan,"PlotOmega_CutVsAntiCut_%i.gif",chan);
    c1->Print(OutCan);
    sprintf(OutCan,"PlotOmega_CutVsAntiCut_%i.eps",chan);
    c1->Print(OutCan);
}

//
// CompOmega_Files - compare omega inv. mass for a specific cut selection between files
//
//                  fAna1 = output 1 from eg2a DMS
//                  fAna2 = output 2 from eg2a DMS
//                  histIndex = 2-D histogram name
//                  tgtIndex = target index
//                  chan = cut bin
//
void CompOmega_Files(char *fAna1, char *fAna2, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0, char *legLine1 = "Leg 1", char *legLine2 = "Leg 2", char *comment="test")
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char legLabel[50];
    
    TFile *fm[2];
    TDirectory *dir[2];
    TH1D *h1D[2];
    TH2D *h2D[2];
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chan);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);

    sprintf(hname,"%s%s",HistName[histIndex],TgtName[tgtIndex]);
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    for(i=0; i<2; i++){
        // data files contain the trees
        switch(i){
            case 0:
                sprintf(legLabel,"%s",legLine1);
                printf("Analyzing file %s\n",fAna1);
                fm[i] = new TFile(fAna1,"READ");
                break;
            case 1:
                sprintf(legLabel,"%s",legLine2);
                printf("Analyzing file %s\n",fAna2);
                fm[i] = new TFile(fAna2,"READ");
                break;
            default: sprintf(legLabel,"Insert File Name"); break;
        }
        
        dir[i] = fm[i]->GetDirectory(TgtName[tgtIndex]);
        
        h2D[i] = (TH2D*)dir[i]->Get(hname);
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D[i]->ProjectionX(strname,chan+1,chan+1,"");
        
        sprintf(title,"Target: %s, Cuts: %s",TgtName[tgtIndex],CutName[chan]);
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(lcol[i]);
        h1D[i]->Draw(fSame[i]);
        
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Files:");
    leg->Draw();
    
    sprintf(OutCan,"CompOmega_%s_%i_%s.gif",hname,chan,comment);
    c1->Print(OutCan);
    sprintf(OutCan,"CompOmega_%s_%i_%s.eps",hname,chan,comment);
    c1->Print(OutCan);
}

void Check_HistIndex(Int_t index){
    if(index<0 || index>=MAX_HIST){
        cout<<"Histogram index "<<index<<" is out of range [0,"<<MAX_HIST-1<<"]"<<endl;
        PrintHistIndex();
        exit(0);
    }
}

void PrintHistIndex()
{
    Int_t i;
    cout<<"Histogram Index:"<<endl;
    for(i=0;i<MAX_HIST;i++){
        cout<<i<<"\t"<<HistName[i]<<endl;
    }
}

void Check_TgtIndex(Int_t index)
{
    if(index<0 || index>=MAX_TGT){
        cout<<"Target index "<<index<<" is out of range [0,"<<MAX_TGT-1<<"]"<<endl;
        PrintTgtIndex();
        exit(0);
    }
}

void PrintTgtIndex()
{
    Int_t i;
    cout<<"Target Index:"<<endl;
    for(i=0;i<MAX_TGT;i++){
        cout<<i<<"\t"<<TgtName[i]<<endl;
    }
}

void Check_CutIndex(Int_t index)
{
    if(index<0 || index>=MAX_CUT){
        cout<<"Cut index "<<index<<" is out of range [0,"<<MAX_CUT-1<<"]"<<endl;
        PrintCutIndex();
        exit(0);
    }
}

void PrintCutIndex()
{
    Int_t i;
    cout<<"Cut Index:"<<endl;
    for(i=0;i<MAX_CUT;i++){
        cout<<i<<"\t"<<CutName[i]<<endl;
    }
}

void Check_CutLoHi(Int_t Lo, Int_t Hi)
{
    if(Hi<Lo){
        cout<<"Lower cut index "<<Lo<<" is greater than upper cut index "<<Hi<<endl;
        exit(0);
    }
}
