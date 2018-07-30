// PlotOmega_dmsOmega.C
//
// macro to project 1D eg2a histograms from 2D histograms
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// gROOT->Reset();   // start from scratch

Int_t lcol[10] = {1,2,4,6,7,8,9,13,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

const Int_t MAX_HIST = 3;
const Int_t MAX_RUN = 4;
const Int_t MAX_TGT = 3;
const Int_t MAX_CUT = 21;

void Check_HistIndex(Int_t index);
void PrintHistIndex();
void Check_TgtIndex(Int_t index);
void PrintTgtIndex();
void Check_CutIndex(Int_t index);
void PrintCutIndex();
void Check_CutLoHi(Int_t Lo, Int_t Hi);

class Cuts
{
    vector<string> LabelCuts;
    vector<string> LegHeader;
    vector<string> HistName;
    vector<string> RunName;
    vector<string> TgtName;
    
public:
    Cuts();
    Int_t Get_nCuts() {return LabelCuts.size();};
    string Get_Cuts(int num) {return LabelCuts[num];};

    Int_t Get_nLegend() {return LegHeader.size();};
    string Get_Legend(int num) {return LegHeader[num];};

    Int_t Get_nHist() {return HistName.size();};
    string Get_Hist(int num) {return HistName[num];};

    Int_t Get_nRun() {return RunName.size();};
    string Get_Run(int num) {return RunName[num];};

    Int_t Get_nTgt() {return TgtName.size();};
    string Get_Tgt(int num) {return TgtName[num];};
};

Cuts::Cuts()
{
    LabelCuts.push_back("None");
    LabelCuts.push_back("All #omega cuts");
    LabelCuts.push_back("M(#pi^{0})");
    LabelCuts.push_back("Q^{2}");
    LabelCuts.push_back("W");
    LabelCuts.push_back("V_{z} Matching");
    LabelCuts.push_back("#theta_{e-,#gamma}");
    LabelCuts.push_back("M(#pi^{+}#pi^{-}");
    LabelCuts.push_back("Part. Topology");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 1");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 2");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 3");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 1");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 2");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 3");
    LabelCuts.push_back("Photon TOF M^{2}");
    LabelCuts.push_back("Dalitz 1");
    LabelCuts.push_back("Proton In Event");
    LabelCuts.push_back("Proton-In-Event/All");
    LabelCuts.push_back("Evt. Particle Comb.");
    LabelCuts.push_back("Evt. Particle Comb./All");

    LegHeader.push_back("Cuts: ");
    LegHeader.push_back("All Cuts Except:");
    LegHeader.push_back("Anti-Cuts:");
    
    HistName.push_back("IMOmega_");
    HistName.push_back("IMOmega_woCut_");
    HistName.push_back("IMOmega_antiCut_");

    RunName.push_back("C12");
    RunName.push_back("Fe56");
    RunName.push_back("Sn");
    RunName.push_back("Pb208");

    TgtName.push_back("NoTarget");
    TgtName.push_back("LD2");
    TgtName.push_back("Nuc");
}

// 
// PlotOmega_CutIndex - plot omega inv. mass for a specific cut selection
//                  
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//                  chanLo = lower bin
//                  chanHi = upper bin
//
void PlotOmega_CutIndex(string fAna, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t runIndex=0, Int_t chanLo = 0, Int_t chanHi=0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char legLabel[50];
    
    Int_t iColor = 0;
    
    TH1D *h1D[MAX_CUT];
    
    Cuts myCuts;
    
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
	printf("Analyzing file %s\n",fAna.c_str());
	TFile *fm = new TFile(fAna.c_str(),"READ");
    TDirectory *tmp = fm->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
	c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D->ProjectionX(strname,2,2,"");
    
    sprintf(title,"Target: %s",myCuts.Get_Tgt(tgtIndex).c_str());
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();

    sprintf(legLabel,"%s",myCuts.Get_Cuts(1).c_str());
    leg->AddEntry(h1D_noCut,legLabel,"l");

    iColor++;
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    
    for(i=chanLo; i<chanHi+1; i++){
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(lcol[iColor]);
        h1D[i]->Draw("same");
        
        sprintf(legLabel,"%s",myCuts.Get_Cuts(i).c_str());
        leg->AddEntry(h1D[i],legLabel,"l");
        
        iColor++;
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(myCuts.Get_Legend(histIndex).c_str());
    leg->Draw();

	sprintf(OutCan,"PlotOmega_%s_%s_%i_%i.gif",hname,myCuts.Get_Run(runIndex).c_str(),chanLo,chanHi);
	c1->Print(OutCan);
	sprintf(OutCan,"PlotOmega_%s_%s_%i_%i.eps",hname,myCuts.Get_Run(runIndex).c_str(),chanLo,chanHi);
	c1->Print(OutCan);
}

//
// OverlayOmega_CutIndex - overlay 2 projections of the omega inv. mass for a specific cut selection
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//                  chan1 = projection 1
//                  chan2 = projection 2
//
void OverlayOmega_CutIndex(string fAna, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan1 = 0, Int_t chan2=0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char legLabel[50];
    
    Int_t iColor = 0;
    
    TH1D *h1D[2];
    
    Cuts myCuts;

    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chan1);
    Check_CutIndex(chan2);
    Check_CutLoHi(chan1,chan2);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    TDirectory *tmp = fm->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,chan1);
    h1D[0] = (TH1D*)h2D->ProjectionX(strname,chan1+1,chan1+1,"");
    h1D[0]->SetLineWidth(2);
    h1D[0]->SetLineColor(lcol[0]);
    h1D[0]->SetTitle(title);
    h1D[0]->GetXaxis()->CenterTitle();
    h1D[0]->GetYaxis()->CenterTitle();
    h1D[0]->GetYaxis()->SetTitle("Counts");
    h1D[0]->GetYaxis()->SetTitleOffset(yoff);
    h1D[0]->Draw();
        
    sprintf(legLabel,"%s",myCuts.Get_Cuts(chan1).c_str());
    leg->AddEntry(h1D[0],legLabel,"l");
    
    sprintf(strname,"%s_%i",hname,chan2);
    h1D[1] = (TH1D*)h2D->ProjectionX(strname,chan2+1,chan2+1,"");
    h1D[1]->SetLineWidth(2);
    h1D[1]->SetLineColor(lcol[1]);
    h1D[1]->Draw("same");
    
    sprintf(legLabel,"%s",myCuts.Get_Cuts(chan2).c_str());
    leg->AddEntry(h1D[1],legLabel,"l");
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(myCuts.Get_Legend(histIndex).c_str());
    leg->Draw();
    
    sprintf(OutCan,"OverlayOmega_%s_%i_%i.gif",hname,chan1,chan2);
    c1->Print(OutCan);
    sprintf(OutCan,"OverlayOmega_%s_%i_%i.eps",hname,chan1,chan2);
    c1->Print(OutCan);
}

//
// PlotOmega_AllGrid - plot omega inv. mass for all cut selections
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//
void PlotOmega_AllGrid(string fAna, Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    
    TH1D *h1D[MAX_CUT];
 
    Cuts myCuts;
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    int nrow = ceil(sqrt(MAX_CUT));
    c1->Divide(nrow,nrow);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    TDirectory *tmp = fm->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    for(i=0; i<MAX_CUT; i++){
        c1->cd(i+1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        
        sprintf(title,"Target: %s, Cut: %s",myCuts.Get_Tgt(tgtIndex).c_str(),myCuts.Get_Cuts(i).c_str());
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
void PlotOmega_CutVsAntiCut(string fAna, Int_t tgtIndex = 0, Int_t chan = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char legLabel[50];
    
    Int_t iColor = 0;
    
    TH1D *h1D_Cut[MAX_CUT];
    TH1D *h1D_AntiCut[MAX_CUT];
    
    Cuts myCuts;
    
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chan);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    TDirectory *tmp = fm->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(0).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D_Cut = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D_Cut->ProjectionX(strname,2,2,"");
    
    sprintf(title,"Target: %s",myCuts.Get_Tgt(tgtIndex).c_str());
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();
    
    sprintf(legLabel,"%s",myCuts.Get_Cuts(1).c_str());
    leg->AddEntry(h1D_noCut,legLabel,"l");
    
    sprintf(strname,"%s_%i",hname,chan);
    h1D_Cut[i] = (TH1D*)h2D_Cut->ProjectionX(strname,chan+1,chan+1,"");
    h1D_Cut[i]->SetLineWidth(2);
    h1D_Cut[i]->SetLineColor(2);
    h1D_Cut[i]->Draw("same");
        
    sprintf(legLabel,"%s",myCuts.Get_Cuts(chan).c_str());
    leg->AddEntry(h1D_Cut[i],legLabel,"l");

    sprintf(hname,"%s%s",myCuts.Get_Hist(2).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D_AntiCut = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,chan);
    h1D_AntiCut[i] = (TH1D*)h2D_AntiCut->ProjectionX(strname,chan+1,chan+1,"");
    h1D_AntiCut[i]->SetLineWidth(2);
    h1D_AntiCut[i]->SetLineColor(4);
    h1D_AntiCut[i]->Draw("same");
    
    sprintf(legLabel,"%s (anti)",myCuts.Get_Cuts(chan).c_str());
    leg->AddEntry(h1D_AntiCut[i],legLabel,"l");
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(myCuts.Get_Legend(0).c_str());
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
void CompOmega_Files(string fAna1, string fAna2, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0, string legLine1 = "Leg 1", string legLine2 = "Leg 2", string comment="test")
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
    
    Cuts myCuts;
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(chan);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);

    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    for(i=0; i<2; i++){
        // data files contain the trees
        switch(i){
            case 0:
                sprintf(legLabel,"%s",legLine1.c_str());
                printf("Analyzing file %s\n",fAna1.c_str());
                fm[i] = new TFile(fAna1.c_str(),"READ");
                break;
            case 1:
                sprintf(legLabel,"%s",legLine2.c_str());
                printf("Analyzing file %s\n",fAna2.c_str());
                fm[i] = new TFile(fAna2.c_str(),"READ");
                break;
            default: sprintf(legLabel,"Insert File Name"); break;
        }
        
        dir[i] = fm[i]->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
        
        h2D[i] = (TH2D*)dir[i]->Get(hname);
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D[i]->ProjectionX(strname,chan+1,chan+1,"");
        
        sprintf(title,"Target: %s, Cuts: %s",myCuts.Get_Tgt(tgtIndex).c_str(),myCuts.Get_Cuts(chan).c_str());
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(lcol[i]);
        (i==0) ? h1D[i]->Draw() : h1D[i]->Draw("same");
        
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Files:");
    leg->Draw();
    
    sprintf(OutCan,"CompOmega_%s_%i_%s.gif",hname,chan,comment.c_str());
    c1->Print(OutCan);
    sprintf(OutCan,"CompOmega_%s_%i_%s.eps",hname,chan,comment.c_str());
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
    Cuts myCuts;
    
    cout<<"Histogram Index:"<<endl;
    for(i=0;i<MAX_HIST;i++){
        cout<<i<<"\t"<<myCuts.Get_Hist(i)<<endl;
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
    Cuts myCuts;
    
    cout<<"Target Index:"<<endl;
    for(i=0;i<MAX_TGT;i++){
        cout<<i<<"\t"<<myCuts.Get_Tgt(i)<<endl;
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

    Cuts myCuts;
    
    cout<<"Cut Index:"<<endl;
    for(i=0;i<MAX_CUT;i++){
        cout<<i<<"\t"<<myCuts.Get_Cuts(i)<<endl;
    }
}

void Check_CutLoHi(Int_t Lo, Int_t Hi)
{
    if(Hi<Lo){
        cout<<"Lower cut index "<<Lo<<" is greater than upper cut index "<<Hi<<endl;
        exit(0);
    }
}
