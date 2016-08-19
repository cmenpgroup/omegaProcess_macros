// 
// Diff_TransMomSq.C - analyze histrograms of transverse momentum squared
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

Int_t const maxRun = 4;
Int_t const maxTgt = 3;
char *RunName[maxRun] = {"C12","Fe56","Sn","Pb208"};
char *TgtName[maxTgt] = {"NoTarget","LD2","Nuc"};
char *HistName[2] = {"PtSq_Omega_AllCuts","PtSq_Omega_AllCuts_IMOmegaCut"};
char *RootPrefix[1] = {"DMS_FirstChannel"};

void OverLay_TransMomSq(char *FilePath, Int_t iRun=0, Int_t iCut=0, Int_t iClose)
{
	Int_t i, j;
    Double_t Mean[2];
    Double_t wMean[2];
    Double_t Diff;
    
    char rootFile[100];
    char title[100];
    char legString[25];
    
    Int_t nTgt = maxTgt -1;
    
    TLine *lMean[2];
    TH1D *hist[2];
    char hname[50];
    
	char plotFilePrefix[100];
    
	// data files contain the trees
    sprintf(rootFile,"%s/%s_%s.root",FilePath, RootPrefix[0],RunName[iRun]);
	printf("Analyzing file %s\n",rootFile);
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	
    sprintf(title,"Overlay of Trans. Mom. Squared");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);

    TLegend *leg = new TLegend(0.6,0.7,1.0,0.875);
    
    for(i=0; i<nTgt; i++){
        sprintf(hname,"%s_%s",HistName[iCut],TgtName[i+1]);
        
        hist[i] = (TH1D*)fd->Get(hname); // get the histogram from the ROOT file

        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        wMean[i] = GetWeightedMean(hist[i]);
        Mean[i] = hist[i]->GetMean();
        cout<< i+1 <<"\t"<<wMean[i] <<"\t"<<Mean[i]<<endl;
        
        hist[i]->SetTitle(0);
        hist[i]->GetYaxis()->SetTitleOffset(yoff);
        hist[i]->SetLineColor(lcol[i]);
        hist[i]->SetLineWidth(2);
        hist[i]->Draw(fSame[i]);
    
        if(i==0){
            sprintf(legString,"%s",TgtName[1]);
        }else{
            sprintf(legString,"%s",RunName[iRun]);
        }
        leg->AddEntry(hist[i],legString,"l");

        lMean[i] = new TLine(wMean[i],0.0,wMean[i],hist[i]->GetMaximum());
        lMean[i]->SetLineWidth(2);
        lMean[i]->SetLineColor(lcol[i+2]);
        lMean[i]->Draw();
    }

    cout << "weighted: "<<wMean[0]-wMean[1]<<endl;
    cout << "GetMean: "<<Mean[0]-Mean[1]<<endl;
    
    leg->SetLineColor(0);
	leg->SetFillStyle(0);
	leg->SetHeader("Targets:");
	leg->Draw();
    
    char plotFilePrefix[100];
    sprintf(plotFilePrefix,"OverLay_%s_%s",HistName[iCut],RunName[iRun]);
    
    char OutCan[100];
    sprintf(OutCan,"%s.gif",plotFilePrefix);
    can1->Print(OutCan);
    sprintf(OutCan,"%s.eps",plotFilePrefix);
    can1->Print(OutCan);

    if(iClose){
        can1->Close();
        fd->Close(); // close the root file
    }
}

void Diff_TransMomSq_AllTargets(char *FilePath, Int_t iCut=1)
{
	Int_t i, j;

    for(i=0; i<maxRun; i++){
        OverLay_TransMomSq(FilePth, i, iCut, 1);
    }
}

Double_t GetWeightedMean(TH1D *hist){

    Int_t i;
    Double_t ret = 0.0;
    Double_t numerator = 0.0;
    Double_t denominator = hist->Integral();
    
    for(i=0; i<hist->GetNbinsX(); i++){
        numerator += hist->GetBinCenter(i)*hist->GetBinContent(i);
    }
    if(denominator) ret = numerator/denominator;
    
    return ret;
}
