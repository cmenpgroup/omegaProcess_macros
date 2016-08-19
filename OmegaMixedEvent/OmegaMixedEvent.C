// 
// OmegaMixedEvent - background subtraction routine by scaling the mixed event spectra and
//                   subtracting it from the actual spectra
//
//            The output is a list of yields.
//
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

char *RunName[4] = {"C12","Fe56","Sn","Pb208"};
char *TgtName[3] = {"NoTarget","LD2","Nuc"};
char *HistName[6] = {"IMOmega","IMOmega_OpAng_Electron_Cut","IMOmega_MassPi0Cut","IMOmega_ZVertCut","IMOmega_QsqCut","IMOmega_AllCuts"};

Float_t FitOmega_MixedEvtBgd(TH1D *hist, TH1D *histME1, TH1D *histME2, char *plotFilePrefix, Int_t iClose)
{
	Float_t PeakLo = 0.7;
	Float_t PeakHi = 0.875;
    
	Int_t i, j, k;
	Float_t yield;
    Double_t sampleData, sampleME;
    
    char hname[50];
    
    char title[100];
	char xtitle[100];
	char ytitle[100];
    
    sprintf(hname,"hAfterBgSub");
    TH1D *hAfterBgSub = (TH1D*)hist->Clone(hname);
    
    sprintf(hname,"hBgdME1");
    TH1D *hBgdME1 = (TH1D*)histME1->Clone(hname);

    sprintf(hname,"hBgdME2");
    TH1D *hBgdME2 = (TH1D*)histME2->Clone(hname);
    
    hBgdME1->Add(hBgdME2,1.0);

    sampleData=hAfterBgSub->Integral(hAfterBgSub->FindBin(1.0),hAfterBgSub->FindBin(2.0));
    sampleME=hBgdME1->Integral(hBgdME1->FindBin(1.0),hBgdME1->FindBin(2.0));
    
    if(sampleME) hBgdME1->Scale(sampleData/sampleME);
    hAfterBgSub->Add(hBgdME1,-1.0);
    
    yield=hAfterBgSub->Integral(hAfterBgSub->FindBin(PeakLo),hAfterBgSub->FindBin(PeakHi));

    sprintf(title,"#omega meson, Mixed Event Background Subtraction, Canvas 1");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);

    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    hist->GetYaxis()->SetTitleOffset(yoff);
    hist->Draw();
    
    hBgdME1->SetLineColor(2);
    hBgdME1->Draw("same");
    
    hAfterBgSub->SetFillColor(4);
    hAfterBgSub->Draw("same");
    
    TLine *lLow = new TLine(PeakLo,0.0,PeakLo,0.5*hAfterBgSub->GetMaximum());
    lLow->SetLineWidth(2);
    lLow->Draw();
    TLine *lHigh = new TLine(PeakHi,0.0,PeakHi,0.5*hAfterBgSub->GetMaximum());
    lHigh->SetLineWidth(2);
    lHigh->Draw();
    
    
    char OutCan[100];
    sprintf(OutCan,"%s.gif",plotFilePrefix);
    can1->Print(OutCan);
    sprintf(OutCan,"%s.eps",plotFilePrefix);
    can1->Print(OutCan);

    if(iClose){
        can1->Close();
        can2->Close();
    }
	return yield;
}

void OmegaMixedEvent_OneTarget(char *rootFile, Int_t iRun=0, Int_t iTgt=0, Int_t iCut=0, Int_t iClose=1)
{
	Int_t i;
	Float_t yield;

	TH1D *hData; // original histogram
    TH2D *hData2D;
    TH2D *hME2D;
    TH1D *hME1;
    TH1D *hME2;
    
    char hnameData[50];
    char hnameData2D[50];
    char hnameME2D[50];
    char hnameME1[50];
    char hnameME2[50];
    
	char plotFilePrefix[100];
    
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	
    sprintf(hnameData2D,"IMOmega_%s",TgtName[iTgt]);
    cout << hnameData2D << endl;
	hData2D = (TH2D*)fd->Get(hnameData2D); // get the histogram from the ROOT file

    sprintf(hnameData,"%s-8",hnameData2D);
    cout << hnameData << endl;
    hData = (TH1D*)hData2D->ProjectionX(hnameData,9,9,"");
    
    sprintf(hnameME2D,"%s_ME_%s",HistName[iCut],TgtName[iTgt]);
    cout << hnameME2D << endl;
	hME2D = (TH2D*)fd->Get(hnameME2D); // get the histogram from the ROOT file
    
    sprintf(hnameME1,"%s-1",hnameME2D);
    cout << hnameME1 << endl;
    hME1 = (TH1D*)hME2D->ProjectionX(hnameME1,1,1,"");
    
    sprintf(hnameME2,"%s-2",hnameME2D);
    cout << hnameME2 << endl;
    hME2 = (TH1D*)hME2D->ProjectionX(hnameME2,2,2,"");
    
    sprintf(plotFilePrefix,"OmegaMixedEvent_%s_%s_P%i",RunName[iRun],hnameData2D);
    
    yield = FitOmega_MixedEvtBgd(hData,hME1,hME2,plotFilePrefix, iClose);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout<<RunName[iRun]<<"\t"<<TgtName[iTgt]<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

	fout.close(); // close the text file
//	fd->Close();  // close ROOT file
}

void FitOmegaMixedEvent_AllTargets(char *rootFile, Int_t iRun=0, Int_t iFit=3)
{
	Int_t i, j;

    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            FitOmega_OneTarget(rootFile, iRun, i, j, iFit, 1);
        }
    }
}



