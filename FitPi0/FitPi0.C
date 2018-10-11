// 
// FitPi0 - background subtraction routine where the peak is fit with a well-known function
//			   like a gaussian.  The backgorund is fit with another well-known function like a
//             polynomial.
//             The output is a list of yields.
//                  
//                  fname = ROOT file containing the histogram
//                  hname = histogram name
//					yieldFile = name of output file without the suffix
//
// M. H. Wood, Canisius College
//--------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MAX_CUT 22
#define MAX_HIST 3
#define MAX_TGT 3
#define MAX_RUN 4

//gROOT->Reset();

Int_t lcol[MAX_CUT] = {1,2,4,6,7,8,9,13,14,15,16,17,18,19};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

Float_t xLo = 0.03;  // lower value of x-axis for drawing histogram
Float_t xHi = 0.5;   // upper value of x-axis for drawing histogram
Float_t PeakLo = 0.11; // lower limit on the peak range
Float_t PeakHi = 0.16; // upper limit on the peak range
Float_t Nsigma = 3.0;  // number of sigma of the peak

void Check_HistIndex(Int_t index);
void PrintHistIndex();
void Check_TgtIndex(Int_t index);
void PrintTgtIndex();
void Check_CutIndex(Int_t index);
void PrintCutIndex();
void Check_CutLoHi(Int_t Lo, Int_t Hi);

Double_t gaussFit(Double_t *x, Double_t *par);
Double_t polFit(Double_t *x, Double_t *par);
Double_t totFit(Double_t *x, Double_t *par);
Double_t breitwigner(Double_t *x, Double_t *par);

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
    LabelCuts.push_back("ChargedPionMom");
    
    LegHeader.push_back("Cuts: ");
    LegHeader.push_back("All Cuts Except:");
    LegHeader.push_back("Anti-Cuts:");
    
    HistName.push_back("IM2Photons_");
    HistName.push_back("IM2Photons_woCut_");
    HistName.push_back("IM2Photons_antiCut_");
    
    RunName.push_back("C12");
    RunName.push_back("Fe56");
    RunName.push_back("Sn");
    RunName.push_back("Pb208");
    
    TgtName.push_back("NoTarget");
    TgtName.push_back("LD2");
    TgtName.push_back("Nuc");
}

Float_t FitPi0(TH1D *hist, char *plotFilePrefix)
{
    Float_t fPeak, fWidth, SumLo, SumHi;
	
	Int_t ibg = 3; // parameter index in the par[] array where the background parameters start

	Int_t i, j, k;
	Float_t x, pval;
	Float_t yield;
	
    cout << "Fitting histogram " << hist->GetName() << endl;
    
	TH1F *hBgFit; // histogram of background
		
	// create canvas
	char title[100];
	char xtitle[100];
	char ytitle[100];
	sprintf(title,"Fitting Analysis"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,600,600); // create the canvas
	
//	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000); 
	
	sprintf(xtitle,"Invariant Mass (GeV)"); // set the x-axis title
	sprintf(ytitle,"Counts"); // set the y-axis title
		
	// histogram of background
	hBgFit = (TH1F*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	hBgFit->SetName("hBgFit");
	hBgFit->SetTitle("Background");

	// fit the histogram
	Double_t par[7]={1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    TF1 *g1 = new TF1("g1",gaussFit,PeakLo,PeakHi,3); // declare fit fcn
    TF1 *pol = new TF1("pol",polFit,xLo,xHi,4);
    TF1 *t1 = new TF1("t1",totFit,xLo,xHi,7);
		
	g1->SetParameters(&par[0]);    // set parameters for initial peak fit
	hist->Fit("g1","R");           // fit the peak
	g1->GetParameters(&par[0]);    // get parameters from initial peak fit

	hist->Fit("pol","R+");         // fit the background
	pol->GetParameters(&par[ibg]); // get parameters fromt background fit

	t1->SetParameters(par);       // set the parameters from initial fits for total fit
	t1->SetLineWidth(5);          // make fit line thicker
	t1->SetLineColor(1);          // set the fit line color
	hist->Fit("t1","R");          // fit spectrum with total fit function
	t1->GetParameters(&par[0]);   //get final parameters
	
	fPeak = t1->GetParameter(1);  // get peak centroid
	fWidth = t1->GetParameter(2); // get peak width
	SumLo = fPeak - Nsigma*fWidth; // calc lower edge at Nsigma away from centroid
	SumHi = fPeak + Nsigma*fWidth; // calc upper edge at Nsigma away from centroid
			
	pol->SetParameters(&par[ibg]); // set the pfinal parameters for the background function
	hBgFit->Add(pol,-1.0); // subtract background function
	for(k=1; k<hBgFit->GetNbinsX(); k++){
		x = hBgFit->GetBinCenter(k); // read bin center value on the x-axis
		if(x>=SumLo && x<=SumHi){ // check that x is in the peak summation region
				pval = hBgFit->GetBinContent(k); // get the number of counts in the bin
			if(pval<0.0) pval = 0.0; // if neg. counts, set to zero
		}else{
			pval =0.0;
		}
		hBgFit->SetBinContent(k,pval); // refill the histogram
	}
    yield = hBgFit->Integral(hBgFit->FindBin(SumLo),hBgFit->FindBin(SumHi)); // sum total counts in peak

	// set up the Pad parameters
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	// draw the original histogram
	hist->SetTitle(0);
	hist->GetXaxis()->SetTitle(xtitle);
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle);
	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleOffset(yoff);
//	hist->GetXaxis()->SetRangeUser(xLo,xHi);  // set the x-axis range for the plot
	hist->SetLineWidth(2);
	hist->SetMinimum(0); // start the y-axis at zero.
	hist->Draw();
	
	// draw the background-subtracted histogram
	hBgFit->SetLineWidth(2);
	hBgFit->SetFillColor(4);
	hBgFit->Draw("same");

	// create the image files
	char OutCan[100];
	sprintf(OutCan,"%s.gif",plotFilePrefix);
	can1->Print(OutCan);

	can1->Close();

    Float_t yield_total = hist->Integral(hBgFit->FindBin(SumLo),hBgFit->FindBin(SumHi)); // signal+background

    cout<<endl;
    cout<<"**************************************************************"<<endl;
    cout<<"Yields:"<<endl;
    cout<<"Signal+Background = "<<yield_total<<endl;
    cout<<"Background = "<<yield_total - (int)yield<<endl;
    cout<<"Signal = "<<(int)yield<<endl;
    cout<<"**************************************************************"<<endl;
    cout<<endl;
    
	return yield;
}

void FitPi0_OneTarget(string fAna, Int_t histIndex=0, Int_t runIndex=0, Int_t tgtIndex = 0, Int_t cutIndex=0)
{
	Int_t i;
	Float_t yield;
    char strname[100];
    char hname[50];
    char plotFilePrefix[100];
    
    Cuts myCuts;
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(cutIndex);
    
	// data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    TDirectory *tmp = fm->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
    sprintf(hname,"%s_%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    cout << hname << endl;
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,i);
    TH1D *hist = (TH1D*)h2D->ProjectionX(strname,cutIndex+1,cutIndex+1,"");
    
    sprintf(plotFilePrefix,"FitPi0_%s_%s",myCuts.Get_Run(runIndex).c_str(),hname);

    yield = FitPi0(hist,plotFilePrefix);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout<<myCuts.Get_Tgt(tgtIndex).c_str()<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

	fout.close(); // close the text file
	fm->Close();  // close ROOT file
}

void FitPi0_AllTargets(char *rootFile, Int_t iRun=0)
{
	Int_t i, j;

    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            FitPi0_OneTarget(rootFile, iRun, i, j);
        }
    }
}

// peak is Gaussian
Double_t gaussFit(Double_t *x, Double_t *par){
	return TMath::Max(1.e-10,par[0]*TMath::Gaus(x[0],par[1],par[2]));
}
// peak is Breit Wigner
Double_t breitwigner(Double_t *x, Double_t *par){
	return TMath::Max(1.e-10,par[0]*TMath::BreitWigner(x[0],par[1],par[2]));
}

// background function is polynomial
Double_t polFit(Double_t *x, Double_t *par){
	Int_t nmax = 3;
	Double_t bck = 0.0;
	for (Int_t i = 0; i<=nmax; i++){
		bck += par[i]*pow(x[0],i);
	}
	return bck;
}

// Sum of background and peak function
// peak uses par[0-2]
// background uses par[3-5]
Double_t totFit(Double_t *x, Double_t *par) {
	return gaussFit(x,par) + polFit(x,&par[3]);
}

//
// PlotPi0_CutIndex - plot Pi0 inv. mass for a specific cut selection
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//                  chanLo = lower bin
//                  chanHi = upper bin
//
void PlotPi0_CutIndex(string fAna, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chanLo = 0, Int_t chanHi=0)
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
    
    TLegend *leg = new TLegend(0.3,0.5,1.0,0.875);
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D = (TH2D*)tmp->Get(hname);

    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D->ProjectionX(strname,1,1,"");
    
    sprintf(title,"Target: %s",myCuts.Get_Tgt(tgtIndex).c_str());
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();
    
    sprintf(legLabel,"%s",myCuts.Get_Cuts(0).c_str());
    leg->AddEntry(h1D_noCut,legLabel,"l");
    
    iColor++;
    
    for(i=chanLo; i<chanHi+1; i++){
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        
        sprintf(title,"Target: %s",myCuts.Get_Tgt(tgtIndex).c_str());
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
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
    
    sprintf(OutCan,"PlotPi0_%s_%i_%i.gif",hname,chanLo,chanHi);
    c1->Print(OutCan);
}

//
// PlotPi0_AllGrid - plot Pi0 inv. mass for all cuts
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//
void PlotPi0_AllGrid(string fAna, Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    
    TH1D *h1D[MAX_CUT];
    
    Cuts myCuts;
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    int nrow = ceil(sqrt(MAX_CUT));
    cout<<"nrow "<<nrow<<endl;
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
    
    sprintf(OutCan,"PlotPi0_AllGrid_%s.gif",hname);
    c1->Print(OutCan);
}

//
// CompPi0_Files - compare pi0 inv. mass for a specific cut selection between files
//
//                  fAna1 = output 1 from eg2a DMS
//                  fAna2 = output 2 from eg2a DMS
//                  histIndex = 2-D histogram name
//                  tgtIndex = target index
//                  chan = cut bin
//
void CompPi0_Files(string fAna1, string fAna2, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0, string legLine1 = "Leg 1", string legLine2 = "Leg 2", string comment="test")
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
    
    TLegend *leg = new TLegend(0.3,0.5,1.0,0.875);
    
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
        if(i==1) h1D[1]->Scale(h1D[0]->Integral()/h1D[1]->Integral());
        (i==0) ? h1D[i]->Draw() : h1D[i]->Draw("same");
        
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Files:");
    leg->Draw();
    
    sprintf(OutCan,"CompPi0_%s_%i_%s.gif",hname,chan,comment.c_str());
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

