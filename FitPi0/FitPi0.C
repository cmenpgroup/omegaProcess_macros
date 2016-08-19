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

#define MAX_CUT 16
#define MAX_HIST 3
#define MAX_TGT 3

gROOT->Reset();

Int_t lcol[MAX_CUT] = {1,2,4,6,7,8,9,13,14,15,16,17,18,19};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[MAX_CUT] = {"","same","same","same","same","same","same","same", "same", "same", "same", "same", "same", "same"};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

Float_t xLo = 0.03;  // lower value of x-axis for drawing histogram
Float_t xHi = 0.5;   // upper value of x-axis for drawing histogram
Float_t PeakLo = 0.11; // lower limit on the peak range
Float_t PeakHi = 0.16; // upper limit on the peak range
Float_t Nsigma = 3.0;  // number of sigma of the peak

char *legHeader[3] = {"Cuts: ","All Cuts Except:","Anti-Cuts:"};
char *RunName[4] = {"C","Fe","Sn","Pb"};
char *TgtName[MAX_TGT] = {"NoTarget","LD2","Nuc"};
char *HistName[MAX_HIST] = {"IM2Photons","IM2Photons_woCut","IM2Photons_OpAng_ElecPhoton_Cut"};
char *CutName[MAX_CUT] = {"None","All #omega cuts","M(#pi^{0})","Q^{2}","W","V_{z} Matching","Part. Topology","#theta_{e-,#gamma}","M(#pi^{+}#pi^{-}","EC 2nd Moment for #gamma's,  Region 1","EC 2nd Moment for #gamma's, Region 2","EC 2nd Moment for #gamma's, Region 3","EC 3rd Moment for #gamma's, Region 1","EC 3rd Moment for #gamma's, Region 2","EC 3rd Moment for #gamma's, Region 3"};

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
	sprintf(OutCan,"%s.eps",plotFilePrefix);
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

void FitPi0_OneTarget(char *fAna, Int_t histIndex=0, Int_t runIndex=0, Int_t tgtIndex = 0, Int_t cutIndex=0)
{
	Int_t i;
	Float_t yield;
    char strname[100];
    char hname[50];
    char plotFilePrefix[100];
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(cutIndex);
    
	// data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(TgtName[tgtIndex]);
    
    sprintf(hname,"%s_%s",HistName[histIndex],TgtName[tgtIndex]);
    cout << hname << endl;
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,i);
    TH1D *hist = (TH1D*)h2D->ProjectionX(strname,cutIndex+1,cutIndex+1,"");
    
    sprintf(plotFilePrefix,"FitPi0_%s_%s",RunName[runIndex],hname);

    yield = FitPi0(hist,plotFilePrefix);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout<<TgtName<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

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
void PlotPi0_CutIndex(char *fAna, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chanLo = 0, Int_t chanHi=0)
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
    
    TLegend *leg = new TLegend(0.3,0.5,1.0,0.875);
    
    sprintf(hname,"%s_%s",HistName[histIndex],TgtName[tgtIndex]);
    TH2D *h2D = (TH2D*)tmp->Get(hname);

    sprintf(strname,"%s_0",hname);
    TH1D *h1D_noCut = (TH1D*)h2D->ProjectionX(strname,1,1,"");
    
    sprintf(title,"Target: %s",TgtName[tgtIndex]);
    h1D_noCut->SetTitle(title);
    h1D_noCut->GetXaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->CenterTitle();
    h1D_noCut->GetYaxis()->SetTitle("Counts");
    h1D_noCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_noCut->SetLineWidth(2);
    h1D_noCut->Draw();
    
    sprintf(legLabel,"%s",CutName[0]);
    leg->AddEntry(h1D_noCut,legLabel,"l");
    
    iColor++;
    
    for(i=chanLo; i<chanHi+1; i++){
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        
        sprintf(title,"Target: %s",TgtName[tgtIndex]);
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        if(chanLo!=chanHi) h1D[i]->SetLineColor(lcol[iColor]);
        h1D[i]->Draw("same");
        
        sprintf(legLabel,"%s",CutName[i]);
        leg->AddEntry(h1D[i],legLabel,"l");
        
        iColor++;
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader(legHeader[histIndex]);
    leg->Draw();
    
    sprintf(OutCan,"PlotPi0_%s_%i_%i.gif",hname,chanLo,chanHi);
    c1->Print(OutCan);
    sprintf(OutCan,"PlotPi0_%s_%i_%i.eps",hname,chanLo,chanHi);
    c1->Print(OutCan);
}

//
// PlotPi0_AllGrid - plot Pi0 inv. mass for all cuts
//
//                  fAna = output from eg2a DMS
//                  tgtIndex = target index
//
void PlotPi0_AllGrid(char *fAna, Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char strname[100];
    
    TH1D *h1D[MAX_CUT];
    
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
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(TgtName[tgtIndex]);
    
    sprintf(hname,"%s_%s",HistName[histIndex],TgtName[tgtIndex]);
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
    
    sprintf(OutCan,"PlotPi0_AllGrid_%s.gif",hname);
    c1->Print(OutCan);
    sprintf(OutCan,"PlotPi0_AllGrid_%s.eps",hname);
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

