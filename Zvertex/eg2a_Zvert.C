// eg2a_Zvert.C
//
// macro to plot eg2a vertex histograms
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

const Int_t NSECTORS= 6;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};
char *hname[2] = {"elecZVertSector","elecZVertSector_Corr"};

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
// PlotZvert - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  inum = histogram number
//                  target = target name
//                  sector = sector number
//
void PlotZvert(char *fAna="Ana.root", int inum, char *target, int sector)
{
    char strname[50];
    
    if(inum<0 || inum>1){
        cout<<"Incorrect histogram number "<<inum<<endl;
        exit(0);
    }
    
    if(sector>0 && sector<=6){
        // Canvas to plot histogram
        TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
        c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        c1->SetBorderSize(5);
        gStyle->SetOptStat(0);
        c1->SetFillStyle(4000);
	
        // data files contain the trees
        printf("Analyzing file %s\n",fAna);
        TFile *fm = new TFile(fAna,"READ");
        TDirectory *tmp = fm->GetDirectory("Kinematics");

        c1->cd();
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
    
        TH2F *hZvert = (TH2F*)tmp->Get(hname[inum]);

        sprintf(strname,"%s_px_%i",hname[inum],sector);
        hZvert_proj = (TH1D*)hZvert->ProjectionX(strname,sector,sector,"");
    
        sprintf(title,"%s, Sector %i",target,sector);
        hZvert_proj->SetTitle(title);
        hZvert_proj->SetXTitle("e^{-} Z Vertex (cm)");
        hZvert_proj->GetXaxis()->CenterTitle();
        hZvert_proj->SetYTitle("Counts");
        hZvert_proj->GetYaxis()->CenterTitle();
        hZvert_proj->GetYaxis()->SetTitleOffset(yoff);
        hZvert_proj->SetLineWidth(2);
        hZvert_proj->Draw();

        sprintf(OutCan,"Plot_%s_%s_S%i.gif",hname[inum],target,sector);
        c1->Print(OutCan);
        sprintf(OutCan,"Plot_%s_%s_S%i.eps",hname[inum],target,sector);
        c1->Print(OutCan);
    }else{
        cout<<"Wrong sector number "<<sector<<endl;
    }
}

//
// OverlayZvertBySector - overlay histogram of e- Z vertex by sector
//
//                  fAna = output from eg2a DMS
//                  inum = histogram number
//                  target = target name
//
void OverlayZvertBySector(char *fAna="Ana.root", int inum, char *target)
{
	Int_t i;
	char legLabel[50];
    char strname[50];
    
    if(inum<0 || inum>1){
        cout<<"Incorrect histogram number "<<inum<<endl;
        exit(0);
    }
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5);
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
	
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");

    c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
 
    TH2F *hZvert = (TH2F*)tmp->Get(hname[inum]);
    
	TH1D *hZSec[NSECTORS];
    
	TLegend *leg = new TLegend(0.60,0.50,1.0,0.85); //declare Legend and give its location
    
    	for(i=0; i<NSECTORS ; i++){
            sprintf(strname,"%s_px_%i",hname[inum],i+1);
            hZSec[i] = (TH1D*)hZvert->ProjectionX(strname,i+1,i+1,"");
        	hZSec[i]->SetTitle(0);
        	hZSec[i]->SetXTitle("e^{-} Z Vertex (cm)");
        	hZSec[i]->GetXaxis()->CenterTitle();
        	hZSec[i]->SetYTitle("Counts");
        	hZSec[i]->GetYaxis()->CenterTitle();
        	hZSec[i]->GetYaxis()->SetTitleOffset(yoff);
        	hZSec[i]->SetLineWidth(2);
        	hZSec[i]->SetLineColor(lcol[i]);
        	hZSec[i]->Draw(fSame[i]);
        
        	sprintf(legLabel,"Sector %i",i+1);
        	leg->AddEntry(hZSec[i],legLabel,"l");
    	}
    	leg->SetLineColor(0);
    	leg->SetFillStyle(0);
    	leg->SetHeader(target);
    	leg->Draw();
    
	sprintf(OutCan,"Overlay_%s_%s.gif",hname[inum],target);
	c1->Print(OutCan);
	sprintf(OutCan,"Overlay_%s_%s.eps",hname[inum],target);
	c1->Print(OutCan);
}

//
// FitZvert - plot histogram with labels
//
//                  fAna = output from eg2a DMS
//                  inum = histogram number
//                  target = target name
//                  sector = sector number
//
void FitZvert(char *fAna="Ana.root", int inum, char *target, int sector)
{
    Float_t fPeak, fWidth, SumLo, SumHi;
    Float_t xLo = -26.5;  // lower value of x-axis for drawing histogram
    Float_t xHi = -22.0;   // upper value of x-axis for drawing histogram
    Float_t PeakLo = -26.0; // lower limit on the peak range
    Float_t PeakHi = -24.0; // upper limit on the peak range
    
    char strname[50];
    
    if(inum<0 || inum>1){
        cout<<"Incorrect histogram number "<<inum<<endl;
        exit(0);
    }
    
    if(sector>0 && sector<=6){
        // Canvas to plot histogram
        TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
        c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        c1->SetBorderSize(5);
        gStyle->SetOptFit(1);
        c1->SetFillStyle(4000);
        
        // data files contain the trees
        printf("Analyzing file %s\n",fAna);
        TFile *fm = new TFile(fAna,"READ");
        TDirectory *tmp = fm->GetDirectory("Kinematics");
        
        c1->cd();
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        TH2F *hZvert = (TH2F*)tmp->Get(hname[inum]);
        
        sprintf(strname,"%s_px_%i",hname[inum],sector);
        hZvert_proj = (TH1D*)hZvert->ProjectionX(strname,sector,sector,"");
        
        sprintf(title,"%s, Sector %i",target,sector);
        hZvert_proj->SetTitle(title);
        hZvert_proj->SetXTitle("e^{-} Z Vertex (cm)");
        hZvert_proj->GetXaxis()->CenterTitle();
        hZvert_proj->SetYTitle("Counts");
        hZvert_proj->GetYaxis()->CenterTitle();
        hZvert_proj->GetYaxis()->SetTitleOffset(yoff);
        hZvert_proj->SetLineWidth(2);
        hZvert_proj->Draw();
        
        // Fit histogram with a combined gaussian for peak and polynomial for background.
        // Set the initial fit parameters.
        // par[0] = gaussian amplitude
        // par[1] = gaussian centroid
        // par[2] = gaussian sigma
        // Polynomial : par[3] + par[4]*x + par[5]*x^2 + par[6]*x^3
        Double_t par[7]={1000.0,-25.0,1.0,1.0,1.0,1.0,1.0};
        TF1 *g1 = new TF1("g1",gaussFit,PeakLo,PeakHi,3); // declare fit fcn
        TF1 *pol = new TF1("pol",polFit,xLo,xHi,4);
        TF1 *t1 = new TF1("t1",totFit,xLo,xHi,7);
        
        g1->SetParameters(&par[0]);    // set parameters for initial peak fit
        hZvert_proj->Fit("g1","R");           // fit the peak
        g1->GetParameters(&par[0]);    // get parameters from initial peak fit
        
        hZvert_proj->Fit("pol","R+");         // fit the background
        pol->GetParameters(&par[3]); // get parameters fromt background fit
        
        t1->SetParNames("Amplitude","Centroid","Sigma","p0","p1","p2","p3");
        t1->SetParameters(par);       // set the parameters from initial fits for total fit
        t1->SetLineWidth(2);          // make fit line thicker
        t1->SetLineColor(2);          // set the fit line color
        hZvert_proj->Fit("t1","R");          // fit spectrum with total fit function
        t1->GetParameters(&par[0]);   //get final parameters
        
        fPeak = t1->GetParameter(1);  // get peak centroid
        fWidth = t1->GetParameter(2); // get peak width
        
        cout<<"***********************************************************"<<endl;
        cout<<"Fit parameters: Centroid = "<<fPeak<<" cm, Sigma = "<<fWidth<<" cm"<<endl;
        cout<<"***********************************************************"<<endl;
        
        sprintf(OutCan,"Fit_%s_%s_S%i.gif",hname[inum],target,sector);
        c1->Print(OutCan);
        sprintf(OutCan,"Fit_%s_%s_S%i.eps",hname[inum],target,sector);
        c1->Print(OutCan);
    }else{
        cout<<"Wrong sector number "<<sector<<endl;
    }
}

// peak is Gaussian
Double_t gaussFit(Double_t *x, Double_t *par){
    return TMath::Max(1.e-10,par[0]*TMath::Gaus(x[0],par[1],par[2]));
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
