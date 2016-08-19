// PlotSCMassSquared.C
//
// macro to analyze TOF mass-squared
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

#define NPART 5
#define MAX_HIST 4
#define MAX_RUN 4

gROOT->Reset();   // start from scratch

Int_t lcol[10] = {1,2,4,6,7,8,9,13,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

Float_t xLo[5] = {-0.4,-0.4,-0.4,-0.4,-0.4};
Float_t xHi[5] = {0.4,0.4,0.4,0.4,0.4};

char *HistName[MAX_HIST] = {"scMassSquared_NC","scMassSquared_EC","scMassSquared_PC","scMassSquared_EPC"};
char *RunName[MAX_RUN] = {"C","Fe","Sn","Pb"};
char *PartName[NPART] = {"Electron","#pi^{-}","#pi^{+}","Photon 1","Photon 2"};

// 
// PlotSCMassSquared_Particle - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  histIndex = histogram index
//                  tgtIndex = target index
//                  chan = particle channel
//
void PlotSCMassSquared_Particle(char *fAna,  Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0)
{
	char OutCan[100];
    char strname[100];
    
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
    
	TH2D *h2D = (TH2D*)tmp->Get(HistName[histIndex]);
    TH1D *h1D;

    sprintf(strname,"%s_px_%i_%i",HistName[histIndex],chan,chan);
    h1D = (TH1D*)h2D->ProjectionX(strname,chan+1,chan+1,"");
    
//    gPad->SetLogy();
    
	h1D->SetTitle(0);
	h1D->GetXaxis()->CenterTitle();
	h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
	h1D->GetYaxis()->SetTitleOffset(yoff);
//    h1D->SetAxisRange(xLo[chan],xHi[chan],"X");
    h1D->SetLineWidth(2);
    h1D->Draw();

	sprintf(OutCan,"Plot_%s_%i_%s.gif",HistName[histIndex],chan,RunName[tgtIndex]);
	c1->Print(OutCan);
	sprintf(OutCan,"Plot_%s_%i_%s.eps",HistName[histIndex],chan,RunName[tgtIndex]);
	c1->Print(OutCan);
}

//
// PlotSCMassSquared_Particle - plot histogram with labels
//
//                  fAna = output from eg2a DMS
//                  histIndex = histogram index
//                  tgtIndex = target index
//                  chan = particle channel
//
void OverlaySCMassSquared_All(char *fAna,  Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    Int_t ymax = 0;
    Int_t ytmp;
    
    char OutCan[100];
    char strname[100];
    char legLabel[50];
    
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
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    TH2D *h2D = (TH2D*)tmp->Get(HistName[histIndex]);
    TH1D *h1D[NPART];
    
    for(i=0; i<NPART; i++){
        sprintf(strname,"%s_px_%i_%i",HistName[histIndex],i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    
        h1D[i]->SetTitle(0);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(i+1);
        ytmp = h1D[i]->GetMaximum();
        if(ytmp > ymax){
            cout<< ymax <<" "<<ytmp<<endl;
            h1D[0]->SetMaximum(ytmp*1.1);
            ymax = ytmp;
        }
        h1D[i]->Draw(fSame[i]);
        
        sprintf(legLabel,"%s",PartName[i]);
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Particles:");
    leg->Draw();
    
    sprintf(OutCan,"OL_%s_%s.gif",HistName[histIndex],RunName[tgtIndex]);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_%s_%s.eps",HistName[histIndex],RunName[tgtIndex]);
    c1->Print(OutCan);
}


//
// FitSCMassSquared_Particle - plot histogram with labels
//
//                  fAna = output from eg2a DMS
//                  histIndex = histogram index
//                  tgtIndex = target index
//                  chan = particle channel
//
void FitSCMassSquared_Particle(char *fAna,  Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0)
{
    char OutCan[100];
    char strname[100];
    
    Float_t xLo = -0.1;  // lower value of x-axis for drawing histogram
    Float_t xHi = 0.1;   // upper value of x-axis for drawing histogram
    Float_t PeakLo = -0.05; // lower limit on the peak range
    Float_t PeakHi = 0.05; // upper limit on the peak range
    Float_t Nsigma = 3.0;  // number of sigma of the peak
    Int_t ibg = 4; // parameter index in the par[] array where the background parameters start
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptFit(1111);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TH2D *h2D = (TH2D*)tmp->Get(HistName[histIndex]);
    TH1D *h1D;
    
    sprintf(strname,"%s_px_%i_%i",HistName[histIndex],chan,chan);
    h1D = (TH1D*)h2D->ProjectionX(strname,chan+1,chan+1,"");

    Double_t par[6]={1.0,0.0,1.0,1.0,1.0,1.0};
    TF1 *v1 = new TF1("v1",voigtFit,PeakLo,PeakHi,4);
    TF1 *pol = new TF1("pol",polFit,xLo,xHi,2);
    TF1 *t1 = new TF1("t1",totVoigtFit,xLo,xHi,6);
    
    t1->SetParName(0,"amplitude");
    t1->SetParName(1,"centroid");
    t1->SetParName(2,"Gsig");
    t1->SetParName(3,"Lgam");
    t1->SetParName(4,"Pol: A_{1}");
    t1->SetParName(5,"Pol: A_{2}");
    
    v1->SetParameters(&par[0]);
    h1D->Fit("v1","R");
    v1->GetParameters(&par[0]);
    
    h1D->Fit("pol","R+");         // fit the background
    pol->GetParameters(&par[ibg]); // get parameters from background fit
    
    t1->SetParameters(par);       // set the parameters from initial fits for total fit
    t1->SetLineWidth(2);          // make fit line thicker
    t1->SetLineColor(2);          // set the fit line color
    h1D->Fit("t1","R");          // fit spectrum with total fit function
    t1->GetParameters(&par[0]);   //get final parameters
    
    h1D->SetTitle(0);
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
    h1D->GetYaxis()->SetTitleOffset(yoff);
    //    h1D->SetAxisRange(xLo[chan],xHi[chan],"X");
    h1D->SetLineWidth(2);
    h1D->Draw();
    
    double ymax = h1D->GetMaximum();
    
    double Gcentroid = t1->GetParameter(1);
    double Gsig = t1->GetParameter(2);
    double Lgam = t1->GetParameter(3);
    
    double fG = 2*Gsig*sqrt(2*log(2)); // FWHM for Gaussian
    double fL = 2*Lgam; // FWHM for Lorentzian
    double fV = 0.5336*fL + sqrt(0.2166*fL*fL + fG*fG); // FWHM for Voigtian

    double xmin = 1.5*(Gcentroid-fV);
    double xmax = 1.5*(Gcentroid+fV);
    char textLeft[50];
    sprintf(textLeft,"%f",xmin);
    TText *tLeft = new TText(0.95*xmin,0.25*ymax,textLeft);
    tLeft->SetTextSize(0.025);
    tLeft->Draw();
    TLine *lLeft = new TLine(xmin,0,xmin,0.5*ymax);
    lLeft->SetLineWidth(2);
    lLeft->Draw();
    char textRight[50];
    sprintf(textRight,"%f",xmax);
    TText *tRight = new TText(1.1*xmax,0.25*ymax,textRight);
    tRight->SetTextSize(0.025);
    tRight->Draw();
    TLine *lRight = new TLine(xmax,0,xmax,0.5*ymax);
    lRight->SetLineWidth(2);
    lRight->Draw();
    
    sprintf(OutCan,"Fit_%s_%i_%s.gif",HistName[histIndex],chan,RunName[tgtIndex]);
    c1->Print(OutCan);
    sprintf(OutCan,"Fit_%s_%i_%s.eps",HistName[histIndex],chan,RunName[tgtIndex]);
    c1->Print(OutCan);
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
    Int_t nmax = 2;
    Double_t bck = 0.0;
    for (Int_t i = 0; i<=nmax; i++){
        bck += par[i]*pow(x[0],i);
    }
    return bck;
}

// Sum of background and peak function
// peak uses par[0-2]
// background uses par[3-4]
Double_t totFit(Double_t *x, Double_t *par) {
    return gaussFit(x,par) + polFit(x,&par[3]);
}

// peak is Gaussian
Double_t voigtFit(Double_t *x, Double_t *par){
    return TMath::Max(1.e-10,par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3]));
}

// Sum of background and peak function
// peak uses par[0-2]
// background uses par[3-4]
Double_t totVoigtFit(Double_t *x, Double_t *par) {
    return voigtFit(x,par) + polFit(x,&par[4]);
}
