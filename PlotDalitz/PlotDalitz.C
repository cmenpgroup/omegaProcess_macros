// 
// PlotDalitz - plot the Dalitz shape from kinematics
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

void PlotDalitz(char *rootFile)
{
	Float_t xMinLo, xMaxLo;  // lower value of x-axis for drawing histogram
	Float_t xMinHi, xMaxHi;   // upper value of x-axis for drawing histogram
	
	Int_t i, j, k;
	
    TH2D *hist; // original histogram
    
    char hname[50];
    
    // data files contain the trees
    printf("Analyzing file %s\n",rootFile);
    TFile *fm = new TFile(rootFile,"READ"); // open up the ROOT file
    
    TDirectory *tmp = fm->GetDirectory("Nuc");
    
    sprintf(hname,"Dalitz_pip_AllCuts_Nuc");
    hist = (TH2D*)tmp->Get(hname); // get the histogram from the ROOT file
    
    // create canvas
	char title[100];
	char xtitle[100];
	char ytitle[100];
	sprintf(title,"#omega Dalitz"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,600,600); // create the canvas
	
//	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000);
		
	// fit the histogram
	Double_t par[4]={0.782,0.138,0.135,0.135};
    xMaxLo = (par[1] + par[2])*(par[1] + par[2]);
    xMaxHi = (par[0] - par[3])*(par[0] - par[3]);
    TF1 *fmax = new TF1("fmax",upperDalitz,xMaxLo,xMaxHi,4); // declare fit fcn
    fmax->SetParameters(&par[0]);    // set parameters for function

    xMinLo = (par[3] + par[2])*(par[3] + par[2]);
    xMinHi = (par[0] - par[1])*(par[0] - par[1]);
    TF1 *fmin = new TF1("fmax",lowerDalitz,xMinLo,xMinHi,4); // declare fit fcn
    fmin->SetParameters(&par[0]);    // set parameters for function
    
    // set up the Pad parameters
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleOffset(yoff);
//	hist->GetXaxis()->SetRangeUser(xLo,xHi);  // set the x-axis range for the plot
//	hist->SetMinimum(0); // start the y-axis at zero.
	hist->Draw("colz");

    fmin->SetLineWidth(4);
    fmin->Draw("same");
    fmax->SetLineWidth(4);
    fmax->Draw("same");
    
	// create the image files
	char OutCan[100];
	sprintf(OutCan,"%s.gif",PlotDalitz);
	can1->Print(OutCan);
	sprintf(OutCan,"%s.eps",PlotDalitz);
	can1->Print(OutCan);

}

// upper curve of Dalitz plot
Double_t upperDalitz(Double_t *x, Double_t *par){

    Double_t Msq = par[0]*par[0];
    Double_t Msq1 = par[1]*par[1];
    Double_t Msq2 = par[2]*par[2];
    Double_t Msq3 = par[3]*par[3];

    Double_t E2 = (x[0] - Msq1 + Msq2)/(2.0*sqrt(x[0]));
    Double_t E3 = (Msq - x[0] - Msq3)/(2.0*sqrt(x[0]));
    
    Double_t Esum = E2 + E3;
    Double_t Ediff = sqrt(E2*E2 - Msq2) -  sqrt(E3*E3 - Msq3);
    
    Double_t ret = Esum*Esum -  Ediff*Ediff;
    
    return ret;
}

// lower curve of Dalitz plot
Double_t lowerDalitz(Double_t *x, Double_t *par){
    
    Double_t Msq = par[0]*par[0];
    Double_t Msq1 = par[1]*par[1];
    Double_t Msq2 = par[2]*par[2];
    Double_t Msq3 = par[3]*par[3];
    
    Double_t E2 = (x[0] - Msq1 + Msq2)/(2.0*sqrt(x[0]));
    Double_t E3 = (Msq - x[0] - Msq3)/(2.0*sqrt(x[0]));
    
    Double_t Esum = E2 + E3;
    Double_t Ediff = sqrt(E2*E2 - Msq2) + sqrt(E3*E3 - Msq3);
    
    Double_t ret = Esum*Esum -  Ediff*Ediff;
    
    return ret;
}

