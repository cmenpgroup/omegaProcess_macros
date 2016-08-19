// PlotEC.C
//
// macro to plot eg2a histograms for EC
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
// PlotEC_xy_local - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  fDir = ROOT directory
//                  HistName = histogram name
//                  target = target name
//
void PlotEC_Sector(char *fAna="Ana.root", char *fDir, char *HistName, char *target)
{
    Int_t i;
    TH2D *hist[6];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,800);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(1111);
	c1->SetFillStyle(4000);
    c1->Divide(3,2);
    
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
	TDirectory *tmp = fm->GetDirectory(fDir);
    
    for(i=0; i<NSECTORS; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
    
        c1->cd(i+1);
        sprintf(hname,"%s%i",HistName,i+1);
        hist[i] = (TH2D*)tmp->Get(hname);
        sprintf(htitle,"Sector %i, %s",i+1,target);
        hist[i]->SetTitle(htitle);
        hist[i]->GetXaxis()->CenterTitle();
        hist[i]->GetYaxis()->CenterTitle();
        hist[i]->GetYaxis()->SetTitleOffset(yoff);
        hist[i]->Draw("colz");

    }
    
	sprintf(OutCan,"Plot_%s_%s.gif",HistName,target);
	c1->Print(OutCan);
	sprintf(OutCan,"Plot_%s_%s.eps",HistName,target);
	c1->Print(OutCan);
}

//
// OverlayEC_xy_local - overlay histogram of EC x vs y
//
//                  fAna = output from eg2a DMS
//                  fDir = ROOT directory
//                  HistName1 = 1st histogram name
//                  HistName2 = 2nd histogram name
//                  target = target name
//
void OverlayEC_xy_local(char *fAna="Ana.root", char *fDir, char *HistName1, char *HistName2, char *target)
{
	Int_t i;
    TH2D *hist1[6];
    TH2D *hist2[6];

    // Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,800);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5);
	gStyle->SetOptStat(1111);
	c1->SetFillStyle(4000);
    c1->Divide(3,2);
    
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(fDir);
    
    for(i=0; i<NSECTORS; i++){
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        c1->cd(i+1);
        sprintf(hname,"EC_XvsY_local_%s%i",HistName1,i+1);
        hist1[i] = (TH2D*)tmp->Get(hname);
        sprintf(htitle,"Sector %i, %s",i+1,target);
        hist1[i]->SetTitle(htitle);
        hist1[i]->GetXaxis()->CenterTitle();
        hist1[i]->GetYaxis()->CenterTitle();
        hist1[i]->GetYaxis()->SetTitleOffset(yoff);
        hist1[i]->Draw("colz");
        
        sprintf(hname,"EC_XvsY_local_%s%i",HistName2,i+1);
        hist2[i] = (TH2D*)tmp->Get(hname);
        hist2[i]->SetMarkerColor(1);
        hist2[i]->Draw("same");

    }
    
    sprintf(OutCan,"OL_EC_XvsY_local_%s_VS_%s_%s.gif",HistName1,HistName2,target);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_EC_XvsY_local_%s_VS_%s_%s.eps",HistName1,HistName2,target);
    c1->Print(OutCan);
}

//
// Analyze_ECvsP - fit the EC/P vs P distributions
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Analyze(char *fAna="Ana.root", char *target)
{
    Int_t i, j;
    
    const Int_t NPARAM = 3;
    TH1D *h1D[NSECTORS][NPARAM];
    TH2D *h2D[NSECTORS];
    
    char HistName2D[50];
    char strname[50];
    char *yname[NPARAM] = {"Amplitude","Mean","Sigma"};
    
    // open text file for the yields
    char OutFile[100];
    sprintf(OutFile,"ECinP_vsECoutP_Fit_%s.dat",target);
    ofstream fout(OutFile);
    
    TCanvas *can[NSECTORS]; // Canvas to plot histogram
    
    Double_t par[3]={1.0,1.0,1.0};
    TF1 *pol;
    TF1 *sig;
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");
    
    for(j=0; j<NSECTORS; j++){
        sprintf(cname,"can%i",j+1);
        sprintf(cname,"Canvas, Sector %i",j+1);
        can[j] = new TCanvas(cname,ctitle,50*j,0,600,600);
        can[j]->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        can[j]->SetBorderSize(5);
        gStyle->SetOptStat(1111);
        can[j]->SetFillStyle(4000);
        can[j]->Divide(2,2);
        
        sprintf(HistName2D,"ECinP_VS_ECoutP_Sector%i",j+1);
        h2D[j] = (TH2D*)tmp->Get(HistName2D);

        can[j]->cd(1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);

        h2D[j]->GetXaxis()->CenterTitle();
        h2D[j]->GetYaxis()->CenterTitle();
        h2D[j]->GetYaxis()->SetTitleOffset(yoff);
    
        h2D[j]->SetAxisRange(0.0,0.3,"X");
        h2D[j]->Draw("colz");

        h2D[j]->FitSlicesY(0,0,-1,20,"QNR");//TODO MAYBE
    
        fout<<j+1<<"\t";
        
        for(i=0; i<NPARAM; i++){
            can[j]->cd(i+2);
            gPad->SetLeftMargin(Lmar);
            gPad->SetRightMargin(Rmar);
            gPad->SetFillColor(0);

            sprintf(strname,"%s_%i",HistName2D,i);
            h1D[j][i] = (TH1D*)gDirectory->Get(strname);
            h1D[j][i]->GetXaxis()->CenterTitle();
            h1D[j][i]->GetYaxis()->CenterTitle();
            h1D[j][i]->GetYaxis()->SetTitle(yname[i]);
            h1D[j][i]->GetYaxis()->SetTitleOffset(yoff);
            h1D[j][i]->SetAxisRange(0.0,0.5,"X");
            if (i==0){
                h1D[j][i]->SetAxisRange(0.0,300,"Y");
            }
            if (i==1) {
                h1D[j][i]->SetAxisRange(0.05,0.3,"Y");
                pol = new TF1("pol","pol1",0.1,0.25); //ADJUST fit range
                h1D[j][i]->Fit("pol","R");
                fout<<pol->GetParameter(0)<<"\t"<<pol->GetParameter(1)<<"\t";//<<pol->GetParameter(2)<<"\t";
            }
            if (i==2) {
                par[0] = 10.0;
                par[1] = 0.1;
                par[2] = 0.1;
                h1D[j][i]->SetAxisRange(0.0,0.3,"Y");
                sig = new TF1("sig","pol1",0.1,0.25); //ADJUST fit range "sig = new TF1("sig",SigmaFit,0.16,0.25,2);"
                //sig->SetParameters(par);
                h1D[j][i]->Fit("sig","R");
                fout<<sig->GetParameter(0)<<"\t"<<sig->GetParameter(1)<<endl;//"\t"<<sig->GetParameter(2)<<endl;
            }
            h1D[j][i]->Draw();
        }
       	sprintf(OutCan,"Analyze_ECinP_vs_ECoutP_S%i_%s.gif",j+1,target);
        can[j]->Print(OutCan);
        sprintf(OutCan,"Analyze_ECinP_vs_ECoutP_S%i_%s.eps",j+1,target);
        can[j]->Print(OutCan);
    }
    fout.close(); // close the text file
}

// function to fit sigma vs momentum distribution
Double_t SigmaFit(Double_t *x, Double_t *par){
    return TMath::Max(1.e-10,sqrt(par[0]*par[0]+par[1]*par[1]/sqrt(x[0])));
}

// background function is polynomial
Double_t polFit(Double_t *x, Double_t *par){
    Int_t nmax = 1;
    Double_t bck = 0.0;
    for (Int_t i = 0; i<=nmax; i++){
        bck += par[i]*pow(x[0],i);
    }
    return bck;
}

// Sum of Mean and Sigma functions
// peak uses par[0-2]
// background uses par[3-5]
Double_t CutBelow(Double_t *x, Double_t *par) {
    //Double_t arr[] = {par[0] - 0.5 * (par[2]), par[1]};
    return polFit(x,par) - 2.0*polFit(x,&par[2]);
}

// Sum of Mean and Sigma functions
// peak uses par[0-2]
// background uses par[3-5]
Double_t CutAbove(Double_t *x, Double_t *par) {
    //Double_t arr[] = {par[0] + 0.5 * (par[2]), par[1]};
    return polFit(x,par) + 2.0*polFit(x,&par[2]);
}

//
// Plot_ECvsP - plot the EC/P vs P distributions with cuts
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot(char *fAna="Ana.root", char *fitParams="ECvsP_Fit.dat", char *target)
{
    Int_t i, j;
    Int_t iSector;
    
    const Int_t NPARAM = 3;
    TF1 *fitMean[NSECTORS];
    TF1 *fitCut[NSECTORS];
    TH2D *h2D[NSECTORS];
    
    char funcName[50];
    
    Double_t a, b, c, d, f, z;
    Double_t parMean[6];
    Double_t parSigma[2];
    
    // open text file for fit parameters
    ifstream fin(fitParams);

    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(3,2);
    
    for(j=0; j<NSECTORS; j++){
        sprintf(hname,"ECinP_VS_ECoutP_Sector%i",j+1);
        h2D[j] = (TH2D*)tmp->Get(hname);
        
        c1->cd(j+1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        h2D[j]->GetXaxis()->CenterTitle();
        h2D[j]->GetYaxis()->CenterTitle();
        h2D[j]->GetYaxis()->SetTitleOffset(yoff);
        
        h2D[j]->SetAxisRange(0.0,3.0,"X");
        h2D[j]->Draw("colz");
        
        fin >> iSector >> a >> b >> c >> d >> f;// >> z;
        parMean[0] = a;
        parMean[1] = b;
        parMean[2] = c;
        parMean[3] = d;
        //parMean[4] = f;
        //parMean[5] = z;
        
        sprintf(funcName,"fitMean%i",j+1);
        fitMean[j] = new TF1(funcName,polFit,0.0,0.26,3);
        fitMean[j]->SetParameters(&parMean[0]);
        fitMean[j]->SetLineWidth(2);
        fitMean[j]->Draw("same");
        
        sprintf(funcName,"fitCut%i",j+1);
        fitCut[j] = new TF1(funcName,CutBelow,0.0,0.26,5);
        fitCut[j]->SetParameters(parMean);
        fitCut[j]->SetLineColor(4);
        fitCut[j]->SetLineWidth(2);
        fitCut[j]->Draw("same");
        
        sprintf(funcName,"fitCut%i",j+1);
        fitCut[j] = new TF1(funcName,CutAbove,0.0,0.26,5);
        fitCut[j]->SetParameters(parMean);
        fitCut[j]->SetLineColor(4);
        fitCut[j]->SetLineWidth(2);
        fitCut[j]->Draw("same");
    }
    sprintf(OutCan,"Plot_ECinPvsECoutP_%s.gif",target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_ECinPvsECoutP_%s.eps",target);
    c1->Print(OutCan);
    
    fin.close();
}

