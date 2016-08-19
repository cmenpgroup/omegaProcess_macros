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

const Int_t NSECTORS = 6;
const Int_t NPHOTONS = 2;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};
char *Plabel[3] = {"#pi^{+}","#pi^{-}","#pi^{+}#pi^{-}"};
char *Detlabel[5] = {"Electron","#pi^{-}","#pi^{+}","photon 1", "photon2"};

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
// PlotEC_Sector - plot histogram with labels for electrons
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
// OverlayEC_xy_local_1Sector - overlay histogram of EC x vs y for one sector
//
//                  fAna = output from eg2a DMS
//                  fDir = ROOT directory
//                  HistName1 = 1st histogram name
//                  HistName2 = 2nd histogram name
//                  target = target name
//
void OverlayEC_xy_local_1Sector(char *fAna="Ana.root", char *fDir, char *HistName1, char *HistName2, int iSector,char *target)
{
    Int_t i;
    TH2D *hist1;
    TH2D *hist2;
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory(fDir);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"EC_XvsY_local_%s%i",HistName1,iSector);
    hist1 = (TH2D*)tmp->Get(hname);
    sprintf(htitle,"Sector %i, %s",iSector,target);
    hist1->SetTitle(htitle);
    hist1->GetXaxis()->CenterTitle();
    hist1->GetYaxis()->CenterTitle();
    hist1->GetYaxis()->SetTitleOffset(yoff);
    hist1->Draw("colz");
        
    sprintf(hname,"EC_XvsY_local_%s%i",HistName2,iSector);
    hist2 = (TH2D*)tmp->Get(hname);
    hist2->SetMarkerColor(1);
    hist2->Draw("same");
    
    sprintf(OutCan,"OL_EC_XvsY_local_%s_VS_%s_S%i_%s.gif",HistName1,HistName2,iSector,target);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_EC_XvsY_local_%s_VS_%s_S%i_%s.eps",HistName1,HistName2,iSector,target);
    c1->Print(OutCan);
}


//
// Analyze_ECvsP - fit the EC/P vs P distributions
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Analyze_ECvsP(char *fAna="Ana.root", char *target)
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
    sprintf(OutFile,"ECvsP_Fit_%s.dat",target);
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
        
        sprintf(HistName2D,"ECtotP_VS_P_Sector%i",j+1);
        h2D[j] = (TH2D*)tmp->Get(HistName2D);

        can[j]->cd(1);
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);

        h2D[j]->GetXaxis()->CenterTitle();
        h2D[j]->GetYaxis()->CenterTitle();
        h2D[j]->GetYaxis()->SetTitleOffset(yoff);
    
        h2D[j]->SetAxisRange(0.0,3.0,"X");
        h2D[j]->Draw("colz");
    
        h2D[j]->FitSlicesY(0,0,-1,20);
    
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
            h1D[j][i]->SetAxisRange(0.0,2.5,"X");
            if (i==0){
                h1D[j][i]->SetAxisRange(0.0,40.0,"Y");
            }
            if (i==1) {
                h1D[j][i]->SetAxisRange(0.15,0.35,"Y");
                pol = new TF1("pol","pol2",0.5,3.0);
                h1D[j][i]->Fit("pol","R");
                fout<<pol->GetParameter(0)<<"\t"<<pol->GetParameter(1)<<"\t"<<pol->GetParameter(2)<<"\t";
            }
            if (i==2) {
                par[0] = 10.0;
                par[1] = 0.1;
                par[2] = 0.1;
                h1D[j][i]->SetAxisRange(0.0,0.3,"Y");
                sig = new TF1("sig",SigmaFit,0.5,3.0,2);
                sig->SetParameters(par);
                h1D[j][i]->Fit("sig","R");
                fout<<sig->GetParameter(0)<<"\t"<<sig->GetParameter(1)<<endl;
            }
            h1D[j][i]->Draw();
        }
       	sprintf(OutCan,"Analyze_ECvsP_S%i_%s.gif",j+1,target);
        can[j]->Print(OutCan);
        sprintf(OutCan,"Analyze_ECvsP_S%i_%s.eps",j+1,target);
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
    Int_t nmax = 2;
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
    return polFit(x,par) - 2.0*SigmaFit(x,&par[3]);
}

// Sum of Mean and Sigma functions
// peak uses par[0-2]
// background uses par[3-5]
Double_t CutAbove(Double_t *x, Double_t *par) {
    return polFit(x,par) + 2.5*SigmaFit(x,&par[3]);
}

//
// Plot_ECvsP - plot the EC/P vs P distributions with cuts
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_ECvsP(char *fAna="Ana.root", char *fitParams="ECvsP_Fit.dat", char *target)
{
    Int_t i, j;
    Int_t iSector;
    
    const Int_t NPARAM = 3;
    TF1 *fitMean[NSECTORS];
    TF1 *fitCut[NSECTORS];
    TH2D *h2D[NSECTORS];
    
    char funcName[50];
    
    Double_t a, b, c, d, f;
    Double_t parMean[5];
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
        sprintf(hname,"ECtotP_VS_P_Sector%i",j+1);
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
        
        fin >> iSector >> a >> b >> c >> d >> f;
        parMean[0] = a;
        parMean[1] = b;
        parMean[2] = c;
        parMean[3] = d;
        parMean[4] = f;
        
        sprintf(funcName,"fitMean%i",j+1);
        fitMean[j] = new TF1(funcName,polFit,0.0,3.0,3);
        fitMean[j]->SetParameters(&parMean[0]);
        fitMean[j]->SetLineWidth(2);
        fitMean[j]->Draw("same");
        
        sprintf(funcName,"fitCut%i",j+1);
        fitCut[j] = new TF1(funcName,CutBelow,0.0,3.0,5);
        fitCut[j]->SetParameters(parMean);
        fitCut[j]->SetLineColor(4);
        fitCut[j]->SetLineWidth(2);
        fitCut[j]->Draw("same");
        
        sprintf(funcName,"fitCut%i",j+1);
        fitCut[j] = new TF1(funcName,CutAbove,0.0,3.0,5);
        fitCut[j]->SetParameters(parMean);
        fitCut[j]->SetLineColor(4);
        fitCut[j]->SetLineWidth(2);
        fitCut[j]->Draw("same");
    }
    sprintf(OutCan,"Plot_ECvsP_%s.gif",target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_ECvsP_%s.eps",target);
    c1->Print(OutCan);
    
    fin.close();
}

//
// Plot_ECvsP_1Sector - plot the EC/P vs P distributions with cuts for one sector
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_ECvsP_1Sector(char *fAna="Ana.root", char *fitParams="ECvsP_Fit.dat", int mySector=1, char *target="C12")
{
    Int_t i, j;
    Int_t iSector;
    
    const Int_t NPARAM = 3;
    TF1 *fitMean;
    TF1 *fitBelow;
    TF1 *fitAbove;
    TH2D *h2D;
    
    char funcName[50];
    
    Double_t a, b, c, d, f;
    Double_t parMean[5];
    Double_t parSigma[2];
    
    // open text file for fit parameters
    ifstream fin(fitParams);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    
    sprintf(hname,"ECtotP_VS_P_Sector%i",mySector);
    h2D = (TH2D*)tmp->Get(hname);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
        
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetYaxis()->SetTitleOffset(yoff);
        
    h2D->SetAxisRange(0.0,3.0,"X");
    h2D->Draw("colz");
    
    for(j=0; j<NSECTORS; j++){
        fin >> iSector >> a >> b >> c >> d >> f;
        if(iSector==mySector){
            parMean[0] = a;
            parMean[1] = b;
            parMean[2] = c;
            parMean[3] = d;
            parMean[4] = f;
        }
    }
    
    fitMean = new TF1("fitMean",polFit,0.0,3.0,3);
    fitMean->SetParameters(&parMean[0]);
    fitMean->SetLineWidth(2);
    fitMean->Draw("same");
    
    fitBelow = new TF1("fitBelow",CutBelow,0.0,3.0,5);
    fitBelow->SetParameters(parMean);
    fitBelow->SetLineColor(1);
    fitBelow->SetLineWidth(2);
    fitBelow->Draw("same");
    
    fitAbove = new TF1("fitAbove",CutAbove,0.0,3.0,5);
    fitAbove->SetParameters(parMean);
    fitAbove->SetLineColor(1);
    fitAbove->SetLineWidth(2);
    fitAbove->Draw("same");

    sprintf(OutCan,"Plot_ECvsP_S%i_%s.gif",mySector,target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_ECvsP_S%i_%s.eps",mySector,target);
    c1->Print(OutCan);
    
    fin.close();
}

//
// PlotECinVsECout - plot histogram with labels for electrons
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void PlotECinVsECout(char *fAna="Ana.root", char *target)
{
    Int_t i;
    TH2D *hist;
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");

    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"ECin_VS_ECout_elecID_00");
    hist = (TH2D*)tmp->Get(hname);
    sprintf(htitle,"%s Runs",target);
    hist->SetTitle(htitle);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleOffset(yoff);
    hist->Draw();
    
    TLine *showCut = new TLine(0.06,0.0,0.06,0.25);
    showCut->SetLineColor(2);
    showCut->SetLineWidth(2);
    showCut->Draw("same");
    
    sprintf(OutCan,"Plot_ECinVsECout_%s.gif",target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_ECinVsECout_%s.eps",target);
    c1->Print(OutCan);
}

//
// Plot_dtECSC - plot histogram with labels for electrons
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_dtECSC(char *fAna="Ana.root", char *target)
{
    Int_t i;
    TH2D *h2D;
    TH1D *h1D_withCut;
    TH1D *h1D_withoutCut;
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"dtime_ECSC_elecID");
    h2D = (TH2D*)tmp->Get(hname);

    h1D_withoutCut = (TH1D*)h2D->ProjectionX("h1D_withoutCut",1,1,"");
    sprintf(htitle,"%s Runs",target);
    h1D_withoutCut->SetTitle(htitle);
    h1D_withoutCut->GetXaxis()->CenterTitle();
    h1D_withoutCut->GetYaxis()->CenterTitle();
    h1D_withoutCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_withoutCut->SetLineWidth(2);
    h1D_withoutCut->Draw();
    
    h1D_withCut = (TH1D*)h2D->ProjectionX("h1D_withCut",9,9,"");
    h1D_withCut->SetFillColor(2);
    h1D_withCut->Draw("same");
    
    sprintf(OutCan,"Plot_dtECSC_%s.gif",target);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_dtECSC_%s.eps",target);
    c1->Print(OutCan);
}

//
// Plot_Mom - plot histogram with labels for electrons
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_Mom(char *fAna="Ana.root", char *target)
{
    Int_t i;
    TH2D *h2D;
    TH1D *h1D_withCut;
    TH1D *h1D_withoutCut;
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("ElectronID");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"Mom_elecID");
    h2D = (TH2D*)tmp->Get(hname);
    
//    h1D_withoutCut = (TH1D*)h2D->ProjectionX("h1D_withoutCut",1,1,"");
    h1D_withoutCut = (TH1D*)h2D->ProjectionX("h1D_withoutCut",3,3,"");
    sprintf(htitle,"%s Runs",target);
    h1D_withoutCut->SetTitle(htitle);
    h1D_withoutCut->GetXaxis()->CenterTitle();
    h1D_withoutCut->GetYaxis()->CenterTitle();
    h1D_withoutCut->GetYaxis()->SetTitleOffset(yoff);
    h1D_withoutCut->SetLineWidth(2);
    h1D_withoutCut->Draw();
    
//    h1D_withCut = (TH1D*)h2D->ProjectionX("h1D_withCut",2,2,"");
//    h1D_withCut->SetFillColor(2);
//    h1D_withCut->Draw("same");
    
    sprintf(OutCan,"Plot_Mom_%s.gif",target);
//    c1->Print(OutCan);
    sprintf(OutCan,"Plot_Mom_%s.eps",target);
//    c1->Print(OutCan);
}

//
// Plot_ECmoments - plot histogram for EC moments
//
//                  fAna = output from eg2a DMS
//                  target = target name
//
void Plot_ECmoments(char *fAna="Ana.root", char *target="C", int iPart=0, int iCuts=0)
{
    Int_t i;
    TH2D *h2D[3];
    TH1D *h1D[3];
    char h1Dname[50];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1800,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(1111);
    c1->SetFillStyle(4000);
    c1->Divide(3,1);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Detectors");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    for (i=0; i<3; i++) {
        c1->cd(i+1);
        switch(iCuts){
            case 0: sprintf(hname,"EChit_M%i",i+2); break;
            case 1: sprintf(hname,"EChit_M%i_cuts",i+2); break;
            default: cout<<"Unknown cut index, "<<iCuts<<endl; exit(0); break;
        }

        h2D[i] = (TH2D*)tmp->Get(hname);
    
        sprintf(h1Dname,"%s_px",hname);
        h1D[i] = (TH1D*)h2D[i]->ProjectionX(h1Dname,iPart+1,iPart+1,"");
        sprintf(htitle,"%s Runs, %s",target,Detlabel[iPart]);
        h1D[i]->SetTitle(htitle);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->SetLineWidth(2);
        h1D[i]->Draw();
    }
    
    sprintf(OutCan,"Plot_ECmoments_%s_%i_%i.gif",target,iPart,iCuts);
    c1->Print(OutCan);
    sprintf(OutCan,"Plot_ECmoments_%s_%i_%i.eps",target,iPart,iCuts);
    c1->Print(OutCan);
}
