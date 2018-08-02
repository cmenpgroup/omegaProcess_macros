// 
// FitOmega - background subtraction routine where the peak is fit with a well-known function
//			   like a gaussian.  The backgorund is fit with another well-known function like a
//             polynomial.
//
//            The output is a list of yields.
//
//
// M. H. Wood, Canisius College
//--------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// gROOT->Reset();

Int_t lcol[10] = {1,2,4,6,7,8,9,12,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

const Int_t MAX_HIST = 3;
const Int_t MAX_RUN = 4;
const Int_t MAX_TGT = 3;
const Int_t MAX_CUT = 21;

Double_t gaussFit(Double_t *x, Double_t *par);
Double_t breitwigner(Double_t *x, Double_t *par);
Double_t polFit(Double_t *x, Double_t *par);
Double_t totFit(Double_t *x, Double_t *par) ;

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

Float_t FitOmega(TH1D *hist, char *plotFilePrefix)
{
	Float_t fPeak, fWidth, SumLo, SumHi;
	Float_t xLo = 0.03;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.5;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.13; // lower limit on the peak range
	Float_t PeakHi = 0.16; // upper limit on the peak range
	Float_t Nsigma = 3.0;  // number of sigma of the peak
	
	Int_t ibg = 3; // parameter index in the par[] array where the background parameters start

	Int_t i, j, k;
	Float_t x, pval;
	Float_t yield;
	
    cout << "Fitting histogram " << hist->GetName() << endl;
    
	TH1D *hBgFit; // histogram of background
		
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
	hBgFit = (TH1D*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
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

	return yield;
}

Float_t FitOmega_justBgd(TH1D *hist, char *plotFilePrefix, Int_t fitIndex, Int_t iClose)
{
    Float_t XrangeLo = 0.0;
	Float_t XrangeHi = 2.0;

	Float_t PeakLo = 0.7;
	Float_t PeakHi = 0.875;

	Float_t fitLo = 0.4; // C12
	Float_t fitHi = 1.15; // C12

//    Float_t fitLo = 0.45; // Pb208
//	Float_t fitHi = 1.25; // Pb208
    
	Int_t i, j, k;
	Float_t x, pval, yield;
	Double_t NDF, ChiSq, ChiSqPerNDF;
    
    char pName[50];
    char pFunc[50];
    char hname[50];
    
    char title[100];
	char xtitle[100];
	char ytitle[100];
    
    if(fitIndex<0){
		cout<<"Invalid polynomial fit of order "<<fitIndex<<endl;
		exit(0);
	}
    
    sprintf(hname,"hAfterBgSub%i",j);
    TH1D *hAfterBgSub = (TH1D*)hist->Clone(hname);
    
    sprintf(hname,"hNoPeak%i",j);
    TH1D *hNoPeak = (TH1D*)hist->Clone(hname);
    
    for(k=1;k<=hNoPeak->GetNbinsX();k++){
        x = hNoPeak->GetBinCenter(k+1);
        if(x>=PeakLo && x<PeakHi){
            hNoPeak->SetBinContent(k,0.0);
            hNoPeak->SetBinError(k,0.0);
        }
    }
    
    // fit the final spectrum
    Double_t par[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    
    sprintf(pName,"pol%i",fitIndex);
    sprintf(pFunc,"pol%i",fitIndex);
    TF1 *pol = new TF1(pName,pFunc,fitLo,fitHi);
    pol->SetLineColor(lcol[1]);
    
    hNoPeak->Fit(pName,"R");
    //    pol->GetParameters(par);
    //    pol->SetParameters(par);
    NDF = pol->GetNDF(); // get number of degrees of freedom
    ChiSq = pol->GetChisquare(); // get chi squared
    if(NDF){
        ChiSqPerNDF = ChiSq/NDF; // reduced chi squared
    }else{
        ChiSqPerNDF = 0.0;
    }
    cout<<fitIndex<<"\t"<<NDF<<"\t"<<ChiSq<<"\t"<<ChiSqPerNDF<<endl;
    
    hAfterBgSub->Add(pol,-1.0);
    for(k=1; k<hAfterBgSub->GetNbinsX(); k++){
        x = hAfterBgSub->GetBinCenter(k);
        if(x>=fitLo && x<=fitHi){
            pval = hAfterBgSub->GetBinContent(k);
            if(pval<0.0) pval = 0.0;
        }else{
            pval =0.0;
        }
        hAfterBgSub->SetBinContent(k,pval);
        if(pval==0.0) hAfterBgSub->SetBinError(k,0.0);
    }
    yield=hAfterBgSub->Integral(hAfterBgSub->FindBin(PeakLo),hAfterBgSub->FindBin(PeakHi));

    sprintf(title,"#omega meson, Just Background Subtraction, Canvas 1");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);

    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    hist->GetXaxis()->SetRangeUser(XrangeLo,XrangeHi);
    hist->GetYaxis()->SetTitleOffset(yoff);
    hist->Draw();
    pol->Draw("same");

    char OutCan[100];
    sprintf(OutCan,"%s-Before.gif",plotFilePrefix);
    can1->Print(OutCan);
    sprintf(OutCan,"%s-Before.eps",plotFilePrefix);
    can1->Print(OutCan);

    sprintf(title,"#omega meson, Just Background Subtraction, Canvas 2");
    TCanvas *can2 = new TCanvas("can2",title,100,0,600,600);

    can2->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can2->SetBorderSize(5);
    can2->SetFillStyle(4000);

    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    hAfterBgSub->GetXaxis()->SetRangeUser(fitLo,fitHi);
    hAfterBgSub->GetYaxis()->SetTitleOffset(yoff);
    hAfterBgSub->Draw();
    
    TLine *lLow = new TLine(PeakLo,0.0,PeakLo,0.5*hAfterBgSub->GetMaximum());
    lLow->SetLineWidth(2);
    lLow->Draw();
    TLine *lHigh = new TLine(PeakHi,0.0,PeakHi,0.5*hAfterBgSub->GetMaximum());
    lHigh->SetLineWidth(2);
    lHigh->Draw();
    
    sprintf(OutCan,"%s-After.gif",plotFilePrefix);
    can2->Print(OutCan);
    sprintf(OutCan,"%s-After.eps",plotFilePrefix);
    can2->Print(OutCan);

    if(iClose){
        can1->Close();
        can2->Close();
    }
	return yield;
}

void FitOmega_OneTarget(char *rootFile, Int_t runIndex=0, Int_t histIndex=0, Int_t tgtIndex=0, Int_t cutIndex=0, Int_t fitIndex=3, Int_t iClose=1)
{
	Int_t i;
	Float_t yield;
    
    Cuts myCuts;
    
    Check_HistIndex(histIndex);
    Check_TgtIndex(tgtIndex);
    Check_CutIndex(cutIndex);
    
    char strname[50];
    char hname[50];
	char plotFilePrefix[100];
    
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	
    TDirectory *tmp = fd->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
    
    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    TH2D *h2D = (TH2D*)tmp->Get(hname);
    
    sprintf(strname,"%s_%i",hname,cutIndex);
    TH1D *h1D = (TH1D*)h2D->ProjectionX(strname,cutIndex+1,cutIndex+1,"");

    sprintf(plotFilePrefix,"FitOmega_%s_%s_P%i",myCuts.Get_Run(runIndex).c_str(),hname,fitIndex);
    
    yield = FitOmega_justBgd(h1D,plotFilePrefix, fitIndex, iClose);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout << myCuts.Get_Run(runIndex).c_str() << "\t" << myCuts.Get_Tgt(tgtIndex).c_str() << "\t" << yield << "\t" << sqrt(yield) << endl;

	fout.close(); // close the text file
//	fd->Close();  // close ROOT file
}

void FitOmega_AllTargets(char *rootFile, Int_t runIndex=0, Int_t histIndex=0, Int_t fitIndex=3)
{
	Int_t i, j;

    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            FitOmega_OneTarget(rootFile, runIndex, histIndex, i, j, fitIndex, 1);
        }
    }
}

void FitOmega_BothTargets(char *rootFile, Int_t runIndex=0, Int_t histIndex=0, Int_t cutIndex=0, Int_t fitIndex=3, Int_t iClose=1)
{
    Int_t i;
    Float_t yield;
    
    Cuts myCuts;
    
    Check_HistIndex(histIndex);
    Check_CutIndex(cutIndex);
    
    char strname[50];
    char hname[50];
    char plotFilePrefix[100];
    
    // data files contain the trees
    printf("Analyzing file %s\n",rootFile);

    TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
    TDirectory *tmp[MAX_TGT];
    
    // open text file for the yields
    char OutFile[100];
    sprintf(OutFile,"FitOmega_%s_%sP%i_both.yld",myCuts.Get_Run(runIndex).c_str(), myCuts.Get_Hist(histIndex).c_str(),fitIndex);
    ofstream fout(OutFile);
    
    for(i=2; i<MAX_TGT; i++){
        tmp[i] = fd->GetDirectory(myCuts.Get_Tgt(i).c_str());
    
        sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(i).c_str());
        TH2D *h2D = (TH2D*)tmp[i]->Get(hname);
    
        sprintf(strname,"%s_%i",hname,cutIndex);
        TH1D *h1D = (TH1D*)h2D->ProjectionX(strname,cutIndex+1,cutIndex+1,"");
 
        sprintf(plotFilePrefix,"FitOmega_%s_%s_P%i_both",myCuts.Get_Run(runIndex).c_str(),hname,fitIndex);

        yield = FitOmega_justBgd(h1D,plotFilePrefix, fitIndex, iClose);
    
        fout<<myCuts.Get_Run(runIndex).c_str()<<"\t"<<myCuts.Get_Tgt(i).c_str()<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;
    }
    
    fout.close(); // close the text file
    //	fd->Close();  // close ROOT file
}

void FitOmega_simBgd(string fAna1, string fAna2, Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0, Int_t fitIndex=3)
{
    TFile *fm[2];
    TDirectory *dir[2];
    TH1D *h1D[2];
    TH2D *h2D[2];

    Cuts myCuts;
    
    Float_t XrangeLo = 0.0;
    Float_t XrangeHi = 2.0;

    Float_t EtaLo = 0.5;
    Float_t EtaHi = 0.625;
    Float_t PeakLo = 0.7;
    Float_t PeakHi = 0.9;
    
    Float_t fitLo = 0.4;
    Float_t fitHi = 2.0;
    
    Int_t xBinScale = 110;
    Int_t yBinScale = 150;

    Int_t i, j, k;
    Float_t x, pval, yield;
    Double_t NDF, ChiSq, ChiSqPerNDF;
    
    char strname[100];
    char pName[50];
    char pFunc[50];
    char hname[50];
    
    char title[100];
    char xtitle[100];
    char ytitle[100];
    char OutCan[100];

    sprintf(hname,"%s%s",myCuts.Get_Hist(histIndex).c_str(),myCuts.Get_Tgt(tgtIndex).c_str());
    
    for(i=0; i<2; i++){
        // data files contain the trees
        switch(i){
            case 0:
                printf("Analyzing file %s\n",fAna1.c_str());
                fm[i] = new TFile(fAna1.c_str(),"READ");
                break;
            case 1:
                printf("Analyzing file %s\n",fAna2.c_str());
                fm[i] = new TFile(fAna2.c_str(),"READ");
                break;
            default: cout << "Analysis out of bounds!" << endl;  break;
        }
        
        dir[i] = fm[i]->GetDirectory(myCuts.Get_Tgt(tgtIndex).c_str());
        
        h2D[i] = (TH2D*)dir[i]->Get(hname);
        sprintf(strname,"%s_%i",hname,i);
        h1D[i] = (TH1D*)h2D[i]->ProjectionX(strname,chan+1,chan+1,"");
    }
    
    sprintf(hname,"hAfterBgSub");
    TH1D *hAfterBgSub = (TH1D*)h1D[0]->Clone(hname);
    
    sprintf(hname,"hNoPeak");
    TH1D *hNoPeak = (TH1D*)h1D[1]->Clone(hname);
    
    for(k=1;k<=hNoPeak->GetNbinsX();k++){
        x = hNoPeak->GetBinCenter(k+1);
        if((x>=EtaLo && x<EtaHi) || (x>=PeakLo && x<PeakHi)){
            hNoPeak->SetBinContent(k,0.0);
            hNoPeak->SetBinError(k,0.0);
        }
    }

    hNoPeak->Scale(h1D[0]->Integral(xBinScale,yBinScale)/hNoPeak->Integral(xBinScale,yBinScale));
    sprintf(title,"#omega meson, Sim Background Subtraction, Canvas 1");
    TCanvas *can1 = new TCanvas("can1",title,0,0,600,600);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    h1D[0]->GetYaxis()->SetTitleOffset(yoff);
    h1D[0]->Draw();
    hNoPeak->SetLineColor(lcol[1]);
    hNoPeak->Draw("same");
    
    sprintf(OutCan,"FitOmega_simBgd-Before.gif");
    can1->Print(OutCan);
    sprintf(OutCan,"FitOmega_simBgd-Before.eps");
    can1->Print(OutCan);
    
    // fit the final spectrum
    Double_t par[10]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    
    sprintf(pName,"pol%i",fitIndex);
    sprintf(pFunc,"pol%i",fitIndex);
    TF1 *pol = new TF1(pName,pFunc,fitLo,fitHi);
    pol->SetLineColor(lcol[1]);
    
    hNoPeak->Fit(pName,"R");
    //    pol->GetParameters(par);
    //    pol->SetParameters(par);
    NDF = pol->GetNDF(); // get number of degrees of freedom
    ChiSq = pol->GetChisquare(); // get chi squared
    if(NDF){
        ChiSqPerNDF = ChiSq/NDF; // reduced chi squared
    }else{
        ChiSqPerNDF = 0.0;
    }
    cout<<fitIndex<<"\t"<<NDF<<"\t"<<ChiSq<<"\t"<<ChiSqPerNDF<<endl;
    
    hAfterBgSub->Add(pol,-1.0);
    for(k=1; k<hAfterBgSub->GetNbinsX(); k++){
        x = hAfterBgSub->GetBinCenter(k);
        if(x>=fitLo && x<=fitHi){
            pval = hAfterBgSub->GetBinContent(k);
            if(pval<0.0) pval = 0.0;
        }else{
            pval =0.0;
        }
        hAfterBgSub->SetBinContent(k,pval);
        if(pval==0.0) hAfterBgSub->SetBinError(k,0.0);
    }
    
    sprintf(title,"#omega meson, Sim Background Subtraction, Canvas 2");
    TCanvas *can2 = new TCanvas("can2",title,100,0,600,600);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    can2->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can2->SetBorderSize(5);
    can2->SetFillStyle(4000);
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    hAfterBgSub->GetYaxis()->SetTitleOffset(yoff);
    hAfterBgSub->Draw();
    
    sprintf(OutCan,"FitOmega_simBgd-After.gif");
    can2->Print(OutCan);
    sprintf(OutCan,"FitOmega_simBgd-After.eps");
    can2->Print(OutCan);
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

