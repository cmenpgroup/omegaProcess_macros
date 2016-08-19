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

Float_t FitOmega_justBgd(TH1D *hist, char *plotFilePrefix, Int_t iFit, Int_t iClose)
{
    Float_t XrangeLo = 0.0;
	Float_t XrangeHi = 2.0;

	Float_t PeakLo = 0.7;
	Float_t PeakHi = 0.875;

//	Float_t fitLo = 0.4; // C12
//	Float_t fitHi = 1.15; // C12

    Float_t fitLo = 0.45; // Pb208
	Float_t fitHi = 1.25; // Pb208
    
	Int_t i, j, k;
	Float_t x, pval, yield;
	Double_t NDF, ChiSq, ChiSqPerNDF;
    
    char pName[50];
    char pFunc[50];
    char hname[50];
    
    char title[100];
	char xtitle[100];
	char ytitle[100];
    
    if(iFit<0){
		cout<<"Invalid polynomial fit of order "<<iFit<<endl;
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
    
    sprintf(pName,"pol%i",iFit);
    sprintf(pFunc,"pol%i",iFit);
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
    cout<<iFit<<"\t"<<NDF<<"\t"<<ChiSq<<"\t"<<ChiSqPerNDF<<endl;
    
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

void FitOmega_OneTarget(char *rootFile, Int_t iRun=0, Int_t iTgt=0, Int_t iCut=0, Int_t iFit=3, Int_t iClose=1)
{
	Int_t i;
	Float_t yield;

	TH1D *hist; // original histogram

    char hname[50];
	char plotFilePrefix[100];
    
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	
    sprintf(hname,"%s_%s",HistName[iCut],TgtName[iTgt]);
    cout << hname << endl;
	hist = (TH1D*)fd->Get(hname); // get the histogram from the ROOT file

    sprintf(plotFilePrefix,"FitOmega_%s_%s_P%i",RunName[iRun],hname,iFit);
    
    yield = FitOmega_justBgd(hist,plotFilePrefix, iFit, iClose);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout<<RunName[iRun]<<"\t"<<TgtName[iTgt]<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

	fout.close(); // close the text file
//	fd->Close();  // close ROOT file
}

void FitOmega_AllTargets(char *rootFile, Int_t iRun=0, Int_t iFit=3)
{
	Int_t i, j;

    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            FitOmega_OneTarget(rootFile, iRun, i, j, iFit, 1);
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

