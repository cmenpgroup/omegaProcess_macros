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
#include <iostream.h>

gROOT->Reset();

char *RunName[1] = {"C12"};//,"Sn50","Fe56","Pb208"};
char *TgtName[1] = {/*"NoTarget","LD2",*/"Nuc"};
char *HistName[1] = {"Recon_Pi0_Mass_PhotID_Cuts"};//"IM2Photons","IM2Photons_OpAng_ElecPhoton_Cut"};

Float_t FitPi0(TH1F *hist, char *plotFilePrefix)
{

	Float_t Lmar = 0.125; // set the left margin
	Float_t Rmar = 0.125; // set the right margin
	Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values
	
	Float_t fPeak, fWidth, SumLo, SumHi;
	Float_t xLo = 0.03;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.5;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.11; // lower limit on the peak range
	Float_t PeakHi = 0.16; // upper limit on the peak range
	Float_t Nsigma = 3.0;  // number of sigma of the peak
	
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
    can1->Divide(2,2);
	
	sprintf(xtitle,"Invariant Mass (GeV)"); // set the x-axis title
	sprintf(ytitle,"Counts"); // set the y-axis title
		
	// histogram of background
	hBgFit = (TH1F*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	hBgFit->SetName("hBgFit");
	hBgFit->SetTitle("Background");

	// fit the histogram
	Double_t par[7]={9000,0.135,1.0,1.0,1.0,1.0,1.0};
    TF1 *g1 = new TF1("g1",gaussFit,PeakLo,PeakHi,3); // declare fit fcn
    TF1 *pol = new TF1("pol",polFit,xLo,xHi,4);
    TF1 *t1 = new TF1("t1",totFit,xLo,xHi,7);
		
    can1->cd(1);
	g1->SetParameters(&par[0]);    // set parameters for initial peak fit
	hist->Fit("g1","WWR");           // fit the peak
	g1->GetParameters(&par[0]);    // get parameters from initial peak fit
    g1->Draw("same");

    can1->cd(2);
	hist->Fit("pol","WWR+");         // fit the background
	pol->GetParameters(&par[ibg]); // get parameters fromt background fit
    pol->Draw("same");

	t1->SetParameters(par);       // set the parameters from initial fits for total fit
	t1->SetLineWidth(5);          // make fit line thicker
	t1->SetLineColor(1);          // set the fit line color
    can1->cd(3);
	hist->Fit("t1","WWR");          // fit spectrum with total fit function
	t1->GetParameters(&par[0]);   //get final parameters
	
    can1->cd();
	fPeak = t1->GetParameter(1);  // get peak centroid
	fWidth = t1->GetParameter(2); // get peak width
	SumLo = fPeak - Nsigma*fWidth; // calc lower edge at Nsigma away from centroid
	SumHi = fPeak + Nsigma*fWidth; // calc upper edge at Nsigma away from centroid
    cout << "fPeak:  " << fPeak << endl;
			
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

    //g1, pol, t1 TODO
    Double_t binWidth = hist->GetXaxis()->GetBinWidth(0);
    cout << "Bin width:  " << binWidth << endl;
    cout << "SumLo:  " << SumLo << "\nSumHi:  " << SumHi << endl;
    TCanvas *can2 = new TCanvas("can2","Fit Functions",0,0,600,600); // create the canvas
    can2->Divide(2);
    can2->cd(1);
    g1->GetXaxis()->SetRangeUser(SumLo,SumHi);
    g1->Draw();
    can2->cd(2);
    pol->GetXaxis()->SetRangeUser(SumLo,SumHi);
    pol->Draw();
    Double_t peakSum = g1->Integral(SumLo,SumHi);
    Double_t bgSum = pol->Integral(SumLo,SumHi);
    cout << "Counts above peak:  " << (peakSum-bgSum)/binWidth << endl;
    cout << "Counts below peak:  " << bgSum/binWidth << endl;

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
	hist->GetXaxis()->SetRangeUser(xLo,xHi);  // set the x-axis range for the plot
	hist->SetLineWidth(2);
	hist->SetMinimum(0); // start the y-axis at zero.
//	hist->Draw();
	
	// draw the background-subtracted histogram
	hBgFit->SetLineWidth(2);
	hBgFit->SetFillColor(4);
//	hBgFit->Draw("same");

	// create the image files
	char OutCan[100];
	sprintf(OutCan,"%s.gif",plotFilePrefix);
	can1->Print(OutCan);
	sprintf(OutCan,"%s.eps",plotFilePrefix);
	can1->Print(OutCan);

	//can1->Close();

	return yield;
}

void FitPi0_OneTarget(char *rootFile, char* fDir="ReconCuts", Int_t iRun=0, Int_t iTgt=0, Int_t iCut=0)
{
	Int_t i;
	Float_t yield;

	TH1F *hist; // original histogram

    char hname[50];
	char plotFilePrefix[100];
    
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
    
    TDirectory *tmp = fd->GetDirectory(fDir);
	
    sprintf(hname,"%s",HistName[iCut]);
    cout << hname << endl;
	hist = (TH1F*)tmp->Get(hname); // get the histogram from the ROOT file

    sprintf(plotFilePrefix,"FitPi0_%s_%s",RunName[iRun],hname);
    
    yield = FitPi0(hist,plotFilePrefix);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",plotFilePrefix);
	ofstream fout(OutFile);
    
	fout<<TgtName<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

	fout.close(); // close the text file
//	fd->Close();  // close ROOT file
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

