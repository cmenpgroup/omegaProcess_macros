// 
// MultiplicityRatio - 
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
Float_t yoff = 1.5;  // set the offset between the y-axis label and the axis values

void Plot_MultiplictyRatio(string yldFile){

    Int_t i, j, k;
    Double_t val1, val2;
    const Double_t num = 3;
    Double_t Tgt[3];
    Double_t errTgt[3];
    Double_t Yield[3][2];
    Double_t errYield[3][2];
    Double_t Ratio[3];
    Double_t errRatio[3];
    Double_t err1, err2;
    
    string run, target;
    
    // open text file for the yields
    ifstream fin(yldFile.c_str());

    while(!fin.eof()){
        fin>>run>>target>>val1>>val2;

        if(target.compare("LD2")==0) j = 0;
        if(target.compare("Nuc")==0) j = 1;
        
        if(run.compare("C12")==0){
            Tgt[0] = pow(12,1./3.);
            errTgt[0] = 0;
            Yield[0][j] = val1;
            errYield[0][j] = val2;
        }elseif(run.compare("Fe56")==0){
            Tgt[1] = pow(56,1./3.);
            errTgt[1] = 0;
            Yield[1][j] = val1;
            errYield[1][j] = val2;
        }elseif(run.compare("Pb208")==0){
            Tgt[2] = pow(208,1./3.);
            errTgt[2] = 0;
            Yield[2][j] = val1;
            errYield[2][j] = val2;
        }else{
            cout<<"Unknown target "<<target<<endl;
        }
    }
    
    for(k=0; k<3; k++){
        if(Yield[k][0]>0 && Yield[k][1]>0){
            Ratio[k] = Yield[k][1]/Yield[k][0];
            err1 = errYield[k][0]/Yield[k][0];
            err2 = errYield[k][1]/Yield[k][1];
            errRatio[k] = Ratio[k]*sqrt(err1*err1 + err2*err2);
        }
    }

    gStyle->SetOptStat(0);
	TCanvas *can1 = new TCanvas("can1","Multiplicity Ratio",0,0,600,600); // create the canvas
    can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
    can1->SetBorderSize(5);
    can1->SetFillStyle(4000);

    // set up the Pad parameters
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TGraphErrors *grRatio = new TGraphErrors(num,Tgt,Ratio,errTgt,errRatio);
    grRatio->Draw("AP");
    grRatio->SetTitle(0);
    grRatio->SetMarkerStyle(21);
    grRatio->SetMarkerColor(4);
    grRatio->GetXaxis()->SetTitle("A^{1/3}");
    grRatio->GetXaxis()->CenterTitle();
    grRatio->GetYaxis()->SetTitle("uncorrected Yield(A)/Yield(D_{2})");
    grRatio->GetYaxis()->CenterTitle();
    grRatio->GetYaxis()->SetTitleOffset(yoff);
    
    // create the image files
    char OutCan[100];
    sprintf(OutCan,"Plot_MultiplicityRatio.gif");
    can1->Print(OutCan);
    sprintf(OutCan,"Plot_MultiplicityRatio.eps");
    can1->Print(OutCan);
    
    fin.close();
    
}
