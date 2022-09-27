#include <TH1.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>

void fitpoly()
{
	const char *files[7] = {"0_01z.txt", "0_05z.txt", "0_1z.txt", "0_2z.txt", "0_3z.txt", "0_4z.txt", "0_5z.txt"};
	double redshifts[7] = {0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50};

	TF1 **fullfits = new TF1*[7];

	TCanvas *c1 = new TCanvas;
	c1->SetLogx();
	c1->SetLogy();
	c1->DrawFrame(0.09,0.001,150,150);

	for (int i = 0; i < 7; i++)
	{
		TGraph *plot = new TGraph(files[i]); 

		TF1 *p91    = new TF1("p91","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) + [6]*pow(x,6) + [7]*pow(x,7) + [8]*pow(x,8) + [9]*pow(x,9)",0.1,0.2);

		TF1 *p92    = new TF1("p92","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) + [6]*pow(x,6) + [7]*pow(x,7) + [8]*pow(x,8) + [9]*pow(x,9)",0.2,1);

		TF1 *p93    = new TF1("p93","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) + [6]*pow(x,6) + [7]*pow(x,7) + [8]*pow(x,8) + [9]*pow(x,9)",1,10);

		TF1 *p94    = new TF1("p94","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) + [6]*pow(x,6) + [7]*pow(x,7) + [8]*pow(x,8) + [9]*pow(x,9)",10,110);

		plot->Fit(p91,"R");
		plot->Fit(p92,"R+");
		plot->Fit(p93,"R+");
		plot->Fit(p94,"R+");

		Double_t par[40];

		p91->GetParameters(&par[0]);
		p92->GetParameters(&par[10]);
		p93->GetParameters(&par[20]);
		p94->GetParameters(&par[30]);

		
                char fname[20];
                sprintf(fname,"fit%d",i);
		
		fullfits[i] = new TF1(fname,"(x>=0.1 && x<=0.2)*([0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) + [6]*pow(x,6) + [7]*pow(x,7) + [8]*pow(x,8) + [9]*pow(x,9)) + (x>0.2  && x<=1)*([10] + [11]*pow(x,1) + [12]*pow(x,2) + [13]*pow(x,3) + [14]*pow(x,4) + [15]*pow(x,5) + [16]*pow(x,6) + [17]*pow(x,7) + [18]*pow(x,8) + [19]*pow(x,9)) + (x>1 && x<=10)*([20] + [21]*pow(x,1) + [22]*pow(x,2) + [23]*pow(x,3) + [24]*pow(x,4) + [25]*pow(x,5) + [26]*pow(x,6) + [27]*pow(x,7) + [28]*pow(x,8) + [29]*pow(x,9)) + (x>10 && x<=110)*([30] + [31]*pow(x,1) + [32]*pow(x,2) + [33]*pow(x,3) + [34]*pow(x,4) + [35]*pow(x,5) + [36]*pow(x,6) + [37]*pow(x,7) + [38]*pow(x,8) + [39]*pow(x,9))" ,0.1 ,110);
			
		fullfits[i]->SetParameters(par);
		
		if (i==0) fullfits[i]->Draw("same");
		else     fullfits[i]->Draw("lsame");

	}
}
