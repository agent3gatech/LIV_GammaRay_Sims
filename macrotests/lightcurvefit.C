
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
#include <TMath>
#include <TRandom3.h>



void lightcurvefit()
	{

	TF1 *SignalSpectrum = new TF1("SignalSpectrum","(pow(x/4,-1.96))/(1 + pow(x/4,-1.96+3.52))",0,1000);
    

	TH1D *LightCurve = new TH1D("LightCurve", "Lightcurve - PKS2155", 500, 0, 1000);

	int k = 0;

        while (k < 10000000)
	{

	 Double_t arrivaltime = TMath::Abs(SignalSpectrum              ->GetRandom());
	
	LightCurve->Fill(arrivaltime);



	k++;
	}

//	for (int i = 0; i < 500; i++)
//      {
//            Double_t bincontent = LightCurve->GetBinContent(i);
//            LightCurve-->SetBinContent(i,bincontent);
//            }





	TCanvas * c1 = new TCanvas("c1");
//        TCanvas * c2 = new TCanvas("c2");

//	c1->cd();
        
//	c2->cd();
	LightCurve->Draw();
                  





























	}
