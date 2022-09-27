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

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TString  sInstrumentResponseFile = "VTS.sensitivity.V6.hard.root";

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\


void pulsarLIV(){ 

	int NumberOfDataSets = 1;                     //Number of datasets to include in OF.

	bool SecondOrder = false;                      //True for second order dispersion term, False for first order dispersion term.
	bool LIV         = false;                       //Include LIV correction or not.
	bool SupLum      = true;                       //Super Luminar assumption if true. Corresponds to a -/+ in deltat.

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	 cout<<"~~~~Reading in the instrument response functions~~~~"<<endl;

        TFile* f = NULL;

        if (!(f = TFile::Open(sInstrumentResponseFile)))
                  {
                  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
                  exit(0);
                  }


        TH2F *hMigrationMatrix = NULL;
        TH1F *hBgRate          = NULL;
        TH1F *hEffAreaTrue     = NULL;
        TH1F *hEffAreaRec      = NULL;

        if (!(hEffAreaRec = (TH1F*)f->Get("EffectiveArea")->Clone("hEffAreaRec")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveArea in the provided Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

  if (!(hEffAreaTrue = (TH1F*)f->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

          if (!(hBgRate = (TH1F*)f->Get("BGRate")->Clone("hBgRate")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram BGRate in the provided root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }

         if (!(hMigrationMatrix = (TH2F*)(f->Get("MigMatrix")->Clone("hMigrationMatrix"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile << endl;
    exit(0);
  }


  cout<<"~~~~Loaded all instrument response function from a CTA file~~~~"<<endl;

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	Double_t ELIV=1e16;
        Double_t Distance=1.511e25;
        Double_t c=2.99792458e8;

	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;
	
	Double_t StartHour      = 15.75;
        Double_t EndHour        = 16;
	Double_t StartObserving = StartHour*60*60;
        Double_t StopObserving  = EndHour*60*60;
        Double_t PulsarPeriod   = 900;

	TF1 *Signal = new TF1("Signal","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) + [3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5]))",0,900);
	
	double par [6] = {9.335, 300, 50, 8.084, 650, 75};

        Signal->SetParameters(par);

	TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);

	TF1 *SignalSpectrum = new TF1("SignalSpectrum","(pow(x/4,-1.96))/(1 + pow(x/4,-1.96+3.52))",0,1000);

	TH1D *h = new TH1D("h","Title",1000,0,24*60*60);
	TH1D *z = new TH1D("z","Title",100,0,1000);

	int meanIterations = 2800;

	char name[14]="pulsardata.txt";

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,meanIterations);

	ofstream file;
        file.open(name);

        for (int n = 0; n < NumberOfDataSets; n++)
        {

	Double_t nonintiterations     = TMath::Abs(sc  ->GetRandom());
	int actualIterations = (int) nonintiterations;

		Double_t *events           = new Double_t [actualIterations];
                Double_t *recenergies      = new Double_t [actualIterations];
                Double_t *thetasquares     = new Double_t [actualIterations];
                int      *ind              = new int      [actualIterations];
                Double_t *rands            = new Double_t [actualIterations];
                Double_t *normareas        = new Double_t [actualIterations];
                Double_t *trueenergies     = new Double_t [actualIterations];
                Double_t *effareas         = new Double_t [actualIterations];

                int RecAreaMaxBin     = hEffAreaRec  -> GetMaximumBin();
                int TrueAreaMaxBin    = hEffAreaTrue -> GetMaximumBin();

                Double_t MaxEffRec  = hEffAreaRec   ->GetBinContent(RecAreaMaxBin);
                Double_t MaxEffTrue = hEffAreaTrue  ->GetBinContent(TrueAreaMaxBin);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\




	int k=0;

	while (k < actualIterations)
		{

		Double_t Rand1    = rand3->Rndm();	
		Double_t Rand2    = rand3->Rndm();

		Double_t Observed = Rand1*24*60*60;
	
		if(Observed > StartObserving && Observed < StopObserving)
		{
	
		int WhatPeriod      = (int) Rand2*((StopObserving - StartObserving)/PulsarPeriod) 
		
//		Double_t arrivaltime = TMath::Abs(Signal->GetRandom());
//		Double_t ArrivalTime = StartObserving + (WhatPeriod*PulsarPeriod); //+ arrivaltime

		Double_t energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
		Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
		z->Fill(energy);

		if(SupLum)
                        {

                          if(SecondOrder)
                                {
                                Double_t deltat = -(Distance/c)*pow(energy/ELIV,2);
                                }
                          else
                                {
                                Double_t deltat = -(Distance/c)*(energy/ELIV);
                                }
                        }

                else
                        {

                        if(SecondOrder)
                                {
                                Double_t deltat = +(Distance/c)*pow(energy/ELIV,2);
                                }
                        else
                                {
                                Double_t deltat = +(Distance/c)*(energy/ELIV);
                                 }
                        }


		Double_t AssignRand = rand3->Rndm();
		
		Double_t LogE = log10(energy);

                int      energybin       = hEffAreaTrue->FindBin(LogE,0,0);
                Double_t energybincenter = hEffAreaTrue->GetBinCenter(energybin);


		if (LogE > energybincenter)
                        {

                        int rightbin = energybin + 1;

			Double_t AreaA = hEffAreaTrue->GetBinContent(energybin);
                        Double_t AreaB = hEffAreaTrue->GetBinContent(rightbin );
                        Double_t LogEA = hEffAreaTrue->GetBinCenter (energybin);
                        Double_t LogEB = hEffAreaTrue->GetBinCenter (rightbin );

                        }

		 else
                        {
                        int leftbin = energybin -1;

                        Double_t AreaA = hEffAreaTrue->GetBinContent(leftbin  );
                        Double_t AreaB = hEffAreaTrue->GetBinContent(energybin);
                        Double_t LogEA = hEffAreaTrue->GetBinCenter (leftbin  );
                        Double_t LogEB = hEffAreaTrue->GetBinCenter (energybin);

                        }

		Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
                Double_t intercept = AreaA - (gradient*LogEA);

                Double_t EffArea           = (intercept + (gradient*LogE));
                Double_t NormalisedEffArea = EffArea/MaxEffTrue;

		 if(AssignRand < NormalisedEffArea)
                        {
			
			int trueenergybin = hMigrationMatrix->GetYaxis()->FindBin(LogE);
			hMigrationMatrix->ProjectionX("xproj",trueenergybin,trueenergybin);

			Double_t LogRecEnergy = xproj->GetRandom();
                        Double_t ERec = pow(10,LogRecEnergy);

		if(LIV)
                                {
                                events[k]       = (ArrivalTime + deltat);
                                }
                        else
                                {
                                events[k]       = (ArrivalTime);
                                }

		 recenergies[k]  = ERec;
                 thetasquares[k] = thetasq;
                 rands[k]        = AssignRand;
                 normareas[k]    = NormalisedEffArea;
                 effareas[k]     = EffArea;
                 trueenergies[k] = energy;

		h->Fill(ArrivalTime);

		++k;
		}
		}	
		}


	TMath::Sort(totalIterations,events,ind,false);

        for(int i=0; i < totalIterations; i++)
                {

                 file << fixed << setprecision(6) << events[ind[i]] << "           " << fixed << setprecision(3) << recenergies[ind[i]] << "           "  << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << fixed << setprecision(3) << effareas[ind[i]] << "           " << fixed << setprecision(3) << normareas[ind[i]] <<  "\n";
                }

        file << "#" << "\n";

        printf("Data Set Complete # %g\n",n+1);

        }

        file.close();




	z->Draw();

	}
