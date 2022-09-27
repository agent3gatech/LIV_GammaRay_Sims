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

	TString  sInstrumentResponseFile = "/nv/hp11/agent3/LIVroot/VTS.sensitivity.V6.hard.root";

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\


void pulsarLIV(){ 

	int NumberOfDataSets = 300;                     //Number of datasets to include in OF.

	bool SecondOrder = false;                      //True for second order dispersion term, False for first order dispersion term.
	bool LIV         = true;                       //Include LIV correction or not.
	bool SupLum      = false;                       //Super Luminar assumption if true. Corresponds to a -/+ in deltat.
	char name[14]="pulsardata.txt";

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
        Double_t Distance=2200*3.0857e16;
        Double_t c=2.99792458e8;

	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;
	
	Double_t StartHour      = 4;
        Double_t EndHour        = 5;
	
	int meanIterations = 2409;	

	Double_t StartObserving = StartHour*60*60;
        Double_t StopObserving  = EndHour*60*60;
        Double_t PulsarPeriod   = 33.5028583/1000;
	
	Double_t HowManyPeriods  = ((StopObserving - StartObserving)/PulsarPeriod);

	Double_t VTSScienceCrabPulseProfile(Double_t *x, Double_t *par) {
     	   par[0]=0;

	
	//Values from published profile

	double amp1=3.41029e+02;
	double mean1=2.97447e-01;
	double sigma1=5.18436e-03;

	double amp2=9.15185e+02;
	double mean2=6.97709e-01;
	double sigma2=1.13830e-02;
	
	double pdf = amp1/(sqrt(2.0 * 3.141592653589793) * sigma1) * exp(-1.0*(x[0]-mean1)*(x[0]-mean1)/(2*sigma1*sigma1));
               pdf+= amp2/(sqrt(2.0 * 3.141592653589793) * sigma2) * exp(-1.0*(x[0]-mean2)*(x[0]-mean2)/(2*sigma2*sigma2));

    return pdf;
}

	TF1 *fPulseProfile = new TF1("fPulseProfile",VTSScienceCrabPulseProfile,0,1,1);
        fPulseProfile->SetNpx(10000);

	TF1 *B = new TF1("B","1",0,1);
        TF1 *H = new TF1("H","1",0,1);

	TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);

	TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"(pow(x/4,-1.96))/(1 + pow(x/4,-1.96+3.52))",0,1000);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-2.5)"                               ,0,1000);
	TF1 *HadronSpectrum     = new TF1("HadronSpectrum"    ,"pow(x,-2.0)"                               ,0,1000);

	TH1D *h = new TH1D("h","Title",1000,0,24*60*60);
	TH1D *z = new TH1D("z","Title",100, 0,1000);
	TH1D *m = new TH1D("m","Title",1000,0,33.5/1000);

	double SignalRatio      = 0.004;
        double BackgroundRatio  = 0.498;
        double HadronRatio      = 0.498;

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

        TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        bc->SetParameter(0,BackgroundRatio*meanIterations);

        TF1 *hc = new TF1("hc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        hc->SetParameter(0,HadronRatio*meanIterations);
	
	ofstream file;
        file.open(name);

        for (int n = 0; n < NumberOfDataSets; n++)
        {

		Double_t nonintsignaliterations     = TMath::Abs(sc  ->GetRandom());
                Double_t nonintbackgrounditerations = TMath::Abs(bc  ->GetRandom());
                Double_t noninthadroniterations     = TMath::Abs(hc  ->GetRandom());

                int actualSignalIterations     = (int) nonintsignaliterations;
                int actualBackgroundIterations = (int) nonintbackgrounditerations;
                int actualHadronIterations     = (int) noninthadroniterations;
                int photonIterations = actualSignalIterations + actualBackgroundIterations;
                int totalIterations  = actualSignalIterations + actualBackgroundIterations + actualHadronIterations;

		Double_t *events           = new Double_t [totalIterations];
                Double_t *recenergies      = new Double_t [totalIterations];
                Double_t *thetasquares     = new Double_t [totalIterations];
                int      *ind              = new int      [totalIterations];
                Double_t *rands            = new Double_t [totalIterations];
                Double_t *normareas        = new Double_t [totalIterations];
                Double_t *trueenergies     = new Double_t [totalIterations];
                Double_t *effareas         = new Double_t [totalIterations];
		Double_t *dts              = new Double_t [totalIterations];
                int      *flags            = new int      [totalIterations];

                int RecAreaMaxBin     = hEffAreaRec  -> GetMaximumBin();
                int TrueAreaMaxBin    = hEffAreaTrue -> GetMaximumBin();

                Double_t MaxEffRec  = hEffAreaRec   ->GetBinContent(RecAreaMaxBin);
                Double_t MaxEffTrue = hEffAreaTrue  ->GetBinContent(TrueAreaMaxBin);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
////=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

		

	int k=0;

	while (k < totalIterations)
		{

		if (k < actualSignalIterations)
                	{
		
		Double_t Rand1    = rand3->Rndm();	
		Double_t Rand2    = rand3->Rndm();
		
		Double_t RandomPeriod    = Rand2*HowManyPeriods;
		int RoundedPeriod        = (int)RandomPeriod;
		
		Double_t arrivalphase = TMath::Abs(fPulseProfile->GetRandom());
		
		Double_t arrivaltime  = arrivalphase*PulsarPeriod;
		
		Double_t daytime = n*24*60*60;

		Double_t ArrivalTime =(daytime + StartObserving + (RoundedPeriod*PulsarPeriod) + arrivaltime);
	
	
		Double_t energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
		Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
		
			}

			 else
                                {
                                if (k < photonIterations)
                                        {
                
					Double_t Rand1    = rand3->Rndm();
			                Double_t Rand2    = rand3->Rndm();

					Double_t RandomPeriod    = Rand2*HowManyPeriods;
			                int RoundedPeriod        = (int)RandomPeriod;
		
					Double_t arrivalphase = TMath::Abs(B                      ->GetRandom());

					Double_t arrivaltime  = arrivalphase*PulsarPeriod;

					Double_t daytime = n*24*60*60;

					Double_t ArrivalTime =(daytime + StartObserving + (RoundedPeriod*PulsarPeriod) + arrivaltime);

		                        Double_t energy      = TMath::Abs(BackgroundSpectrum     ->GetRandom());
                                        Double_t thetasq     = TMath::Abs(Th1                    ->GetRandom());
					}
		
			else
                                        {
                                        
					Double_t Rand1    = rand3->Rndm();
                                        Double_t Rand2    = rand3->Rndm();

                                        Double_t RandomPeriod    = Rand2*HowManyPeriods;
                                        int RoundedPeriod        = (int)RandomPeriod;

                                        Double_t arrivalphase = TMath::Abs(H                      ->GetRandom());

                                        Double_t arrivaltime  = arrivalphase*PulsarPeriod;

                                        Double_t daytime = n*24*60*60;

                                        Double_t ArrivalTime =(daytime + StartObserving + (RoundedPeriod*PulsarPeriod) + arrivaltime);

					Double_t energy      = TMath::Abs(HadronSpectrum ->GetRandom());
                                        Double_t thetasq     = TMath::Abs(Th1            ->GetRandom());
					}
			}


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
        	 dts[k]         = deltat;
                 trueenergies[k] = energy;

		h->Fill(ArrivalTime);

		  if (k < actualSignalIterations)
                                {
                                flags[k] = 0;
                                }

                         else
                                {

                                if (k < photonIterations)
                                        {
                                        flags[k] = 1;
                                        }
                                else
                                        {
                                        flags[k] = 2;
                                        }
                                }


		++k;
		}
			
		}


	TMath::Sort(totalIterations,events,ind,false);

        file << "#" << "      " << totalIterations << "\n";

	for(int i=0; i < totalIterations; i++)
                {

                  file << setw(11) << fixed << setprecision(6) << events[ind[i]] << "          " << setw(9) << fixed << setprecision(6) << dts[ind[i]] << "           " << setw(6)<< fixed << setprecision(3) << recenergies[ind[i]] << "           " << setw(6) << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << setw(9) << fixed << setprecision(3) << effareas[ind[i]] << "           " << setw(5)<< fixed << setprecision(3) << normareas[ind[i]] <<  "          " << flags[ind[i]] <<   "\n";
                }

        printf("Data Set Complete # %g\n",n+1);

        }

        file.close();



//	m->Draw();
	//z->Draw();

	}
