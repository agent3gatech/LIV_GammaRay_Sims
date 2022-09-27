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

//=================================================================================================================================================

	 TString  sInstrumentResponseFile1 = "/gpfs/pace1/project/phy-otte/agent3/LIV/IRFs/HESS.PG1553.1.root";
	 TString  sInstrumentResponseFile2 = "/gpfs/pace1/project/phy-otte/agent3/LIV/IRFs/HESS.PG1553.2.root";
	 TString  sInstrumentResponseFile3 = "/gpfs/pace1/project/phy-otte/agent3/LIV/IRFs/HESS.PG1553.3.root";
	 TString  sInstrumentResponseFile4 = "/gpfs/pace1/project/phy-otte/agent3/LIV/IRFs/HESS.PG1553.4.root";

//=================================================================================================================================================

void SIMPG1553() {

	int NumberOfDataSets = 500;		       //Number of datasets to include in OF.

        bool SecondOrder = false;       	       //True for second order dispersion term, False for first order dispersion term.
	bool LIV         = false;                       //Include LIV correction or not.
        bool SupLum      = false;                       //Super Luminal assumption if true. Corresponds to a -/+ in deltat.

	char name[256] = "PG1553NOLIV6.txt";

//=================================================================================================================================================

	cout<<"~~~~Reading in the instrument response functions~~~~"<<endl;

	TFile* f = NULL;

	if (!(f = TFile::Open(sInstrumentResponseFile1)))
		  {
    		  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile1 << endl;
    		  exit(0);
  		  }

	TFile* g = NULL;

        if (!(g = TFile::Open(sInstrumentResponseFile2)))
                  {
                  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile2 << endl;
                  exit(0);
                  }

	TFile* o = NULL;

        if (!(o = TFile::Open(sInstrumentResponseFile3)))
                  {
                  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile3 << endl;
                  exit(0);
                  }

	TFile* p = NULL;

        if (!(p = TFile::Open(sInstrumentResponseFile4)))
                  {
                  cout << "!!!!!!!ERROR: Could not open Root performance file!!!!!!!" << sInstrumentResponseFile4 << endl;
                  exit(0);
                  }


	TH1F *hEffAreaTrue1     = NULL;
	TH1F *hEffAreaTrue2     = NULL;
	TH1F *hEffAreaTrue3     = NULL;
	TH1F *hEffAreaTrue4     = NULL;

	TH2F *hMigrationMatrix1 = NULL;
	TH2F *hMigrationMatrix2 = NULL;
	TH2F *hMigrationMatrix3 = NULL;
	TH2F *hMigrationMatrix4 = NULL;


 if (!(hEffAreaTrue1 = (TH1F*)f->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue1")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile1 << endl;
    exit(0);
  }

	 if (!(hMigrationMatrix1 = (TH2F*)(f->Get("MigMatrix")->Clone("hMigrationMatrix1"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile1 << endl;
    exit(0);
  }

 if (!(hEffAreaTrue2 = (TH1F*)g->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue2")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile2 << endl;
    exit(0);
  }

         if (!(hMigrationMatrix2 = (TH2F*)(g->Get("MigMatrix")->Clone("hMigrationMatrix2"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile2 << endl;
    exit(0);
  }

 if (!(hEffAreaTrue3 = (TH1F*)o->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue3")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile3 << endl;
    exit(0);
  }

         if (!(hMigrationMatrix3 = (TH2F*)(o->Get("MigMatrix")->Clone("hMigrationMatrix3"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile3 << endl;
    exit(0);
  }

 if (!(hEffAreaTrue4 = (TH1F*)p->Get("EffectiveAreaEtrue")->Clone("hEffAreaTrue4")))
  {
    cout << "!!!!!!!ERROR: Did not find histogram EffectiveAreaEtrue in the provided Root performance file!!!!!!!" << sInstrumentResponseFile4 << endl;
    exit(0);
  }

         if (!(hMigrationMatrix4 = (TH2F*)(p->Get("MigMatrix")->Clone("hMigrationMatrix4"))))
  {
    cout << "!!!!!!!ERROR: Did not find migration matrix in the provided root performance file!!!!!!!" << sInstrumentResponseFile4 << endl;
    exit(0);
  }

  cout<<"~~~~Loaded all instrument response function from a CTA file~~~~"<<endl;

//=================================================================================================================================================

	Double_t ELIV      = 1e16;			//LIV energy parameter in TeV.
	Double_t c         = 2.99792458e8;
//	Double_t Distance  = 4.429e24;			//Initialize Distance Parameter for each source.

	Double_t z = 0.5;
        Double_t tH0 = 4.55e17;
        Double_t omegaM = 0.3;
        Double_t omegaG = 0.7;
        Double_t hz = sqrt(omegaG + (omegaM*pow(1+z,3)));
        Double_t D1 = z*(1+z)/hz;
        Double_t D2 = z*(1+z)*(1+z)/hz;

	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;

//=================================================================================================================================================
	int meanIterations = 154;              //Mean number of events, actual value randomly generated from Poissonian distribution.

	TF1 *S = new TF1("S","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) + [3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5]))",0,8000);   /* PG1553();*/
        TF1 *B = new TF1("B","1",0,8000);
	TF1 *H = new TF1("H","1",0,8000);

	double par [6] = {9.335,2685,1088,8.084,6296,1002};
	S->SetParameters(par);

        TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"pow(x,-4.8)",0.3,0.8);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-4.8)",0.3,0.8);
	TF1 *HadronSpectrum     = new TF1("HadronSpectrum"    ,"pow(x,-2.7)",0.3,0.8);

	double SignalRatio      = 0.56;
	double BackgroundRatio  = 0.3;
	double HadronRatio      = 0.14;

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

        TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        bc->SetParameter(0,BackgroundRatio*meanIterations);

	TF1 *hc = new TF1("hc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        hc->SetParameter(0,HadronRatio*meanIterations);

//=================================================================================================================================================

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
        	int      *ind              = new int      [totalIterations];
		Double_t *normareas        = new Double_t [totalIterations];
		Double_t *trueenergies     = new Double_t [totalIterations];
		Double_t *effareas         = new Double_t [totalIterations];
		Double_t *dts              = new Double_t [totalIterations];
		int      *flags            = new int      [totalIterations];

   		 Double_t energy;
   		 Double_t arrivaltime;
		 Double_t deltat;
		 Double_t AreaB;
		 Double_t AreaA;
		 Double_t LogEB;
    		 Double_t LogEA;
    	   	 Double_t NormalisedEffArea;
    		 Double_t ERec;
    		 Double_t EffArea;
    		 TH1D* xproj;
    		 Double_t MaxEffArea = -100000;

		int TrueAreaMaxBin1     = hEffAreaTrue1  ->GetMaximumBin();
            	Double_t MaxEffTrue1    = hEffAreaTrue1  ->GetBinContent(TrueAreaMaxBin1);
                if(MaxEffArea < MaxEffTrue1) MaxEffArea = MaxEffTrue1;

		int TrueAreaMaxBin2     = hEffAreaTrue2  ->GetMaximumBin();
                Double_t MaxEffTrue2    = hEffAreaTrue2  ->GetBinContent(TrueAreaMaxBin2);
                if(MaxEffArea < MaxEffTrue2) MaxEffArea = MaxEffTrue2;

		int TrueAreaMaxBin3     = hEffAreaTrue3  ->GetMaximumBin();
                Double_t MaxEffTrue3    = hEffAreaTrue3  ->GetBinContent(TrueAreaMaxBin3);
                if(MaxEffArea < MaxEffTrue3) MaxEffArea = MaxEffTrue3;

		int TrueAreaMaxBin4     = hEffAreaTrue4  ->GetMaximumBin();
                Double_t MaxEffTrue4    = hEffAreaTrue4  ->GetBinContent(TrueAreaMaxBin4);
                if(MaxEffArea < MaxEffTrue4) MaxEffArea = MaxEffTrue4;

//=================================================================================================================================================
		int k = 0;

		while (k < totalIterations)
			{

			if (k < actualSignalIterations)
				{
				arrivaltime = TMath::Abs(S                  ->GetRandom());
       	        		energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
				}

			else
				{
				if (k < photonIterations)
					{
					arrivaltime = TMath::Abs(B                      ->GetRandom());
              				energy      = TMath::Abs(BackgroundSpectrum     ->GetRandom());
					}
				else
					{
					arrivaltime = TMath::Abs(H                  ->GetRandom());
            			 	energy      = TMath::Abs(HadronSpectrum ->GetRandom());
					}
				}


		if(SupLum)
			{

              		  if(SecondOrder)
				{
         		        deltat = -(tH0*D2)*pow(energy/ELIV,2);
                		}
               		  else
				{
                		deltat = -(tH0*D1)*(energy/ELIV);
                		}
              		}

                else
			{

              		if(SecondOrder)
				{
                		deltat = +(tH0*D2)*pow(energy/ELIV,2);
               			}
                	else
				{
				deltat = +(tH0*D1)*(energy/ELIV);
               			 }
               		}


		if(LIV)
                                {
                                arrivaltime       = (arrivaltime + deltat);
                                }
                        else
                                {
                                deltat = 0;
                                arrivaltime       = (arrivaltime + deltat);
                                }



		Double_t AssignRand = rand3->Rndm();

		Double_t LogE = log10(energy);

////=================================================================================================================================================

	if (arrivaltime <= 1800)
		{
		int      energybin       = hEffAreaTrue1->FindBin(LogE,0,0);
		Double_t energybincenter = hEffAreaTrue1->GetBinCenter(energybin);

		if (LogE > energybincenter)
			{

			int rightbin = energybin + 1;

			AreaA = hEffAreaTrue1->GetBinContent(energybin);
			AreaB = hEffAreaTrue1->GetBinContent(rightbin );
			LogEA = hEffAreaTrue1->GetBinCenter (energybin);
			LogEB = hEffAreaTrue1->GetBinCenter (rightbin );

			}

		else
			{
			int leftbin = energybin -1;

			AreaA = hEffAreaTrue1->GetBinContent(leftbin  );
		  AreaB = hEffAreaTrue1->GetBinContent(energybin);
		  LogEA = hEffAreaTrue1->GetBinCenter (leftbin  );
		  LogEB = hEffAreaTrue1->GetBinCenter (energybin);

			}

		Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
		Double_t intercept = AreaA - (gradient*LogEA);

		EffArea           = (intercept + (gradient*LogE));
		NormalisedEffArea = EffArea/MaxEffArea;

		int trueenergybin = hMigrationMatrix1->GetYaxis()->FindBin(LogE);

                xproj = hMigrationMatrix1->ProjectionX("xproj",trueenergybin,trueenergybin);

		Double_t LogRecEnergy = xproj->GetRandom();
                ERec = pow(10,LogRecEnergy);
		}


	if (arrivaltime > 1800 && arrivaltime <= 3600)
                {
                int      energybin       = hEffAreaTrue2->FindBin(LogE,0,0);
                Double_t energybincenter = hEffAreaTrue2->GetBinCenter(energybin);

                if (LogE > energybincenter)
                        {

                        int rightbin = energybin + 1;

                        AreaA = hEffAreaTrue2->GetBinContent(energybin);
                        AreaB = hEffAreaTrue2->GetBinContent(rightbin );
                        LogEA = hEffAreaTrue2->GetBinCenter (energybin);
                        LogEB = hEffAreaTrue2->GetBinCenter (rightbin );

                        }

                else
                        {
                        int leftbin = energybin -1;

                        AreaA = hEffAreaTrue2->GetBinContent(leftbin  );
                        AreaB = hEffAreaTrue2->GetBinContent(energybin);
                        LogEA = hEffAreaTrue2->GetBinCenter (leftbin  );
                        LogEB = hEffAreaTrue2->GetBinCenter (energybin);

                        }

                Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
                Double_t intercept = AreaA - (gradient*LogEA);

                EffArea           = (intercept + (gradient*LogE));
                NormalisedEffArea = EffArea/MaxEffArea;

                int trueenergybin = hMigrationMatrix2->GetYaxis()->FindBin(LogE);

                xproj = hMigrationMatrix2->ProjectionX("xproj",trueenergybin,trueenergybin);

		Double_t LogRecEnergy = xproj->GetRandom();
                ERec = pow(10,LogRecEnergy);
                }


	 if (arrivaltime > 3600 && arrivaltime <= 5400)
                {
                int      energybin       = hEffAreaTrue3->FindBin(LogE,0,0);
                Double_t energybincenter = hEffAreaTrue3->GetBinCenter(energybin);

                if (LogE > energybincenter)
                        {

                        int rightbin = energybin + 1;

                        AreaA = hEffAreaTrue3->GetBinContent(energybin);
                        AreaB = hEffAreaTrue3->GetBinContent(rightbin );
                        LogEA = hEffAreaTrue3->GetBinCenter (energybin);
                        LogEB = hEffAreaTrue3->GetBinCenter (rightbin );

                        }

                else
                        {
                        int leftbin = energybin -1;

                        AreaA = hEffAreaTrue3->GetBinContent(leftbin  );
                        AreaB = hEffAreaTrue3->GetBinContent(energybin);
                        LogEA = hEffAreaTrue3->GetBinCenter (leftbin  );
                        LogEB = hEffAreaTrue3->GetBinCenter (energybin);

                        }

                Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
                Double_t intercept = AreaA - (gradient*LogEA);

                EffArea           = (intercept + (gradient*LogE));
                NormalisedEffArea = EffArea/MaxEffArea;

                int trueenergybin = hMigrationMatrix3->GetYaxis()->FindBin(LogE);

                xproj = hMigrationMatrix3->ProjectionX("xproj",trueenergybin,trueenergybin);

		Double_t LogRecEnergy = xproj->GetRandom();
                ERec = pow(10,LogRecEnergy);
                }


	if (arrivaltime > 5400 && arrivaltime <= 7200)
                {
                int      energybin       = hEffAreaTrue2->FindBin(LogE,0,0);
                Double_t energybincenter = hEffAreaTrue2->GetBinCenter(energybin);

                if (LogE > energybincenter)
                        {

                        int rightbin = energybin + 1;

                        AreaA = hEffAreaTrue2->GetBinContent(energybin);
                        AreaB = hEffAreaTrue2->GetBinContent(rightbin );
                        LogEA = hEffAreaTrue2->GetBinCenter (energybin);
                        LogEB = hEffAreaTrue2->GetBinCenter (rightbin );

                        }

                else
                        {
                        int leftbin = energybin -1;

                        AreaA = hEffAreaTrue2->GetBinContent(leftbin  );
                        AreaB = hEffAreaTrue2->GetBinContent(energybin);
                        LogEA = hEffAreaTrue2->GetBinCenter (leftbin  );
                        LogEB = hEffAreaTrue2->GetBinCenter (energybin);

                        }

                Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
                Double_t intercept = AreaA - (gradient*LogEA);

                EffArea           = (intercept + (gradient*LogE));
                NormalisedEffArea = EffArea/MaxEffArea;

                int trueenergybin = hMigrationMatrix2->GetYaxis()->FindBin(LogE);

                xproj = hMigrationMatrix2->ProjectionX("xproj",trueenergybin,trueenergybin);

		Double_t LogRecEnergy = xproj->GetRandom();
                ERec = pow(10,LogRecEnergy);
                }


	if (arrivaltime > 7200 && arrivaltime <= 9000)
                {
                int      energybin       = hEffAreaTrue4->FindBin(LogE,0,0);
                Double_t energybincenter = hEffAreaTrue4->GetBinCenter(energybin);

                if (LogE > energybincenter)
                        {

                        int rightbin = energybin + 1;

                        AreaA = hEffAreaTrue4->GetBinContent(energybin);
                        AreaB = hEffAreaTrue4->GetBinContent(rightbin );
                        LogEA = hEffAreaTrue4->GetBinCenter (energybin);
                        LogEB = hEffAreaTrue4->GetBinCenter (rightbin );

                        }

                else
                        {
                        int leftbin = energybin -1;

                        AreaA = hEffAreaTrue4->GetBinContent(leftbin  );
                        AreaB = hEffAreaTrue4->GetBinContent(energybin);
                        LogEA = hEffAreaTrue4->GetBinCenter (leftbin  );
                        LogEB = hEffAreaTrue4->GetBinCenter (energybin);

                        }

                Double_t gradient = (AreaB - AreaA)/(LogEB - LogEA);
                Double_t intercept = AreaA - (gradient*LogEA);

                EffArea           = (intercept + (gradient*LogE));
                NormalisedEffArea = EffArea/MaxEffArea;

                int trueenergybin = hMigrationMatrix4->GetYaxis()->FindBin(LogE);

                xproj = hMigrationMatrix4->ProjectionX("xproj",trueenergybin,trueenergybin);

		Double_t LogRecEnergy = xproj->GetRandom();
                ERec = pow(10,LogRecEnergy);
                }

//////=================================================================================================================================================

		if(AssignRand <= NormalisedEffArea )
			{

			events[k]       = arrivaltime;
			recenergies[k]  = ERec;
			normareas[k]    = NormalisedEffArea;
			effareas[k]     = EffArea;
			dts[k]         = deltat;
			trueenergies[k] = energy;

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

			k++;
			}
		}

//=================================================================================================================================================
//=================================================================================================================================================
//=================================================================================================================================================
	TMath::Sort(totalIterations,events,ind,false);

	file << "#" << "      " << totalIterations << "\n";

	for(int i=0; i < totalIterations; i++)
		{

 		 file << setw(11) << fixed << setprecision(6) << events[ind[i]] << "          " << setw(9) << fixed << setprecision(6) << dts[ind[i]] << "           " << setw(6)<< fixed << setprecision(3) << recenergies[ind[i]] << "           " << setw(6) << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << setw(9) << fixed << setprecision(3) << effareas[ind[i]] << "           " << setw(5)<< fixed << setprecision(3) << normareas[ind[i]] <<  "          " << flags[ind[i]] <<   "\n";
 		}


	printf("Data Set Complete # %i\n",n+1);

	}

	file.close();

//=================================================================================================================================================
	printf("All done.\n");
	}


//=================================================================================================================================================
