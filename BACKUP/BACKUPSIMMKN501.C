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
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TString  sInstrumentResponseFile = "VTS.sensitivity.V6.hard.root";

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

void SIMMKN501() {

	int NumberOfDataSets = 1;		       //Number of datasets to include in OF. 

        bool SecondOrder = false;       	       //True for second order dispersion term, False for first order dispersion term.
	bool LIV         = true;                       //Include LIV correction or not.
        bool SupLum      = false;                       //Super Luminar assumption if true. Corresponds to a -/+ in deltat.
        
	char name[16] = "LIVMKN501.txt";

	bool graphs      = false;

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

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
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
	Double_t ELIV      = 1e16;			//LIV energy parameter in TeV.
	Double_t c         = 2.99792458e8;	        //Initialize Distance Parameter for each source.
//	Double_t Distance  = 4.429e24;

	Double_t z = 0.034;
        Double_t tH0 = 4.55e17;
        Double_t omegaM = 0.3;
        Double_t omegaG = 0.7;
        Double_t hz = sqrt(omegaG + (omegaM*pow(1+z,3)));
        Double_t D1 = z*(1+z)/hz;
        Double_t D2 = z*(1+z)*(1+z)/hz;

	TRandom3 *rand3 = new TRandom3();
        gRandom = rand3;
 
//	TF1 *sc = new TF1("sc","Gaus(SignalRatio*meanIterations    ,pow(SignalRatio*meanIterations    ,0.5))",0,10*SignalRatio*meanIterations);
//	TF1 *bc = new TF1("bc","Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5))",0,10*BackgroundRatio*meanIterations);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	int meanIterations = 1491;              //Mean number of events, actual value randomly generated from Poissonian distribution. 


	TF1 *S = new TF1("S","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",0,2750);   /* PG1553();*/

        TF1 *B = new TF1("B","1",0,2750);
	TF1 *H = new TF1("H","1",0,2750);

	double par [3] = {100, 805, 219};

	S->SetParameters(par);

        TF1 *Th1 = new TF1("Th1","0.0633/cosh(4.3116*sqrt(x))",0,20);
        TF1 *Th2 = new TF1("Th2","0.0633/10",0,20);

        TF1 *SignalSpectrum     = new TF1("SignalSpectrum"    ,"pow(x,-2.4)",0.15,12);
        TF1 *BackgroundSpectrum = new TF1("BackgroundSpectrum","pow(x,-2.7)",0.15,12);

        TH1D *h = new TH1D("h", "", 1225, 0, 2750);
        TH1D *a = new TH1D("a", "", 1225, 0.15,12);
        TH1D *m = new TH1D("m", "", 1225, 0.15,12);

	double SignalRatio      = 0.369;
	double BackgroundRatio  = 0.631;
	double HadronRatio      = 0;

//      double SignalCounts     = S->Integral(0,8000);
//      double BackgroundCounts = B->Integral(0,8000);
//      double SignalRatio      = SignalCounts/(SignalCounts+BackgroundCounts);
//	double BackgroundRatio  = BackgroundCounts/(SignalCounts+BackgroundCounts);	

	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        sc->SetParameter(0,SignalRatio*meanIterations);

        TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        bc->SetParameter(0,BackgroundRatio*meanIterations);

	TF1 *hc = new TF1("hc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
        hc->SetParameter(0,HadronRatio*meanIterations);


//	printf("###################################%g\n");
//      printf("###################################%g\n");
//      printf("SignalCounts                     = %g\n",SignalCounts);
//      printf("BackgroundCounts                 = %g\n",BackgroundCounts);
//      printf("SignalRatio                      = %g\n",SignalRatio);
//	printf("BackgroundRatio                  = %g\n",BackgroundRatio);
//      printf("###################################%g\n");

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
//	TF1 *FlatProb = new TF1("FlatProb","1",0,1);
//	TF1 *Acceptance = new TF1("Acceptance","sqrt(x)/2",0,4);

//	TF1 *sc = new TF1("sc","Gaus(SignalRatio*meanIterations    ,pow(SignalRatio*meanIterations    ,0.5))",0,10*SignalRatio*meanIterations);
//      TF1 *bc = new TF1("bc","Gaus(BackgroundRatio*meanIterations,pow(BackgroundRatio*meanIterations,0.5))",0,10*BackgroundRatio*meanIterations);
	
//	Double_t rootmeanIterations = pow(meanIterations,0.5);
//	Double_t rootsignalratiomeanIterations = pow(SignalRatio*meanIterations,0.5);
//	Double_t rootbackratiomeanIterations = pow(BackgroundRatio*meanIterations,0.5);

//	TF1 *sc = new TF1("sc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
//	sc->SetParameter(0,SignalRatio*meanIterations);
//	sc->SetParameter(1,meanIterations);
//	sc->SetParameter(2,2);

//	TF1 *bc = new TF1("bc","exp(-(x-[0])*(x-[0])/(2*[0]))",0,2*meanIterations);
//      bc->SetParameter(0,BackgroundRatio*meanIterations);
//      bc->SetParameter(1,BackgroundRatio*meanIterations);
//      bc->SetParameter(2,rootbackratiomeanIterations);

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	ofstream file;
        file.open(name);

	for (int n = 0; n < NumberOfDataSets; n++)
	{
  
		Double_t nonintsignaliterations     = TMath::Abs(sc  ->GetRandom());
		Double_t nonintbackgrounditerations = TMath::Abs(bc  ->GetRandom());
		Double_t noninthadroniterations     = 0;

		int actualSignalIterations     = (int) nonintsignaliterations;
        	int actualBackgroundIterations = (int) nonintbackgrounditerations;
		int actualHadronIterations     = (int) noninthadroniterations;
		int photonIterations = actualSignalIterations + actualBackgroundIterations;
		int totalIterations  = actualSignalIterations + actualBackgroundIterations + actualHadronIterations;

//      printf("actualSignalIterations           = %g\n",actualSignalIterations);
//	printf("actualBackgroundIterations       = %g\n",actualBackgroundIterations);
//     	printf("###################################%g\n");
//      printf("###################################%g\n");

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
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

		int k = 0;
	
		while (k < totalIterations) 
			{
	
			if (k < actualSignalIterations)
				{
				Double_t arrivaltime = TMath::Abs(S                  ->GetRandom());
       	        		Double_t energy      = TMath::Abs(SignalSpectrum     ->GetRandom());
       	        		Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
//				h->Fill(arrivaltime);
//				z->Fill(energy);
				}

			else
				{		
				if (k < photonIterations)
					{
					Double_t arrivaltime = TMath::Abs(B                      ->GetRandom());
              				Double_t energy      = TMath::Abs(BackgroundSpectrum     ->GetRandom());
             		   		Double_t thetasq     = TMath::Abs(Th1                    ->GetRandom());		
				//	z->Fill(arrivaltime);
				//	z->Fill(energy);
					}
				else
					{
					Double_t arrivaltime = TMath::Abs(H                  ->GetRandom());
            			 	Double_t energy      = TMath::Abs(BackgroundSpectrum ->GetRandom());
 			                Double_t thetasq     = TMath::Abs(Th1                ->GetRandom());
				//	m->Fill(arrivaltime);
				//	z->Fill(energy);
					}


				}



		if(SupLum)
			{

              		  if(SecondOrder)
				{
         		        Double_t deltat = -(tH0*D2)*pow(energy/ELIV,2);
                		}
               		  else
				{
                		Double_t deltat = -(tH0*D1)*(energy/ELIV);
                		}
              		}

                else
			{

              		if(SecondOrder)
				{
                		Double_t deltat = +(tH0*D2)*pow(energy/ELIV,2);
               			}
                	else
				{
				Double_t deltat = +(tH0*D1)*(energy/ELIV);
               			 }
               		}	

	
		Double_t AssignRand = rand3->Rndm();
	
		Double_t LogE = log10(energy);

		int      energybin       = hEffAreaTrue->FindBin(LogE,0,0);
		Double_t energybincenter = hEffAreaTrue->GetBinCenter(energybin);
	
		if (LogE > energybincenter)
			{

			int rightbin = energybin + 1;
//			int leftstart = energybin - 1;

//			int y = rightstart;
//			int x = leftstart;

//			while (hEffAreaRec->GetBinContent(y)==0){
//			++y;
//			}
//			int rightbin = y;
//			Double_t righteffarea = hEffAreaRec->GetBinContent(rightbin);	

//			while (hEffAreaRec->GetBinContent(x)==0){
//			--x;
//			}
//			int leftbin = x;
//			Double_t lefteffarea = hEffAreaRec->GetBinContent(leftbin);

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

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
	
		if(AssignRand < NormalisedEffArea)
			{
//			int numberxbins   = hMigrationMatrix->GetNbinsX();
//			int numberybins   = hMigrationMatrix->GetNbinsY();

			int trueenergybin = hMigrationMatrix->GetYaxis()->FindBin(LogE);
//			TH1D *MigSlice = new TH1D("MigSlice","Title",numberxbins,0,numberxbins);
	
//			for (int i = 0; i < numberxbins; i++)
//				{
//				Double_t bincontent = hMigrationMatrix->GetBinContent(i,trueenergybin);
//				MigSlice->SetBinContent(i,bincontent);
//				}

			hMigrationMatrix->ProjectionX("xproj",trueenergybin,trueenergybin);
//			hMigrationMatrix->ProjectionY("yproj",trueenergybin,trueenergybin);	

//			int recenergybin = MigSlice->GetRandom();  
//			Double_t LogRecEnergy = hMigrationMatrix->GetXaxis()->GetBinCenter(recenergybin);

			Double_t LogRecEnergy = xproj->GetRandom();	
			Double_t ERec = pow(10,LogRecEnergy);

			if(LIV)
				{
				events[k]       = (arrivaltime + deltat);
		//		h->Fill(arrivaltime + deltat);
				}
			else
				{
				events[k]       = (arrivaltime);
		//		h->Fill(arrivaltime);
				}
		
			recenergies[k]  = ERec;
			thetasquares[k] = thetasq;	
			rands[k]        = AssignRand;
			normareas[k]    = NormalisedEffArea;	
			effareas[k]     = EffArea;
			dts[k]         = deltat;
			trueenergies[k] = energy;	
//			h->Fill(arrivaltime);
//     			z->Fill(energy);
//    		        m->Fill(thetasq);
//                      m->Fill(ERec);

			if (k < actualSignalIterations)
                                {
                        	flags[k] = 0;	
			        }

			 else
                                {
                        
			        if (k < photonIterations)
                                        {
                        		flags[k] = 12;
		                        }
                                else
                                        {
					flags[k] = 2;
                                        }
                                }




	
//			delete MigSlice;
//			MigSlice->Delete();




//			printf("energy= %g\n",energy);
//			printf("LogE= %g\n",LogE);
//			printf("energybin= %g\n",energybin);
//			printf("rightbin= %g\n",rightbin);
//			printf("AreaA= %g\n",AreaA);
//			printf("AreaB= %g\n",AreaB);
//			printf("LogEA= %g\n",LogEA);
//			printf("LogEB= %g\n",LogEB);
//			printf("m= %g\n",gradient);
//			printf("c= %g\n",intercept);
//			printf("EffArea= %g\n",EffArea);
//			printf("AssignRand= %g\n",AssignRand);
//			printf("MaxEff= %g\n",MaxEff);
//			printf("NormalisedEffArea= %g\n",NormalisedEffArea);

//			printf("trueenergybin= %g\n",trueenergybin);
//			printf("recenergybin= %g\n",recenergybin);
//			printf("LogRecEnergy= %g\n",LogRecEnergy);
//			printf("LogEnergy= %g\n",LogE);
//			printf("ERec= %g\n",ERec);
//			printf("energy= %g\n",energy);	

			k++;
			}
//		else{
//		printf("energy= %g\n",energy);
//		printf("LogE= %g\n",LogE);
//		printf("AssignRand= %g\n",AssignRand);
//		printf("NormalisedEffArea= %g\n",NormalisedEffArea);
//		}
		}

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	TMath::Sort(totalIterations,events,ind,false);

	file << "#" << "      " << totalIterations << "\n";

	for(int i=0; i < totalIterations; i++) 
		{

 		 file << fixed << setprecision(6) << events[ind[i]] << "          " << fixed << setprecision(6) << dts[ind[i]] << "          " << fixed << setprecision(3) << recenergies[ind[i]] << "           "  << fixed << setprecision(3) << trueenergies[ind[i]] << "           " << fixed << setprecision(3) << effareas[ind[i]] << "           " << fixed << setprecision(3) << normareas[ind[i]] <<  "          " << flags[ind[i]] <<   "\n";
 		}
	

	printf("Data Set Complete # %g\n",n+1);	 
	
	}	

	file.close();

//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\

	gStyle->SetOptStat(0);
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetTitleAlign(13);

        gStyle->SetCanvasColor(10);

        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameFillColor(0);

        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(0);
        gStyle->SetPadTopMargin(0.07);
        gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.11);
        gStyle->SetPadBottomMargin(0.1);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        gStyle->SetHistFillStyle(0);
        gStyle->SetOptTitle(0);

        gStyle->SetTitleSize(0.22);
        gStyle->SetTitleFontSize(2);
        gStyle->SetTitleFont(42);
        gStyle->SetTitleFont(62,"xyz");
        gStyle->SetTitleYOffset(1.0);
        gStyle->SetTitleXOffset(1.0);
        gStyle->SetTitleXSize(0.07);
        gStyle->SetTitleYSize(0.07);
        gStyle->SetTitleX(.15);
        gStyle->SetTitleY(.98);
        gStyle->SetTitleW(.70);
        gStyle->SetTitleH(.05);

        gStyle->SetStatFont(42);
        gStyle->SetStatX(.91);
        gStyle->SetStatY(.90);
        gStyle->SetStatW(.15);
        gStyle->SetStatH(.15);

        gStyle->SetLabelFont(42,"xyz");
        gStyle->SetLabelSize(0.035,"xyz");
        gStyle->SetGridColor(16);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetOptStat(111111);
        gStyle->SetOptFit(111111);

//	TCanvas * c1 = new TCanvas("c1");
//        TCanvas * c2 = new TCanvas("c2");
//      TCanvas * c3 = new TCanvas("c3");
      
//	TF1 *f1 = new TF1("f1","[1]*pow(x,-[0])",0.15,12);

//	c1->cd();
        h->GetXaxis()->SetTitle("Event Time (s)");
        h->GetYaxis()->SetTitle("Counts");
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.06);
//        h->Draw();




//        c2->cd();
//        f1->SetParNames("Spectral Index","Arb. Intercept");
//        a->Fit("f1","R");
        a->GetXaxis()->SetTitle("Energy (TeV)");
        a->GetYaxis()->SetTitle("Counts");
        a->GetXaxis()->CenterTitle();
        a->GetYaxis()->CenterTitle();
        a->GetXaxis()->SetTitleSize(0.05);
        a->GetYaxis()->SetTitleSize(0.06);
        a->GetXaxis()->SetLabelSize(0.05);
        a->GetYaxis()->SetLabelSize(0.06);
//        a->Draw();

//	c1->SaveAs("MKNIC1.C");
//        c2->SaveAs("MKNIC2.C");


	printf("Distrust even Einstein %g\n");
	}


//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\
//=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\\






























































