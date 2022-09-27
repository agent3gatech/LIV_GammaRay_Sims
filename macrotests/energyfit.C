

	void energyfit(){


	TF1 *SignalSpectrum     = new TF1("SignalSpectrum","pow(x,-3.46)",0.3,4);

	TH1D *SpectrumHistogram = new TH1D("SpectrumHistogram", "Energy Spectrum", 2000, 0.3,4);

	
	int k = 0;
	while (k < 10000000)
	{

	Double_t energy      = TMath::Abs(SignalSpectrum ->GetRandom());

	SpectrumHistogram->Fill(energy);

	k++;
	}

	TF1 *f1 = new TF1("f1","[0]*pow(x,[1])",0.3,4);
        SpectrumHistogram->Fit("f1","R");
        SpectrumHistogram->Draw();

	gPad->SetLogx(); //Set LogLog axis for Energy Histogram
        gPad->SetLogy(); //Set LogLog axis for Energy Histogram









	















































}
