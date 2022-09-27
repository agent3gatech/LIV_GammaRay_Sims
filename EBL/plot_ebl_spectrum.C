
// MeV
double average_ebl_polynomial(double _e)
{

	// _x en µm
	double _x = (6.582118e-22 * 299792458. * 1e6)/_e;
	double res = 1.2 +
	  0.533686 * pow(log10(_x),1) +
	  -1.81116 * pow(log10(_x),2) +
	  -1.28135 * pow(log10(_x),3) +
	   1.94804 * pow(log10(_x),4) +
	  0.655142 * pow(log10(_x),5) +
	  -1.19371 * pow(log10(_x),6) +
	  0.405981 * pow(log10(_x),7) +
	-0.0437552 * pow(log10(_x),8);
	return pow(10,res);

}

// MeV
// _mu = cosinus de l'angle entre les deux photons
double sigma_gg(double _eg, double _eb, double _mu)
{

	double eth  = 2 * pow(0.511,2)/(_eg * (1 - _mu));
	cout << "eth  = " << eth << endl;
	double beta = sqrt(1 - eth/_eb);
	cout << "beta = " << beta << endl;
	double sigt = 6.65e-25; // en cm^2
	double sig  = (3*sigt/16.) * (1 - pow(beta,2)) *
						(2*beta*(pow(beta,2)-2)+(3-pow(beta,4))*log((1+beta)/(1-beta)));
	if(TMath::IsNaN(sig)) return 9999.;
	else return sig;

}


// Eb is set to 1e-7 MeV
double sigma_gg_10m7(double _eg, double _mu)
{

	return sigma_gg(_eg, 1e-7, _mu);

}



void plot_ebl_spectrum()
{

	TCanvas *c = new TCanvas("plot_ebl_spectrum_mev", "plot_ebl_spectrum_mev", 20, 20, 600, 600);
	c->SetLogy();
	c->SetLogx();
	c->SetTickx();
	c->SetTicky();	
	c->SetRightMargin(0.06);
	c->SetLeftMargin(0.15);
	c->SetTopMargin(0.06);
	
	TH1F *hframe = c->DrawFrame(1.97e-10, 1, 1.97e-6, 100);
	hframe->SetXTitle("#epsilon (MeV)");
	hframe->GetXaxis()->CenterTitle();
	hframe->SetYTitle("#nuI_{#nu} (nW m^{-2} sr^{-1})");
	hframe->GetYaxis()->SetTitleOffset(1.5);
	hframe->GetYaxis()->CenterTitle();

	TF1 *f_average_polynomial = new TF1("f_average_ebl_polynomial", "average_ebl_polynomial(x)", 1.97e-10, 1.97e-6);
	f_average_polynomial->SetNpx(1000);
	f_average_polynomial->Draw("SAME");

	TCanvas *c = new TCanvas("plot_sigma", "plot_sigma", 620, 20, 600, 600);
	c->SetLogy();
	c->SetLogx();
	c->SetTickx();
	c->SetTicky();	
	c->SetRightMargin(0.06);
	c->SetLeftMargin(0.15);
	c->SetTopMargin(0.06);
	
/*	TH1F *hframe = c->DrawFrame(100, 0, 1e5, 1e-26);
	hframe->SetXTitle("E_{#gamma} (GeV)");
	hframe->GetXaxis()->CenterTitle();
	hframe->SetYTitle("#sigma_{#gamma#gamma} (cm^{2})");
	hframe->GetYaxis()->SetTitleOffset(1.5);
	hframe->GetYaxis()->CenterTitle();*/

/*	TF2 *f_sigma = new TF2("f_sigma", "sigma_gg(x,y,-0.5)", 100, 1000000, 1.97e-10, 1.97e-6);
	f_sigma->SetNpx(100);
	f_sigma->SetNpy(100);
	f_sigma->Draw();*/

	TF1 *f_average_polynomial2 = new TF1("f_sigma2", "sigma_gg_10m7(x,-0.5)", 100, 100000);
	f_average_polynomial2->SetLineColor(kRed);
	f_average_polynomial2->SetNpx(1000);
	f_average_polynomial2->Draw();



}