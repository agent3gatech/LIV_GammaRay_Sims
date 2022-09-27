
/*
par[0] = index 1
par[1] = normalisation
par[2] = index 2
par[3] = z cut
*/
double BrokenPowerLaw(double *x, double *par)
{

	double z = x[0];
	if(z <= par[3]) return par[1] * TMath::Power(1 + z, par[0])/TMath::Power(1 + par[3], par[0]);
	else return par[1] * TMath::Power(1 + par[3], - par[2]) * TMath::Power(1 + z, par[2]);

}


double SFR_z(double *x, double *par)
{

	return BrokenPowerLaw(x, par) + BrokenPowerLaw(x, par + 4);

}




plot_rho_z()
{

	TCanvas *c = new TCanvas("plot_rho_z", "plot_rho_z", 600, 600);
	c->SetLogy();
	c->SetTickx();
	c->SetTicky();	
	c->SetRightMargin(0.06);
	c->SetLeftMargin(0.15);
	c->SetTopMargin(0.06);
	
	TH1F *hframe = c->DrawFrame(0, 1e-3, 4, 1);
	hframe->SetXTitle("z");
	hframe->GetXaxis()->CenterTitle();
	hframe->SetYTitle("SFR (M_{#odot} yr^{-1} Mpc^{-3})");
	hframe->GetYaxis()->SetTitleOffset(1.5);
	hframe->GetYaxis()->CenterTitle();

	TF1 *f_SFR_z = new TF1("f_SFR_z", SFR_z, 0, 4, 8);
	f_SFR_z->SetParameters(3.5, 0.1, -1.2, 1.2,
									4.5, 0.1, 0, 1.0);
	f_SFR_z->SetNpx(1000);
	
	TF1 *f_SFR_OPT_z = new TF1("f_SFR_OPT_z", BrokenPowerLaw, 0, 4, 4);
	f_SFR_OPT_z->SetParameters(3.5, 0.1, -1.2, 1.2);
	f_SFR_OPT_z->SetNpx(1000);
	f_SFR_OPT_z->SetLineStyle(2);

	TF1 *f_SFR_LIG_z = new TF1("f_SFR_LIG_z", BrokenPowerLaw, 0, 4, 4);
	f_SFR_LIG_z->SetParameters(4.5, 0.1, 0, 1.0);
	f_SFR_LIG_z->SetNpx(1000);
	f_SFR_LIG_z->SetLineStyle(5);
	
	f_SFR_z->Draw("SAME");
	f_SFR_OPT_z->Draw("SAME");
	f_SFR_LIG_z->Draw("SAME");


	TLegend *leg = new TLegend(0.6, 0.15, 0.9, 0.3);
	leg->SetBorderSize(1);
	leg->SetFillColor(kWhite);
	leg->AddEntry(f_SFR_z, "SFR Total", "l");
	leg->AddEntry(f_SFR_OPT_z, "SFR OPT", "l");
	leg->AddEntry(f_SFR_LIG_z, "SFR LIG", "l");
	leg->Draw();


}