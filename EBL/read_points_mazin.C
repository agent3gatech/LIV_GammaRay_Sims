
// Axe x x1=89.083 x2=539.334
//          0.01      10
float convert_x(float _x)
{

	float res = (log10(10) - log10(0.01))/(539.334 - 89.083)* _x - 5.93555595766175781e-01 - 2;
	return pow(10,res);

}

// Axe y y1=601.14 y2=217.64
//          0.01      2
float convert_y(float _y)
{

	float res = (log10(2) - log10(0.01))/(217.64 - 601.14)* _y + 1.60736104105584876e+00;
	return pow(10, res);

}

double pol_log(double *x, double *par)
{

	return pow(10, 	par[0] + 
							cosh(x[0]* par[1] + par[2])
							);

}







/*
TGraph *interpolate(float z, TGraph *gr1, float z1, TGraph *gr2, flaot z2)
{

	const int Min_N = TMath::Min(gr1->GetN(), gr2->GetN());
	float X[Min_N], Y[Min_N];
	for(int k = 0; k < Min_N; k++) {
		double x1, y1; x2, y2;
		gr1->GetPoint(k, x1, y1);
		gr2->GetPoint(k, x2, y2);
		X[k] = x1;
		Y[k] = 
		
	
	}



}
*/

double find_E_from_expmtau(TGraph *gr, double expmtau, double xmin, double xmax)
{

	double step = 1e-4;
	double x_res, lim = 1e9;
	for(double x = xmin; x < xmax; x += step) {
	
		if(fabs(gr->Eval(x, 0, "S") - expmtau) <= lim) {
			lim = fabs(gr->Eval(x, 0, "S") - expmtau);
			x_res = x;
//			printf("%12.7lf %15.7lf %15.7lf %15.7lf\n", x, gr->Eval(x, 0, "S"), expmtau, fabs(gr->Eval(x, 0, "S") - expmtau));
		}
	
	}

//	cout << endl;
	return x_res;

}


TGraph *get_graph_for_fixed_expmtau(TGraph **gr_tab, double *z_tab, int n_graph, double expmtau, int style)
{

	double *X = new double[n_graph];
	double *Y = new double[n_graph];
	for(int k = 0; k < n_graph; k++) {
	
		X[k] = z_tab[k];
		Y[k] = find_E_from_expmtau(gr_tab[k], expmtau, 0.01, 10);
//		printf("%d %12.4lf %12.4lf\n", k, X[k], Y[k]);
	
	}

//	cout << endl;

	TGraph *gr = new TGraph(n_graph, X, Y);
	gr->SetMarkerStyle(style);
	return gr;

}






void read_points_mazin()
{

	float xr, yr;
	float x[50], y[50];
	TString line;
	
	TFile *fr = new TFile("z_expmtau.root", "recreate");
	
	TCanvas *c = new TCanvas("plot_exp_tau", "plot_exp_tau", 20, 20, 600, 600);
	c->SetLogy();
	c->SetLogx();
//	c->SetGridy();
//	c->SetGridx();
	c->SetTickx();
	c->SetTicky();	
	c->SetRightMargin(0.06);
	c->SetLeftMargin(0.15);
	c->SetTopMargin(0.06);
	c->cd();

	TH1F *hframe = c->DrawFrame(0.01, 0.01, 10, 2);
	hframe->SetXTitle("E (TeV)");
	hframe->GetXaxis()->CenterTitle();
	hframe->SetYTitle("exp(-#tau)");
	hframe->GetYaxis()->SetTitleOffset(1.5);
	hframe->GetYaxis()->CenterTitle();

	int N = 0, I = 0;
	TGraph *gr[5];

	FILE *fp = fopen("points_fig2_mazin.txt", "r");
	while(!feof(fp)) {
	
		line.Gets(fp);
		if(line.Length() == 0) {
			//cout << endl;
			gr[I] = new TGraph(N, x, y);
			N = 0;
			I++;
		}
		int res = sscanf(line.Data(), "%f,%f", &xr, &yr); 
		if (res == 2) {
			//cout << line << "    " << convert_x(xr) << "    " << convert_y(yr) << "    " << N << endl;
			x[N] = convert_x(xr); y[N] = convert_y(yr);
			N++;
		}

	}

	fclose(fp);
	
	double Z[5] = {0.03, 0.06, 0.1, 0.2, 0.5};

	TLegend *leg = new TLegend(0.2, 0.2, 0.4, 0.4);
	leg->SetBorderSize(0);
	leg->SetFillColor(kWhite);
	for(int g = 0; g < 5; g++) {
	
		gr[g]->SetLineStyle(g+1);
		gr[g]->SetLineWidth(2);
		TString leg_t, gr_name0;
		gr_name0.Form("gr_epmtau_vs_e_z_%.2f",Z[g]);
		gr_name0.ReplaceAll(".", "_");
		leg_t.Form("z = %.2lf", Z[g]);
		leg->AddEntry(gr[g], leg_t, "l");
		gr[g]->SetName(gr_name0.Data());
		gr[g]->Draw("LSAME"); 
		gr[g]->Write();
	
	}
	leg->Draw();

	TLine *limit = new TLine(0.01, exp(-1), 10, exp(-1));
	limit->SetLineStyle(7);
	limit->Draw();
	
	TLatex *txt_limit = new TLatex(1.5e-2, 0.4,"#tau = 1");
	txt_limit->Draw();
	
	c->Update();
	

	TCanvas *c2 = new TCanvas("plot_e_z", "plot_e_z", 620, 20, 600, 600);
	c2->SetLogy();
	c2->SetTickx();
	c2->SetTicky();	
	c2->SetRightMargin(0.06);
	c2->SetLeftMargin(0.15);
	c2->SetTopMargin(0.06);
	c2->cd();

	TH1F *hframe2 = c2->DrawFrame(0, 0.01, 0.52, 10);
	hframe2->SetXTitle("z");
	hframe2->GetXaxis()->CenterTitle();
	hframe2->SetYTitle("E (TeV)");
	hframe2->GetYaxis()->SetTitleOffset(1.5);
	hframe2->GetYaxis()->CenterTitle();

//	TF1 *fit_func = new TF1("fit_func", pol_log, 0, 0.6, 3);
	
	TGraph *t_gr_e_vs_epmtau;
	double expmtau_min = 0.3, expmtau_max = 1.0, expmtau_step = 0.01;
	
	for(double expmtau = expmtau_min; expmtau <= expmtau_max; expmtau += expmtau_step) {
	
		cout << "> " << expmtau << endl;
		TString gr_name; gr_name.Form("gr_e_vs_z_epmtau_%.3f", expmtau);
		gr_name.ReplaceAll(".", "_");
		t_gr_e_vs_epmtau = get_graph_for_fixed_expmtau(gr, Z, 5, expmtau, 1);
		t_gr_e_vs_epmtau->SetName(gr_name.Data());
		t_gr_e_vs_epmtau->Write();
		t_gr_e_vs_epmtau->Draw("PLSAME");
		c2->Update();
	
	}
	
		cout << "> " << 0.995 << endl;
		TString gr_name; gr_name.Form("gr_e_vs_z_epmtau_%.3f", 0.995);
		gr_name.ReplaceAll(".", "_");
		t_gr_e_vs_epmtau = get_graph_for_fixed_expmtau(gr, Z, 5, 0.995, 1);
		t_gr_e_vs_epmtau->SetName(gr_name.Data());
		t_gr_e_vs_epmtau->Write();
		t_gr_e_vs_epmtau->Draw("PLSAME");
		c2->Update();

		cout << "> " << 0.999 << endl;
		TString gr_name; gr_name.Form("gr_e_vs_z_epmtau_%.3f", 0.999);
		gr_name.ReplaceAll(".", "_");
		t_gr_e_vs_epmtau = get_graph_for_fixed_expmtau(gr, Z, 5, 0.999, 1);
		t_gr_e_vs_epmtau->SetName(gr_name.Data());
		t_gr_e_vs_epmtau->Write();
		t_gr_e_vs_epmtau->Draw("PLSAME");
		c2->Update();
	
	
/*
	TGraph *gr03 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.3, 29);
	gr03->Draw("PLSAME"); c2->Update();
	TGraph *gr04 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.4, 2);
	gr04->Draw("PLSAME"); c2->Update();
	TGraph *gr05 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.5, 3);
	gr05->Draw("PLSAME"); c2->Update();
	TGraph *gr06 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.6, 4);
	gr06->Draw("PLSAME"); c2->Update();
	TGraph *gr07 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.7, 5);
	gr07->Draw("PLSAME"); c2->Update();
	TGraph *gr08 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.8, 20);
	gr08->Draw("PLSAME"); c2->Update();
	TGraph *gr09 = get_graph_for_fixed_expmtau(gr, Z, 5, 0.9, 21);
	gr09->Draw("PLSAME"); c2->Update();

	TLegend *leg2 = new TLegend(0.7, 0.6, 0.9, 0.9);
	leg2->SetBorderSize(0);
	leg2->SetFillColor(kWhite);
	leg2->SetHeader("exp(-#tau) = ");
	leg2->AddEntry(gr03, "0.3", "p");
	leg2->AddEntry(gr04, "0.4", "p");
	leg2->AddEntry(gr05, "0.5", "p");
	leg2->AddEntry(gr06, "0.6", "p");
	leg2->AddEntry(gr07, "0.7", "p");
	leg2->AddEntry(gr08, "0.8", "p");
	leg2->AddEntry(gr09, "0.9", "p");
	leg2->Draw();
*/
//	fr->Write();
	fr->Close();


}