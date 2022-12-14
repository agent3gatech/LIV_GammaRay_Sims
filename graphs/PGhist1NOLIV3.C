{
//=========Macro generated from canvas: c1/
//=========  (Wed May 31 11:16:00 2017) by ROOT version5.34/36
   TCanvas *c1 = new TCanvas("c1", "",11,67,700,500);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   
   TH1D *h = new TH1D("h","  Lightcurve - PG1553",100,0,8000);
   h->SetBinContent(2,1);
   h->SetBinContent(3,2);
   h->SetBinContent(4,1);
   h->SetBinContent(5,2);
   h->SetBinContent(6,1);
   h->SetBinContent(7,2);
   h->SetBinContent(8,2);
   h->SetBinContent(9,2);
   h->SetBinContent(10,1);
   h->SetBinContent(11,1);
   h->SetBinContent(12,3);
   h->SetBinContent(13,4);
   h->SetBinContent(14,3);
   h->SetBinContent(15,3);
   h->SetBinContent(16,5);
   h->SetBinContent(17,5);
   h->SetBinContent(18,3);
   h->SetBinContent(19,5);
   h->SetBinContent(20,5);
   h->SetBinContent(21,6);
   h->SetBinContent(22,8);
   h->SetBinContent(23,3);
   h->SetBinContent(24,9);
   h->SetBinContent(25,9);
   h->SetBinContent(26,5);
   h->SetBinContent(27,2);
   h->SetBinContent(28,5);
   h->SetBinContent(29,5);
   h->SetBinContent(30,3);
   h->SetBinContent(31,5);
   h->SetBinContent(32,7);
   h->SetBinContent(33,4);
   h->SetBinContent(34,4);
   h->SetBinContent(35,6);
   h->SetBinContent(36,5);
   h->SetBinContent(37,7);
   h->SetBinContent(38,2);
   h->SetBinContent(39,6);
   h->SetBinContent(40,4);
   h->SetBinContent(41,6);
   h->SetBinContent(42,9);
   h->SetBinContent(43,4);
   h->SetBinContent(44,5);
   h->SetBinContent(45,4);
   h->SetBinContent(46,7);
   h->SetBinContent(47,4);
   h->SetBinContent(48,5);
   h->SetBinContent(49,4);
   h->SetBinContent(50,5);
   h->SetBinContent(52,5);
   h->SetBinContent(53,1);
   h->SetBinContent(54,3);
   h->SetBinContent(55,4);
   h->SetBinContent(56,2);
   h->SetBinContent(57,5);
   h->SetBinContent(58,3);
   h->SetBinContent(59,4);
   h->SetBinContent(60,5);
   h->SetBinContent(61,1);
   h->SetBinContent(62,3);
   h->SetBinContent(63,6);
   h->SetBinContent(65,5);
   h->SetBinContent(66,3);
   h->SetBinContent(67,6);
   h->SetBinContent(68,6);
   h->SetBinContent(69,7);
   h->SetBinContent(70,8);
   h->SetBinContent(71,5);
   h->SetBinContent(72,8);
   h->SetBinContent(73,5);
   h->SetBinContent(74,3);
   h->SetBinContent(75,4);
   h->SetBinContent(76,12);
   h->SetBinContent(77,7);
   h->SetBinContent(78,5);
   h->SetBinContent(79,1);
   h->SetBinContent(80,6);
   h->SetBinContent(81,5);
   h->SetBinContent(82,6);
   h->SetBinContent(83,8);
   h->SetBinContent(84,1);
   h->SetBinContent(85,3);
   h->SetBinContent(86,7);
   h->SetBinContent(87,4);
   h->SetBinContent(88,6);
   h->SetBinContent(89,5);
   h->SetBinContent(90,2);
   h->SetBinContent(91,4);
   h->SetBinContent(92,3);
   h->SetBinContent(93,5);
   h->SetBinContent(94,2);
   h->SetBinContent(95,3);
   h->SetBinContent(96,1);
   h->SetBinContent(97,2);
   h->SetBinContent(98,4);
   h->SetBinContent(99,4);
   h->SetBinContent(100,2);
   h->SetEntries(419);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h->SetLineColor(ci);
   h->GetXaxis()->SetLabelFont(42);
   h->GetXaxis()->SetLabelSize(0.035);
   h->GetXaxis()->SetTitleSize(0.035);
   h->GetXaxis()->SetTitleFont(42);
   h->GetYaxis()->SetLabelFont(42);
   h->GetYaxis()->SetLabelSize(0.035);
   h->GetYaxis()->SetTitleSize(0.035);
   h->GetYaxis()->SetTitleFont(42);
   h->GetZaxis()->SetLabelFont(42);
   h->GetZaxis()->SetLabelSize(0.035);
   h->GetZaxis()->SetTitleSize(0.035);
   h->GetZaxis()->SetTitleFont(42);
   h->Draw("");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
