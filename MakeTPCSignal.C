
// Save the figure in pdf, eps and png formats:
Bool_t bSaveFigures = kTRUE;




//const char* file = "MergedOutput.root";


void MakeTPCSignal(){

gStyle->SetLegendBorderSize(0);

     TFile* file = new TFile("MergedOutput.root"); //



    TH1F* h1  = (TH1F*) file->Get("Output_highMult")->FindObject("Protons")->FindObject("fHistTPCSignal")->Clone("h1"); //

    //TH1F* h2 = (TH1F*) file->Get("Output_highMult")->FindObject("AProtons")->FindObject("fHistTPCSignal")->Clone("h2"); //

  TH1F* h3  = (TH1F*) file->Get("Output_highMult")->FindObject("Deuterons")->FindObject("fHistTPCSignal")->Clone("h3"); //

    h1->SetTitle();
    h1->SetLineColor(kBlue);
    //h2->SetLineColor(kRed);
    h1->SetStats(kFALSE);
   // h2->SetStats(kFALSE);
    h1->GetXaxis()->SetTitleOffset(0.98844);
    h1->GetYaxis()->SetTitleOffset(1.3); //h->GetYaxis()->SetTitleOffset(0.744);

    TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
    TLegend *legend = new TLegend(0.9,0.6,0.7,0.7);// (x1,y1,x2,y2)x1->> 1st parameter

    legend->SetFillStyle(0);
   // legend->SetHeader("Raw Spectrum");
    legend->AddEntry(h1,"protons");
   // legend->AddEntry(h2,"Anti-protons");
    legend->AddEntry(h3,"deuterons");

     h1->Draw(); //"e1"
     cfig->Update();
    // h2->Draw("same"); //"e1 same"
    h3->Draw();

    TFile* fileout = new TFile("RawSpectrum_antip_p.root","RECREATE");
    h1->Write();
   // h2->Write();
    legend->Draw();
    h3->Write();

  if(bSaveFigures)
  {
  //cfig->SaveAs("ToyMC_SC42_const_2018-10-23.pdf");
  //cfig->SaveAs("ToyMC_SC42_const_2018-10-23.eps");
  cfig->SaveAs("TPCSignal_antip_p.png");
 } // if(bSaveFigures)
    fileout->Close();
}


