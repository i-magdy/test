// fit range
float xmin = 0.2;
float xmax = 2.0;
// min and max for signal region
double xmin_peak = 0.6;
double xmax_peak = 1.6;
// DCA cut
float DCAlow  = -0.1;
float DCAhigh =  0.1;
// momentum bins to extract raw signal
int binlow   = 4;  // lowest momentum bin for analysis
int binmid   = 8;  // raw signal directly from histP up to this bin
int binmid2  = 26; // raw signal from TOFm2 fit (fitting background only) up to this bin
int binhigh  = 40; // raw signal from TOFm2 fit up to this bin
float pstep = 0.1; // momentum bin size, same for protons

// input settings
const char* file = "High_Mult_TOF_PID_LHC16_k_l_o_p.root";
const char* list = "Protons";
const char* part = "prot";
const char* outp = "Output_highMult_TOFpid";
//const char* outp_tof = "Output_TOFPid";
const char* hist = "fHistTOFmass2_DCAxy_p";
const char* outd = "TOFm2fits";
const char* fout = "RawSpectrum";

#include "TPaveStats.h"

//declare functions
void FillSpectrumLowP(TH3F* hist, TH1F* rawSpectrum);
void FitTOFm2Background(TH3F* hist, TH1F* rawSpectrum, TH1F* purity);
void FitTOFm2Norm(TH3F* hist, TH1F* rawSpectrum);
void WriteHistos(TH1F* spectrum, TH1F* purity);
Double_t TOFbackgroundExcludePeak(Double_t* x, Double_t *par);

void MakeRawSignalProt(){

    // hist raw spectrum
    TH1F* hRawSpectrum = new TH1F("hRawSpectrum", Form("Raw spectrum %s",list), 45, 0, 4.5);
    hRawSpectrum->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    hRawSpectrum->GetYaxis()->SetTitle("counts");

    // hist raw spectrum with TOF3sigma cut
    TH1F* hRawSpectrum_tofcut = new TH1F("hRawSpectrum_tofcut", Form("Raw spectrum with TOF cut %s",list), 45, 0, 4.5);
    hRawSpectrum_tofcut->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    hRawSpectrum_tofcut->GetYaxis()->SetTitle("counts");

    TH1F* purity = new TH1F("hPurity",Form("Purity %s",list),40,0,4.0);
    purity->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    purity->GetYaxis()->SetTitle("purity");

    TFile* file = new TFile(file);
    // 3D hist
    TH3F* hTOFm2_DCAxy_p = (TH3F*)file->Get(outp)->FindObject(list)->FindObject(hist)->Clone("hTOFm2_DCAxy_p");
    hTOFm2_DCAxy_p->Sumw2();

    TH1F* hist_tofcut_p = (TH1F*)file->Get(outp)->FindObject(list)->FindObject("fHistP")->Clone("hist_tofcut_p");

    gSystem->Exec(Form("mkdir -p %s_%s", outd, part));
    gStyle->SetOptTitle(0);

    // fill low mometum part using counts from 3d hist (with cut on DCAxy < 0.1 cm!)
    FillSpectrumLowP(hTOFm2_DCAxy_p,hRawSpectrum);

    // fill intermediate part (0.7 - 2.5 GeV) using fit to TOFm2 by fitting background only
    FitTOFm2Background(hTOFm2_DCAxy_p, hRawSpectrum, purity);

    // fill high momentum part using fit to TOFm2 with normalised signal function
    FitTOFm2Norm(hTOFm2_DCAxy_p, hRawSpectrum);

    // correct raw spectrum for PID purity
    for (int i=8;i<=24;i++){
        hRawSpectrum->SetBinContent(i, (hist_tofcut_p->GetBinContent(i))*(purity->GetBinContent(i)));
        hRawSpectrum->SetBinError(i, hist_tofcut_p->GetBinError(i));
    }

    WriteHistos(hRawSpectrum, purity);
}

void FillSpectrumLowP(TH3F* hist, TH1F* rawSpectrum){

    TH2F* hDCAxy_p = (TH2F*)hist->Project3D("zx")->Clone("hDCAxy_p");
    // cut on DCAxy
    hDCAxy_p->GetYaxis()->SetRangeUser(DCAlow,DCAhigh);
    TH1F* h_p = (TH1F*)hDCAxy_p->ProjectionX("oe")->Clone("h_p");

    for (Int_t i = binlow; i< binmid; i++){
        rawSpectrum->SetBinContent(i, h_p->GetBinContent(i));
        rawSpectrum->SetBinError(i, h_p->GetBinError(i));
    }

}

void FitTOFm2Background(TH3F* hist, TH1F* rawSpectrum, TH1F* purity){

    // cut on DCA
    hist->GetZaxis()->SetRangeUser(DCAlow,DCAhigh);
    // 2D hist
    TH2F* hTOFm2_p = (TH2F*)hist->Project3D("yx")->Clone("hTOFm2_p");
    // 1D hist
    count(hTOFm2_p);
    TH1F* hTOFm2;
    // fit functions
    TF1 fTOFbackground("fTOFbackground", "[0] + [1]*x + [2]*TMath::Exp([3]*x) + (x <= ([5] + [6]*[7])) * [4] * TMath::Gaus(x, [5], [6]) + (x > ([5]+[6]*[7])) * [4] * TMath::Exp(-(x-[5]-[7]*[6]*0.5)*[7]/[6])", xmin, xmax);

    // [0] + [1]*x + [2]*TMath::Exp([3]*x)
    // TMath::Exp([0]+[1]*x+[2]*x*x)

    // TOFm2 outside peak region
    TF1 *fitFcnBckgrndExc = new TF1("fitFcnBckgrndExc",TOFbackgroundExcludePeak,xmin,xmax,8);
    // start values for bin 15
    Double_t StartValuesBackground[8] = {1166, -424, 79148, -7, 343566, 0.24, 0.02, 1.0};

    // fixed values for protons
    Double_t FixedValuesBackground8Prot[8] = {1.33952e+02, -4.72519e+01, 4.80282e+05, -7.40715e+00, -3.96490e+05, 2.58284e-02, 3.34786e-03, 2.48370e-02};

    // fixed values for anti-protons
    Double_t FixedValuesBackground8[8] = {2.13561e+02, -4.36305e+01, 6.74672e+05, -1.61308e+01, -3.57301e+05, 4.00920e-02, 2.92436e-03, 4.74999e-02};
    Double_t FixedValuesBackground9[8] = {307, -81.8, 655590, -24, -310654, 0.037, 0.00244, 0.0616};
    Double_t FixedValuesBackground10[8] = {315, -94.2, 2276870, -30.6, -514170, 0.0747, 0.0019, 0.071};

    // yield and yield error
    Double_t yield = 1.0;
    Double_t yieldErr = 0.0;

    // purity
    Double_t pur = 1.0;

    // help histogram to subtract background from data
    Int_t nBins = 200;
    TH1F* fHistHelp = new TH1F("fHistHelp","fHistHelp",nBins,0,xmax);

    // canvas
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->SetLogy();

    // ---------------------------------------------------------------

    // fit in momentum bins starting from bin 1.4-1.5 to 2.5
    for (Int_t ibin = 15; ibin <= binmid2-1; ibin++){

        hTOFm2 = (TH1F*)hTOFm2_p->ProjectionY("_py",ibin,ibin,"oe")->Clone(Form("hTOFm2_ptbin_%i",ibin));

        // set (start) values
        for (Int_t i = 0; i <= 7; i++){
            fitFcnBckgrndExc->SetParameter(i,StartValuesBackground[i]);
        }

        hTOFm2->Fit("fitFcnBckgrndExc","","",xmin,xmax);

        // set values to background in full range
        for (Int_t i = 0; i <= 7; i++){
            fTOFbackground.SetParameter(i,fitFcnBckgrndExc->GetParameter(i));
        }

        //remove fitted function from data
        hTOFm2->GetFunction("fitFcnBckgrndExc"); //->Delete();

        // draw options
        hTOFm2->SetMarkerStyle(kOpenCircle);
        hTOFm2->SetMarkerColor(kBlue);
        hTOFm2->SetLineColor(kBlue);
        hTOFm2->GetXaxis()->SetRangeUser(xmin,xmax);
        hTOFm2->Draw("e1");

        fTOFbackground.SetLineColor(kRed);
        fTOFbackground.DrawCopy("same");

        for(Int_t ival =0; ival <= 7; ival++) {
            StartValuesBackground[ival] = fitFcnBckgrndExc->GetParameter(ival);
            // cout << StartValues[ival] << endl;
        }

        // extract raw signal as data histogram minus background
        Double_t xValue      = 0.595;
        Double_t xStep       = 0.01;
        Double_t yBackground = 0.0;
        Double_t yData       = 0.0;

        for (int i=61; i<=160; i++){
            xValue += xStep;
            yBackground = fTOFbackground.Eval(xValue);
            yData       = hTOFm2->GetBinContent(i);
            fHistHelp->SetBinContent(i, yData - yBackground);
            fHistHelp->SetBinError(i,hTOFm2->GetBinError(i));
        }

        fHistHelp->SetLineColor(kGreen);
        fHistHelp->SetMarkerColor(kGreen);
        fHistHelp->SetLineWidth(2);
        fHistHelp->DrawCopy("e1,same");

        TLegend * leg = new TLegend(0.1, 0.85, 0.4, 0.99);
        leg->AddEntry(hTOFm2, Form("Data %s",part == "prot"? "proton" : "anti-proton"),   "P");
        leg->AddEntry(fHistHelp, "Signal", "LPE");
        leg->AddEntry("fTOFbackground", "Background", "L");
        leg->SetFillColor(0);
        leg->SetTextSize(gStyle->GetTextSize()*0.8);
        leg->Draw();

        TLatex *momentumbin = new TLatex(0.7, 0.6, TString::Format("%g < #it{p} < %g GeV/#it{c}",(ibin-1)*pstep,ibin*pstep));
        momentumbin->SetNDC();
        momentumbin->SetTextFont(42);
        momentumbin->SetTextSize(0.035);
        momentumbin->Draw();

        TLatex * text0 = new TLatex (1.85,  1.2, "ALICE Simulation");
        text0->Draw();

        // yield and error
        yield = fHistHelp->IntegralAndError(61,160,yieldErr);
        pur = (fHistHelp->Integral(81,120))/(hTOFm2->Integral(81,120));

        fHistHelp->Reset();

        rawSpectrum->SetBinContent(ibin, yield);
        rawSpectrum->SetBinError(ibin, yieldErr);

        purity->SetBinContent(ibin,pur);

        // Create and plot a box with stats
        c1->Modified();
        c1->Update();

        // Retrieve the stat box
        TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
        ps->SetName("mystats");

        ps->SetX1NDC(0.65);
        ps->SetX2NDC(0.99);
        ps->SetY1NDC(0.7);
        ps->SetY2NDC(0.9);

        TList *list = ps->GetListOfLines();

        // Remove Entries, mean and RMS lines
        TText *tEntries = ps->GetLineWith("Entries");
        TText *tMean    = ps->GetLineWith("Mean");
        TText *tRMS     = ps->GetLineWith("RMS");
        list->Remove(tEntries);
        list->Remove(tMean);
        list->Remove(tRMS);

        int textFont = 42;
        float textSize = 0.03;
        int npointsBckgrnd = 80;

        // Add lines in the stat box
        TLatex *myt = new TLatex(0, 0, TString::Format("#chi^{2} / ndf (bckgnd only) = %g / %i", fitFcnBckgrndExc.GetChisquare(), npointsBckgrnd ));
        myt->SetTextFont(textFont);
        myt->SetTextSize(textSize);

        TLatex *myt0 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", "yield p", yield, yieldErr));
        myt0->SetTextFont(textFont);
        myt0->SetTextSize(textSize);

        list->Add(myt);
        list->Add(myt0);

        // the following line is needed to avoid that the automatic redrawing of stats
        hTOFm2->SetStats(0);

        c1->Modified();
        c1->Update();

        c1->Print(Form("%s_%s/TOFm2fit_%s_bin_%i.pdf",outd, part, part,ibin));

    } // loop over momentum bins

    // ---------------------------------------------------------------

    // fit in momentum bins starting from bin 1.4-1.5 to 0.7
    for (Int_t ibin = 14; ibin >= binmid; ibin--){

        hTOFm2 = (TH1F*)hTOFm2_p->ProjectionY("_py",ibin,ibin,"oe")->Clone(Form("hTOFm2_ptbin_%i",ibin));

        // set (start) values
        for (Int_t i = 0; i <= 7; i++){
            fitFcnBckgrndExc->SetParameter(i,StartValuesBackground[i]);
        }

        if (ibin == 8 && part == "prot"){
            for (Int_t i = 0; i <= 7; i++){
                fitFcnBckgrndExc->FixParameter(i,FixedValuesBackground8Prot[i]);
            }
        }

        if (ibin == 8 && part == "aprot"){
            for (Int_t i = 0; i <= 7; i++){
                fitFcnBckgrndExc->FixParameter(i,FixedValuesBackground8[i]);
            }
        }

        if (ibin == 9 && part == "aprot"){
            for (Int_t i = 0; i <= 7; i++){
                fitFcnBckgrndExc->FixParameter(i,FixedValuesBackground9[i]);
            }
        }

        if (ibin == 10 && part == "aprot"){
            for (Int_t i = 0; i <= 7; i++){
                fitFcnBckgrndExc->FixParameter(i,FixedValuesBackground10[i]);
            }
        }

        hTOFm2->Fit("fitFcnBckgrndExc","","",xmin,xmax);

        // set values to background in full range
        for (Int_t i = 0; i <= 7; i++){
            fTOFbackground.SetParameter(i,fitFcnBckgrndExc->GetParameter(i));
        }

        //remove fitted function from data
       hTOFm2->GetFunction("fitFcnBckgrndExc"); //->Delete();

        // draw options
        hTOFm2->SetMarkerStyle(kOpenCircle);
        hTOFm2->SetMarkerColor(kBlue);
        hTOFm2->SetLineColor(kBlue);
        hTOFm2->GetXaxis()->SetRangeUser(xmin,xmax);
        hTOFm2->Draw("e1");

        fTOFbackground.SetLineColor(kRed);
        fTOFbackground.DrawCopy("same");

        for(Int_t ival =0; ival <= 7; ival++) {
            StartValuesBackground[ival] = fitFcnBckgrndExc->GetParameter(ival);
            // cout << StartValues[ival] << endl;
        }

        // extract raw signal as data histogram minus background
        Double_t xValue      = 0.595;
        Double_t xStep       = 0.01;
        Double_t yBackground = 0.0;
        Double_t yData       = 0.0;

        int lowbin = 61;
        int highbin = 160;
        if ((ibin == 8 || ibin == 9 || ibin == 10) && part == "aprot") highbin = 180;

        for (int i=lowbin; i<=highbin; i++){
            xValue += xStep;
            yBackground = fTOFbackground.Eval(xValue);
            yData       = hTOFm2->GetBinContent(i);
            fHistHelp->SetBinContent(i, yData - yBackground);
            fHistHelp->SetBinError(i,hTOFm2->GetBinError(i));
        }

        fHistHelp->SetLineColor(kGreen);
        fHistHelp->SetMarkerColor(kGreen);
        fHistHelp->SetLineWidth(2);
        fHistHelp->DrawCopy("e1,same");

        TLegend * leg = new TLegend(0.1, 0.85, 0.4, 0.99);
        leg->AddEntry(hTOFm2, Form("Data %s",part == "prot"? "proton" : "anti-proton"),   "P");
        leg->AddEntry(fHistHelp, "Signal", "LPE");
        leg->AddEntry("fTOFbackground", "Background", "L");
        leg->SetFillColor(0);
        leg->SetTextSize(gStyle->GetTextSize()*0.8);
        leg->Draw();

        TLatex *momentumbin = new TLatex(0.7, 0.6, TString::Format("%g < #it{p} < %g GeV/#it{c}",(ibin-1)*pstep,ibin*pstep));
        momentumbin->SetNDC();
        momentumbin->SetTextFont(42);
        momentumbin->SetTextSize(0.035);
        momentumbin->Draw();

        // yield and error
        yield = fHistHelp->IntegralAndError(lowbin,highbin,yieldErr);
        pur = (fHistHelp->Integral(81,120))/(hTOFm2->Integral(81,120));

        fHistHelp->Reset();

        rawSpectrum->SetBinContent(ibin, yield);
        rawSpectrum->SetBinError(ibin, yieldErr);
        purity->SetBinContent(ibin, pur);

        // Create and plot a box with stats
        c1->Modified();
        c1->Update();

        // Retrieve the stat box
        TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
        ps->SetName("mystats");

        ps->SetX1NDC(0.65);
        ps->SetX2NDC(0.99);
        ps->SetY1NDC(0.7);
        ps->SetY2NDC(0.9);

        TList *list = ps->GetListOfLines();

        // Remove Entries, mean and RMS lines
        TText *tEntries = ps->GetLineWith("Entries");
        TText *tMean    = ps->GetLineWith("Mean");
        TText *tRMS     = ps->GetLineWith("RMS");
        list->Remove(tEntries);
        list->Remove(tMean);
        list->Remove(tRMS);

        int textFont = 42;
        float textSize = 0.03;
        int npointsBckgrnd = 80;

        if (ibin == 8 && part == "prot") {
            fitFcnBckgrndExc->SetChisquare(62.4762);
        }

        if (ibin == 8 && part == "aprot") {
            fitFcnBckgrndExc->SetChisquare(136.615);
            npointsBckgrnd = 60;
        }

        if (ibin == 9 && part == "aprot") {
            fitFcnBckgrndExc->SetChisquare(309.181);
            npointsBckgrnd = 60;
        }

        if (ibin == 10 && part == "aprot") {
            fitFcnBckgrndExc->SetChisquare(201.486);
            npointsBckgrnd = 60;
        }

        // Add lines in the stat box
        TLatex *myt = new TLatex(0, 0, TString::Format("#chi^{2} / ndf (bckgnd only) = %g / %i", fitFcnBckgrndExc->GetChisquare(), npointsBckgrnd ));
        myt->SetTextFont(textFont);
        myt->SetTextSize(textSize);

        TLatex *myt0 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", "yield p", yield, yieldErr));
        myt0->SetTextFont(textFont);
        myt0->SetTextSize(textSize);

        list->Add(myt);
        list->Add(myt0);

        // the following line is needed to avoid that the automatic redrawing of stats
        hTOFm2->SetStats(0);

        c1->Modified();
        c1->Update();

        c1->Print(Form("%s_%s/TOFm2fit_%s_bin_%i.pdf",outd, part, part,ibin));

    } // loop over momentum bins

    // ---------------------------------------------------------------
}

void FitTOFm2Norm(TH3F* hist, TH1F* rawSpectrum){

    // cut on DCA
    hist->GetZaxis()->SetRangeUser(DCAlow,DCAhigh);
    // 2D hist
    TH2F* hTOFm2_p = (TH2F*)hist->Project3D("yx")->Clone("hTOFm2_p");
    // 1D hist
    TH1F* hTOFm2;
    // fit functions
    TF1 fTOFcombined("fTOFcombined", "[0]*1/(TMath::Sqrt(0.5*TMath::Pi())*([2]+[2]*TMath::Erf([3]/TMath::Sqrt(2)))+[2]/[3]*TMath::Exp(-[3]*[3]/2))*TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])* (x < [1]+[3]*[2]) + (x > [1]+[3]*[2])*[0]*1/(TMath::Sqrt(0.5*TMath::Pi())*([2]+[2]*TMath::Erf([3]/TMath::Sqrt(2)))+[2]/[3]*TMath::Exp(-[3]*[3]/2))*TMath::Exp(-(x-[1]-[3]*[2]*0.5)*[3]/[2]) + TMath::Exp([4]+[5]*x+[6]*x*x) + (x <= ([8] + [9]*[10])) * [7] * TMath::Gaus(x, [8], [9]) + (x > ([8]+[9]*[10])) * [7] * TMath::Exp(-(x-[8]-[10]*[9]*0.5)*[10]/[9])", xmin, xmax);

    TF1 fTOFsignal("fTOFsignal","[0]*1/(TMath::Sqrt(0.5*TMath::Pi())*([2]+[2]*TMath::Erf([3]/TMath::Sqrt(2)))+[2]/[3]*TMath::Exp(-[3]*[3]/2))*TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])* (x < [1]+[3]*[2]) + (x > [1]+[3]*[2])*[0]*1/(TMath::Sqrt(0.5*TMath::Pi())*([2]+[2]*TMath::Erf([3]/TMath::Sqrt(2)))+[2]/[3]*TMath::Exp(-[3]*[3]/2))*TMath::Exp(-(x-[1]-[3]*[2]*0.5)*[3]/[2])",xmin,xmax);

    TF1 fTOFbackground("fTOFbackground", "TMath::Exp([0]+[1]*x+[2]*x*x) + (x <= ([4] + [5]*[6])) * [3] * TMath::Gaus(x, [4], [5]) + (x > ([4]+[5]*[6])) * [3] * TMath::Exp(-(x-[4]-[6]*[5]*0.5)*[6]/[5])", xmin, xmax);

    // start values for bin 26
    Double_t StartValues[11] = {9888, 0.88, 0.082, 0.95, 8.6, -1.9, 0.0002, 67700, 0.22, 0.09, 1.3};
    // set par limits
    fTOFcombined.SetParLimits(1,0.8,1.2); // mean proton
    fTOFcombined.SetParLimits(2,0.0,0.2); // width proton
    fTOFcombined.SetParLimits(3,0.,2.);   // tail proton

    fTOFcombined.SetParLimits(4,6.,10.);    // bckgrnd constant
    fTOFcombined.SetParLimits(5,-2.0,-1.5); // bckgrnd lin. term
    fTOFcombined.SetParLimits(6,-0.,0.3);   // bckgrnd quad. term

    fTOFcombined.SetParLimits(8,0.15,0.4); // mean kaon
    fTOFcombined.SetParLimits(9,0.0,0.2); // width kaon
    fTOFcombined.SetParLimits(10,0.,2.);   // tail kaon

    // set par names
    fTOFcombined.SetParName(0,"yield p"); // yield proton (norm)
    fTOFcombined.SetParName(1,"M^{2}_{p}"); // mean proton
    fTOFcombined.SetParName(2,"#sigma_{p}"); // width proton
    fTOFcombined.SetParName(3,"#tau_{p}");   // tail proton

    fTOFcombined.SetParName(4,"Bkg const.");    // bckgrnd constant
    fTOFcombined.SetParName(5,"Bkg lin. term"); // bckgrnd lin. term
    fTOFcombined.SetParName(6,"Bkg quad. term");   // bckgrnd quad. term

    fTOFcombined.SetParName(7,"yield K"); // yield kaon
    fTOFcombined.SetParName(8,"M^{2}_{K}"); // mean kaon
    fTOFcombined.SetParName(9,"#sigma_{K}"); // width kaon
    fTOFcombined.SetParName(10,"#tau_{K}");   // tail kaon

    // yield and yield error
    Double_t yield = 1.0;
    Double_t yieldErr = 0.0;

    // canvas
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->SetLogy();

    // fit in momentum bins starting from bin 2.5-2.6
    for (Int_t ibin = binmid2; ibin <= binhigh; ibin++){

        hTOFm2 = (TH1F*)hTOFm2_p->ProjectionY("_py",ibin,ibin,"oe")->Clone(Form("hTOFm2_ptbin_%i",ibin));

        // set (start) values
        for (Int_t i = 0; i <= 10; i++){
            fTOFcombined.SetParameter(i,StartValues[i]);
        }

        hTOFm2->Fit(&fTOFcombined,"QNR");

        // draw options
        hTOFm2->SetMarkerStyle(kOpenCircle);
        hTOFm2->SetMarkerColor(kBlue);
        hTOFm2->SetLineColor(kBlue);
        hTOFm2->GetXaxis()->SetRangeUser(xmin,xmax);
        hTOFm2->Draw("e1");
        // total fit
        fTOFcombined.SetLineColor(kBlack);
        fTOFcombined.DrawCopy("same");
        // signal and background parameters
        for(Int_t ival =0; ival <= 10; ival++) {
            if (ival <=3){
                fTOFsignal.SetParameter(ival,fTOFcombined.GetParameter(ival));
            } else {
                fTOFbackground.SetParameter(ival-4,fTOFcombined.GetParameter(ival));
            }
        }
        fTOFsignal.SetLineColor(kGreen);
        fTOFsignal.DrawCopy("same");
        fTOFbackground.SetLineColor(kRed);
        fTOFbackground.DrawCopy("same");

        TLegend * leg = new TLegend(0.1, 0.85, 0.4, 0.999);
        leg->AddEntry(hTOFm2, Form("Data %s",part == "prot"? "proton" : "anti-proton"),   "P");
        leg->AddEntry("fTOFsignal", "Signal", "L");
        leg->AddEntry("fTOFbackground", "Background", "L");
        leg->AddEntry("fTOFcombined", "Total fit", "L");
        leg->SetFillColor(0);
        leg->SetTextSize(gStyle->GetTextSize()*0.8);
        leg->Draw();

        TLatex *momentumbin = new TLatex(0.7, 0.45, TString::Format("%g < #it{p} < %g GeV/#it{c}",(ibin-1)*pstep,ibin*pstep));
        momentumbin->SetNDC();
        momentumbin->SetTextFont(42);
        momentumbin->SetTextSize(0.035);
        momentumbin->Draw();

        for(Int_t ival =0; ival <= 10; ival++) {
            StartValues[ival] = fTOFcombined.GetParameter(ival);
            // cout << StartValues[ival] << endl;
        }

        // Create and plot a box with stats
        c1->Modified();
        c1->Update();

        // Retrieve the stat box
        TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
        ps->SetName("mystats");

        ps->SetX1NDC(0.7);
        ps->SetX2NDC(0.99);
        ps->SetY1NDC(0.5);
        ps->SetY2NDC(0.9);

        TList *list = ps->GetListOfLines();

        // Remove Entries, mean and RMS lines
        TText *tEntries = ps->GetLineWith("Entries");
        TText *tMean    = ps->GetLineWith("Mean");
        TText *tRMS     = ps->GetLineWith("RMS");
        list->Remove(tEntries);
        list->Remove(tMean);
        list->Remove(tRMS);

        int textFont = 42;
        float textSize = 0.03;
        int npoints = 180;

        // Add lines in the stat box
        TLatex *myt = new TLatex(0, 0, TString::Format("#chi^{2} / ndf = %g / %i", fTOFcombined.GetChisquare(), npoints ));
        myt->SetTextFont(textFont);
        myt->SetTextSize(textSize);

        TLatex *myt0 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", fTOFcombined.GetParName(0), (fTOFcombined.GetParameter(0))/(hTOFm2->GetBinWidth(1)), (fTOFcombined.GetParError(0))/((hTOFm2->GetBinWidth(1)))));
        myt0->SetTextFont(textFont);
        myt0->SetTextSize(textSize);

        TLatex *myt1 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", fTOFcombined.GetParName(1), fTOFcombined.GetParameter(1), fTOFcombined.GetParError(1)));
        myt1->SetTextFont(textFont);
        myt1->SetTextSize(textSize);

        TLatex *myt2 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", fTOFcombined.GetParName(2), fTOFcombined.GetParameter(2), fTOFcombined.GetParError(2)));
        myt2->SetTextFont(textFont);
        myt2->SetTextSize(textSize);

        TLatex *myt3 = new TLatex(0, 0, TString::Format("%s = %g #pm %g", fTOFcombined.GetParName(3), fTOFcombined.GetParameter(3), fTOFcombined.GetParError(3)));
        myt3->SetTextFont(textFont);
        myt3->SetTextSize(textSize);

        list->Add(myt);
        list->Add(myt0);
        list->Add(myt1);
        list->Add(myt2);
        list->Add(myt3);
void count(TH2F* hData_DCAxy_p){

  for (size_t bin = 0; bin < hData_DCAxy_p->GetNbinsY(); bin++) {
    /* code */
    TH1F* hData_DCAxy;
    hData_DCAxy = (TH1F*)hData_DCAxy_p->ProjectionY("_py",bin,bin,"oe")->Clone(Form("hData_DCAxy_ptbin_%i",bin));

  //hData_DCAxy->Rebin(2);
  cout<<"********************************************************"<<endl;
  cout<<"******************     "<<bin<<"      *******************"<<endl;
  for (size_t i = 0; i < hData_DCAxy->GetNbinsX(); i++) {
    double d = hData_DCAxy->GetBinContent(i);
    if (d != 0.0){
    cout<<"****** Bin Num :"<<i<<"   ******** is : "<<d<<"**********"<<endl;
    }
  }
}
}
        // the following line is needed to avoid that the automatic redrawing of stats
        hTOFm2->SetStats(0);

        c1->Modified();
        c1->Update();

        c1->Print(Form("%s_%s/TOFm2fit_%s_bin_%i.pdf",outd, part, part,ibin));

//        // yield and stat error
//        yield    = (fTOFsignal.Integral(0.,2.0))/(hTOFm2->GetBinWidth(1));
//        yieldErr = TMath::Sqrt(yield);
//        rawSpectrum->SetBinContent(ibin, yield);
//        rawSpectrum->SetBinError(ibin, yieldErr);

        rawSpectrum->SetBinContent(ibin, (fTOFcombined.GetParameter(0))/(hTOFm2->GetBinWidth(1)));
        rawSpectrum->SetBinError(ibin, (fTOFcombined.GetParError(0))/(hTOFm2->GetBinWidth(1)));

    } // loop over momentum bins
}

void WriteHistos(TH1F* spectrum, TH1F* purity){

    TFile* fileout = new TFile(Form("%s_%s.root",fout,part),"RECREATE");
    spectrum->Write();
    purity->Write();
    fileout->Close();

}

Double_t TOFbackgroundExcludePeak(Double_t* x, Double_t *par){

    if (x[0] > xmin_peak && x[0] < xmax_peak) {
        TF1::RejectPoint();
        return 0;
    }

    return par[0] + par[1]*x[0] + par[2]*TMath::Exp(par[3]*x[0]) + (x[0] <= (par[5] + par[6]*par[7])) * par[4] * TMath::Gaus(x[0], par[5], par[6]) + (x[0] > (par[5]+par[6]*par[7])) * par[4] * TMath::Exp(-(x[0]-par[5]-par[7]*par[6]*0.5)*par[7]/par[6]);

}

void count(TH2F* hData_DCAxy_p){

  for (size_t bin = 0; bin < hData_DCAxy_p->GetNbinsY(); bin++) {
    /* code */
    TH1F* hData_DCAxy;
    hData_DCAxy = (TH1F*)hData_DCAxy_p->ProjectionY("_py",bin,bin,"oe")->Clone(Form("hData_DCAxy_ptbin_%i",bin));

  //hData_DCAxy->Rebin(2);
  cout<<"********************************************************"<<endl;
  cout<<"******************     "<<bin<<"      *******************"<<endl;
  for (size_t i = 0; i < hData_DCAxy->GetNbinsX(); i++) {
    double d = hData_DCAxy->GetBinContent(i);
    if (d != 0.0){
    cout<<"****** Bin Num :"<<i<<"   ******** is : "<<d<<"**********"<<endl;
    }
  }
}
}
