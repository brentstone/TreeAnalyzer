//Fit HHRes
{
  TFile * f = new TFile("hSolTrees_HPTRes.root");
  std::vector<float> ptBins {300,400,500,600,800,1000,1200,1400,1600,1800,2000};
  
  TGraphErrors * gF = new TGraphErrors();
  TGraphErrors * gC = new TGraphErrors();
  TGraphErrors * gCor = new TGraphErrors();
  int iGP = 0;
  for(unsigned int iP = 0; iP < ptBins.size()-1; ++iP){
    TString hName = TString::Format("signal_pt_%.0f_%.0f",ptBins[iP],ptBins[iP+1]);    
    TH1 * h = 0;
    TH1 * hc = 0;
    TH1 * hPT = 0;
    
    f->GetObject(hName+"_HPtRes",h);         
    f->GetObject(hName+"_true_Hpt",hPT);
    f->GetObject(hName+"_corrHPtRes",hc);
    std::cout <<hName<<" "<< h <<" "<< hPT << std::endl;
    if(h==0||hPT==0) continue;
    gF->SetPoint(iGP,hPT->GetMean(),h->GetMean() );
    gF->SetPointError(iGP,0,h->GetMeanError());
    h->GetXaxis()->SetRange(31,71);
    gC->SetPoint(iGP,hPT->GetMean(),h->GetMean() );
    gC->SetPointError(iGP,0,h->GetMeanError());
    
    if(hc){
      hc->GetXaxis()->SetRange(31,71);
      gCor->SetPoint(iGP,hPT->GetMean(),hc->GetMean() );
      gCor->SetPointError(iGP,0,hc->GetMeanError());              
    }

    iGP++;    
  }
  
  Plotter * p = new Plotter();
  p->addGraph(gF,"Full range");
  p->addGraph(gC,"0.6 #rightarrow 1.4");
  p->addGraph(gCor,"Corrected 0.6 #rightarrow 1.4");
  
  p->draw();
}


///Look at dists
{
  TFile * f = new TFile("hhtrain_plots.root");
  std::vector<TString> sigs = {"radion","graviton"};
  std::vector<int> masses = {800,1000,3500};
  std::vector<TString> vars = {"hwwMag","true_hwwMag","hwwMagRes","extraMetPerp","extraMetParRelhwwMag","wqqPTRes","virtualSD_approxW","onshellSD_approxW","virtualSD_approxH","onshellSD_approxH","virtualSD_qqSDMassCoarse","onshellSD_qqSDMassCoarse"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    for(unsigned int iM = 0; iM < masses.size(); ++iM){
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = TString::Format("%s_m%i",sigs[iS].Data(),masses[iM]) +"_"+vars[iV];
      f->GetObject(hname,h);      
      p->addHistLine(h,TString::Format("%s: %.1f",sigs[iS].Data(),masses[iM]/1000.));
    }
  }
    p->normalize();
    p->drawSplitRatio(0,"stack",false,false,vars[iV]);
  }
  
}


///Figure out binning (fine)
{
  TFile * f = new TFile("hSolTrees_temp.root");
  std::vector<TString> sigs = {"signal_low","signal_test","signal_high"};
  std::vector<TString> shells = {"","_vqq","_osqq"};
  std::vector<double> avgs(3);
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    TH1 * h = 0;
    f->GetObject(sigs[iS]+"_true_hwwMag",h);
    if(h) avgs[iS] = h->GetMean();
  }      
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","approxW","approxH","approxJ","approxJ2","qqSDMassCoarse","hWW","Wlnu"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[1] +"_"+vars[iV];
      f->GetObject(hname,h);      
      if(h==0) continue;
      nH++;    
      p->addHistLine(h,sigs[iS]);
  }
  if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    p->draw(false,vars[iV]);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;    
      for(unsigned int iS = 0; iS < shells.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[0]+shells[iS] +"_"+vars[iV];
      f->GetObject(hname,h); 
      if(h==0) continue;
      nH++;         
      p->addHistLine(h,sigs[0]+shells[iS]);
  }
    if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]+"_shells");
    p->draw(false,vars[iV]+"_shells");
    // objs.push_back(c);
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}
///Figure out binning (final binning)
{
  TFile * f = new TFile("hSolTrees_temp.root");
  std::vector<TString> sigs = {"signal_low","signal_test","signal_high"};
  std::vector<TString> shells = {"_vqq","_osqq"};
    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","approxW","approxH","qqSDMassCoarse"};
  std::vector<int> vBinV      = {30  ,50  ,50  ,15,25,1};
  std::vector<double> vBinLEV = {-150,-0.6,-0.5,60,90,1};
  std::vector<double> vBinHEV = {150 ,0.4 ,0.5 ,90,190,1};
  
  std::vector<int> vBinO      = {30  ,50  ,50  ,15,20,1};
  std::vector<double> vBinLEO = {-150,-0.6,-0.5,0 ,110,1};
  std::vector<double> vBinHEO = {150 ,0.4 ,0.5 ,60,150,1};
  
  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
            Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[iSh] +"_"+vars[iV];                  
      f->GetObject(hname,h);      
      if(!h) continue;
      
      int nNB; double nL; double nH;
      if(iSh == 0){
        nNB = vBinV[iV];
        nL = vBinLEV[iV];
        nH = vBinHEV[iV];
      } else {
        nNB = vBinO[iV];
        nL = vBinLEO[iV];
        nH = vBinHEO[iV];
      }
      
      if(nL==nH) continue;
      TH1 * hR = new TH1F(TString(h->GetName()) +"_reb",h->GetTitle(),nNB,nL,nH);
      hR->SetBinContent(0,h->GetBinContent(0));
      hR->SetBinContent(nNB+1,h->GetBinContent(h->GetNbinsX()+1));      
      for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB ){
        int hRB = hR->FindFixBin(h->GetBinCenter(iB));
        
        hR->SetBinContent(hRB, hR->GetBinContent(hRB) + h->GetBinContent(iB) );
      }
      
      
      p->addHistLine(hR,sigs[iS]);
  
}
    p->normalize();
    p->setBotMinMax(0,2);
    TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]);
    objs.push_back(c);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}


TEST RESULTS
  //Mass shape comparison
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"radion_m1000","radion_m2000","ttbar","wjets"};
  std::vector<TString> vars = {"simpleHH","chi2HH","likeliHH"};

  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * p = new Plotter();
    for(unsigned int iV = 0; iV < vars.size(); ++iV){
      TH1* h = 0;
      f->GetObject(sigs[iS]+"_"+vars[iV],h);
      if(h==0) continue;
      p->addHistLine(h,vars[iV]);
    }
    p->setBotMinMax(0,2);
    p->drawSplitRatio(0,"stack",false,false,sigs[iS]);
  }


}

//Efficiency curves of the likli/chi
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"radion_m1000","radion_m2000","ttbar_m1000","wjets_m1000","ttbar_m2000","wjets_m2000"};  
  std::vector<TString> vars = {"chi2","likeli","likeliNoSD"};


    for(unsigned int iV = 0; iV < vars.size(); ++iV){  
    Plotter * p = new Plotter();
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1* h = 0;
      f->GetObject(sigs[iS]+"_"+vars[iV],h);
      if(h==0) continue;
      h = PlotTools::getIntegral(h,false,true);
      p->addHistLine(h,sigs[iS]);
    }
    // p->setBotMinMax(0,2);
    // p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // p->normalize();
    p->draw(false,vars[iV]);
  }
  
  
}


///Straight eff check
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets"};
  std::vector<TString> vars = {"chi2HH","likeliHH"};
  std::vector<TString> cuts = {"chilt20","chilt9","chilt5","lllt37","lllt32","lllt29p5"};


  //start With sigal
  Plotter * p = new Plotter();
  TH1* hd = 0;
  f->GetObject("signal_sampParam",hd);
  if(hd!= 0){
    hd = (TH1*) hd->Clone();
    p->addHistLine(hd,"noCut");
    for(unsigned int iC = 0; iC < cuts.size(); ++iC){          
      TH1* h = 0;
      f->GetObject("signal_"+cuts[iC]+"_sampParam",h);      
      if(h==0) continue;  
      p->addHistLine(h,cuts[iC]);
    }
    p->setMinMax(0.5,1);
    p->drawRatio(0,"stack",false,false,"signal_sampParam");
  }
  
  //Now BKG
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * p = new Plotter();
    TH1* hd = 0;
    f->GetObject(sigs[iS]+"_simpleHH",hd);      
    if(hd==0) continue;
    hd = (TH1*) hd->Clone();
    p->addHistLine(hd,"noC");
  for(unsigned int iV = 0; iV < vars.size(); ++iV){  
      for(unsigned int iC = 0; iC < cuts.size(); ++iC){            
    TH1* h = 0;
    f->GetObject(sigs[iS]+"_"+cuts[iC]+"_"+vars[iV],h);      
      if(h==0) continue;      
      p->addHistLine(h,cuts[iC]);
    }

  }
  p->setMinMax(0,1.5);
  p->drawRatio(0,"stack",false,false,sigs[iS]+"_cutComps");  
}
 
}


///BKG Shape
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets"};
  std::vector<TString> vars = {"chi2HH","likeliHH"};
  std::vector<TString> dens = {"simpleHH","simpleHH"};

  std::vector<TString> cuts = {"chiltinf","chilt20","chilt9","chilt5","llltinf","lllt37","lllt32","lllt29p5"};


    for(unsigned int iV = 0; iV < vars.size(); ++iV){  
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      Plotter * p = new Plotter();
      TH1* hd = 0;
      f->GetObject(sigs[iS]+"_"+dens[iV],hd);
      
      if(hd==0) continue;
      hd = (TH1*) hd->Clone();
      hd->Scale(1./hd->Integral(hd->FindFixBin(1000),hd->FindFixBin(4000)));
              p->addHistLine(hd,"noC");
  for(unsigned int iC = 0; iC < cuts.size(); ++iC){
    
    
    
          
    TH1* h = 0;
    f->GetObject(sigs[iS]+"_"+cuts[iC]+"_"+vars[iV],h);      
      if(h==0) continue;      
            h->Scale(1./h->Integral(hd->FindFixBin(1000),hd->FindFixBin(4000)));
      p->addHistLine(h,cuts[iC]);
    }
    p->setBotMinMax(0,2);
    // p->normalize();
    p->drawSplitRatio(0,"stack",false,false,sigs[iS]+"_"+vars[iV]);

    // p->draw(false,vars[iV]);
  }
  
}
  
  
}


{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets"};
  std::vector<TString> vars = {"chi2HH","likeliHH"};
  std::vector<TString> dens = {"simpleHH","simpleHH"};

  std::vector<TString> cuts1 = {"chiltinf","chilt20","chilt9","chilt5"};
  std::vector<TString> cuts2 = {"llltinf","lllt37","lllt32","lllt29p5"};


    for(unsigned int iC = 0; iC < cuts1.size(); ++iC){  
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      Plotter * p = new Plotter();
      TH1* hd = 0;
      f->GetObject(sigs[iS]+"_"+dens[0],hd);
      
      if(hd==0){
        std::cout << sigs[iS]+"_"+dens[iC] << std::endl;
        continue;
      } 
      hd = (TH1*) hd->Clone();
              // p->addHistLine(hd,"noC");
              
              TH1* h1 = 0;
              f->GetObject(sigs[iS]+"_"+cuts1[iC]+"_"+vars[0],h1);      
                if(h1==0) {
                  std::cout << sigs[iS]+"_"+cuts1[iC]+"_"+vars[0]<<std::endl;
                  continue;      
                }
                p->addHistLine(h1,cuts1[iC]+"_"+vars[0]);        
                
                TH1* h2 = 0;
                  f->GetObject(sigs[iS]+"_"+cuts2[iC]+"_"+vars[1],h2);      
                  if(h2==0){
                                      std::cout << sigs[iS]+"_"+cuts2[iC]+"_"+vars[1]<<std::endl;
                                      continue;      
                  } 
                  p->addHistLine(h2,cuts2[iC]+"_"+vars[1]);     
                  
                  
                  p->setBotMinMax(0,2);
                  p->drawSplitRatio(0,"stack",false,false,sigs[iS]+"_"+cuts1[iC]);   
              
              

  }
  
}
  
  
}



///////. BACKGROUND LIKLI


///Figure out binning (fine)
{
  TFile * f = new TFile("hSolTrees_bkgTempStudy_bkg.root");
  std::vector<TString> sigs = {"ttbarPW","ttbar","wjets"};
  std::vector<TString> shells = {"","_low","_test","_high"};

  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","approxW","approxH2"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[0] +"_"+vars[iV];
      f->GetObject(hname,h);      
      if(h==0) continue;
      nH++;    
      p->addHistLine(h,sigs[iS]);
  }
  if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    p->draw(false,vars[iV]);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;    
      for(unsigned int iS = 0; iS < shells.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[0]+shells[iS] +"_"+vars[iV];
      f->GetObject(hname,h); 
      if(h==0) continue;
      nH++;         
      p->addHistLine(h,sigs[0]+shells[iS]);
  }
    if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]+"_shells");
    p->draw(false,vars[iV]+"_shells");
    // objs.push_back(c);
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}
///Figure out binning (final binning)
{
  TFile * f = new TFile("hSolTrees_bkgTempStudy_bkg.root");
  std::vector<TString> sigs = {"ttbarPW"};
  std::vector<TString> shells = {"_low","_high"};
    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","approxW","approxH"};
  std::vector<int> vBinV      = {30  ,50  ,20,50};
  std::vector<double> vBinLEV = {-150,-0.6,60,100};
  std::vector<double> vBinHEV = {150 ,0.4 ,100,600};
  
  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
            Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[iSh] +"_"+vars[iV];                  
      f->GetObject(hname,h);      
      if(!h) continue;
      
      int nNB; double nL; double nH;

        nNB = vBinV[iV];
        nL = vBinLEV[iV];
        nH = vBinHEV[iV];
  
  
      
      if(nL==nH) continue;
      TH1 * hR = new TH1F(TString(h->GetName()) +"_reb",h->GetTitle(),nNB,nL,nH);
      hR->SetBinContent(0,h->GetBinContent(0));
      hR->SetBinContent(nNB+1,h->GetBinContent(h->GetNbinsX()+1));      
      for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB ){
        int hRB = hR->FindFixBin(h->GetBinCenter(iB));
        
        hR->SetBinContent(hRB, hR->GetBinContent(hRB) + h->GetBinContent(iB) );
      }
      
      
      p->addHistLine(hR,sigs[iS]);
  
}
    p->normalize();
    p->setBotMinMax(0,2);
    TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]);
    objs.push_back(c);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}

///ROC CURVES
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"radion_m1000","radion_m2000"};
  std::vector<TString> bkgCs = {"m1000","m2000"};
  
  // std::vector<TString> sigs = {"radion_m1000_lllt32"};
  // std::vector<TString> bkgCs = {"lllt32"};
  std::vector<TString> bkgs = {"ttbar","wjets","qcd"};
  
  
  // std::vector<TString> vars = {"chi2","likeli","Blikeli","SoBlikeli"};
  std::vector<TString> vars = {"likeli","Blikeli","SoBlikeli"};
    // std::vector<TString> vars = {"likeliBHWW","likeliHWW"};
  int cutGT = 1;
  for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
    Plotter * p = new Plotter();
    for(unsigned int iV = 0; iV < vars.size(); ++iV){
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
        TH1* hs = 0;
        f->GetObject(sigs[iS]+"_"+vars[iV],hs);     
        if(hs==0)continue;
        hs=(TH1*)hs->Clone();
        TH1* hb = 0;
        f->GetObject(bkgs[iB]+"_"+bkgCs[iS]+"_"+vars[iV],hb);     
        if(hb==0)continue;
        hb=(TH1*)hb->Clone();
        
        auto * roc = PlotTools::getRocCurve(hs,hb,iV==cutGT, "signal eff","bkg. eff");
        p->addGraph(roc,sigs[iS]+"_"+vars[iV]);                
      }                  
    }
    p->draw(false,bkgs[iB])      ;
  }
}


////QCD TESTING
{
  TFile * f = new TFile("hSolTrees_qcdTempStudy.root");
  std::vector<TString> sigs = {"qcd","radion_m1000","radion_m2000"};
  std::vector<TString> vars = {"extraMetParRelhwwMag","approxH","wPTRelhwwMag","extraMetParRelhwwMag2","approxH3"};
    // std::vector<TString> vars = {"deltaPhi","deltaPhi2","deltaTheta","deltaTheta2"};
    // std::vector<TString> vars = {"approxW","approxH","ptDR","approxW2","approxH2","ptDR2","DR","DR2"};

  std::vector<TString> cuts = {""};


    for(unsigned int iC = 0; iC < cuts.size(); ++iC){  
          for(unsigned int iV = 0; iV < vars.size(); ++iV){  
                  Plotter * p = new Plotter();
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){

              
              TH1* h1 = 0;
              f->GetObject(sigs[iS]+"_"+cuts[iC]+vars[iV],h1);      
                if(h1==0) {
                  std::cout << sigs[iS]+"_"+cuts[iC]+vars[iV]<<std::endl;
                  continue;      
                }
                p->addHistLine(h1,sigs[iS]);                                                    

  }
  p->rebin(4);
  p->normalize();
  p->draw(false,vars[iV]+"_"+cuts[iC]);   
  
}
  
  
}
}