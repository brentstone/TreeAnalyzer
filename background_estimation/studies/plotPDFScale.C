#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include <vector>
#include "HistoPlotting/include/Plotter.h"

using namespace std;

double getComb(TFile *f, TString prefix, TString suffix) {
	double tot = 0.0;
	for(TString s : {"1l","2l"}) {
		TH1 *h = (TH1*)f->Get(prefix+"_"+s+"_"+suffix);
		if(!h)
			cout<<prefix+"_"+s+"_"+suffix<<endl;
		tot += h->GetBinContent(1);
	}
	return tot;
}



void doOneSignal(TString ss, TFile *f) {

//    TGraph * gn = new TGraph;
//    TGraph * gsu = new TGraph;
//    TGraph * gsd = new TGraph;
//    TGraph * gpu = new TGraph;
//    TGraph * gpd = new TGraph;
//    TGraph * gtu = new TGraph;
//    TGraph * gtd = new TGraph;

	std::vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};

	double avgUnc = 0;

    auto doOneMassPoint = [&](int nP, int sM, TString nameWithMass) {
    	double nomAcc = 0;

    	printf("\nMX = %d\n",sM);
    	// acceptance, no pdf or scale weights applied
    	TH1 *nom_incl = (TH1*)f->Get(nameWithMass+"_nom_yield");
    	nomAcc = getComb(f,nameWithMass,"nom_yield") / nom_incl->GetBinContent(1);

//    	printf("nomacc = %f\n\nScale accs / nom acc:\n",nomAcc);
    	std::vector<double> scales, pdfs;

    	// get scale acceptances
    	for(int s=0; s<8; s++) {
    		TH1 *scale_incl = (TH1*)f->Get(nameWithMass+TString::Format("_s_%d_yield",s));
    		double acc = getComb(f,nameWithMass,TString::Format("s_%d_yield",s)) / scale_incl->GetBinContent(1);
    		scales.push_back(acc);
//    		std::cout<<acc/nomAcc<<std::endl;
    	}

    	// get pdf acceptances
        double x = 0;
        double x2 = 0;

        double pdfVariance = 0;

        double acc0 = 0;
    	for(int p=0; p<101; p++) {
    		TH1 *pdf_incl = (TH1*)f->Get(nameWithMass+TString::Format("_p_%d_yield",p));
    		double acc = getComb(f,nameWithMass,TString::Format("p_%d_yield",p)) / pdf_incl->GetBinContent(1);

    		if(p == 0) {
    			acc0 = acc;
    		} else {
        		pdfs.push_back(acc);
        		pdfVariance += (acc - acc0)*(acc - acc0);
    		}
    	}

    	std::sort(scales.begin(),scales.end(), [](double a, double b) {return a < b;} );
    	std::sort(pdfs.begin()  ,pdfs.end()  , [](double a, double b) {return a < b;} );

    	double pdfUnc = std::sqrt(pdfVariance);
    	double scaleUnc = (scales.back()-scales.front())/2;
    	double unc = std::sqrt(pdfVariance+scaleUnc*scaleUnc);

//    	printf("\nPDF uncertainty = %f\n",pdfUnc);
//    	printf("Generous scale uncertainty = %f\n",scaleUnc);
    	printf("Total uncertainty = %f\n",unc);

    	avgUnc += unc/float(masses.size());

//        gn->SetPoint(nP,float(sM),nomAcc);
//        gsu->SetPoint(nP,float(sM),scales.back() /nomAcc -1.0 );
//        gsd->SetPoint(nP,float(sM),scales.front()/nomAcc-1.0);

//        gpu->SetPoint(nP,float(sM),(mean+stdDev)/nomAcc-1.0);
//        gpd->SetPoint(nP,float(sM),(mean-stdDev)/nomAcc-1.0);
//        gtu->SetPoint(nP,float(sM),std::sqrt((mean+stdDev)*(mean+stdDev)+scales.back() *scales.back() )/nomAcc);
//        gtd->SetPoint(nP,float(sM),std::sqrt((mean-stdDev)*(mean-stdDev)+scales.front()*scales.front())/nomAcc);
    };

	for(unsigned im = 0; im < masses.size(); ++im)
		doOneMassPoint(im,masses[im],TString::Format("%s_m%d",ss.Data(),masses[im]));

	std::cout<<std::endl<<"avgUnc = "<<avgUnc<<std::endl;
	std::cout<<"-------------------------------------------------------------------------------------------"<<std::endl;
//    Plotter * p = new Plotter;
//
//    p->addGraph(gn, "nominal");
//    p->addGraph(gsd,"scale: Down");
//    p->addGraph(gsu, "scale: Up");
//    p->addGraph(gpd,"pdf: Down");
//    p->addGraph(gpu, "pdf: Up");
//    p->addGraph(gtd,"tot: Down");
//    p->addGraph(gtu, "tot: Up");
//
////    p->setMinMax(-0.5,0.5);
//    p->setYTitle(ss);
//
//    p->draw(false,ss);

}

void plotPDFScale(TString year) {

	TFile *fr = TFile::Open("/Users/brentstone/Dropbox/Physics/HHbbWW/pdfscale/pdfScaleDist"+year+"_spin0.root");
	TFile *fb = TFile::Open("/Users/brentstone/Dropbox/Physics/HHbbWW/pdfscale/pdfScaleDist"+year+"_spin2.root");

	doOneSignal("radion",fr);
	doOneSignal("blkgrav",fb);

}






