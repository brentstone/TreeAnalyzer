
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "TRandom3.h"
#include <string.h>



class StatTesterAnalyzer {
public:

    StatTesterAnalyzer(const TH1D* inputModel, const TH1D* inputData, const unsigned int nToys, const bool saveToys,std::string hName, std::string outName = "")
{

        const auto * model = inputModel->GetArray();
        const auto * data  = inputData->GetArray();
        const int nBins    = inputModel->GetSize();

        ts_nom_sa = calcSA(nBins,model,data);
        ts_nom_ks = calcKS(nBins,model,data);
        std::vector<double> ts_sa; ts_sa.reserve(nToys);
        std::vector<double> ts_ks; ts_ks.reserve(nToys);

        std::unique_ptr<TRandom3> rand(new TRandom3(1234));

        for(unsigned int iT = 0; iT < nToys; ++iT){
            auto toyData = generatePseudoData(nBins,model,&*rand);
            ts_sa.push_back(calcSA(nBins,model,&toyData[0]));
            ts_ks.push_back(calcKS(nBins,model,&toyData[0]));
            if(saveToys){
                auto toyH = makeTH1(std::string("toy_")+ASTypes::int2Str(iT),nBins,&toyData[0],inputData );
                plotter.add1D(toyH);
            }
        }
        processTS(ts_nom_sa,ts_sa,hName +"_saturated",true);
        processTS(ts_nom_ks,ts_ks,hName+"_ks",false);
//        std::cout << outName<<std::endl;
        if(outName.size())plotter.write(outName);
}

    typedef std::pair<TH1D*,TH1D*> ModelAndData;
    StatTesterAnalyzer(const ModelAndData& nominal, const std::vector<ModelAndData>& toys,std::string hName,  std::string outName = "")
{
        const int nBins    = nominal.first->GetSize();
        const unsigned int nToys = toys.size();
        ts_nom_sa = calcSA(nBins,nominal);
        ts_nom_ks = calcKS(nBins,nominal);
        std::vector<double> ts_sa; ts_sa.reserve(nToys);
        std::vector<double> ts_ks; ts_ks.reserve(nToys);

        for(unsigned int iT = 0; iT < nToys; ++iT){
            ts_sa.push_back(calcSA(nBins,toys[iT]));
            ts_ks.push_back(calcKS(nBins,toys[iT]));
        }

        processTS(ts_nom_sa,ts_sa,hName+"_saturated",true);
        processTS(ts_nom_ks,ts_ks,hName+"_ks",false);
//        std::cout << outName<<std::endl;
        if(outName.size())plotter.write(outName);
}

    typedef std::pair<TH2D*,TH2D*> ModelAndData2D;
    StatTesterAnalyzer(const ModelAndData2D& nominal, const std::vector<ModelAndData2D>& toys,std::string hName, std::string outName = "")
{
        const int nBinsX  = nominal.first->GetNbinsX();
        const int nBinsY  = nominal.first->GetNbinsY();
        const unsigned int nToys = toys.size();
        ts_nom_sa = calcSA2D(nBinsX,nBinsY,nominal);
        std::vector<double> ts_sa; ts_sa.reserve(nToys);

        for(unsigned int iT = 0; iT < nToys; ++iT){
            ts_sa.push_back(calcSA2D(nBinsX,nBinsY,toys[iT]));
        }

        processTS(ts_nom_sa,ts_sa,hName+"_saturated2D",true);
        if(outName.size())plotter.write(outName);
}

    StatTesterAnalyzer(const std::vector<ModelAndData2D> allNoms, const std::vector<std::vector<ModelAndData2D>> allToys, std::string outName = "")
{
    	// loop over the categories (dif selections):
    	if(allNoms.size() != allToys.size())
    		throw std::invalid_argument("Sizes of nominal and toys must be the same!");

        std::cout<<allNoms.size()<<std::endl;
        std::cout<<allToys.size()<<std::endl;
        const unsigned int nToys = allToys.front().size();
        std::cout<<"ntoys = "<<nToys<<std::endl;
        std::vector<double> ts_sa(nToys,0.0);
        std::cout<<ts_sa.size()<<std::endl;


    	for(unsigned i=0; i<allNoms.size(); ++i) {
    		std::cout<<"in loop"<<std::endl;
    		auto nominal = allNoms[i];
    		auto toys = allToys[i];
//            const int nBinsX  = nominal.first->GetNbinsX();
//            const int nBinsY  = nominal.first->GetNbinsY();
//            const int nBinsX  = nominal.first->GetNbinsX();
//            const int nBinsY  = nominal.first->GetNbinsY();
            ts_nom_sa += calcSA2D(30,82,nominal);
            ts_nom_sa += 1.0;

            for(unsigned int iT = 0; iT < nToys; ++iT){
            	std::cout<<"toy "<<iT<<" "<<toys[iT].first<<" "<<toys[iT].second<<std::endl;
//            	std::cout<<"name "<<toys[iT].first->GetName()<<toys[iT].second->GetName()<<std::endl;
//
//            	std::cout<<toys[iT].first->GetNbinsX()<<toys[iT].first->GetNbinsY()<<std::endl;
//            	ts_sa[iT] += calcSA2D(30,82,toys[iT]);
            	ts_sa[iT] += 1.0;

            	std::cout<<ts_sa[iT]<<std::endl;
            }
    	}
		std::cout<<"exit loop"<<std::endl;

//        processTS(ts_nom_sa,ts_sa,hName+"_saturated2D",true);
        if(outName.size())plotter.write(outName);
}

    void processTS(const double nominal_ts,std::vector<double>& toy_ts, const std::string hName, bool isSA){
        std::sort(toy_ts.begin(), toy_ts.end(), [](const double a, const double b){return a < b;});

        double nToys = toy_ts.size();

        auto getQuantile = [&](const double quantile) -> double{
            return toy_ts[int( nToys*quantile  )];
        };

        double quant_5   = getQuantile(0.05);
        double quant_159 = getQuantile(0.15865);
        double quant_50  = getQuantile(0.50);
        double quant_841 = getQuantile(0.84135);
        double quant_95  = getQuantile(0.95);

        double width = (toy_ts.back() - toy_ts.front())/nToys;
        int firstAbove = -1;
        TH1 * h = new TH1D(hName.c_str(),";test statistic",toy_ts.size(), toy_ts[0] -width , toy_ts[toy_ts.size()-1]+width);
        double total     = 0;

        for(unsigned int iT = 0; iT < nToys; ++iT){
            const double ts = toy_ts[iT];
            if(firstAbove <0 && ts >= nominal_ts) firstAbove = iT;
            h->Fill(ts);
            total += ts;
        }

        total /= nToys;
        double pV = firstAbove >= 0 ? 1 - (firstAbove)/nToys :0;
        std::cout << std::endl<< hName <<" -> Data value: "<< nominal_ts <<" ";
        std::cout <<"Toy values mean(5%,15.9%,50%,84.1%,95%) p>x: : "<< total <<"("<<  quant_5<<","<<quant_159<<","<<quant_50<<","<<quant_841<<","<<quant_95<<") "<<pV<<std::endl;

        if(isSA){
            ts_avg_sa  =total;
            ts_up_sa   =quant_841;
            ts_down_sa =quant_159;
        } else {
            ts_avg_ks  =total;
            ts_up_ks   =quant_841;
            ts_down_ks =quant_159;
        }

        plotter.add1D(h);
    }

    TH1D* makeTH1(const std::string& name, const unsigned int nBins, const double * input, const TH1D * refHist){
        TH1D* h((TH1D*)(refHist->Clone(name.c_str())));
        h->Reset("M");
        h->SetDirectory(0);
        for(unsigned int iB = 1; iB + 1  < nBins; ++iB){
            h->SetBinContent(iB,input[iB]);
        }
        return h;
    }

    std::unique_ptr<double[]> generatePseudoData(const unsigned int nBins, const double * input, TRandom* rand ){
        std::unique_ptr<double[]> pd (new double[nBins]() );
        for(unsigned int iB = 1; iB + 1  < nBins; ++iB){
            pd[iB] = rand->Poisson(input[iB]);
        }
        return pd;
    }

    double getIntegral(const unsigned int nBins, const double *input){
        double total = 0;
        for(unsigned int iB =1 ; iB + 1 < nBins; ++iB){
            total += input[iB];
        }
        return total;
    }

    double calcSA(const unsigned int nBins, const ModelAndData& modAData) {return calcSA(nBins,modAData.first->GetArray(),modAData.second->GetArray());}
    double calcSA(const unsigned int nBins, const double *inputModel, const double *inputData){
        double total = 0;
        for(unsigned int iB =1 ; iB + 1  < nBins; ++iB){
            //            std::cout << total <<" :: "<< iB <<"  "<<inputModel[iB] <<"  "<< inputData[iB] <<std::endl;
            total += inputModel[iB];
            if(inputData[iB]) total += inputData[iB]*std::log(inputData[iB]/inputModel[iB]) - inputData[iB];
        }

//        std::cout << total*2<<std::endl;
        return total*2;
    }

    double calcSA2D(const unsigned int nBinsX, const unsigned int nBinsY, const ModelAndData2D& modAData) {return calcSA2D(nBinsX,nBinsY,modAData.first,modAData.second);}
    double calcSA2D(const unsigned int nBinsX, const unsigned int nBinsY, const TH2 *inputModel, const TH2 *inputData){
//        std::cout<<"calculating"<<std::endl;
//        std::cout<<inputModel<<std::endl;
//        std::cout<<inputData<<std::endl;
//        std::cout<<"nBins = "<<nBinsX<<", "<<nBinsY<<std::endl;


    	double total = 0;
        for(unsigned int iBX = 1 ; iBX <= nBinsX; ++iBX) for(unsigned int iBY = 1 ; iBY <= nBinsY; ++iBY){
//            std::cout<<"iBins = "<<iBX<<", "<<iBY<<std::endl;
        	double nModel = inputModel->GetBinContent(iBX,iBY);
//            std::cout<<"nModel = "<<nModel<<std::endl;

            total += nModel;
            double nData = inputData->GetBinContent(iBX,iBY);
//            std::cout<<"nData = "<<nData<<std::endl;
            if(nData)
            	total += nData*std::log(nData/nModel) - nData;
        }

//        std::cout << "total = " << total*2<<std::endl;
        return total*2;
    }

    double calcKS(const unsigned int nBins, const ModelAndData& modAData) {return calcKS(nBins,modAData.first->GetArray(),modAData.second->GetArray());}
    double calcKS(const unsigned int nBins, const double *inputModel, const double *inputData){
        auto modInt = getIntegral(nBins,inputModel);
        auto datInt = getIntegral(nBins,inputData);

        double runningModel = 0;
        double runningData = 0;
        double maxTS = 0;
        for(unsigned int iB =1 ; iB+1 < nBins; ++iB){
            runningModel += inputModel[iB]/modInt;
            runningData += inputData[iB]/datInt;
            double result = std::fabs(runningModel-runningData);
            if(result>maxTS) maxTS = result;
        }
        return maxTS;
    }


    HistGetter plotter;
    double ts_nom_sa =0;
    double ts_nom_ks =0;
    double ts_avg_sa =0;
    double ts_avg_ks =0;
    double ts_up_sa =0;
    double ts_up_ks =0;
    double ts_down_sa =0;
    double ts_down_ks =0;

};

#endif

void StatTester(const TH1D* inputModel, const TH1D* inputData, const int nToys, const bool saveToys, std::string hName, std::string outName = ""){
    StatTesterAnalyzer(inputModel, inputData, nToys, saveToys,hName, outName);
}
