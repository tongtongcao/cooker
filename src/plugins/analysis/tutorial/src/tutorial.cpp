#include <tutorial.h>
#include <iostream>
#include <cmath>
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "actree.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
#include <map>
using namespace std;
// Constructor/destructor
tutorial::tutorial(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
	//treeCali=0;
	tree = new TTree();
	//tree = in_;
};

tutorial::~tutorial()
{
};
Long_t tutorial::histos(){
	for(int i = 0; i < 12; i++){
		string name("ACu_ch");
		name+=(i+1)/10+'0';
		name+=(i+1)%10+'0';
		name+="_ac";
		string title("adc Ch for ac u");
		title+=(i+1)/10+'0';
		title+=(i+1)%10+'0';
		title+="Ch_ADC";
		h1adc_ACu[i] = dH1(name.c_str(),title.c_str(),200,0,4500);

		name=string("ACd_ch");
		name+=(i+1)/10+'0';
		name+=(i+1)%10+'0';
		name+="_ac";
		title=string("adc Ch for ac d");
		title+=(i+1)/10+'0';
		title+=(i+1)%10+'0';
		title+="Ch_ADC";
		h1adc_ACd[i] = dH1(name.c_str(),title.c_str(),200,0,4500);

	}
	return 0;
}

Long_t tutorial::startup()
{
	getBranchObject("RawAcInfo",(TObject **) &treeRaw);
	getBranchObject("RawBeamInfo",(TObject **) &treeBeam);

	return 0;
};

Long_t tutorial::process()
{
	for(int i = 0; i < 12; i++){
		h1adc_ACu[i]->Fill(treeRaw->ADC_ACU[i]);
		h1adc_ACd[i]->Fill(treeRaw->ADC_ACD[i]);
	}
	return 0;
};

// Cooker command line option function
Long_t tutorial::cmdline(char *cmd)
{

	return 0;
};

Long_t tutorial::done(){
	TCanvas* c3 = new TCanvas("c3", " AC ADC Upstream ", 100, 100, 800, 650);
	c3->Divide(3,4);
	TCanvas* c4 = new TCanvas("c4", " AC ADC Downstream ", 100, 100, 800, 650);
	c4->Divide(3,4);

	for(int i = 0; i < 12; i++){

    string htitle("AC ADC gap: ");
    htitle+=(i+1)/10+'0';
    htitle+=(i+1)%10+'0';
    htitle+="Gap";

		c3->cd(i+1);
                h1adc_ACu[i]->SetTitle(htitle.c_str());
		h1adc_ACu[i]->SetLineColor(kMagenta-7);
		h1adc_ACu[i]->GetXaxis()->SetTitle("ADC Ch");
		//		gPad->SetLogy();
		h1adc_ACu[i]->Draw("HIST");
		c4->cd(i+1);
                h1adc_ACd[i]->SetTitle(htitle.c_str());
		h1adc_ACd[i]->SetLineColor(kMagenta-7);
		h1adc_ACd[i]->GetXaxis()->SetTitle("ADC Ch");
		//		gPad->SetLogy();
		h1adc_ACd[i]->Draw("HIST");

	}

	c3->Write("AC_adcU");
	c4->Write("AC_adcD");

}
extern "C"
{
	Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
	{
		return (Plugin *) new tutorial(in,out,inf_,outf_,p);
	}
}

ClassImp(tutorial);
