#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include "mle.hpp"
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <limits>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
	traj tr;

	mleparams mlep;

	int extralabels = 5;
	double initbeta = 1.0/90.0;
	double mumult = 0.1;
	double maxinitW = initbeta/2;
	double maxT = std::numeric_limits<double>::infinity();
	double minT = 0.0;
	int allinitW = 1; // actually boolean
	unsigned int initseed = 0;
	string outputfilename = "crimedata/testout.txt";

	po::options_description odesc("options");
	odesc.add_options()
		("help","write help message and exit")
		("show","write parameters and exit")
		("maxT",po::value<double>(&maxT)->default_value(std::numeric_limits<double>::infinity()),"max ending time")
		("minT",po::value<double>(&minT)->default_value(0.0),"min starting time")
		("nhidden",po::value<int>(&extralabels)->default_value(5),
			"number hidden labels")
		("nem",po::value<int>(&mlep.nsteps)->default_value(100),
			"number of EM iterations")
		("nsamp",po::value<int>(&mlep.nsamp)->default_value(50000),
			"number of MC sampler / EM iteration")
		("nburnin0",po::value<int>(&mlep.nburnin0)->default_value(100000),
			"number of initial burnin steps (per thread)")
		("nburnin",po::value<int>(&mlep.nburnin)->default_value(0),
			"number of burnin steps (per thread) after 1st EM iteration")
		("nskip",po::value<int>(&mlep.nskip)->default_value(100),
			"number of samples to skip before recording")
		("kappa",po::value<double>(&mlep.kappa)->default_value(5.0),
			"kappa")
		("minW",po::value<double>(&mlep.minW)->default_value(1e-6),
			"min value for any element of W")
		("maxW",po::value<double>(&mlep.maxW)->default_value(10),
			"min value for any element of W")
		("minbeta",po::value<double>(&mlep.minbeta)->default_value(1e-4),
			"min value for beta")
		("nthread",po::value<int>(&mlep.nthread)->default_value(8),
			"number of threads")
		("lambda",po::value<double>(&mlep.lambda)->default_value(0.0),
			"lambda (L1-regularizing strength on Ws)")
		("initbeta",po::value<double>(&initbeta)->default_value(1.0/90.0),
			"starting beta value")
		("mumult",po::value<double>(&mumult)->default_value(0.1),
			"multiplier for starting mu value (multiplied by emperial rate for label)") 
		("maxinitW",po::value<double>(&maxinitW)->default_value(initbeta/2.0),
			"maximum initial W")
		("allinitW",po::value<int>(&allinitW)->default_value(false),
			"whether to set non-zero Ws between observed labels")
		("clampWitt",po::value<int>(&mlep.clampWitt)->default_value(0),
			"if not allinitW, number of iterations to keep clamped")
		("out",po::value<string>(&outputfilename)->default_value(string("crimedata/testout.txt")),
			"output model filename")
		("config",po::value<string>()->default_value(string("")),"configuration filename")
		("initseed",po::value<unsigned int>(&initseed)->default_value(0),"initial seed (0=use random device)");

	po::variables_map omap;
	po::store(po::parse_command_line(argc,argv,odesc),omap);
	string cfname = omap["config"].as<string>();
	if (cfname!="") {
		cout << "parsing " << cfname << endl;
		po::store(po::parse_config_file<char>(cfname.c_str(),odesc),omap);
	}
	po::notify(omap);
	if (omap.count("help")) {
		cout << odesc << endl;
		return 1;
	}
	if (omap.count("show")) {
		cout << "nhidden=" << extralabels << endl
			<< "nem=" << mlep.nsteps << endl
			<< "maxT=" << maxT << endl
			<< "minT=" << minT << endl
			<< "nsamp=" << mlep.nsamp << endl
			<< "nburnin0=" << mlep.nburnin0 << endl
			<< "nburnin=" << mlep.nburnin << endl
			<< "nskip=" << mlep.nskip << endl
			<< "kappa=" << mlep.kappa << endl
			<< "minW=" << mlep.minW << endl
			<< "maxW=" << mlep.maxW << endl
			<< "minbeta=" << mlep.minbeta << endl
			<< "nthread=" << mlep.nthread << endl
			<< "lambda =" << mlep.lambda << endl
			<< "initbeta=" << initbeta << endl
			<< "mumult=" << mumult << endl
			<< "maxinitW=" << maxinitW << endl
			<< "allinitW=" << allinitW << endl
			<< "clampWitt=" << mlep.clampWitt << endl
			<< "out=" << outputfilename << endl
			<< "config=" << cfname << endl;
		return 1;
	}

	int nlabels = 77;

	std::vector<traj> data;

	ifstream cdata("icpsr.txt");
	if (!cdata.good()) {
		cout << "could not open icpsr.txt" << endl;
		return 1;
	}

	data.push_back(traj{});
	double T;
	cdata >> T;
	T = std::min(T,maxT);
	T -= minT;
	data.back().tend = T;
	data.back().events.resize(nlabels);
	data.back().unobs.resize(nlabels);
	vector<int> nevents(nlabels,0);
	int nev = 0;
	while(true) {
		int l;
		double t;
		cdata >> l >> t;
		t -= minT;
		if (cdata.eof() || t>T) break;
		if (t<0.0) continue;
		data.back().events[l].insert(t);
		nevents[l]++;
		nev++;
	}
	std::cout << "read " << nev << " events" << std::endl;
	vector<double> initmu;
	int maxc = *(max_element(nevents.begin(),nevents.end()));
	for(auto &c : nevents) initmu.push_back(c/T * mumult);


	std::random_device rd;
	if (initseed == 0) initseed = rd();
	std::mt19937 rand(initseed);
	std::uniform_real_distribution<> wdist(0.0,maxinitW);

	vector<pair<double,double>> noobs(1,pair<double,double>{0.0,T});
	for(int i=0;i<extralabels;i++) {
		initmu.push_back(maxc/T * 0.1);
		data.back().events.emplace_back();
		data.back().unobs.emplace_back(noobs);
	}

	vector<vector<double>> initW(nlabels+extralabels,vector<double>(nlabels+extralabels,0.0));
	for(int i=0;i<nlabels+extralabels;i++) 
		for(int j=0;j<nlabels+extralabels;j++)
			initW[i][j] = !allinitW && i<nlabels && j<nlabels ? 0 : wdist(rand);

	hp<multikernel<singleexpkernel>> process(initmu,initW,1,initbeta);

	ofstream ofs;
	if (outputfilename!="") {
		ofs.open(outputfilename.c_str());
		if (!ofs.good()) {
			cout << "could not open " << outputfilename << " for writing" << endl;
			return 2;
		}
	}

	mle(process,data,rand,mlep);

	if (outputfilename!="") {
		ofs << "configfile = " << cfname << endl;
		ofs << "initseed = " << initseed << endl;
		process.kernel.save(ofs);
	}
}

