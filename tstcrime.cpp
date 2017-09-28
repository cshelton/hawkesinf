#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include "mle.hpp"
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <boost/program_options.hpp>
#include "tstpred.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
	traj tr;

	mleparams mlep;

	double initbeta = 1.0/90.0;
	double mumult = 0.1;
	double maxinitW = initbeta/2;
	double maxT = std::numeric_limits<double>::infinity();
	double tstT = 10956.0;
	double minT = tstT-365.0;
	int allinitW = 1; // actually boolean
	unsigned int initseed = 0;
	string inputfilename = "crimedata/h5-l5-c50-0";
	string outputfilename = "crimedata/results5.txt";

	po::options_description odesc("options");
	odesc.add_options()
		("help","write help message and exit")
		("show","write parameters and exit")
		("maxT",po::value<double>(&maxT)->default_value(std::numeric_limits<double>::infinity()),"max ending time")
		("minT",po::value<double>(&minT)->default_value(10956-365.0),"starting time")
		("tstT",po::value<double>(&tstT)->default_value(10956.0),"start testing time")
		("nsamp",po::value<int>(&mlep.nsamp)->default_value(10000000),
			"number of MC sampler / pred")
		("nburnin0",po::value<int>(&mlep.nburnin0)->default_value(1000000),
			"number of initial burnin steps (per thread)")
		("nburnin",po::value<int>(&mlep.nburnin)->default_value(1000),
			"number of burnin steps (per thread) after 1st pred")
		("kappa",po::value<double>(&mlep.kappa)->default_value(5.0),
			"kappa")
		("nthread",po::value<int>(&mlep.nthread)->default_value(8),
			"number of threads")
		("in",po::value<string>(&inputfilename)->default_value(string("crimedata/h5-l5-c50-0")),
			"input model filename")
		("out",po::value<string>(&outputfilename)->default_value(string("crimedata/results5.txt")),
			"output results filename")
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
		cout << "maxT=" << maxT << endl
			<< "minT=" << minT << endl
			<< "tstT=" << tstT << endl
			<< "nsamp=" << mlep.nsamp << endl
			<< "nburnin0=" << mlep.nburnin0 << endl
			<< "nburnin=" << mlep.nburnin << endl
			<< "kappa=" << mlep.kappa << endl
			<< "nthread=" << mlep.nthread << endl
			<< "in=" << inputfilename << endl
			<< "out=" << outputfilename << endl
			<< "config=" << cfname << endl;
		return 1;
	}

	int nlabels = 77;

	traj data;

	ifstream cdata("icpsr.txt");
	if (!cdata.good()) {
		cout << "could not open icpsr.txt" << endl;
		return 1;
	}

	double T;
	cdata >> T;
	T = std::min(std::min(T,maxT),tstT);
	T -= minT;
	tstT -= minT;
	data.tend = T;
	data.events.resize(nlabels);
	data.unobs.resize(nlabels);
	map<double,int> toadd;
	int nev = 0, nevtst=0;
	
	while(true) {
		int l;
		double t;
		cdata >> l >> t;
		t -= minT;
		if (cdata.eof() || t>maxT) break;
		if (t<0.0) continue;
		if (t<tstT) {
			data.events[l].insert(t);
			nev++;
		} else {
			toadd[t] = l;
			nevtst++;
		}
	}
	std::cout << "read " << nev << " events and " << nevtst << " after for testing " << std::endl;

	std::random_device rd;
	if (initseed == 0) initseed = rd();
	std::mt19937 rand(initseed);

	ifstream ifs;
	ifs.open(inputfilename.c_str());
	if (!ifs.good()) {
		cout << "could not open " << inputfilename << " for reading" << endl;
		return 2;
	}
	multikernel<singleexpkernel> k
			(vector<double>(),vector<vector<double>>(),1,1.0);
	string str;
	getline(ifs,str);
	getline(ifs,str);
	k.load(ifs);

	hp<multikernel<singleexpkernel>> process(k);

	ofstream ofs(outputfilename.c_str());

	testpred(process,data,toadd,mlep,77,ofs,rand);

}

