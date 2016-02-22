#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <numeric>
#include <time.h>
#include <stdexcept>
#include <signal.h>
#include <unistd.h>
#include <string>
#include <sys/stat.h>
#include <fstream>

#define EXPK
//#define POWK

//using namespace std::string_literals;

using namespace std;

template<typename T>
using initl = std::initializer_list<T>;

void pause(int millisec) {
auto start = std::chrono::high_resolution_clock::now();
	std::this_thread::sleep_for(std::chrono::milliseconds{millisec});
}

void waitforkey() {
	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail());
	std::cin.get();
}

pair<vector<double>,vector<vector<double>>> subsamp(const vector<double> x,
		const vector<vector<double>> y) {
	int nstep = x.size()/1000;
	pair<vector<double>,vector<vector<double>>> ans;
	for(int i=0;i<y.size();i++) ans.second.emplace_back();
	for(int i=0;i<x.size();i+=nstep) {
		ans.first.emplace_back(x[i]);
		for(int j=0;j<y.size();j++)
			ans.second[j].emplace_back(y[j][i]);
	}
	return ans;
}

static void timerhandler(int sig, siginfo_t *si, void *uc) {
	*((int *)(si->si_value.sival_ptr)) = 1;
}

template<typename F>
int runfortime(F f, double tsec) {
	volatile int timerdone = 0;

	struct sigaction sa;
	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = timerhandler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGUSR1,&sa,NULL)==-1) 
		throw runtime_error("sigaction failed");

	timer_t timerid;
	struct sigevent sev;
	sev.sigev_notify = SIGEV_SIGNAL;
	sev.sigev_signo = SIGUSR1;
	sev.sigev_value.sival_ptr = (void *)(&timerdone);
	// or CLOCK_PROCESS_CPUTIME_ID
	if (timer_create(CLOCK_THREAD_CPUTIME_ID,&sev,&timerid) == -1)
		throw runtime_error("timer_create failed");

	struct itimerspec its;
	its.it_interval.tv_sec = 0;	
	its.it_interval.tv_nsec = 0;	
	its.it_value.tv_sec = (int)floor(tsec);
	its.it_value.tv_nsec = (long long)floor((tsec-its.it_value.tv_sec)
										*1000000000);
	//cout << its.it_value.tv_sec << " & " << its.it_value.tv_nsec << endl;

	if (timer_settime(timerid,0,&its,NULL)==-1)
		throw runtime_error("timer_settime failed");

	//cout << "runfortime loc: " << (void *)(&timerdone);

	//kill(getpid(),SIGUSR1);
	struct itimerspec its2;

	int count = 0;
	while(!timerdone) {
		f(); count++;
	}
	return count;
}

#ifdef EXPK
static string kname("exp");
using hprocess = hp<multikernel<singleexpkernel>>;
#endif

#ifdef POWK
static string kname("pow");
using hprocess = hp<multikernel<singlepowerkernel>>;
#endif

struct problem {
	traj evid;
	hprocess process;
	int pnum;

	static hprocess getproblemprocess(int num) {
		switch(num) {
			case 1: return {vector<double>{0.1,0.0001},
				vector<vector<double>>{vector<double>{1/4.0,2/4.0},
					vector<double>{1/4.0,1/4.0}},
#ifdef EXPK
					1,1};
#endif
#ifdef POWK
					1,-2,1};
#endif
			case 2: return {vector<double>{0.0001,0.0001},
				vector<vector<double>>{vector<double>{3/9.0,1/9.0},
					vector<double>{1/9.0,3/9.0}},
#ifdef EXPK
					1,0.5};
#endif
#ifdef POWK
					1,-1.5,1};
#endif
			case 3: return {vector<double>{0.01,0.000001,0.000001},
				vector<vector<double>>{vector<double>{0.0, 1.0, 0.0},
						vector<double>{0.0, 0.0, 1.0},
					vector<double>{0.0, 0.0, 0.0}},
#ifdef EXPK
					1.0,2};
#endif
#ifdef POWK
					1,-2,1};
#endif
			case 4:
			case 5:
				return {vector<double>{0.001,0.000001,0.000001,0.000001,0.000001},
               vector<vector<double>>{vector<double>{0.0, 1.0, 1.0, 0.0, 0.0},
                    vector<double>{0.0, 0.0, 1.0, 1.0, 0.0},
                    vector<double>{0.0, 0.0, 0.0, 1.0, 1.0},
                    vector<double>{0.0, 0.0, 0.0, 0.0, 1.0},
				vector<double>{0.0, 0.0, 0.0, 0.0, 0.0}},
#ifdef EXPK
					1.0,0.5};
#endif
#ifdef POWK
					1.0,-1.5,1};
#endif
		}
	}

	problem(int num) : pnum(num), process(getproblemprocess(num)) {
		switch(num) {
			case 1:
				evid.tend = 5;

				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();

				evid.unobs[0].emplace_back(1,3);

				evid.events[1].emplace(4);
			break;
			case 2:
				evid.tend = 10;

				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();

				evid.unobs[0].emplace_back(1,3);
				evid.unobs[0].emplace_back(6.5,7.5);
				evid.unobs[1].emplace_back(5,7);

				evid.events[0].emplace(6.1);
				evid.events[0].emplace(6.3);
				evid.events[0].emplace(8.0);
				evid.events[0].emplace(8.2);

				evid.events[1].emplace(0.7);
			break;
			case 3:
				evid.tend = 3;

				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();

				evid.unobs[0].emplace_back(0,3);
				evid.unobs[1].emplace_back(0,3);
				evid.unobs[2].emplace_back(0,2);

				evid.events[2].emplace(2.5);
			break;
			case 4:
				evid.tend = 3;

				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.events.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();
				evid.unobs.emplace_back();

				evid.unobs[0].emplace_back(0,3);
				evid.unobs[1].emplace_back(0,3);
				evid.unobs[2].emplace_back(0,3);
				evid.unobs[3].emplace_back(0,3);
				evid.unobs[4].emplace_back(0,2);

				evid.events[4].emplace(2.5);
			break;
			case 5:
#ifdef EXPK
				loaddata("p5data.txt");
#endif
#ifdef POWK
				loaddata("p5data2.txt");
#endif
				evid.unobs[0].emplace_back(0,evid.tend);
				evid.unobs[1].emplace_back(1,evid.tend);
				evid.unobs[3].emplace_back(3,evid.tend);
				evid.events[0].clear();
				evid.events[1].clear();
				evid.events[3].clear();
			break;
		}
	}
	
	void loaddata(const string &fname) {
		ifstream data(fname.c_str());
		if (!data.good())
			throw runtime_error("failed to open "+fname);
		data >> evid.tend;
		int me = -1;
		while(1) {
			int ii,ne;
			data >> ii >> ne;
			if (data.eof()) break;
			if (ii!=me+1) 
				throw runtime_error("bad file format: "+fname);
			evid.events.emplace_back();
			evid.unobs.emplace_back();
			for(int j=0;j<ne;j++) {
				double t;
				data >> t;
				evid.events[ii].emplace(t);
			}
			me = std::max(me,ii);
		}
	}

	int netypes() const {
		return evid.events.size();
	}

	double stat(const traj &tr) const {
		return pnum<4 ? tr.events[0].size() : tr.events[3].size();
	}
};

//------

template<typename RAND>
struct gsampler {
	const problem &p;
	hprocess::gibbsstate s;
	int bin,itt;
	RAND &r;
	double v;
	int c;

	gsampler(const problem &pr, double kappa,
				int burnin, RAND &rand) :
			p(pr),
			s(pr.process.initgibbs(p.evid,kappa,rand)),
			bin(burnin),itt(0),v(0),c(0),
			r(rand) {
	}

	void step() {
		p.process.gibbsstep(s,r);
		itt++;
		if (itt>bin) {
			v += p.stat(s.trajectory());
			c++;
		}
	}

	double estat() const {
		return v/c;
	}
	
	double nsamp() const {
		return c;
	}
};

template<typename RAND>
auto makegsampler(const problem &pr, double kappa,
			int burnin, RAND &rand) {
	return gsampler<RAND>(pr,kappa,burnin,rand);
}

//------

template<typename RAND>
struct issampler {
	const problem &p;
	RAND &r;
	double v;
	double wt;
	double wt2;

	issampler(const problem &pr, RAND &rand) :
			p(pr),v(0),wt(0),r(rand) {
	}

	void step() {
		auto samp = p.process.isample(p.evid,r);
		double w = exp(samp.second);
		wt += w;
		wt2 += w*w;
		v += p.stat(samp.first)*w;
	}

	double estat() const {
		return v/wt;
	}

	double nsamp() const {
		return wt*wt/wt2;
	}
};

template<typename RAND>
auto makeissampler(const problem &pr, RAND &rand) {
	return issampler<RAND>(pr,rand);
}

//------

template<typename SAMPLER>
tuple<double,double,int> runsampler(SAMPLER s, double time) {
	auto step = [&s](){ s.step(); };
	int c = runfortime(step,time);
	return tuple<double,double,int>(s.estat(),s.nsamp(),c);
}

string streamname(int pnum, int burnin, int kappa) {
	string stem = string("data/var-")+kname+"-"+to_string(pnum)+"-"+to_string(burnin)+"-"+to_string(kappa);
	string ext = ".dat";
	struct stat buffer;
	if (stat((stem+ext).c_str(),&buffer)!=0) return stem+ext;
	for(int i=0;i<100;i++) {
		string totry = stem+"+"+to_string(i)+ext;
		if (stat(totry.c_str(),&buffer)!=0) return totry;
	}
	throw runtime_error("failed to create output file");
}

template<typename R>
void sampleproblem(int pnum, double T, R &rand) {
	problem pr(pnum);
	auto tr = pr.process.sample(pr.netypes(),T,rand);
	cout << T << endl;
	for(int i=0;i<tr.events.size();i++) {
		cout << i << ' ' << tr.events[i].size() << endl;
		for(auto &e : tr.events[i]) 
			cout << '\t' << e << endl;
	}
}

int main(int argc, char **argv) {

	std::random_device rd;
     std::mt19937_64 rand(rd());


	int pnum = argc>1 ? atoi(argv[1]) : 1;
	double maxt = argc>2 ? atof(argv[2]) : 10;
	double mint = argc>3 ? atof(argv[3]) : 0.01;
	int nrep = argc>4 ? atoi(argv[4]) : 1000;
	int npts = argc>5 ? atoi(argv[5]) : 20;
	int burnin = argc>6 ? atoi(argv[6]) : 1000;
	double kappa = argc>7 ? atof(argv[7]) : 2;
	int only = argc>8 ? atoi(argv[8]) : 0;

	if (pnum<0) {
		sampleproblem(-pnum,maxt,rand);
		exit(0);
	}

	string sname = streamname(pnum,burnin,kappa);

	problem pr(pnum);

	auto makesamplerA = [&pr,kappa,burnin,&rand]() {
				return makegsampler(pr,kappa,burnin,rand); };
	const char *snameA = "gibbs";

	auto makesamplerB = [&pr,&rand]() {
				return makeissampler(pr,rand); };
	const char *snameB = "IS";

	Gnuplot plot;
	plot.set_title(to_string(pnum)+" -- "+string(snameA)+"/"+snameB+": "+to_string(burnin)+" "+to_string(kappa)+" (itt #0/"+to_string(nrep)+")  "+sname);
	plot.set_style("linespoints");

	vector<double> times;
	for(double l=log(mint); l<log(maxt);l+=(log(maxt)-log(mint))/npts) {
		times.emplace_back(exp(l));
	}
	array<vector<double>,2> valsum,valsum2,nsamps,cs;
	for(int k=0;k<2;k++) {
		valsum[k] = vector<double>(times.size(),0.0);
		valsum2[k] = vector<double>(times.size(),0.0);
		nsamps[k] = vector<double>(times.size(),0.0);
		cs[k] = vector<double>(times.size(),0.0);
	}
	for(int i=0;i<nrep;i++) {
		for(int j=0;j<times.size();j++) {
			double v,nsamp;
			int c;
			if (!only || only==1) {
				std::tie(v,nsamp,c) = runsampler(makesamplerA(),times[j]);
				valsum[0][j] += v;
				valsum2[0][j] += v*v;
				nsamps[0][j] += nsamp;
				cs[0][j] += c;
			}
			if (!only || only==2) {
				std::tie(v,nsamp,c) = runsampler(makesamplerB(),times[j]);
				valsum[1][j] += v;
				valsum2[1][j] += v*v;
				nsamps[1][j] += nsamp;
				cs[1][j] += c;
			}
			cout << '.';
			cout.flush();
		}

		array<vector<double>,2> valvar;
		valvar[0] = vector<double>(times.size(),0);
		valvar[1] = vector<double>(times.size(),0);
		array<double,2> trueest;
		trueest[0] = valsum[0].back()/(i+1);
		cout << " (" << trueest[0] << ") ";
		trueest[1] = valsum[1].back()/(i+1);
		cout << " (" << trueest[1] << ")" << endl;
		for(int k=0;k<2;k++) 
			for(int j=0;j<times.size();j++)
				valvar[k][j] = (valsum2[k][j] - 2*trueest[k]*valsum[k][j])/(i+1) + trueest[k]*trueest[k];

		plot.reset_plot();
		plot.set_title(to_string(pnum)+" -- "+string(snameA)+"/"+snameB+": "+to_string(burnin)+" "+to_string(kappa)+" (itt #"+to_string(i+1)+"/"+to_string(nrep)+")  "+sname);
		plot.set_xlogscale(10);
		plot.set_ylogscale(10);
		plot.plot_xys(times,valvar);

		{
			ofstream outdata(sname.c_str(),std::ios::trunc);
			outdata << "# repnum maxrepnum" << endl;
			outdata << i+1 << ' ' << nrep << endl;
			outdata << "# trueest_gibbs trueest_isamp" << endl;
			for(int k=0;k<2;k++)
				outdata << ' ' << trueest[k];
			outdata << endl;
			outdata << "# runtime var_gibbs var_isamp mean_gibbs mean_isamp mean_nsamp_gibbs mean_effsamp_isamp mean_nitt_gibbs mean_nitt_isamp" << endl;
			for(int j=0;j<times.size();j++) {
				outdata << times[j];
				for(int k=0;k<2;k++)
					outdata << ' ' << valvar[k][j];
				for(int k=0;k<2;k++)
					outdata << ' ' << valsum[k][j]/(i+1);
				for(int k=0;k<2;k++)
					outdata << ' ' << nsamps[k][j]/(i+1);
				for(int k=0;k<2;k++)
					outdata << ' ' << cs[k][j]/(i+1);
				outdata << endl;
			}
		}
	}
	waitforkey();
}