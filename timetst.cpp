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
#include <sstream>
#include <algorithm>
#include <fenv.h>

//#define EXPK
#define POWK

#define USESPARSE

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
	if (timer_create(CLOCK_THREAD_CPUTIME_ID,&sev,&timerid) == -1) {
		switch (errno) {
			case EAGAIN: throw runtime_error("timer_create failed: temp error");
			case EINVAL: throw runtime_error("timer_create failed: invalid clock ID");
			case ENOMEM: throw runtime_error("timer_create failed: mem alloc");
			default: throw runtime_error("timer_create failed: unknown");
		}
	}
			

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
	if (count==1)
		cout << "warning: " << tsec << " produced only 1 iteration" << endl;
	return count;
}

#ifdef EXPK
	#ifdef USESPARSE
		static string kname("exp(s)");
		using hprocess = hp<sparsemultikernel<singleexpkernel>>;
	#else
		static string kname("exp");
		using hprocess = hp<multikernel<singleexpkernel>>;
	#endif
#endif

#ifdef POWK
	#ifdef USESPARSE
		static string kname("pow(s)");
		using hprocess = hp<sparsemultikernel<singlepowerkernel>>;
	#else
		static string kname("pow");
		using hprocess = hp<multikernel<singlepowerkernel>>;
	#endif
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
					1,1
#endif
#ifdef POWK
					1,-2,1
#endif
				};
			case 2: return {vector<double>{0.0001,0.0001},
				vector<vector<double>>{vector<double>{3/9.0,1/9.0},
					vector<double>{1/9.0,3/9.0}},
#ifdef EXPK
					1,0.5
#endif
#ifdef POWK
					1,-1.5,1
#endif
				};
			case 3: return {vector<double>{0.01,0.000001,0.000001},
				vector<vector<double>>{vector<double>{0.0, 1.0, 0.0},
						vector<double>{0.0, 0.0, 1.0},
					vector<double>{0.0, 0.0, 0.0}},
#ifdef EXPK
					1.0,2
#endif
#ifdef POWK
					1,-2,1
#endif
				};
			case 4:
			case 5:
				return {vector<double>{0.001,0.000001,0.000001,0.000001,0.000001},
               vector<vector<double>>{vector<double>{0.0, 1.0, 1.0, 0.0, 0.0},
                    vector<double>{0.0, 0.0, 1.0, 1.0, 0.0},
                    vector<double>{0.0, 0.0, 0.0, 1.0, 1.0},
                    vector<double>{0.0, 0.0, 0.0, 0.0, 1.0},
				vector<double>{0.0, 0.0, 0.0, 0.0, 0.0}},
#ifdef EXPK
					1.0,0.5
#endif
#ifdef POWK
					1.0,-1.5,1
#endif
				};
			case 6:
				return loadgraph(100,"graph100",0.05,0.25,0.125);
			case 7:
				return loadgraph(500,"graph500",0.05,0.125,0.0625);
			case 8:
				return loadgraph(500,"graph500",0.1,0.125,0.0625);

		}
	}

	static hprocess loadgraph(int n,string fname, double mu, double selfalpha, double linkalpha) {
		ifstream f(fname.c_str());
		if (!f.good()) throw runtime_error("failed to open "+fname);
#ifdef USESPARSE
		std::vector<std::vector<std::pair<int,double>>> W(n);
#else
		std::vector<std::vector<double>> W(n,std::vector<double>(n,0));
#endif

		while(1) {
			int i;
			f >> i;
			if (f.eof()) break;
			string str;
			getline(f,str);
			istringstream ss(str);
			int j;
			vector<int> ind(1,i);
			while(ss>>j) {
#ifdef USESPARSE
				ind.emplace_back(j);
#else
				W[i][j] = linkalpha;
#endif
			}
#ifdef USESPARSE
			sort(ind.begin(),ind.end());
			for(int j : ind)
				W[i].emplace_back(j,j==i ? selfalpha : linkalpha);
#else
			W[i][i] = selfalpha;
#endif
		}
		return {vector<double>(n,mu),W,
#ifdef EXPK
					1,1
#endif
#ifdef POWK
					1,-2,1
#endif
			};
	}

	problem(int num) : pnum(num), process(getproblemprocess(num)) {
		for(int i=0;i<process.kernel.baserates.size();i++) {
			evid.events.emplace_back();
			evid.unobs.emplace_back();
		}
		switch(num) {
			case 1:
				evid.tend = 5;

				evid.unobs[0].emplace_back(1,3);

				evid.events[1].emplace(4);
			break;
			case 2:
				evid.tend = 10;

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

				evid.unobs[0].emplace_back(0,3);
				evid.unobs[1].emplace_back(0,3);
				evid.unobs[2].emplace_back(0,2);

				evid.events[2].emplace(2.5);
			break;
			case 4:
				evid.tend = 3;

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
			case 6:
				{
#ifdef EXPK
				loaddata("graph100data.txt");
#endif
#ifdef POWK
				loaddata("graph100data2.txt");
#endif
			
				constexpr double mint=0,maxt=10;
				for(int i=0;i<evid.unobs.size();i+=2) {
					evid.unobs[i].emplace_back(mint,maxt);
					set<double> e;
					swap(e,evid.events[i]);
					for(double t : e)
						if (t<mint || t>=maxt) evid.events[i].emplace(t);
				}
				}
			break;
			case 7:
				{
#ifdef EXPK
				loaddata("graph500data.txt");
#endif
#ifdef POWK
				loaddata("graph500data2.txt");
#endif
				constexpr double mint=0,maxt=10;
				for(int i=0;i<evid.unobs.size();i+=2) {
					evid.unobs[i].emplace_back(mint,maxt);
					set<double> e;
					swap(e,evid.events[i]);
					for(double t : e)
						if (t<mint || t>=maxt) evid.events[i].emplace(t);
				}
				}
			break;
			case 8:
				{
#ifdef EXPK
				loaddata("graph500double.txt");
#endif
#ifdef POWK
				loaddata("graph500double2.txt");
#endif
				constexpr double mint=0,maxt=10;
				for(int i=0;i<evid.unobs.size();i+=2) {
					evid.unobs[i].emplace_back(mint,maxt);
					set<double> e;
					swap(e,evid.events[i]);
					for(double t : e)
						if (t<mint || t>=maxt) evid.events[i].emplace(t);
				}
				}
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
			//evid.events.emplace_back();
			//evid.unobs.emplace_back();
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
		switch(pnum) {
			case 0: case 1: case 2: case 3:
					return tr.events[0].size();
			case 4: case 5:
					return tr.events[3].size();
			case 6:
			case 7:
			case 8: {
				int ret = 0;
				for(int i=0;i<tr.events.size();i+=2)
					ret += tr.events[i].size();
				return ret;
				}
			default:
				return 0;
		}
	}
};

struct logsp {
	logsp(double d, double logexp=0.0) : x(d), y(logexp) {}
	operator double() const { return x*exp(y); }
	
	// x * exp(y)
	double x,y;
};

logsp operator+(const logsp &a, const logsp &b) {
	return a.x==0.0 ? b :
		(b.x==0.0 ? a :
		 (a.y<b.y ? logsp(b.x + a.x*exp(a.y-b.y),b.y) :
		  logsp(a.x + b.x*exp(b.y-a.y),a.y)
		 )
		);
}

logsp &operator+=(logsp &a, const logsp &b) {
	return (a = a+b);
}

logsp operator-(const logsp &a, const logsp &b) {
	return a+logsp(-b.x,b.y);
}

logsp operator*(const logsp &a, const logsp &b) {
	return {a.x*b.x,a.y+b.y};
}

logsp operator*(const logsp &a, const double &d) {
	return {a.x*d,a.y};
}

logsp operator/(const logsp &a, const logsp &b) {
	return {a.x/b.x,a.y-b.y};
}

logsp operator/(double a, const logsp &b) {
	return {a/b.x,-b.y};
}
	

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

	double ttl() const {
		return v;
	}
	
	double nsamp() const {
		return c;
	}

	double ttlwt() const {
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
	logsp v;
	logsp wt;
	logsp wt2;

	issampler(const problem &pr, RAND &rand) :
			p(pr),v(0),wt(0),wt2(0),r(rand) {
	}

	void step() {
		auto samp = p.process.isample(p.evid,r);
		logsp w(1,samp.second);
		wt += w;
		wt2 += w*w;
		v += w*p.stat(samp.first);
	}

	logsp estat() const {
		return v/wt;
	}

	logsp ttl() const {
		return v;
	}

	double nsamp() const {
		return wt*wt/wt2;
	}

	logsp ttlwt() const {
		return wt;
	}
};

template<typename RAND>
auto makeissampler(const problem &pr, RAND &rand) {
	return issampler<RAND>(pr,rand);
}

//------

template<typename SAMPLER>
tuple<logsp,double,int,logsp> runsampler(SAMPLER s, double time) {
	auto step = [&s](){ s.step(); };
	int c = runfortime(step,time);
	double num = s.ttl();
	double den = s.ttlwt();
	double ans = s.estat();
	//cout << time << ": " << num << '/' << den << "=" << num/den << " =?= " << ans << endl;
	return tuple<logsp,double,int,logsp>(s.ttl(),s.nsamp(),c,s.ttlwt());
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

//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	std::random_device rd;
     std::mt19937_64 rand(rd());
     //std::ranlux48 rand(rd());


	int pnum = argc>1 ? atoi(argv[1]) : 1;
	double maxt = argc>2 ? atof(argv[2]) : 10;
	double mint = argc>3 ? atof(argv[3]) : 0.01;
	double kappa = argc>4 ? atof(argv[4]) : 2;
	int burnin = argc>5 ? atoi(argv[5]) : 1000;
	int nrep = argc>6 ? atoi(argv[6]) : 100;
	int npts = argc>7 ? atoi(argv[7]) : 16;
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
	double lmint = log(mint);
	double lmaxt = log(maxt);
	double dt = (lmaxt-lmint)/npts;
	for(double l=lmint; l<lmaxt+dt/2;l+=dt) {
		times.emplace_back(exp(l));
	}
	array<vector<double>,2> valsum,valsum2,nsamps,cs;
	array<vector<logsp>,2> wts,sums;
	for(int k=0;k<2;k++) {
		valsum[k] = vector<double>(times.size(),0.0);
		valsum2[k] = vector<double>(times.size(),0.0);
		nsamps[k] = vector<double>(times.size(),0.0);
		cs[k] = vector<double>(times.size(),0.0);
		wts[k] = vector<logsp>(times.size(),0.0);
		sums[k] = vector<logsp>(times.size(),0.0);
	}
	for(int i=0;i<nrep;i++) {
		for(int j=0;j<times.size();j++) {
			logsp ttl(0);
			double nsamp;
			logsp ttlwt(0);
			int c;
			if (!only || only==1) {
				std::tie(ttl,nsamp,c,ttlwt) = runsampler(makesamplerA(),times[j]);
				double v = ttl/ttlwt;
				valsum[0][j] += v;
				valsum2[0][j] += v*v;
				nsamps[0][j] += nsamp;
				cs[0][j] += c;
				sums[0][j] += ttl;
				wts[0][j] += ttlwt;
			}
			if (!only || only==2) {
				std::tie(ttl,nsamp,c,ttlwt) = runsampler(makesamplerB(),times[j]);
				double v = ttl/ttlwt;
				valsum[1][j] += v;
				valsum2[1][j] += v*v;
				nsamps[1][j] += nsamp;
				cs[1][j] += c;
				sums[1][j] += ttl;
				wts[1][j] += ttlwt;
			}
			cout << '.';
			cout.flush();
		}

		array<vector<double>,2> valvar;
		valvar[0] = vector<double>(times.size(),0);
		valvar[1] = vector<double>(times.size(),0);
		array<double,2> trueest;
		for(int k=0;k<2;k++) {
			logsp num(0.0),den(0.0);
			for(int j=0;j<times.size();j++) {
				num += sums[k][j];
				den += wts[k][j];
			}
			trueest[k] = num/den;
		}
		//trueest[0] = valsum[0].back()/(i+1);
		cout << " (" << trueest[0] << ") ";
		//trueest[1] = valsum[1].back()/(i+1);
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
			outdata << "# runtime var_G var_I mean_G mean_I mean_nsamp_G mean_effsamp_I mean_nitt_G mean_nitt_I sum_G sum_I den_G den_I" << endl;
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
				for(int k=0;k<2;k++)
					outdata << ' ' << sums[k][j];
				for(int k=0;k<2;k++)
					outdata << ' ' << wts[k][j];
				for(int k=0;k<2;k++)
					outdata << ' ' << valsum[k][j];
				for(int k=0;k<2;k++)
					outdata << ' ' << valsum2[k][j];
				outdata << endl;
			}
		}
	}
	waitforkey();
}
