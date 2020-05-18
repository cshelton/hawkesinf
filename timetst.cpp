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
#include <cmath>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <fenv.h>
#include <iomanip>
#include <boost/program_options.hpp>

//#define OLDGRAPHLOAD
//using namespace std::string_literals;

using namespace std;
namespace po = boost::program_options;


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
	timer_delete(timerid);
	return count;
}




template<typename K, typename WT>
typename enable_if<is_same<typename K::basekernel,singleexpkernel>::value,hp<K>>::type
makemultiK(vector<double> &&mus, WT &&ws,
				double expa, double expb, 
				double powa, double powb, double powg) {
	return {std::move(mus), std::forward<WT>(ws), expa, expb};
}

template<typename K, typename WT>
typename enable_if<is_same<typename K::basekernel,singlepowerkernel>::value,hp<K>>::type
makemultiK(vector<double> &&mus, WT &&ws,
				double expa, double expb, 
				double powa, double powb, double powg) {
	return {std::move(mus), std::forward<WT>(ws), powa, powb, powg};
}

template<typename K>
typename enable_if<!K::issparse,hp<K>>::type
loadgraph(int n,string fname, double mu, double selfalpha, double linkalpha) {
	ifstream f(fname.c_str());
	if (!f.good()) throw runtime_error("failed to open "+fname);
	std::vector<std::vector<double>> W(n,std::vector<double>(n,0));

	std::vector<int> indegrees(n,0);
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
			W[i][j] = linkalpha;
			indegrees[j]++;
		}
		W[i][i] = selfalpha;
	}
#ifndef OLDGRAPHLOAD
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++) W[i][j] /= indegrees[j];
#endif
	return makemultiK<K>(vector<double>(n,mu),W,
				1,1,
				1,-2,1
		);
}

template<typename K>
typename enable_if<K::issparse,hp<K>>::type
loadgraph(int n,string fname, double mu, double selfalpha, double linkalpha) {
	ifstream f(fname.c_str());
	if (!f.good()) throw runtime_error("failed to open "+fname);
	std::vector<std::vector<std::pair<int,double>>> W(n);

	std::vector<int> indegrees(n,0);
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
			ind.emplace_back(j);
			indegrees[j]++;
		}
		sort(ind.begin(),ind.end());
		for(int j : ind)
			W[i].emplace_back(j,j==i ? selfalpha : linkalpha);
	}
#ifndef OLDGRAPHLOAD
	for(int i=0;i<n;i++)
		for(auto &p : W[i]) p.second /= indegrees[p.first];
#endif
	return makemultiK<K>(vector<double>(n,mu),W,
				1,1,
				1,-2,1
		);
}

template<typename K,typename T>
typename enable_if<is_same<typename K::basekernel, singleexpkernel>::value,T>::type
pickarg(T &&fnameexp, T &&fnamepow) {
	return fnameexp;
}

template<typename K,typename T>
typename enable_if<is_same<typename K::basekernel, singlepowerkernel>::value,T>::type
pickarg(T &&fnameexp, T &&fnamepow) {
	return fnamepow;
}

template<typename K>
struct problem {
	traj evid;
	hp<K> process;
	int pnum;

	static hp<K> getproblemprocess(int num) {
		switch(num) {
			case 1: return makemultiK<K>(vector<double>{0.1,0.0001},
				vector<vector<double>>{vector<double>{1/4.0,2/4.0},
					vector<double>{1/4.0,1/4.0}},
					1,1,
					1,-2,1
				);
			case 2: return makemultiK<K>(vector<double>{0.0001,0.0001},
				vector<vector<double>>{vector<double>{3/9.0,1/9.0},
					vector<double>{1/9.0,3/9.0}},
					1,0.5,
					1,-1.5,1
				);
			case 3: return makemultiK<K>(vector<double>{0.01,0.000001,0.000001},
				vector<vector<double>>{vector<double>{0.0, 1.0, 0.0},
						vector<double>{0.0, 0.0, 1.0},
					vector<double>{0.0, 0.0, 0.0}},
					1.0,2,
					1,-2,1
				);
			case 4:
			case 5:
				return makemultiK<K>(vector<double>{0.001,0.000001,0.000001,0.000001,0.000001},
               vector<vector<double>>{vector<double>{0.0, 1.0, 1.0, 0.0, 0.0},
                    vector<double>{0.0, 0.0, 1.0, 1.0, 0.0},
                    vector<double>{0.0, 0.0, 0.0, 1.0, 1.0},
                    vector<double>{0.0, 0.0, 0.0, 0.0, 1.0},
				vector<double>{0.0, 0.0, 0.0, 0.0, 0.0}},
					1.0,0.5,
					1.0,-1.5,1
				);
			case 12:
				return makemultiK<K>(vector<double>{0.001,0.000001,0.000001,0.000001,0.000001},
               vector<vector<double>>{vector<double>{0.0, 1.0, 1.0, 0.0, 0.0},
                    vector<double>{0.0, 0.0, 1.0, 1.0, 0.0},
                    vector<double>{0.0, 0.0, 0.0, 1.0, 1.0},
                    vector<double>{0.0, 0.0, 0.0, 0.0, 1.0},
				vector<double>{0.0, 0.0, 0.0, 0.0, 0.0}},
					1.0,0.5,
					1.0,-2,1
				);
			case 6:
#ifdef OLDGRAPHLOAD
				return loadgraph<K>(100,"graph100",0.05,0.25,0.125);
#else
				return loadgraph<K>(100,"graph100",0.1,0.5,0.125);
#endif
			case 14:
			case 7:
#ifdef OLDGRAPHLOAD
				return loadgraph<K>(500,"graph500",0.05,0.125,0.0625);
#else
				return loadgraph<K>(500,"graph500",0.1,0.5,0.125);
#endif
			case 8:
#ifdef OLDGRAPHLOAD
				return loadgraph<K>(500,"graph500",0.1,0.125,0.0625);
#else
				return loadgraph<K>(500,"graph500",0.2,0.5,0.25);
#endif
			case 13:
			case 9:
				return loadgraph<K>(100,"graph100",0.1,0.5,0.125);
			case 10:
				return loadgraph<K>(100,"graph100",0.1,0.5,0.125);
			case 11:
				return loadgraph<K>(100,"graph100",0.1,0.5,0.125);

		}
		throw std::runtime_error("invalid problem number");
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
			case 12:
				loaddata("p5data-corr.txt", "p5data2-corr.txt");
				evid.unobs[0].emplace_back(0,evid.tend);
				evid.unobs[1].emplace_back(1,evid.tend);
				evid.unobs[3].emplace_back(3,evid.tend);
				evid.events[0].clear();
				evid.events[1].clear();
				evid.events[3].clear();
			break;
			case 5:
				loaddata("p5data.txt", "p5data2.txt");
				evid.unobs[0].emplace_back(0,evid.tend);
				evid.unobs[1].emplace_back(1,evid.tend);
				evid.unobs[3].emplace_back(3,evid.tend);
				evid.events[0].clear();
				evid.events[1].clear();
				evid.events[3].clear();
			break;
			case 6:
				{
				loaddata("graph100data.txt","graph100data2.txt");
				graphremove(0,10);
				}
			break;
			case 14:
			case 7:
				{
				loaddata("graph500data.txt","graph500data2.txt");
				graphremove(0,10);
				}
			break;
			case 8:
				{
				loaddata("graph500double.txt","graph500double2.txt");
				graphremove(0,10);
				}
			break;
			case 13:
			case 9:
				{
				loaddata("graph100-10.dat","graph100-10-2.dat");
				graphremove(0,10);
			
				}
			break;
			case 10:
				{
				loaddata("graph100-100.dat","graph100-100-2.dat");
				graphremove(0,100);
			
				}
			break;
			case 11:
				{
				loaddata("graph100-1000.dat","graph100-1000-2.dat");
				graphremove(0,1000);
			
				}
			break;
		}
	}

	void graphremove(double mint, double maxt) {
		for(int i=0;i<evid.unobs.size();i+=2) {
			evid.unobs[i].emplace_back(mint,maxt);
			set<double> e;
			swap(e,evid.events[i]);
			for(double t : e)
				if (t<mint || t>=maxt) evid.events[i].emplace(t);
		}
	}
	

	void loaddata(const string &fnameexp, const string &fnamepow) {
		const string &fname = pickarg<K>(fnameexp,fnamepow);
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
			case 4: case 5: case 12:
					return tr.events[3].size();
			case 6: case 7: case 8: case 9: case 10: case 11: {
				int ret = 0;
				for(int i=0;i<tr.events.size();i+=2)
					ret += tr.events[i].size();
				return ret;
				}
			case 13: case 14: {
				int ret = 0;
				for(int i=0;i<tr.events.size();i+=2) {
					vector<decltype(tr.events[i].begin())> oi,oe;
					for(int j : process.kernel.fromW(i)) {
						if (j==i) continue;
						oi.push_back(tr.events[j].begin());
						oe.push_back(tr.events[j].end());
					}
					bool tocount = false;
					for(auto &e : tr.events[i]) {
						for(int j=0;j<oi.size();j++)
							while(oi[j]!=oe[j] && *(oi[j]) < e) {
								oi[j]++;
								tocount = false;
							}
						if (tocount) ret++;
						tocount = true;
					}
				}
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


template<typename RAND, typename K>
struct gsampler {
	const problem<K> &p;
	typename hp<K>::gibbsstate s;
	long long bin,itt,c;
	RAND &r;
	double v;

	gsampler(const problem<K> &pr, double kappa,
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

template<typename RAND,typename K>
auto makegsampler(const problem<K> &pr, double kappa,
			int burnin, RAND &rand) {
	return gsampler<RAND,K>(pr,kappa,burnin,rand);
}

//------

template<typename RAND, typename K>
struct issampler {
	const problem<K> &p;
	RAND &r;
	logsp v;
	logsp wt;
	logsp wt2;

	issampler(const problem<K> &pr, RAND &rand) :
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

template<typename RAND, typename K>
auto makeissampler(const problem<K> &pr, RAND &rand) {
	return issampler<RAND,K>(pr,rand);
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

template<typename K>
string kname();

template<>
string kname<multikernel<singleexpkernel>>() { return {"exp"}; }
template<>
string kname<multikernel<singlepowerkernel>>() { return {"pow"}; }
template<>
string kname<sparsemultikernel<singleexpkernel>>() { return {"exp(s)"}; }
template<>
string kname<sparsemultikernel<singlepowerkernel>>() { return {"pow(s)"}; }

template<typename K>
string streamname(int pnum, int burnin, int kappa) {
	string stem = string("data/var-")+kname<K>()+"-"+to_string(pnum)+"-"+to_string(burnin)+"-"+to_string(kappa);
	string ext = ".dat";
	struct stat buffer;
	if (stat((stem+ext).c_str(),&buffer)!=0) return stem+ext;
	for(int i=0;i<100;i++) {
		string totry = stem+"+"+to_string(i)+ext;
		if (stat(totry.c_str(),&buffer)!=0) return totry;
	}
	throw runtime_error("failed to create output file");
}

template<typename K, typename R>
void sampleproblem(int pnum, double T, R &rand) {
	problem<K> pr(pnum);
	auto tr = pr.process.sample(pr.netypes(),T,rand);
	cout << T << endl;
	for(int i=0;i<tr.events.size();i++) {
		cout << i << ' ' << tr.events[i].size() << endl;
		for(auto &e : tr.events[i]) 
			cout << '\t' << e << endl;
	}
}

struct params {
	int pnum;
	double maxt;
	double mint;
	double kappa;
	int burnin;
	int nrep;
	int npts;
	int only;
	int usesparse;
	int kernel;
	bool sampleonly;
	double value;
	bool plot;
	bool outdata;
};

template<typename K, typename R>
void runit(const params &p, R &rand) {

	string sname = streamname<K>(p.pnum,p.burnin,p.kappa);

	if (p.sampleonly) {
		sampleproblem<K>(p.pnum,p.maxt,rand);
		return;
	}

	problem<K> pr(p.pnum);

	auto makesamplerA = [&pr,p,&rand]() {
				return makegsampler(pr,p.kappa,p.burnin,rand); };
	const char *snameA = "gibbs";

	auto makesamplerB = [&pr,&rand]() {
				return makeissampler(pr,rand); };
	const char *snameB = "IS";

	Gnuplot plot;
	if (p.plot) {
		plot.set_title(to_string(p.pnum)+" -- "+string(snameA)+"/"+snameB+": "+to_string(p.burnin)+" "+to_string(p.kappa)+" (itt #0/"+to_string(p.nrep)+")  "+sname);
		plot.set_style("linespoints");
	}

	vector<double> times;
	double lmint = log(p.mint);
	double lmaxt = p.npts==1 ? lmint*2 : log(p.maxt);
	double dt = p.npts==1 ? lmint*5 : (lmaxt-lmint)/(p.npts-1);
	for(double l=lmint; l<lmaxt+dt/2;l+=dt) {
		times.emplace_back(exp(l));
	}
	double tt = 0;
	for(auto t : times) tt += t*2;
	cout << "total time estimated at " << tt*p.nrep << " seconds" << endl;
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
	for(int i=0;i<p.nrep;i++) {
		for(int j=0;j<times.size();j++) {
			logsp ttl(0);
			double nsamp;
			logsp ttlwt(0);
			int c;
			if (!p.only || p.only==1) {
				std::tie(ttl,nsamp,c,ttlwt) = runsampler(makesamplerA(),times[j]);
				double v = ttl/ttlwt;
				valsum[0][j] += v;
				valsum2[0][j] += v*v;
				nsamps[0][j] += nsamp;
				cs[0][j] += c;
				sums[0][j] += ttl;
				wts[0][j] += ttlwt;
			}
			if (!p.only || p.only==2) {
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
		if (!std::isnan(p.value))
			trueest[0] = trueest[1] = p.value;
		for(int k=0;k<2;k++) 
			for(int j=0;j<times.size();j++)
				valvar[k][j] = (valsum2[k][j] - 2*trueest[k]*valsum[k][j])/(i+1) + trueest[k]*trueest[k];

		if (p.plot) {
			plot.reset_plot();
			plot.set_title(to_string(p.pnum)+" -- "+string(snameA)+"/"+snameB+": "+to_string(p.burnin)+" "+to_string(p.kappa)+" (itt #"+to_string(i+1)+"/"+to_string(p.nrep)+")  "+sname);
			plot.set_xlogscale(10);
			plot.set_ylogscale(10);
			plot.plot_xys(times,valvar);
		}

		if (p.outdata) {
			ofstream outdata(sname.c_str(),std::ios::trunc);
			outdata << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
			outdata << "# repnum maxrepnum" << endl;
			outdata << i+1 << ' ' << p.nrep << endl;
			if (std::isnan(p.value))
				outdata << "# trueest_gibbs trueest_isamp" << endl;
			else 
				outdata << "# true_gibbs true_isamp" << endl;
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
}

int main(int argc, char **argv) {

	long rseed;
	params p;

	po::options_description odesc("options");
	odesc.add_options()
		("help","write help message and exit")
		("show","write parameters and exit")
		("problemnum,p",po::value<int>(&p.pnum)->default_value(1),
				"problem number")
		("maxtime,T",po::value<double>(&p.maxt)->default_value(10.0),
				"maximum runtime (sec)")
		("mintime,t",po::value<double>(&p.mint)->default_value(0.01),
				"minimum runtime (sec)")
		("numpts,n",po::value<int>(&p.npts)->default_value(16),
				"number of runtime points")
		("kappa,k",po::value<double>(&p.kappa)->default_value(2),
				"Gibbs sampler kappa")
		("burnin,b",po::value<int>(&p.burnin)->default_value(1000),
				"number burn-in iterations")
		("numrep,r",po::value<int>(&p.nrep)->default_value(100),
				"number of experiment repetitions")
		("alg,a",po::value<int>(&p.only)->default_value(0),
				"algorithm to be tested (0=both, 1=Gibbs, 2=IS)")
		("kernel",po::value<int>(&p.kernel)->default_value(1),
				"kernel type (1=exp, 2=pow)")
		("sparse,s",po::value<int>(&p.usesparse)->default_value(true),
				"sparse kernel? (1=yes, 0=no)")
		("randseed",po::value<long>(&rseed)->default_value(-1),
				"random seed (<0 => randomly selected)")
		("sampleonly",po::value<bool>(&p.sampleonly)->default_value(false),
				"true if only to sample trajectory")
		("config",po::value<string>()->default_value(string("")),
				"configuration file name")
		("truevalue",po::value<double>(&p.value)->default_value(std::numeric_limits<double>::quiet_NaN()),
				"true value of calculation")
		("nopause","do not wait for keystroke when finished")
		("noplot","do not plot with gnuplot")
		("nowritedata","do not write to data log")
		; 

	po::variables_map omap;
	po::store(po::parse_command_line(argc, argv, odesc),omap);
	string cfname = omap["config"].as<string>();
	if (cfname!="") {
		cout << "parsing " << cfname << endl;
		po::store(po::parse_config_file<char>(cfname.c_str(),odesc),omap);
	} else {
		cout << "cfname = " << cfname << endl;
	}
	po::notify(omap);
	p.plot = !omap.count("noplot");
	p.outdata = !omap.count("nowritedata");
	if (omap.count("help")) {
		cout << odesc << endl;
		return 1;
	}
	if (omap.count("show")) {
		cout << "pnum=" << p.pnum << endl
			<< "maxt=" << p.maxt << endl
			<< "mint=" << p.mint << endl
			<< "npts=" << p.npts << endl
			<< "kappa=" << p.kappa << endl
			<< "burnin=" << p.burnin << endl
			<< "nrep=" << p.nrep << endl
			<< "only=" << p.only << endl
			<< "kernel=" << p.kernel << endl
			<< "sparse=" << p.usesparse << endl
			<< "randseed=" << rseed << endl
			<< "sampleonly=" << p.sampleonly << endl
			<< "config=" << cfname << endl;
		return 1;
	}
//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);


	std::random_device rd;
     std::mt19937_64 rand(rseed<0 ? rd() : rseed);
     //std::ranlux48 rand(rd());

	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

	switch(p.kernel) {
		case 1: if (p.usesparse) runit<sparsemultikernel<singleexpkernel>>(p,rand);
			else runit<multikernel<singleexpkernel>>(p,rand);
		break;
		case 2: if (p.usesparse) runit<sparsemultikernel<singlepowerkernel>>(p,rand);
			else runit<multikernel<singlepowerkernel>>(p,rand);
		break;
		default:
			cout << "illegal kernel parameters" << endl;
	}
	if (!omap.count("nopause")) waitforkey();
}

