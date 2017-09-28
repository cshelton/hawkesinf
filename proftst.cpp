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

//using namespace std::string_literals;

using namespace std;

using hprocess = hp<multikernel<singleexpkernel>>;


struct problem {
	traj evid;
	hprocess process;

	problem() :
		process(vector<double>{0.1,0.0001},
               vector<vector<double>>{vector<double>{1/4.0,2/4.0},
                    vector<double>{1/4.0,1/4.0}},
               1,1) {
		evid.tend = 5;

		evid.events.emplace_back();
		evid.events.emplace_back();
		evid.unobs.emplace_back();
		evid.unobs.emplace_back();

		evid.unobs[0].emplace_back(1,3);

		evid.events[1].emplace(4);
	};

	double stat(const traj &tr) const {
		return tr.events[0].size();
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

	issampler(const problem &pr, RAND &rand) :
			p(pr),v(0),wt(0),r(rand) {
	}

	void step() {
		auto samp = p.process.isample(p.evid,r);
		double w = exp(samp.second);
		wt += w;
		v += p.stat(samp.first)*w;
	}

	double estat() const {
		return v/wt;
	}
};

template<typename RAND>
auto makeissampler(const problem &pr, RAND &rand) {
	return issampler<RAND>(pr,rand);
}

//------

template<typename SAMPLER>
double runsampler(SAMPLER s, double time) {
	auto step = [&s](){ s.step(); };
	runfortime(step,time);
	return s.estat();
}

int main(int argc, char **argv) {

	std::random_device rd;
     std::mt19937_64 rand(rd());


	problem pr;

/*
	auto makesampler = [&pr,&rand]() {
				return makegsampler(pr,2,1000,rand); };
*/
	auto makesampler = [&pr,&rand]() {
				return makeissampler(pr,rand); };

	auto sampler = makesampler();
	for(int i=0;i<10000;i++) {
		//cout << "step " << i << endl;
		sampler.step();
	}
}
