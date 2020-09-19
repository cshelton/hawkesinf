#ifndef MLE_HPP
#define MLE_HPP

#include "hp.hpp"
#include "kernels.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <future>

struct mless {
	mless(int nlabels) : nl(nlabels), Nplus(0),
		Nl(nlabels,0), Nzero(nlabels,0),
		N(nlabels,std::vector<int>(nlabels,0)),
		R(0), TT(0), Ti(nlabels,std::pair<double,double>(0.0,0.0)) {}
	std::size_t nl;
	// total number of events that have non-root parents
	int Nplus;
	// total number of events spawned from each label
	std::vector<int> Nl;
	// total number of events that have root parents, for each label
	std::vector<int> Nzero;
	// N[i][j] is number of events of label j spawned from label i
	std::vector<std::vector<int>> N;
	// sum of times from event to child
	double R;
	// total amount of time
	double TT;
	// for each label, for each event, T - t_i and t_{par} - t_i 
	//  (last is negative if no parent)
	std::vector<std::pair<double,double>> Ti;
	

	template<typename K>
	void addstate(const typename hp<multikernel<K>>::gibbsstate &samp, const multikernel<K> &k) {
		using etype= typename hp<multikernel<K>>::gibbsstate::etype;
		double T = samp.curr.tend;
		TT += T;
		for(auto &e : samp.events) {
			
			if (e.second.e==etype::root || e.second.e==etype::virt) 
				continue;

			Ti[e.first.label].first += k.skernel.intphi(0.0, T-e.first.t);
			Ti[e.first.label].second += e.second.par->second.e!=etype::root ? e.first.t-e.second.par->first.t : -1.0; 

			if (e.second.par->second.e==etype::root) {
				Nzero[e.first.label]++;
			} else {
				R += e.first.t - e.second.par->first.t;
				N[e.second.par->first.label][e.first.label]++;
				Nl[e.second.par->first.label]++;
				Nplus++;
			}
		}
	}
};

// doubling to bracket, following by golden section line search
template<typename F>
double fnmin(F &f, double x, double res=1e-3) {
	double x0 = res, x1 = x, x2 = x*2;
	double f0 = f(x0), f1 = f(x1), f2 = f(x2);
	while(f0<f1 || f1>f2) {
		x0 = x1; f0 = f1;
		x1 = x2; f1 = f2;
		x2 = x2*2; f2 = f(x2);
	}
	constexpr double ratio=0.618;
	double xb = x0 + (1.0-ratio)*(x2-x0);
	double xc = x0 + ratio*(x2-x0);
	double fb = f(xb);
	double fc = f(xc);
	while(x2-x0 > x2*res) {
		if (fb>fc) {
			x0 = xb; f0 = fb;
			xb = xc; fb = fc;
			xc = x0 + ratio*(x2-x0);
			fc = f(xc);
		} else {
			x2 = xc; f2 = fc;
			xc = xb; fc = fb;
			xb = x0 + (1.0-ratio)*(x2-x0);
			fb = f(xb);
		}
	}
	return x0;
}


// for the moment, only for 1-parameter base kernel!!
template<typename K>
void mleopt(multikernel<K> &k, const std::vector<mless> &ss,
		int nitt=1, double minW =1e-3, double maxW = 1e3, double minbeta=1e-3,
		double lambda=0.0, bool clamp=false) {
	//k.skernel.alpha = 1.0;
	//k.skernel.beta = fnmin(negllh,k.skernel.beta,minbeta);
	k.skernel.setparam(0,1.0);

	/*
	auto negllh1 = [&ss](double beta) {
		double ret = ss.Nplus*std::log(beta)-beta*ss.R;
		for(int l=0;l<ss.nl;l++) {
			double sum = 0.0;
			for(auto &dt : ss.Ti[l]) {
				double Gi = std::exp(-beta*dt.first);
				sum += 1-Gi;
			}
			ret -= ss.Nl[l]*std::log(sum);
		}
		return -ret;
	};
	*/

	auto oneval = [&ss,&k](int i, int l) {
		return ss[i].Ti[l];
	};


	auto negllh = [&ss,&k,&oneval,&lambda](double beta) {
		double ret = 0.0;
		k.skernel.setparam(1,beta);
		for(int l=0;l<ss[0].nl;l++) {
			double sum = 0.0;
			int Nl = 0;
			std::vector<std::future<std::pair<double,double>>> futs(ss.size());
			for(int i=0;i<futs.size();i++)
				futs[i] = std::async(std::launch::async,oneval,i,l);
			for(int i=0;i<futs.size();i++) {
				auto res = futs[i].get();
				sum += res.first;
				ret += res.second;
				Nl += ss[i].Nl[l];
			}
			if (Nl>0)
				ret -= Nl*std::log(lambda+sum);
		}
		return -ret;
	};

	double bestbeta = fnmin(negllh,k.skernel.getparam(1),minbeta);
	k.skernel.setparam(1, bestbeta);

	std::cout << "W = " << std::endl;
	for(int l=0;l<ss[0].nl;l++) {
		//std::cout << l << std::endl;
		double den = lambda;
		for(auto &ssi : ss)
			den += ssi.Ti[l].first;
		for(int lp=0;lp<ss[0].nl;lp++) {
			int Nllp = 0;
			for(auto &ssi : ss)
				Nllp += ssi.N[l][lp];
			if (clamp && k.W[l][lp]==0)
				k.W[l][lp] = Nllp/den;
			else
				k.W[l][lp] = std::min(maxW,std::max(minW,Nllp/den));
			if (k.W[l][lp]>minW*1000)
				std::cout << '(' << l << ',' << lp << ") = " << k.W[l][lp] << std::endl;
		}
		//std::cout << std::endl;
	}
	std::cout << "mu = " << std::endl;
	for(int l=0;l<ss[0].nl;l++) {
		size_t TT = 0;
		int Nzero = 0;
		for(auto &ssi : ss) {
			TT += ssi.TT;
			Nzero += ssi.Nzero[l];
		}
		//std::cout << "for mu " << l << " N = " << Nzero << " and T = " << TT << std::endl;
		k.baserates[l] = std::min(maxW,std::max(minW,(double)Nzero/TT));
		std::cout << k.baserates[l] << ' ';
	}
	//std::cout << std::endl << "end beta = " << k.skernel.beta << std::endl;
	k.setWstats();
}

template<typename K, typename R>
void mlestep(const hp<multikernel<K>> &p,
		std::vector<typename hp<multikernel<K>>::gibbsstate> &states,
		mless &ss, int nsamp, int nburnin, int nskip,
		R &rand) {
	for(auto &s : states) {
		for(int i=0;i<nburnin;i++)
			p.gibbsstep(s,rand);
		for(int i=0;i<nsamp;i++) {
			p.gibbsstep(s,rand);
			for(int j=0;j<nskip;j++)
				p.gibbsstep(s,rand);
			ss.addstate<K>(s, p.kernel);
		}
	}
}

template<typename K>
double llh(hp<multikernel<K>> &p, const traj &data) {
	//for every event in data
	//take the log of the mu of that event
	std::vector<std::set<double>::const_iterator> it(data.events.size());
	int cnt;
	for(int i=0;i<data.events.size();++i) {
		it.at(i) = data.events.at(i).begin();
		cnt += data.events.at(i).size();
	}
	std::vector<std::pair<int,double>> order;
	order.reserve(cnt);

	double T = data.tend;
	//std::cout << T << std::endl;

	int ndone = 0;
	double musum = 0, phisum = 0, mutsum = 0, phitsum = 0;
	double oldphisum, newphisum;
	while(ndone<it.size()) {
		//find the soonest event
		int l=-1;
		double t=std::numeric_limits<double>::infinity();
		for(int j=0;j<data.events.size();++j) {
			if(it.at(j) != data.events.at(j).end() 
				and *(it.at(j)) < t) {
				l = j;
				t = *it.at(j);
			}
		}
		if(++it.at(l) == data.events.at(l).end()) ++ndone;

		musum += log(p.kernel.mu(l));	
		
		newphisum = p.kernel.phi(order.back().first, l, order.back().second-t);
		phisum += newphisum;
		oldphisum = newphisum;
/*
		for(const std::pair<int,double>& e : order){
			phisum += p.kernel.phi(e.first, l, e.second-t);
		}
*/

		phitsum += p.kernel.intphi(l, 0, T - t);

		order.emplace_back(l,t);

		if(order.size() % 10000 == 0) {
			std::cout << "seen " << order.size() << " events" << std::endl
				<< "musum " << musum 
				<< "\tphisum " << phisum
				<< "\tphitsum " << phitsum << std::endl;
		}	
	}

	for(int i=0;i<data.events.size();++i) {
		mutsum += p.kernel.mu(i) * T;
	}
	

	return musum + phisum - (mutsum + phitsum);
}


struct mleparams {
	double lambda = 0.0;

	int nsteps=10;
	int nsamp=100;
	int nburnin0=1000;
	int nburnin=0;
	int nskip=0;
	double kappa=2;
	double minW=0.001;
	double maxW=1000;
	double minbeta=0.001;
	int nthread=4;
	int clampWitt=0;
};


template<typename K,typename R>
void mle(hp<multikernel<K>> &p, const std::vector<traj> &data, R &rand,
		const mleparams &params) {
	std::vector<std::vector<typename hp<multikernel<K>>::gibbsstate>>
				states(params.nthread);
	for(int i=0;i<params.nthread;i++)
		for(auto &x : data)
			states[i].emplace_back(p.initgibbs(x,params.kappa,rand));
	std::vector<R> rs;
	for(int i=0;i<params.nthread;i++) rs.emplace_back(rand());

	std::cout << "perfomring mle" << std::endl;
	for(int step=0;step<params.nsteps;step++) {
		std::cout << "step " << step << std::endl;
		std::vector<mless> ss(params.nthread,p.kernel.baserates.size());
		std::vector<std::future<void>> futs(params.nthread);
		
		for(int i=0;i<params.nthread;i++) {
			auto &rr = rs[i];
			auto &si = states[i];
			auto &ssi = ss[i];
			int nburn = step==0 ? params.nburnin0 : params.nburnin;
			int ns = params.nsamp*(i+1)/params.nthread - params.nsamp*i/params.nthread;
			futs[i] = std::async(std::launch::async,
					[&rr,&si,&ssi,&p,nburn,&params,ns]() {
						mlestep(p,si,ssi,ns,nburn,params.nskip,rr);
					});
		}
		for(auto &f : futs)
			f.wait();
		mleopt(p.kernel,ss,1,params.minW,params.maxW,params.minbeta,params.lambda,
				step<params.clampWitt);
		/*
		for(auto &s : states) {
			for(int i=0;i<(step ? params.nburnin : params.nburnin0);i++)
				p.gibbsstep(s,rand);
			for(int i=0;i<params.nsamp;i++) {
				p.gibbsstep(s,rand);
				for(int j=0;j<params.nskip;j++)
					p.gibbsstep(s,rand);
				ss.addstate<K>(s);
			}
		}
		mleopt(p.kernel,ss,1,params.minW,params.maxW,params.minbeta,params.lambda);
		*/

		std::cout << "end iteration " << step << std::endl;
	}
}

#endif
