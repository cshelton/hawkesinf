#include "mle.hpp"
#include "kernels.hpp"
#include "hp.hpp"
#include <random>
#include <limits>
#include <iostream>

template<typename K, typename R>
void testpred(const hp<multikernel<K>> &p, traj data,
		std::map<double,int> toadd, const mleparams &params,
		int nobs, std::ostream &os, R &rand) {

	data.unobs.resize(p.kernel.baserates.size());
	data.events.resize(data.unobs.size());

	std::vector<R> rs;
	for(int i=0;i<params.nthread;i++) rs.emplace_back(rand());

	std::cout << "looping for " << toadd.size() << " new events" << std::endl;
	for(auto &next : toadd) {

		double t0 = next.first;
		double wint = data.tend+(t0-data.tend+1)*10;
		std::cout << "predicting from " << data.tend << " to " << t0 << " for event of label " << next.second << " with window until " << wint << std::endl;
		for(int i=0;i<nobs;i++) {
			data.unobs[i].clear();
			data.unobs[i].emplace_back(data.tend,wint);
		}
		for(int i=nobs;i<p.kernel.baserates.size();i++) {
			data.unobs[i].clear();
			data.unobs[i].emplace_back(0,wint);
		}
		data.tend = wint;

		std::vector<typename hp<multikernel<K>>::gibbsstate> states;
		for(int i=0;i<params.nthread;i++)
			states.emplace_back(p.initgibbs(data,params.kappa,rand));

		std::vector<std::vector<int>> ans;
		for(int i=0;i<params.nthread;i++)
			ans.emplace_back(nobs,0);

		std::vector<std::future<void>> futs(params.nthread);
		for(int i=0;i<params.nthread;i++) {
			auto &rr = rs[i];
			auto &si = states[i];
			auto &a = ans[i];
			int ns = params.nsamp*(i+1)/params.nthread
					- params.nsamp*i/params.nthread;
			int nburn = params.nburnin0;
			futs[i] = std::async(std::launch::async,
					[&rr,&si,&p,&a,t0,nburn,ns]() {
						pred1(p,si,a,t0,nburn,ns,rr);
						});
		}
		for(auto &f : futs)
			f.wait();

		os << t0 << ' ' << next.second;

		int tot = 0;
		for(int i=0;i<nobs;i++)
			for(auto &v : ans) tot += v[i];

		int sright = 0;
		for(auto &v : ans) sright += v[next.second];

		os << ' ' << (double)sright/tot;

		os << ' ' << sright << ' ' << tot;

		for(int i=0;i<nobs;i++) {
			int s = 0;
			for(auto &v : ans) s += v[i];
			os << ' ' << s;
		}
		os << std::endl;

		data.events[next.second].emplace(t0);
		data.tend = next.first;
	}
}

template<typename K, typename R>
void pred1(const hp<multikernel<K>> &p,
		typename hp<multikernel<K>>::gibbsstate &state,
		std::vector<int> &ansdist, double t0, 
		int nburn, int nsamp, R &rand) {
	for(int i=0;i<nburn;i++)
		p.gibbsstep(state,rand);
	for(int i=0;i<nsamp;i++) {
		p.gibbsstep(state,rand);
		int lnext=-1;
		double lt = std::numeric_limits<double>::infinity();
		for(int l=0;l<ansdist.size();l++) {
			auto lb = state.curr.events[l].upper_bound(t0);
			if (lb != state.curr.events[l].end() 
					&& *lb<lt) {
				lnext = l;
				lt = *lb;
			}
		}
		if (lnext!=-1) ansdist[lnext]++;
	}
}
