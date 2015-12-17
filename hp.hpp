#ifndef HP_HPP
#define HP_HPP

#include <utility>
#include <queue>
#include <random>
#include "traj.hpp"
#include <cmath>

// A Hawkes process
// type K (the kernel) must implement the following member functions:
//    phi(int i, int j, double t)
//      which is the rate of event j due to 
//      an event of type i happening t time units ago
// and
//    intphi(int i, int j, double t0, double t1)
//       which is the integral of phi(i,j,t) over t from t0 to t1
//       that is
//          if intphi(i,j,t0,t1) = r then
//          \int_{t=t0}^{t1} phi(i,j,t) dt = r;
//       as an example,
//          if phi(i,j,t)=k (that is, phi is a constant),
//          then intphi(i,j,t0,t1) = (t1-t0)*k;
// and
//    invintphi(int i, int j, double s, double t0)
//      which is the inverse of the integral of phi(i,j,t) over t
//      that is,
//          if invintphi(i,j,s,t0) = r, then
//          \int_{t=t0}^{r} phi(i,j,t) dt = s;
//      as an example,
//          if phi(i,j,t)=k (that is, phi is a constant),
//          then invintphi(i,j,s,t0) = s/k+t0;
// and
//    mu(int i)
//      which is the base rate of event i
template<typename K>
struct hp {
	K kernel;

	template<typename... T>
	constexpr hp(T &&...x) : kernel(std::forward<T>(x)...) {}

	struct eventtype {
		int i;
		double t;
		int origi;
		double origt;
		eventtype(int ii, double tt, int oii, double ott)
			: i(ii), t(tt), origi(oii), origt(ott) {}
		bool operator<(const eventtype &e) const { return t<e.t; }
		bool operator<=(const eventtype &e) const { return t<=e.t; }
		bool operator>(const eventtype &e) const { return t>e.t; }
		bool operator>=(const eventtype &e) const { return t>=e.t; }
		bool operator==(const eventtype &e) const { return t==e.t; }
		bool operator!=(const eventtype &e) const { return t!=e.t; }
	};

	template<typename R>
	eventtype sampleevent(int fromi, int toi, double fromt, double currt,
							R &rand) const {
		std::exponential_distribution<> expdist(1);
		eventtype ret;
		ret.origi = fromi;
		ret.origt = fromt;
		ret.i = toi;
		if (std::isfinite(ret.origt)) {
			ret.t = kernel.invintphi(fromi,toi,expdist(rand),currt-fromt);
		else ret.t = currt+expdist(rand)/kernel.mu(toi);
		return ret;
	}

	template<typename R>
	traj sample(int neventtypes, double T, R &rand) const {
		std::priority_queue<eventtype,
				std::vector<eventtype>,
				std::greater<eventtype>
			> events;

		auto addevent = [&events,&rand,T]
				(int origi, int newi, double origt, double currt) {
					eventtype e = sampleevent(origi,newi,origt,currt,rand);
					if (e.t < T) events.emplace(std::move(nne));
				}

		traj ret;
		ret.tend = T;
		for(int i=0;i<neventtypes;i++) ret.events.emplace_back();
		for(int i=0;i<neventtypes;i++) 
			addevent(i,i,std::numeric_limits<double>::infinity(),0);
		while(!events.empty()) {
			eventtype e = events.top();
			events.pop();
			ret.events[e.i].push_back(e.t);
			eventtype ne = sampleevent(e.origi,e.i,e.origt,e.t,rand);
			if (ne.t < T) events.emplace(std::move(ne));
			for(int i=0;i<neventtypes;i++)
				addevent(e.i,i,e.t,e.t);
		}
		return ret;
	}

	template<typename R>
	std::pair<traj,double> isample(traj tr, R &rand) const {
		std::priority_queue<eventtype,
				std::vector<eventtype>,
				std::greater<eventtype>
			> events;
		std::vector<int> unobsi;
		std::vector<bool> isunobs;

		auto timenextunobs = [unobsi,tr](int i, double t) { 
			return (unobsi[i]==tr.unobs[i].size()
						|| t>tr.unobs[i][unobsi[i]].first) ? t
					tr.unobs[i][unobsi[i]].first;
		}

		auto nexttimeunobs = [unobsi,tr](int i, double t, int j) {
			// replace with binary search?
			while(j<tr.unobs[i].size() && t>tr.unobs[i][j].second) ++j;
			if (j==tr.unobs[i].size())
				 return std::make_pair(std::numeric_limits<double>::infinity(),-1);
			if (t>tr.unobs[i][j].first) return std::make_pair(t,j);
			return std::make_pair(tr.unobs[i][j].first,j);
		}

		const double T = tr.tend;
		auto addevent = [unobsi,tr,&events,T,&rand]
				(int origi, int newi, double oldt, double currt) {
					double lowert = timenextunobs(newi,currt);
					int j= unobsi[newi];
					while(1) {
						eventtype e = sampleevent(origi,newi,oldt,lowert,rand);
						std::tie(lowert,j) = nexttimeunobs(newi,e.t,j);
						if (!std::isfinite(e.t)) break;
						if (lowert==e.t) {
							if (e.t < T) events.emplace(std::move(e));
							break;
						}
					}
				}
				
		for(int i=0;i<tr.events.size();i++) {
			unobsi.emplace_back(0);
			for(auto &t : tr.events[i])
				events.emplace(i,t,-1,std::numeric_limits<double>::infinity());
			for(int j=0;j<tr.unobs[i].size();j++) {
				auto &u = tr.unobs[i][j].
				events.emplace(i,u.first,-2,std::numeric_limits<double>::infinity());
				events.emplace(i,u.second,-3,std::numeric_limits<double>::infinity());
			}
			addevent(i,i,std::numeric_limits<double>::infinity(),0);
			isunobs.emplace_back(timenextunobs(i,0)==0);
		}
		double logwt = 0;
		auto state = kernel.basestate(0);
		while(!events.empty()) {
			eventtype e = events.top();
			events.pop();
			logwt += kernel.advstate(e.t,state,isunobs,true)
			if (e.origi==-3) {
				isunobs[e.i]=false;
				unobsi[e.i]++;
				continue;
			} 
			if (e.origi==-2) {
				isunobs[e.i]=true;
				continue;
			}
			if (e.origi==-1) logwt += kernel.eventrate(e.i,state,true);
			else addevent(e.origi,e.i,e.origt,e.t);
			ret.events[e.i].push_back(e.t);
			kernel.eventtostate(e.i,state);
			for(int i=0;i<neventtypes;i++)
				addevent(e.i,i,e.t,e.t);
		}
		return logwt;
	}
};

#endif
