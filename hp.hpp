#ifndef HP_HPP
#define HP_HPP

#include <utility>
#include <queue>
#include <random>
#include "traj.hpp"
#include <cmath>
#include <tuple>
#include <map>
#include <cassert>
#include "missingstd.hpp"

//#define DEBUG
//#define DEBUGPCV
#define RANDADV
#define PARENTCHANGEV 

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
							R &rand, double ratemult=1.0) const {
		std::exponential_distribution<> expdist(1);
		return {toi,
			(fromi>=0 && std::isfinite(fromt))
				? fromt+kernel.invintphi(fromi,toi,expdist(rand)/ratemult,currt-fromt)
				: currt+expdist(rand)/(kernel.mu(toi)*ratemult),
			fromi,
			fromt};
	}

	template<typename R>
	traj sample(int neventtypes, double T, R &rand) const {
		std::priority_queue<eventtype,
				std::vector<eventtype>,
				std::greater<eventtype>
			> events;

		auto addevent = [&events,&rand,T,this]
				(int origi, int newi, double origt, double currt) {
					eventtype e = sampleevent(origi,newi,origt,currt,rand);
					if (e.t < T) events.emplace(std::move(e));
				};

		traj ret;
		ret.tend = T;
		for(int i=0;i<neventtypes;i++) ret.events.emplace_back();
		for(int i=0;i<neventtypes;i++) 
			addevent(i,i,std::numeric_limits<double>::infinity(),0);
		while(!events.empty()) {
			eventtype e = events.top();
			events.pop();
			ret.events[e.i].emplace_hint(ret.events[e.i].end(),e.t);
			addevent(e.origi,e.i,e.origt,e.t);
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
		const double T = tr.tend;

		auto timenextunobs = [&unobsi,&tr](int i, double t) { 
			return (unobsi[i]==tr.unobs[i].size()
					? std::numeric_limits<double>::infinity()
					: (t>tr.unobs[i][unobsi[i]].first
						? t
						: tr.unobs[i][unobsi[i]].first)
					);
		};

		auto nexttimeunobs = [&unobsi,&tr](int i, double t, int j) {
			// replace with binary search?
			while(j<tr.unobs[i].size() && t>=tr.unobs[i][j].second) ++j;
			if (j==tr.unobs[i].size())
				 return std::make_pair(std::numeric_limits<double>::infinity(),-1);
			if (t>tr.unobs[i][j].first) return std::make_pair(t,j);
			return std::make_pair(tr.unobs[i][j].first,j);
		};

		auto addevent = [&unobsi,&tr,&events,T,&rand,&nexttimeunobs,&timenextunobs,this]
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
				};
				
		for(int i=0;i<tr.events.size();i++) {
			unobsi.emplace_back(0);
			for(auto &t : tr.events[i])
				events.emplace(i,t,-1,std::numeric_limits<double>::infinity());
			for(auto &u : tr.unobs[i]) {
				events.emplace(i,u.first,-2,std::numeric_limits<double>::infinity());
				events.emplace(i,u.second,-3,std::numeric_limits<double>::infinity());
			}
			addevent(i,i,std::numeric_limits<double>::infinity(),0);
			isunobs.push_back(timenextunobs(i,0)==0);
		}
		double logwt = 0;
		auto state = kernel.basestate(0);
		while(!events.empty()) {
			eventtype e = events.top();
			events.pop();
			auto dwt = kernel.advstate(e.t,state,isunobs,true);
			logwt += dwt;
			if (e.origi==-3) {
				isunobs[e.i]=false;
				unobsi[e.i]++;
				continue;
			} 
			if (e.origi==-2) {
				isunobs[e.i]=true;
				continue;
			}
			if (e.origi==-1) {
				auto ww = kernel.eventrate(e.i,state,true);
				logwt += ww;
			} else {
				addevent(e.origi,e.i,e.origt,e.t);
				tr.events[e.i].emplace_hint(tr.events[e.i].end(),e.t);
			}
			kernel.eventtostate(e.i,state);
			for(auto &i : kernel.fromW(e.i))
				addevent(e.i,i,e.t,e.t);
		}
		auto ddwt = kernel.advstate(T,state,isunobs,true);
		logwt += ddwt;
		return {tr,logwt};
	}

	struct gibbsstate {
		double kappa;
		traj orig,curr;
		enum class etype { norm, virt, evid, root };
		struct eventtime;
		struct eventinfo;
		using eiterator = typename std::map<eventtime,eventinfo>::iterator;
		struct eventtime {
			eventtime(double tt, int ll) : t(tt), label(ll) {}
			double t;
			int label;
			bool operator<(const eventtime &e) const {
				return label==e.label ? t<e.t : label<e.label;
				//return t==e.t ? label<e.label : t<e.t;
			}
		};
		struct eventinfo {
			eventinfo(etype ee, eiterator pp
#ifdef RANDADV
				, int ii
#endif
				) :
				e(ee), par(std::move(pp)), numrealchildren(0)
#ifdef RANDADV
				, i(ii)
#endif
				{}
			etype e;
			eiterator par;
			// perhaps should be a set to make removal faster?
			std::vector<eiterator> vchildren;
			int numrealchildren;
#ifdef RANDADV
			int i;
#endif
		};

		static const char *type2str(etype e) {
			switch(e) {
				case etype::norm: return "N";
				case etype::virt: return "V";
				case etype::evid: return "E";
				case etype::root: return "R";
				default: return "?";
			}
		}

		void print(std::ostream &os) const {
			std::map<eventtime,int> num;
			int i=0;
			for(auto &e : events) num[e.first] = i++;
			i = num[ce->first];
			for(auto &e : events) {
				os << (num[e.first]==i ? '*' : ' ');
				os << '[' << num[e.first] << "] " << e.first.t << '/' << e.first.label << " (type=" << type2str(e.second.e) << ", par=" << num[e.second.par->first] << ", #realch=" << e.second.numrealchildren << ", vchildren={";
				for(int j=0;j<e.second.vchildren.size();j++) {
					if (j>0) os << ',';
					os << num[e.second.vchildren[j]->first];
				}
				os << "})" << std::endl;
			}
			os << "curr:" << std::endl;
			for(int i=0;i<curr.events.size();i++) {
				os << "   " << i << ":";
				for(auto t : curr.events[i]) 
					os << ' ' << t;
				os << std::endl;
			}
			os << std::endl;
		}

		std::map<eventtime,eventinfo> events;
#ifdef RANDADV
		std::vector<eiterator> elist;
#endif
		eiterator ce;
		int nimmobile;

		gibbsstate(traj tr, double kk) : orig(std::move(tr)), kappa(kk) {
			curr.tend = orig.tend;
			curr.events.resize(orig.events.size());
			curr.unobs.resize(orig.unobs.size());
			nimmobile = 1;
			auto r = addevent(0,-1,etype::root,events.end());
			for(int l=0;l<orig.events.size();l++) {
				nimmobile += orig.events[l].size();
				for(auto t : orig.events[l])
					addevent(t,l,etype::evid,r);
			}
			ce = events.begin();
		}

		eiterator addevent(double t, int l, etype e, const eiterator &p) {
			assert(p==events.end() || t>p->first.t);
#ifdef RANDADV
			auto place = events.emplace(eventtime{t,l},eventinfo{e,p,(int)(elist.size())});
			elist.push_back(place.first);
#else
			auto place = events.emplace(eventtime{t,l},eventinfo{e,p});
#endif
			if (e==etype::virt)
				p->second.vchildren.emplace_back(place.first);
			else {
				p->second.numrealchildren++;
				if (e!=etype::root)
					curr.events[l].emplace(t);
			}
			return place.first;
		}
		// virtpar: whether to bother to erase from parent's
		// vchildren list if virtual (in some cases calling fn will do it)
		void delevent(const eiterator &p, bool virtpar=true) {
			if (p->second.e==etype::norm)
				curr.events[p->first.label].erase(p->first.t);
			if (virtpar && p->second.e==etype::virt) {
				auto &vc = p->second.par->second.vchildren;
				for(int i=0;i<vc.size();i++)
					if (vc[i]==p) {
						vc[i] = vc.back();
						vc.resize(vc.size()-1);
						break;
					}
			}
#ifdef RANDADV
			if (p->second.i < elist.size()-1) {
				elist.back()->second.i = p->second.i;
				elist[p->second.i] = elist.back();
			}
			elist.resize(elist.size()-1);
#endif
			if (p->second.e!=etype::virt) {
				assert(p->second.par->second.numrealchildren>0);
				p->second.par->second.numrealchildren--;
			}
			auto np = events.erase(p);
			if (ce==p) {
				ce=np;
				if (ce==events.end()) ce = events.begin();
			}
			
		}

		void makeeventnorm(const eiterator &e) {
			auto &vc = e->second.par->second.vchildren;
			for(int i=0;i<vc.size();i++)
				if (vc[i]==e) {
					vc[i] = vc.back();
					vc.resize(vc.size()-1);
					break;
				}
			e->second.par->second.numrealchildren++;
			e->second.e=etype::norm;
			curr.events[e->first.label].emplace(e->first.t);
		}

		void makeeventvirt(eiterator &e) {
			assert(e->second.numrealchildren==0);
			assert(e->second.par->second.numrealchildren>0);
			for(auto c : e->second.vchildren) delevent(c,false);
			e->second.par->second.numrealchildren--;
			e->second.par->second.vchildren.emplace_back(e);
			e->second.vchildren.clear();
			e->second.e = etype::virt;
			curr.events[e->first.label].erase(e->first.t);
		}

		const traj &trajectory() const { return curr; }

		template<typename R>
		bool advance(R &rand, bool initinc=true) {
#ifdef RANDADV
			std::uniform_int_distribution<> samp(0,elist.size()-1);
			ce = elist[samp(rand)];
#else
			if (initinc) {
				++ce;
			}
			if (ce==events.end()) ce=events.begin();
#endif
			return ce==events.begin();
		}
	};

	template<typename R>
	gibbsstate initgibbs(traj tr, double kappa, R &rand) const {
		return {tr,kappa};
	}

	template<typename R>
	bool gibbsstep(gibbsstate &state, R &rand) const {

		using etype=typename gibbsstate::etype;
		using eiterator=typename gibbsstate::eiterator;
		eiterator ce = state.ce;

		auto rangecmp = [](const std::pair<double,double> &range, double t) {
					return  range.second<t;
				};

		auto nexttimeunobs =
			[&rangecmp](const std::vector<std::pair<double,double>> &unobs, double t,
				std::vector<std::pair<double,double>>::const_iterator j) {
			j = mystd::lower_bound(j,unobs.end(),t,rangecmp);
			if (j==unobs.end())
				 return std::make_pair(std::numeric_limits<double>::infinity(),unobs.end());
			if (t>j->first) return std::make_pair(t,j);
			return std::make_pair(j->first,j);
		};

		auto sampvirtevent = [&state,&rand,&rangecmp,&nexttimeunobs,this]
				(int origi, int newi, double oldt, double currt) {
					auto &unobs = state.orig.unobs[newi];
					std::vector<std::pair<double,double>>::const_iterator j
						= mystd::lower_bound(
							unobs.begin(),
							unobs.end(),
							currt,
							rangecmp);
					if (j==unobs.end()) return state.orig.tend;
					if (currt<j->first) currt = j->first;
					while(1) {
						eventtype e = sampleevent(origi,newi,oldt,currt,
										rand,state.kappa-1.0);
						if (!std::isfinite(e.t)) return state.orig.tend;
						std::tie(currt,j) = nexttimeunobs(unobs,e.t,j);
						if (!std::isfinite(currt)) return state.orig.tend;
						if (currt==e.t && e.t < state.orig.tend)
							return e.t;
					}
				};

		auto sampvirt =
			[&state,&rand,&sampvirtevent,this](double t0,int origi, int newi,
						const eiterator &p) {
				double t = t0;
				std::vector<double> ret;
				while((t=sampvirtevent(origi,newi,t0,t))<state.orig.tend)
					ret.push_back(t);
					//state.addevent(t,newi,etype::virt,p);
				return ret;
			};

		auto resampvchildren1 = [&state,&sampvirt,&rand,this](eiterator ev) {
			std::vector<std::pair<int,double>> vetimes;
			for(auto &l : kernel.fromW(ev->first.label))
				if (!state.orig.unobs[l].empty()) {
					auto ts = sampvirt(ev->first.t,ev->first.label,l,ev);
					for(const auto &t : ts) vetimes.emplace_back(l,t);
				}
			return vetimes;
			};

		auto resampvchildrenrate = [&state](eiterator ev, const std::vector<std::pair<int,double>> &vetimes) {
			int m=vetimes.size();
			int n = state.events.size(); //*2;// - state.nimmobile;
			int deln = m-ev->second.vchildren.size();
			return n/(double)(n+deln);
			};

		auto expratio = [](double lambda, int N, double frac) {
			double ret = 0.0;
			double nfact = 1;
			double expnl = std::exp(-lambda);
			double ln = 1.0;
			double totalpr = 0.0;
			int n = 0;
			while(totalpr<frac) {
				double pr = expnl*ln/nfact;
				ret += N*pr/(N+n);
				totalpr += pr;
				ln *= lambda;
				n++;
				nfact *= n;
			}
			return ret/totalpr;
			};
				

		auto averesampvchildrenrate = [&expratio,&state,&rangecmp,this](eiterator ev) {
				double deln = -ev->second.vchildren.size();
				double t0 = ev->first.t;
				for(int l : kernel.fromW(ev->first.label)) {
					auto &unobs = state.orig.unobs[l];
					for(auto j = mystd::lower_bound(
							unobs.begin(), unobs.end(),
							t0, rangecmp); j!=unobs.end();++j)
						deln += (state.kappa-1)*kernel.intphi(
							ev->first.label, l,
							std::max(j->first-t0,0.0),
							j->second-t0);
				}
				int n = state.events.size(); //*2;// - state.nimmobile;
				//return n/(double)(n+deln);
				return expratio(deln,n,0.999);
			};

		auto unsampvchildrenrate = [&state](eiterator ev) {
			int n = state.events.size();
			int deln = ev->second.vchildren.size();
			return n/(double)(n-deln);
			};

		auto resampvchildren2 = [&state,this](eiterator ev, const std::vector<std::pair<int,double>> &vetimes) {
			for(auto e : ev->second.vchildren)
				state.delevent(e,false);
			ev->second.vchildren.clear();
			for(auto &p : vetimes)
				state.addevent(p.second,p.first,etype::virt,ev);
		};

		auto resampvchildren = [this,&resampvchildren1,&resampvchildrenrate,&resampvchildren2,&rand](eiterator ev) {
			auto vetimes = resampvchildren1(ev);
			double accpr = resampvchildrenrate(ev,vetimes);
			std::uniform_real_distribution<> samp(0.0,1.0);
			if (samp(rand)<=accpr) resampvchildren2(ev,vetimes);
		};

		std::uniform_real_distribution<> sampunif(0,1);
	
		// resample virtualness
		if (ce->second.e==etype::virt) {
			double wvirt = (state.kappa-1);
			double wnorm = exp(-kernel.intphi(ce->first.label,0.0,state.orig.tend-ce->first.t));
			//auto vetimes = resampvchildren1(ce);
			//wnorm *= resampvchildrenrate(ce,vetimes);
			wnorm *= averesampvchildrenrate(ce);
			std::uniform_real_distribution<> samp(0,wvirt+wnorm);
			if (samp(rand)>=wvirt) {
				state.makeeventnorm(ce);
				auto vetimes = resampvchildren1(ce);
				resampvchildren2(ce,resampvchildren1(ce));
			} else
				return state.advance(rand);
		} else if (ce->second.e==etype::norm) {
			if (ce->second.numrealchildren==0) {
				double wnorm = exp(-kernel.intphi(ce->first.label,0.0,state.orig.tend-ce->first.t));
				double wvirt = unsampvchildrenrate(ce)*(state.kappa-1);
				std::uniform_real_distribution<> samp(0,wvirt+wnorm);
				if (samp(rand)<wvirt) {
					state.makeeventvirt(ce);
					return state.advance(rand);
				}
			}
		}
		
		resampvchildren(ce);

#ifdef PARENTCHANGEV
		auto isvirtable = [](eiterator e, eiterator ce) {
			return (e->second.e==etype::virt || (e->second.e==etype::norm
				    && (e->second.numrealchildren==0
				        || (e==ce->second.par && e->second.numrealchildren==1))));
		};
#endif

		// resample parent:
		if (ce->second.e==etype::evid || ce->second.e==etype::norm) {
			std::vector<std::pair<eiterator,double>> poss;
			double wtsum = 0.0;
			eiterator e;

			for(e = state.events.begin();e!=state.events.end() && e->first.label==-1;++e) {
				double wt = kernel.phi(-1,ce->first.label,
									ce->first.t-e->first.t);
				if (wt>0.0) {
					poss.emplace_back(e,wt);
					wtsum += wt;
				}
			}
			for(auto &l : kernel.toW(ce->first.label)) {
				typename gibbsstate::eventtime levent(-1,l);
				e = state.events.lower_bound(levent);
				for(e = state.events.lower_bound(levent);e!=state.events.end() 
						&& e->first.label==l && e->first.t<ce->first.t;++e) {
#ifndef PARENTCHANGEV
					if (e->second.e==etype::virt) continue;
#endif
					double wt = kernel.phi(l,ce->first.label,
										ce->first.t-e->first.t);
					if (wt>0.0) {
						poss.emplace_back(e,wt);
						wtsum += wt;
					}
				}
			}

			std::uniform_real_distribution<> samp(0,wtsum);
			double s = samp(rand);
			for(auto &p : poss) {
				s -= p.second;
				if (s<=0) {
					eiterator &e = p.first;
					eiterator &eold = ce->second.par;
					if (e==eold) break;
#ifdef PARENTCHANGEV
					double acc = 1.0;
					if (!isvirtable(e,ce))
							acc *= state.kappa;
					if (!isvirtable(eold,ce)
								|| (e->second.e==etype::virt
									&& e->second.par == eold))
						acc /= state.kappa;

					double u = 1.0;
					double up = 1.0;
					double Wden = wtsum;
					int r2den = state.events.size();

					bool makev = false;
					std::uniform_real_distribution<> samp(0,1.0);
					std::vector<std::pair<int,double>> vetimes;
					if (e->second.e==etype::virt) {
						up = exp(-kernel.intphi(e->first.label,0.0,
									state.orig.tend-e->first.t));
						vetimes = resampvchildren1(e);
						r2den += vetimes.size();
						for(auto &p : vetimes)
							if (p.second<ce->first.t) 
								Wden += kernel.phi(p.first,ce->first.label,
										ce->first.t-p.second);
					}
					if (eold->second.e==etype::norm
					 && eold->second.numrealchildren==1
					 && (e->second.e!=etype::virt
							|| e->second.par!=eold)
					 && samp(rand)>=1.0/state.kappa) {
						makev = true;
						u = exp(-kernel.intphi(eold->first.label,0.0,
									state.orig.tend-eold->first.t));
						r2den -= eold->second.vchildren.size();
						for(auto c : eold->second.vchildren) 
							if (c->first.t<ce->first.t)
								Wden -= kernel.phi(c->first.label,ce->first.label, ce->first.t-c->first.t);
					}
					acc *= (up/u)*(wtsum/Wden)
						*(state.events.size()/(double)r2den);
					if (acc<1.0 && samp(rand)>=acc) {
						break; // reject change
					}
#endif
#ifdef PARENTCHANGEV
					if (e->second.e==etype::virt) {
						state.makeeventnorm(e);
						resampvchildren2(e,vetimes);
					}
#endif
					eold->second.numrealchildren--;
#ifdef PARENTCHANGEV
					if (makev) state.makeeventvirt(eold);
#endif
					ce->second.par = e;
					e->second.numrealchildren++;
					break;
				}
			}
		}
		return state.advance(rand);
	}
};

#endif
