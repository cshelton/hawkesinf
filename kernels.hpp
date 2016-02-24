#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <cmath>
#include <vector>
#include "missingstd.hpp"

struct singleexpkernel {
	double alpha,beta; // phi = alpha*exp(-beta*t)
	constexpr singleexpkernel(double a, double b) : alpha(a), beta(b) {}

	constexpr double phi(double t) const { return alpha*exp(-beta*t); }
	constexpr double intphi(double t0, double t1) const {
		return alpha*exp(-beta*t0)*(-std::expm1(-beta*(t1-t0)))/beta;
	}
	//constexpr double invintphi(double s, double t0) const {
	double invintphi(double s, double t0) const {
		const double lim = alpha/beta/std::exp(beta*t0);
		if (s<lim) {
			//std::cout << "t0 = " << t0 << std::endl;
			//std::cout << "s = " << s << std::endl;
			//std::cout << "lim = " << lim << std::endl;
			//std::cout << "beta = " << beta << std::endl;
			//std::cout << "ans = " << t0-std::log1p(-s/lim)/beta << std::endl;
		}
		return s >= lim
			? std::numeric_limits<double>::infinity()
			: t0-std::log1p(-s/lim)/beta;
	}

	struct state {
		double t;
		double lambda; // perhaps should keep log(lambda)?
	};

	state basestate(double t) const { return {t,0}; }
	double advstate(double t, state &s, bool uselog=true) const {
		double delt = t-s.t;
		double ret = -s.lambda*(-std::expm1(-beta*delt))/beta;
		s.t = t;
		//s.lambda *= std::exp(-delt*s.lambda/beta);
		s.lambda *= std::exp(-delt*beta);
		return uselog ? ret : std::exp(ret);
	}
	void eventtostate(state &s) const {
		s.lambda += alpha;
	}
	double eventrate(const state &s, bool uselog=true) const {
		return uselog ? std::log(s.lambda) : s.lambda;
	}
};

struct singlepowerkernel {
	double alpha,beta,gamma; // phi = alpha*(t+gamma)^beta
		// beta (or b) should be less than -1
	constexpr singlepowerkernel(double a, double b, double g) :
			alpha(a), beta(b), gamma(g) {}

	constexpr double phi(double t) const {
		return alpha*std::pow(t+gamma,beta);
	}

	constexpr double intphi(double t0, double t1) const {
		return alpha*(std::pow(t1+gamma,beta+1)-std::pow(t0+gamma,beta+1))/(beta+1);
	}
	//constexpr double invintphi(double s, double t0) const {
	double invintphi(double s, double t0) const {
		// it seems like this should be doable in a more
		// stable fashion...
		double ret = std::pow((beta+1)*s/alpha + std::pow(gamma+t0,beta+1),
					1.0/(beta+1))
				-gamma;
		return ret>t0 ? ret : std::numeric_limits<double>::infinity();
	}

	struct state {
		double t;
		std::vector<double> ts;
	};

	state basestate(double t) const { return {t,std::vector<double>{}}; }
	double advstate(double t, state &s, bool uselog=true) const {
		double ret = 0.0;
		//for(auto &tt : s.ts) ret -= std::log(intphi(s.t-tt,t-tt));
		for(auto &tt : s.ts) ret -= intphi(s.t-tt,t-tt);
		if (!uselog) ret = std::exp(ret);
		s.t = t;
		return ret;
	}
	double eventtostate(state &s) const {
		s.ts.emplace_back(s.t);
	}
	double eventrate(const state &s, bool uselog=true) const {
		double ret = 0.0;
		for(auto &&tt : s.ts) ret += phi(s.t-tt);
		return uselog ? std::log(ret) : ret;
	}
		
};

template<typename SK>
struct multikernel {
	std::vector<double> baserates;
	double baseratesum;
	// W[i][j] is the multiplier for events from i generating new events in j
	std::vector<std::vector<double>> W,Wtrans;
	std::vector<double> Wsum;
	SK skernel;

	template<typename... T>
	multikernel(std::vector<double> mus,
			std::vector<std::vector<double>> ws,
			T &&...skparams) 
		: baserates(std::move(mus)), W(std::move(ws)), Wsum(W.size(),0),
			Wtrans(W.size(),std::vector<double>(W.size(),0)), skernel(std::forward<T>(skparams)...) {
		baseratesum = 0;
		for(int i=0;i<W.size();i++) {
			double s = 0;
			for(int j=0;j<W.size();j++) {
				Wtrans[j][i] = W[i][j];
				s += W[i][j];
			}
			Wsum[i] = s;
			baseratesum += baserates[i];
		}
	}

	struct toitt {
		std::vector<double>::const_iterator i,e;
		int ind;
		toitt(std::vector<double>::const_iterator b, std::vector<double>::const_iterator en, bool end=false) {
			e=en;
			if (end) {
				i=e;
				ind = 0;
			} else {
				i = b;
				for(ind = 0;i!=e && *i==0;++i,++ind)
					;
			}
		}

		const int &operator*() const { return ind; }
		const int *operator->() const { return &ind; }
		toitt &operator++() { 
			for(++i,++ind;i!=e && *i==0;++i,++ind)
				;
			return *this;
		}
		toitt &operator++(int) { 
			toitt ret(*this);
			++(*this);
			return ret;
		}
		bool operator==(const toitt &itt) const {
			return i==itt.i;
		}
		bool operator!=(const toitt &itt) const {
			return i!=itt.i;
		}
	};

	struct toittrange {
		std::vector<double>::const_iterator b,e;
		toittrange(const multikernel &kernel, int i, bool isto) {
			if (isto) {
				b = kernel.Wtrans[i].begin();
				e = kernel.Wtrans[i].end();
			} else if (i<0) {
				b = kernel.baserates.begin();
				e = kernel.baserates.end();
			} else {
				b = kernel.W[i].begin();
				e = kernel.W[i].end();
			}
		}

		toitt begin() const { return toitt(b,e); }
		toitt end() const { return toitt(b,e,true); }
	};

	constexpr toittrange fromW(int from) const {
		return {*this,from,false};
	}
	
	constexpr toittrange toW(int to) const {
		return {*this,to,true};
	}

	constexpr double phi(int i, int j, double t) const {
		return i<0 ? baserates[j] : W[i][j]*skernel.phi(t);
	}
	// same as above, but summed over all j
	constexpr double phi(int i, double t) const {
		return i<0 ? baseratesum : Wsum[i]*skernel.phi(t);
	}

	constexpr double intphi(int i, int j, double t0, double t1) const {
		return i<0 ? baserates[j]*(t1-t0) : skernel.intphi(t0,t1)*W[i][j];
	}
	// same as above, but summed over all j
	constexpr double intphi(int i, double t0, double t1) const {
		return i<0 ? baseratesum*(t1-t0) : skernel.intphi(t0,t1)*Wsum[i];
	}

	constexpr double invintphi(int i, int j, double s, double t0) const {
		if (i<0) return t0+s/baserates[j];
		if (W[i][j]<=0.0) return std::numeric_limits<double>::infinity();
		return skernel.invintphi(s/W[i][j],t0);
	}

	constexpr double invintphi(int i, double s, double t0) const {
		if (i<0) return t0+s/baseratesum;
		if (Wsum[i]<=0.0) return std::numeric_limits<double>::infinity();
		return skernel.invintphi(s/Wsum[i],t0);
	}

	constexpr double mu(int i) const {
		return baserates[i];
	}

	struct state {
		double t;
		std::vector<typename SK::state> states;
	};

	state basestate(double t) const {
		state ret{};
		for(int i=0;i<baserates.size();i++)
			ret.states.emplace_back(skernel.basestate(t));
		ret.t = t;
		return ret;
	}
	double advstate(double t, state &s, const std::vector<bool> &omit, bool uselog=true) const {
		double ret = 0.0;
		for(int i=0;i<W.size();i++) {
			double mult = 0.0;
			for(int j=0;j<W.size();j++)
				if (!omit[j]) mult+=W[i][j];
			ret += skernel.advstate(t,s.states[i],true)*mult;
		}
		for(int i=0;i<omit.size();i++)
			if (!omit[i]) ret += baserates[i]*(s.t-t);
		s.t = t;
		return uselog ? ret : std::exp(ret);
	}
	void eventtostate(int i, state &s) const {
		skernel.eventtostate(s.states[i]);
	}
	double eventrate(int i, const state &s, bool uselog=true) const {
		double ret = baserates[i];
		for(int j=0;j<W.size();j++)
			if (W[j][i]>0)
				ret += skernel.eventrate(s.states[j],false)*W[j][i];
		//std::cout << "event " << i << " with ";
		//std::cout << "br=" << baserates[i] << " and ";
		//for(auto z : s.states)
		//	std::cout << z.t << ',' << z.lambda << " (" << skernel.eventrate(z,false) << ") ";
		//std::cout << " => " << ret << std::endl;
		return uselog ? std::log(ret) : ret;
	}
};
		
template<typename SK>
struct sparsemultikernel {
	std::vector<double> baserates;
	double baseratesum;
	// W[i][j] is the multiplier for events from i generating new events in j
	std::vector<std::vector<std::pair<int,double>>> W;
	std::vector<std::vector<std::pair<int,double>>> Wtrans;
	std::vector<double> Wsum;
	SK skernel;

	template<typename... T>
	sparsemultikernel(std::vector<double> mus,
			std::vector<std::vector<std::pair<int,double>>> ws,
			T &&...skparams) 
		: baserates(std::move(mus)), W(std::move(ws)), Wtrans(W.size()), skernel(std::forward<T>(skparams)...) {
		baseratesum = 0;
		for(int i=0;i<W.size();i++) {
			double s = 0;
			for(auto &j : W[i]) {
				s += j.second;
				Wtrans[j.first].emplace_back(i,j.second);
			}
			Wsum.emplace_back(s);
			baseratesum += baserates[i];
		}
	}

	template<typename... T>
	sparsemultikernel(std::vector<double> mus,
			const std::vector<std::vector<double>> &ws,
			T &&...skparams) 
		: baserates(std::move(mus)), Wtrans(ws.size()), W(ws.size()),
			skernel(std::forward<T>(skparams)...) {
		baseratesum = 0;
		for(int i=0;i<ws.size();i++) {
			double s = 0;
			for(int j=0;j<ws[i].size();j++)
				if (ws[i][j]>0.0) {
					s += ws[i][j];
					Wtrans[j].emplace_back(i,ws[i][j]);
					W[i].emplace_back(j,ws[i][j]);
				}
			Wsum.emplace_back(s);
			baseratesum += baserates[i];
		}
	}

	struct toitt {
		bool isbase;
		union {
			struct {
				std::vector<std::pair<int,double>>::const_iterator i,e;
			} w;
			struct {
				std::vector<double>::const_iterator i,e;
				int ind;
			} b;
		};
		toitt(std::vector<double>::const_iterator be, std::vector<double>::const_iterator en, bool end=false) {
			isbase = true;
			b.e=en;
			if (end) {
				b.i=b.e;
				b.ind = 0;
			} else {
				b.i = be;
				for(b.ind = 0;b.i!=b.e && *b.i==0;++b.i,++b.ind)
					;
			}
		}
		toitt(std::vector<std::pair<int,double>>::const_iterator b, std::vector<std::pair<int,double>>::const_iterator en, bool end=false) {
			isbase = false;
			w.e = en;
			if (end) {
				w.i=w.e;
			} else {
				w.i=b;
			}
		}

		const int &operator*() const { return isbase ? b.ind : w.i->first; }
		toitt &operator++() { 
			if (isbase) {
				for(++b.i,++b.ind;b.i!=b.e && *b.i==0;++b.i,++b.ind)
					;
			} else {
				++w.i;
			}
			return *this;
		}
		toitt &operator++(int) { 
			toitt ret(*this);
			++(*this);
			return ret;
		}
		bool operator==(const toitt &itt) const {
			return isbase ? (b.i==itt.b.i) : (w.i==itt.w.i);
		}
		bool operator!=(const toitt &itt) const {
			return !(*this==itt);
		}
	};

	struct toittrange {
		bool isbase;
		union {
			struct {
				std::vector<double>::const_iterator b,e;
			} b;
			struct {
				std::vector<std::pair<int,double>>::const_iterator b,e;
			} w;
		};
		toittrange(const sparsemultikernel &kernel, int i, bool isto=false) {
			if (isto) {
				isbase = false;
				w.b = kernel.Wtrans[i].begin();
				w.e = kernel.Wtrans[i].end();
			} else if (i<0) {
				isbase = true;
				b.b = kernel.baserates.begin();
				b.e = kernel.baserates.end();
			} else {
				isbase = false;
				w.b = kernel.W[i].begin();
				w.e = kernel.W[i].end();
			}
		}

		toitt begin() const { return isbase ? toitt(b.b,b.e)
									: toitt(w.b,w.e); }
		toitt end() const { return isbase ? toitt(b.b,b.e,true)
									: toitt(w.b,w.e,true); }
	};

	constexpr toittrange fromW(int from) const {
		return {*this,from,false};
	}

	constexpr toittrange toW(int to) const {
		return {*this,to,true};
	}

	constexpr double getW(int i, int j) const {
		auto l = mystd::lower_bound(W[i].begin(),W[i].end(),j,
				[](const std::pair<int,double> &p, int j) {
					return p.first < j; });
		return (l!=W[i].end() && l->first==j) ? l->second : 0.0;
	}

	constexpr double phi(int i, int j, double t) const {
		return i<0 ? baserates[j] : getW(i,j)*skernel.phi(t);
	}
	// same as above, but summed over all j
	constexpr double phi(int i, double t) const {
		return i<0 ? baseratesum : Wsum[i]*skernel.phi(t);
	}

	constexpr double intphi(int i, int j, double t0, double t1) const {
		return i<0 ? baserates[j]*(t1-t0) : skernel.intphi(t0,t1)*getW(i,j);
	}
	// same as above, but summed over all j
	constexpr double intphi(int i, double t0, double t1) const {
		return i<0 ? baseratesum*(t1-t0) : skernel.intphi(t0,t1)*Wsum[i];
	}

	constexpr double invintphi(int i, int j, double s, double t0) const {
		if (i<0) return t0+s/baserates[j];
		if (getW(i,j)<=0.0) return std::numeric_limits<double>::infinity();
		return skernel.invintphi(s/getW(i,j),t0);
	}

	constexpr double invintphi(int i, double s, double t0) const {
		if (i<0) return t0+s/baseratesum;
		if (Wsum[i]<=0.0) return std::numeric_limits<double>::infinity();
		return skernel.invintphi(s/Wsum[i],t0);
	}

	constexpr double mu(int i) const {
		return baserates[i];
	}

	struct state {
		double t;
		std::vector<typename SK::state> states;
	};

	state basestate(double t) const {
		state ret{};
		for(int i=0;i<baserates.size();i++)
			ret.states.emplace_back(skernel.basestate(t));
		ret.t = t;
		return ret;
	}
	double advstate(double t, state &s, const std::vector<bool> &omit, bool uselog=true) const {
		double ret = 0.0;
		for(int i=0;i<W.size();i++) {
			double mult = 0.0;
			for(auto &j : W[i])
				if (!omit[j.first]) mult+=j.second;
			ret += skernel.advstate(t,s.states[i],true)*mult;
		}
		for(int i=0;i<omit.size();i++)
			if (!omit[i]) ret += baserates[i]*(s.t-t);
		s.t = t;
		return uselog ? ret : std::exp(ret);
	}
	void eventtostate(int i, state &s) const {
		skernel.eventtostate(s.states[i]);
	}
	double eventrate(int i, const state &s, bool uselog=true) const {
		double ret = baserates[i];
		for(auto &j : Wtrans[i])
			ret += skernel.eventrate(s.states[j.first],false)*j.second;
		return uselog ? std::log(ret) : ret;
	}
};
#endif
