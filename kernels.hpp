#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <cmath>
#include <vector>

struct singleexpkernel {
	double alpha,beta; // phi = alpha*exp(-beta*t)
	constexpr singleexpkernel(double a, double b) : alpha(a), beta(b) {}

	constexpr double phi(double t) const { return alpha*exp(-beta*t); }
	constexpr double intphi(double t0, double t1) const {
		return alpha*(-std::expm1(-beta*(t1-t0)))/beta;
	}
	constexpr double invintphi(double s, double t0) const {
		if (t>= alpha/beta) return std::numeric_limits<double>::infinity();
		return t0-std::log1p(-beta*s/alpha)/beta;
	}

	struct state {
		double t;
		double lambda; // perhaps should keep log(lambda)?
	};

	state basestate(double t) { return {t,0}; }
	double advstate(double t, state &s, bool uselog=true) {
		double delt = t-s.t;
		double ret = -s.lambda*(-std::expm1(-beta*delt))/beta
		s.t = t;
		s.lambda *= std::exp(-delt*s.lambda/beta);
		return uselog ? ret : std::exp(ret);
	}
	void eventtostate(state &s) {
		s.lambda += alpha;
	}
	double eventrate(const state &s, bool uselog=true) {
		return uselog ? std::log(lambda) : lambda;
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
	constexpr double invintphi(double s, double t0) const {
		// it seems like this should be doable in a more
		// stable fashion...
		return std::pow((beta+1)*s/alpha + std::pow(gamma+t0,beta+1),
					1.0/(beta+1))
				-gamma;
	}

	struct state {
		double t;
		std::vector<double> ts;
	};

	state basestate(double t) { return {t,std::vector<double>{}}; }
	double advstate(double t, state &s, bool uselog=true) {
		double ret = 0.0;
		for(auto &tt : s.ts) ret -= std::log(intphi(s.t,t));
		if (!uselog) ret = std::exp(ret);
		s.t = t;
		return ret;
	}
	double eventtostate(state &s) {
		s.ts.emplace_back(s.t);
	}
	double eventrate(const state &s, bool uselog=true) {
		double ret = 0.0;
		for(auto &&tt : s.ts) ret += phi(s.t-tt);
		return uselog ? std::log(ret) : ret;
	}
		
};

template<typename SK>
struct multikernel {
	std::vector<double> baserates;
	// W[i][j] is the multiplier for events from i generating new events in j
	std::vector<std::vector<double>> W;
	SK skernel;

	template<typename... T>
	multikernel(std::vector<double> mus,
			std::vector<std::vector<double>> ws,
			T &&...skparams) 
		: baserates(mus), W(ws), skernel(std::forward<T>(skparams)...) {
	}

	constexpr double phi(int i, int j, double t) const {
		return W[i][j]*skernel.phi(t);
	}

	constexpr double intphi(int i, int j, double t0, double t1) const {
		return skernel.intphi(t0,t1)*W[i][j];
	}

	constexpr double invintphi(int i, int j, double s, double t0) const {
		if (W[i][j]<=0.0) return std::numeric_limits<double>::infinity();
		return skernel.invintphi(s/W[i][j],t0);
	}

	constexpr double mu(int i) const {
		return baserates[i];
	}

	struct state {
		std::vector<typename SK::state> states;
	};

	state basestate(double t) {
		state ret{};
		for(int i=0;i<baserates.size();i++)
			ret.emplace_back(skernel.basestate(t));
		return ret;
	}
	double advstate(double t, state &s, const std::vector<bool> &inc, bool uselog=true) {
		double ret = 0.0;
		for(int i=0;i<W.size();i++) {
			double mult = 0.0;
			for(int j=0;j<W.size();j++)
				if (inc[j]) mult+=W[i][j];
			ret += skernel.advstate(t,s.states[i],true)*mult;
		}
		return uselog ? ret : std::exp(ret);
	}
	void eventtostate(int i, state &s) {
		skernel.eventtostate(s.states[i]);
	}
	double eventrate(int i, const state &s, bool uselog=true) {
		double ret = 1.0;
		for(int j=0;j<W.size();j++)
			if (W[j][i]>0)
				ret += skernel.eventrate(s.states[i],false)*W[j][i];
		return uselog ? std::log(ret) : ret;
	}
};
		
#endif
