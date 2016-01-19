#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <numeric>

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

int main(int argc, char **argv) {
	traj tr;
	hp<multikernel<singleexpkernel>> process
		(vector<double>{0.1,0.0001},
			vector<vector<double>>{vector<double>{1/4.0,2/4.0},
				vector<double>{1/4.0,1/4.0}},
			1,1);

	std::random_device rd;
	std::mt19937_64 rand(rd());
	
	tr.tend = 5;
	tr.unobs.emplace_back();
	tr.unobs.emplace_back();
	tr.unobs[0].emplace_back(1,3);
	tr.events.emplace_back();
	tr.events.emplace_back();
	tr.events[1].emplace(4);


	int n = 1000000;
	int m = 1000;
	int nval = 2;
	int npts = 1000;
	if (argc>1) n = atoi(argv[1]);
	if (argc>2) m = atoi(argv[2]);

	vector<double> xs(n,0);
	vector<vector<double>> ys;
	vector<vector<double>> y2s;
	vector<Gnuplot> plots(2);
	for(int i=0;i<nval;i++) {
		ys.emplace_back(n,0);
		y2s.emplace_back(n,0);
		plots[i].set_style("linespoints");
	}
	iota(xs.begin(),xs.end(),1);
	vector<double> ssx(npts,0),ssy(npts,0);
	
	for(int j=0;j<m;j++) {
		vector<double> s(nval,0.0);
		vector<double> c(nval,0.0);
		auto state = process.initgibbs(tr,2.0,rand);
		for(int i=0;i<n;i++) {
			while(!process.gibbsstep(state,rand))
				;
			auto samp = state.trajectory();
			c[0] = 0;
			c[1] = 0;
			for(auto t : samp.events[0])
				if (t>=1 && t<2) c[0]++;
				else if (t>=2 && t<3) c[1]++;

			for(int k=0;k<nval;k++) {
				s[k] += c[k];
				ys[k][i]+=s[k]/(i+1);
				y2s[k][i]+=(s[k]*s[k])/(i+1);
			}
		}
		cout << j << '/' << m << ":";
		for(int k=0;k<nval;k++) cout << ' ' << s[k]/n;
		cout << endl;
		if (j%10==0) {
		//{
			int nm = j+1;
			for(int k=0;k<nval;k++) {
				double val = ys[k].back()/nm;
				plots[k].reset_plot();
				int ii = 0;
				for(double l=log(1);l<log(n);l+=(log(n)-log(1))/npts) {
					int li = floor(exp(l));
					ssx[ii] = xs[li];
					ssy[ii] = (y2s[k][li] - 2*val*ys[k][li])/nm + val*val;
					++ii;
				}
				plots[k].plot_xy(ssx,ssy);
				plots[k].set_xlogscale(10);
				plots[k].set_ylogscale(10);
				pause(1);
			}
		}
	}
	waitforkey();
}

