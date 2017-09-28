#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include "mle.hpp"
#include <string>
#include <random>

//using namespace std::string_literals;

using namespace std;

template<typename T>
using initl = std::initializer_list<T>;

int main(int argc, char **argv) {
	traj tr;
/*

	tr.events.emplace_back(initl<double>{1.2, 1.3, 3.4, 5.0});
	tr.events.emplace_back(initl<double>{0.0, 2.3, 2.35});
	tr.events.emplace_back(initl<double>{1.0});

	tr.unobs.emplace_back(initl<std::pair<double,double>>{std::pair<double,double>{1.5,3.0}});
	tr.unobs.emplace_back();
	tr.unobs.emplace_back(initl<std::pair<double,double>>{std::pair<double,double>{0.0,0.5},
					std::pair<double,double>{1.5,2.5}});

	trajplot plot("test"s);
	plot.plot(tr);

*/

	hp<multikernel<singleexpkernel>> process
		(vector<double>{0.1,0,0},
			vector<vector<double>>{vector<double>{0.1,2,0.1},
				vector<double>{0.1,0.1,2},
				vector<double>{2,0.1,0.1}},
			1,2.5);

	std::random_device rd;
	std::mt19937 rand(rd());

	std::vector<traj> data;
	for(int i=0;i<1000;i++) {
		data.emplace_back(process.sample(3,10.0,rand));
		removewindow(data.back(),0,1.0,5.0);
		removewindow(data.back(),2,3.0,7.0);
	}

	
	decltype(process) lp
		(vector<double>{0.5,0.5,0.5},
			vector<vector<double>>(3,vector<double>{0.5,0.5,0.5}),
			1,1);

	int nthread = 1;
	if (argc>1)
		nthread = atoi(argv[1]);

	mle(lp,data,rand,10,1000,1000,0,10,2,0.001,0.001,nthread);
}

