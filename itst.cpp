#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
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
	
	tr.tend = 10;
	tr.unobs.emplace_back();
	tr.unobs.emplace_back();
	tr.unobs.emplace_back();
	tr.unobs[0].emplace_back(5,8);
	tr.unobs[1].emplace_back(4,9);
	tr.unobs[2].emplace_back(1,5);
	tr.events.emplace_back();
	tr.events.emplace_back();
	tr.events.emplace_back();
	tr.events[0].emplace(0.5);
	tr.events[0].emplace(1.25);
	tr.events[0].emplace(9.5);
	tr.events[1].emplace(9.75);
	tr.events[2].emplace(6.1);
	tr.events[2].emplace(6.3);

	trajplot plot(string("expkernel"));
	plot.plot(tr);

	while(true) {
		std::cin.clear();
		std::cin.ignore(std::cin.rdbuf()->in_avail());
		std::cin.get();
		auto samp = process.isample(tr,rand);
		std::cout << "weight = " << samp.second << std::endl;
		plot.plot(samp.first);
	}
}

