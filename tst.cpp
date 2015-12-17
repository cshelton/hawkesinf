#include "plottraj.hpp"
#include "hp.hpp"
#include "kernels.hpp"
#include <string>
#include <random>

using namespace std::string_literals;

using namespace std;

template<typename T>
using initl = std::initializer_list<T>;

int main(int argc, char **argv) {
/*
	traj tr;

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
		(vector<double>{0.5,0,0},
			vector<vector<double>>{vector<double>{0,1,0},
				vector<double>{0,0,1},
				vector<double>{10,0,5}},
			1,2.5);

	std::random_device rd;
	std::mt19937 rand(rd());
	

	auto tr = process.sample(3,20,rand);

	trajplot plot("expkernel"s);
	plot.plot(tr);

	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail());
	std::cin.get();
}

