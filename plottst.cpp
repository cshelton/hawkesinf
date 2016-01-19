#include "gnuplot-cpp/gnuplot_i.hpp"
#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <numeric>

using namespace std;

void waitforkey() {
	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail());
	std::cin.get();
}

int main(int argc, char **argv) {
	Gnuplot plot("lines");
	vector<int> a{2,3,4};
	vector<int> b{4,1,2};
	plot.plot_xy(a,b);
	waitforkey();
}

