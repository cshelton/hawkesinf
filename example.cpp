#include "hp.hpp"
#include "kernels.hpp"
#include "mle.hpp"
#include "traj.hpp"
#include <random>
#include <iostream>

using namespace std;

std::vector<traj> getdata() {
	// here you would load the data and set the label-time-windows
	// where it was unobserved
	//
	// In this case, we'll just sample from a process and remove 
	// some time windows
	hp<multikernel<singleexpkernel>> process
		(vector<double>{0.5,0.5,0.5}, // mu values -- base rates
		 vector<vector<double>>{
			vector<double>{1.0,0.0,0.5},
			vector<double>{0.0,1.0,0.5},
			vector<double>{1.0,0.5,0.5}}, // wts among labels
		 1.0,2.1 // base kernel parameters: alpha & beta
		);

	std::random_device rd;
	std::mt19937 rand(rd());

	std::vector<traj> ret;
	for (int i=0;i<100;i++) { // sample 10 trajectories
		ret.emplace_back(process.sample(3,20,rd));
			// 3 labels (consistent with process above), 20 time units long
		removewindow(ret.back(),1,5.0,15.0);
			// on each trajectory, erase the events from t=5 to t=15
			// for label 1 (and mark those as unobserved)
			// could be multiple windows, could be different for each traj
	}
	return ret;
}

int main(int argc, char **argv) {
	
	auto data = getdata();
	
	// starting point for MLE (just arbitrarily set at all 1s
	// but could be set randomly)
	hp<multikernel<singleexpkernel>> process
		(vector<double>{1.0,1.0,1.0}, // mu values -- base rates
		 vector<vector<double>>{
			vector<double>{1.0,1.0,1.0},
			vector<double>{1.0,1.0,1.0},
			vector<double>{1.0,1.0,1.0}}, // wts among labels
		 1.0,1.0 // base kernel parameters: alpha & beta
		);

	std::random_device rd;
	std::mt19937 rand(rd());

	// run maximum likelihood estimation (see mle.hpp)
	// currently this only works for a 1-parameter base kernel
	// (beta in the case of exp kernel, b/c alpha is redundant with
	//  the weights and therefore set to 1.0)

	mleparams optparams;
	optparams.nburnin0 = 1000; // burn in first EM iteration
	optparams.nburnin = 100; // burn in after that
	optparams.nsamp = 100; // 100 samples per example
	optparams.nsteps = 100; // run for 100 EM iterations
			// see mle.cpp for other options (search for "struct mleparams")
	mle(process,data,rand,optparams);

	cout << "final parameters:" << endl;
	process.kernel.save(cout); // dump final process parameters
}
				
