#ifndef TRAJ_HPP
#define TRAJ_HPP

#include <vector>
#include <set>

struct traj {
	double tend; // t0 always =0
	// events[i] are the times of events of type i
	std::vector<std::set<double>> events;

	// unobs[i][j] is the jth interval that is not observed
	//  for event type i (.first is the time of start of the interval
	//  and .second is the time of the end of the interval)
	std::vector<std::vector<std::pair<double,double>>> unobs;
};
	
#endif
