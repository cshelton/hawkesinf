#ifndef TRAJ_HPP
#define TRAJ_HPP

#include <vector>

struct traj {
	double tend; // t0 always =0
	// events[i][j] is the time of the jth event for event type i
	std::vector<std::vector<double>> events;

	// unobs[i][j] is the jth interval that is not observed
	//  for event type i (.first is the time of start of the interval
	//  and .second is the time of the end of the interval)
	std::vector<std::vector<std::pair<double,double>>> unobs;
};
	
#endif
