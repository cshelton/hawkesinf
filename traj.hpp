#ifndef TRAJ_HPP
#define TRAJ_HPP

#include <vector>
#include <set>
#include <algorithm>

struct traj {
	double tend; // t0 always =0
	// events[i] are the times of events of type i
	std::vector<std::set<double>> events;

	// unobs[i][j] is the jth interval that is not observed
	//  for event type i (.first is the time of start of the interval
	//  and .second is the time of the end of the interval)
	std::vector<std::vector<std::pair<double,double>>> unobs;
};

void removewindow(traj &tr, int label, double t0, double t1) {
	for(auto i = tr.events[label].begin();i!=tr.events[label].end();)
		if (*i < t1 && *i>=t0)
			i= tr.events[label].erase(i);
		else ++i;

	for(auto i=tr.unobs[label].begin();i!=tr.unobs[label].end();)
		if (!(t1<i->first || t0>i->second)) {
			t0 = std::min(t0,i->first);
			t1 = std::max(t1,i->second);
			i = tr.unobs[label].erase(i);
		} else ++i;

	tr.unobs[label].emplace_back(t0,t1);
}
	
#endif
