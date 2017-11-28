#ifndef PLOTTRAJ_HPP
#define PLOTTRAJ_HPP

#include "gnuplot-cpp/gnuplot_i.hpp"
#include <string>
#include <array>
#include "traj.hpp"
#include <iostream>

class trajplot {
public:
	trajplot(const std::string &title, const traj &tr) : gplot(title) {
		plot(tr);
	}
	trajplot(std::string title) : gplot(title) {
	}

	void settitle(const std::string &title) {
		gplot.set_title(title);
	}

	void plot(const traj &tr) {
		gplot.reset_plot();
		gplot << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb \"#666666\" behind";
		int np = tr.events.size();
		gplot.set_yrange(-0.5,np-0.5);
		gplot.set_xrange(0,tr.tend);
		gplot.set_pointsize(1.6).set_style("points");
		for(int i=0;i<np;i++)
			if (!tr.events[i].empty())
				gplot.plot_xy(tr.events[i],std::vector<int>(tr.events[i].size(),i));
		gplot.set_style("lines");
		for(int i=0;i<np && i<tr.unobs.size();i++) 
			for(auto &p : tr.unobs[i])
				gplot.plot_xy(std::array<double,2>{p.first,p.second},
						std::array<int,2>{i,i});
	}
private:
	Gnuplot gplot;
};

#endif
