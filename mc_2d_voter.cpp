#include <boost/random.hpp>
#include <numeric>
#include <iostream>
#include <map>
#include <utility>
#include <cstdlib>
#include <boost/multi_array.hpp>

extern "C" {
#include <fcntl.h>
#include <unistd.h>
}

boost::mt19937 rng;

struct grid_lattice {
	typedef char cell_t;

	explicit grid_lattice(int N_): time(0.0), N(N_), g(boost::extents[N][N]) {
		active[std::make_pair(N/2,N/2)] = 0;
		set(N/2, N/2, 1);
	}

	cell_t cell(int i, int j) const {sanitise(i,j); return g[i][j];}
	void set(int i, int j, cell_t o) {
		// precondition: active.find(make_pair(i,j)) != active.end();
		// also, if cell(i,j) then its neighbours are in active

		// asymmetry in flipping
		if(cell(i,j) == o) return;
		else if(o) { // up
			cell_(i,j) = o;
			inc_neighbours(i+1,j);
			inc_neighbours(i-1,j);
			inc_neighbours(i,j+1);
			inc_neighbours(i,j-1);
		} else { // down
			cell_(i,j) = o;
			if(labelled_neighbours(i,j) == 0)
				active.erase(std::make_pair(i,j));
			dec_neighbours(i+1,j);
			dec_neighbours(i-1,j);
			dec_neighbours(i,j+1);
			dec_neighbours(i,j-1);
		}
	}

	int labelled_neighbours(int i, int j) const {
		return cell(i+1,j) + cell(i-1,j) + cell(i,j+1) + cell(i,j-1);
	}

	int size() const {
		int s = 0;
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j)
				s += g[i][j];
		}
		return s;
	}

	bool empty() const {return active.empty();}

	double next_event() const {
		using namespace boost;
		exponential_distribution<> time_distribution(active.size());
		variate_generator<mt19937 &, exponential_distribution<> >
			event(rng, time_distribution);
		return event();
	}

	void flip() {
		// precondition: !empty();
		using namespace boost;
		uniform_int<> active_cell(0, active.size()-1);
		variate_generator<mt19937 &, uniform_int<> >
			active_cell_chooser(rng, active_cell);
		active_list_t::iterator c = active.begin();
		std::advance(c, active_cell_chooser());

		uniform_int<> four(1, 4);
		variate_generator<mt19937 &, uniform_int<> > d4(rng, four);
		set(c->first.first, c->first.second, 
			(cell_t)(d4()) > c->second ? 0 : 1);
	}

	void restart() {
		// precondition: empty();
		time = 0.0;
		active[std::make_pair(N/2,N/2)] = 0;
		set(N/2, N/2, 1);
	}

	double time;
	const int N; // grid size

private:
	boost::multi_array<cell_t, 2> g;
	typedef std::map<std::pair<int,int>, cell_t> active_list_t;
	active_list_t active;

	void sanitise(int & i, int & j) const {i = (i+N)%N; j = (j+N)%N;} // positive dividend
	cell_t & cell_(int i, int j) {sanitise(i,j); return g[i][j];}
	void inc_neighbours(int i, int j) {
		sanitise(i,j);
		active_list_t::iterator c = active.find(std::make_pair(i,j));
		if(c == active.end())
			active.insert(std::make_pair(std::make_pair(i,j), 1));
		else
			c->second++;
	}
	void dec_neighbours(int i, int j) {
		sanitise(i,j);
		active_list_t::iterator c = active.find(std::make_pair(i,j));
		if(--(c->second) == 0 && cell(i,j) == 0)
			active.erase(c);
	}
};

std::ostream & operator<<(std::ostream & os, const grid_lattice & g)
{
	for(int i = 0; i < g.N; ++i) {
		for(int j = 0; j < g.N; ++j)
			os << static_cast<int>(g.cell(i,j)) << " ";
		os << std::endl;
	}
	return os;
}

int main(int argc, char ** argv)
{
	using namespace std;
	using namespace boost;

	if(argc <= 3) {
		cout << "usage: " << argv[0] << " <grid size> <time> <runs>" << endl;
		return -1;
	}

	int system_random = open("/dev/random", O_RDONLY);
	uint32_t seed;
	read(system_random, &seed, 4);
	close(system_random);

	rng.seed(seed);

	for(int count = 0; count < atoi(argv[3]);) {
		grid_lattice grid(atoi(argv[1]));
		while(true) {
			double dt = grid.next_event();
			grid.time += dt;
			if(grid.time > atof(argv[2])) break;
			grid.flip();
			if(grid.empty()) grid.restart();
		}

		if(!grid.empty()) {
			cout << grid.size() << endl;
//			cout << grid << endl;
			count++;
		}
	}

	return 0;
}
 
