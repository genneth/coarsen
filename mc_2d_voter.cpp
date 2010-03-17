#include <boost/random.hpp>
#include <numeric>
#include <iostream>
#include <map>
#include <utility>

boost::mt19937 rng;

template <size_t N>
struct grid_lattice {
	explicit grid_lattice(): time(0.0), g({{0}}) {
		active[std::make_pair(N/2,N/2)] = 0;
		set(N/2, N/2, 1);
	}

	unsigned int cell(size_t i, size_t j) const {return g[i%N][j%N];}
	void set(size_t i, size_t j, unsigned int o) {
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

	size_t labelled_neighbours(size_t i, size_t j) const {
		return cell(i+1,j) + cell(i-1,j) + cell(i,j+1) + cell(i,j-1);
	}

	size_t size() const {
		size_t s = 0;
		for(size_t i = 0; i < N; ++i) {
			for(size_t j = 0; j < N; ++j)
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
			(unsigned int)(d4()) > c->second ? 0 : 1);
	}

	double time;

private:
	unsigned int g[N][N];
	typedef std::map<std::pair<size_t,size_t>, unsigned int> active_list_t;
	active_list_t active;

	unsigned int & cell_(size_t i, size_t j) {return g[i%N][j%N];}
	void inc_neighbours(size_t i, size_t j) {
		active_list_t::iterator c = active.find(std::make_pair(i,j));
		if(c == active.end())
			active[std::make_pair(i,j)] = 1;
		else
			c->second++;
	}
	void dec_neighbours(size_t i, size_t j) {
		active_list_t::iterator c = active.find(std::make_pair(i,j));
		if(--c->second == 0 && cell(i,j) == 0)
			active.erase(c);
	}
};

template <size_t N>
std::ostream & operator<<(std::ostream & os, const grid_lattice<N> & g)
{
	for(unsigned int i = 0; i < N; ++i) {
		for(unsigned int j = 0; j < N; ++j)
			os << g.cell(i,j) << " ";
		os << std::endl;
	}
	return os;
}

int main()
{
	using namespace std;
	using namespace boost;

	const size_t N = 101;

	for(int count = 0; count < 1;) {
		grid_lattice<N> grid;
		do {
			double dt = grid.next_event();
			grid.flip();
			grid.time += dt;
		} while(grid.time < 50.0 && !grid.empty());

		if(!grid.empty()) {
			cout << grid.size() << endl;
			cout << grid << endl;
			count++;
		}
	}

	return 0;
}
 
