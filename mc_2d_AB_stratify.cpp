#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <cstdlib>
#include <map>
#include <boost/assert.hpp>

extern "C" {
#include <fcntl.h>
#include <unistd.h>
}

boost::mt19937 rng; // returns uint32_t

struct grid_lattice {

public:
	double time;
	const size_t N;
	uint32_t nB;
private:
	// the lowest bit is used to indicate A(1) or B(0); the upper 31 are a label
	boost::multi_array<uint32_t, 2> g;

public:
	explicit grid_lattice(size_t N_): time(0.0), N(N_), nB(N*N), g(boost::extents[N][N]) {
		using namespace boost;
		uniform_real<> type_real(0.0, 1.0);
		variate_generator<mt19937 &, uniform_real<> >
			type(rng, type_real);
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j) {
				if(type() < 0.36) { // A
					g[i][j] = 0x1;
					nB--;
				} else { // B
					g[i][j] = 0;
				}
			}
		}
		relabel();
	}

	void relabel() {
		uint32_t label = 1;
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j) {
				if(g[i][j] & 0x1) { // A
					g[i][j] = (label << 1) | 0x1;
					label++;
				} else { // B
					g[i][j] = 0;
				}
			}
		}
	}

	std::pair<int16_t, int16_t> pick_stratifying_cell() const
	{
		BOOST_ASSERT(nB > 0);
		using namespace boost;
		uniform_int<> coord_int(0, N-1);
		variate_generator<mt19937 &, uniform_int<> >
			coord(rng, coord_int);
		int16_t i, j;
		do {
			i = coord(); j = coord();
		} while(g[i][j] & 0x1);
		return std::make_pair(i,j);
	}

	double next_event() const {
		BOOST_ASSERT(nB > 0);
		using namespace boost;
		exponential_distribution<> time_distribution(nB);
		variate_generator<mt19937 &, exponential_distribution<> >
			event(rng, time_distribution);
		return event();
	}

	void migrate_vacancy(std::pair<int16_t, int16_t> location) {
		using namespace boost;
		using namespace std;

		BOOST_ASSERT(nB < N*N);
		uniform_int<> neighbour_int(0, 3);
		variate_generator<mt19937 &, uniform_int<> >
			neighbour(rng, neighbour_int);
		int i = location.first, j = location.second;
//		std::cerr << "at " << i << ", " << j << ", ";
		int dir = neighbour();
		int next_i = (dir & 0x1) ? i : i + dir - 1;
		int next_j = (dir & 0x1) ? j + (dir&(~0x1)) - 1 : j;
		sanitise(next_i, next_j);

		if(g[next_i][next_j] & 0x1) { // divide A 
//			std::cerr << "divide " << next_i << ", " << next_j << std::endl;
			uniform_real<double> choice_real(0, 1.0);
			variate_generator<mt19937 &, uniform_real<double> >
				choice(rng, choice_real);
			double r = 0.20;
			double c = choice();
			if(c < r) { // AA
				g[i][j] = g[next_i][next_j];
				nB--;
			} else if(c < 0.5) { // AB
				g[i][j] = g[next_i][next_j] & (~0x1);
			} else if(c < (1-r)) { // BA
				g[i][j] = g[next_i][next_j];
				g[next_i][next_j] &= ~0x1;
			} else { // BB
				g[next_i][next_j] &= ~0x1;
				g[i][j] = g[next_i][next_j];
				nB++;
			}
		} else { // move B into vacancy
//			std::cerr << "migrate " << next_i << ", " << next_j << std::endl;
			g[i][j] = g[next_i][next_j];
			migrate_vacancy(make_pair(next_i, next_j));
		}
	}

	std::map<uint32_t, uint32_t> histogram() const
	{
		std::map<uint32_t, uint32_t> hist;
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j) {
				if(g[i][j] != 0) { // not an unlabelled B
					std::map<uint32_t, uint32_t>::iterator
						k = hist.find(g[i][j] >> 1);
					if(k == hist.end())
						hist[g[i][j] >> 1] = 1;
					else
						k->second++;
				}
			}
		}
		return hist;
	}
	
	friend std::ostream & operator<<(std::ostream &, const grid_lattice &);

private:
	void sanitise(int & i, int & j) const {
		// precondition: i,j >= -1
		i = (i+N)%N;
		j = (j+N)%N;
	}

};

std::ostream & operator<<(std::ostream & os, const grid_lattice & g)
{
	for(int i = 0; i < g.N; ++i) {
		for(int j = 0; j < g.N; ++j)
			os << g.g[i][j] << " ";
		os << std::endl;
	}
	os << g.nB << std::endl;
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

	for(int i = 0; i < atoi(argv[3]); ++i) {
		grid_lattice grid(atoi(argv[1]));

		grid.time += grid.next_event();
		while(grid.time < 10.0 && grid.nB > 0 && grid.nB < grid.N*grid.N) {
			std::pair<int16_t, int16_t> vac = grid.pick_stratifying_cell();
			grid.migrate_vacancy(vac);
			grid.time += grid.next_event();
		}

		grid.time = 0;
		grid.relabel();

		grid.time += grid.next_event();
		while(grid.time < atof(argv[2]) && grid.nB > 0 && grid.nB < grid.N*grid.N) {
			std::pair<int16_t, int16_t> vac = grid.pick_stratifying_cell();
			grid.migrate_vacancy(vac);
			grid.time += grid.next_event();
		}

		std::map<uint32_t, uint32_t> hist = grid.histogram();
		for(std::map<uint32_t, uint32_t>::iterator i = hist.begin();
				i != hist.end(); ++i) {
			std::cout << i->second << std::endl;
		}
	}

	return 0;

}
