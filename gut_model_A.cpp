/* Square lattice model, with progenitors on a quarter-size sub-lattice.
 * Differentiated cells chose to stratify/apoptose, causing one of the 
 * neighbouring progenitors (either two or four) to differentiate and migrate
 * into vacancy. The vacancy at the progenitor site is then filled by one of
 * four nearest progenitors.
*/

#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <boost/multi_array.hpp>
#include <cstdlib>
#include <map>
#include <boost/assert.hpp>
#include <fstream>
#include <vector>
#include "hilbert.hpp"

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
		BOOST_ASSERT(N%2 == 0);

		using namespace boost;
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j) {
				if(i % 2 == 0 && j % 2 == 0) { // A
					std::vector<uint32_t> p;
					p.push_back(i);
					p.push_back(j);
					g[i][j] = (hilbert_index(2, ceil(log(N) / log(2)), p) << 1) | 0x1;
					nB--;
				} else { // B
					g[i][j] = 0;
				}
			}
		}
	}

	void stratify_cell()
	{
		using namespace boost;

		// pick B cell
		BOOST_ASSERT(nB > 0);
		uniform_int<> coord_int(0, N-1);
		variate_generator<mt19937 &, uniform_int<> >
			coord(rng, coord_int);
		int i, j;
		do {
			i = coord(); j = coord();
		} while(g[i][j] & 0x1);
		
		// replace with neighbouring A
		BOOST_ASSERT(i % 2 == 1 || j % 2 == 1);
		uniform_int<> two(0, 1);
		variate_generator<mt19937 &, uniform_int<> >
			coin(rng, two);
		int ai = i, aj = j;
		if(i % 2 == 0 && j % 2 == 1)
			// A's should be above and below
			aj += coin() * 2 - 1;
		else if(i % 2 == 1 && j % 2 == 0)
			// A's should be left and right
			ai += coin() * 2 - 1;
		else {
			// four corners
			ai += coin() * 2 - 1;
			aj += coin() * 2 - 1;
		}
		sanitise(ai, aj);
		g[i][j] = g[ai][aj] & (~0x1);
		
		// replace A with neighbouring A
		// continuing using coin() from above
		int ci = ai, cj = aj;
		if(coin()) ci += 2*(coin() * 2 - 1);
		else       cj += 2*(coin() * 2 - 1);
		sanitise(ci, cj);
		g[ai][aj] = g[ci][cj];
	}

	double next_event() const {
		BOOST_ASSERT(nB > 0);
		using namespace boost;
		exponential_distribution<> time_distribution(nB);
		variate_generator<mt19937 &, exponential_distribution<> >
			event(rng, time_distribution);
		return event();
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

	void save_as_P3(const char * filename) {
		using namespace std;
		ofstream file(filename);
		file << "P3" << endl;
		file << N << " " << N << endl;
		file << "255" << endl;
		for(int i = 0; i < N; ++i) {
			for(int j = 0; j < N; ++j) {
				if(g[i][j]) {
					std::vector<uint32_t> rgb = hilbert_point(3, 8, g[i][j] + 0x2ffff);
					file << rgb[0] << " " << rgb[1] << " " << rgb[2] << " ";
				} else
					file << "0 0 0 ";
			}
			file << endl;
		}
	}

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

	if(argc <= 4) {
		cout << "usage: " << argv[0] << " <grid size> <time> <runs> <picture>" << endl;
		return -1;
	}

	int system_random = open("/dev/random", O_RDONLY);
	uint32_t seed;
	read(system_random, &seed, 4);
	close(system_random);
	
// doesn't actually work yet. need 1.43
//	random_device system_seed("/dev/random");

	rng.seed(seed);

	for(int i = 0; i < atoi(argv[3]); ++i) {
		grid_lattice grid(atoi(argv[1]));

		grid.time += grid.next_event();
		while(grid.time < atoi(argv[2])) {
			grid.stratify_cell();
			grid.time += grid.next_event();
		}

		std::map<uint32_t, uint32_t> hist = grid.histogram();
		for(std::map<uint32_t, uint32_t>::iterator i = hist.begin();
				i != hist.end(); ++i) {
			std::cout << i->second << std::endl;
		}

		if(i == 0) grid.save_as_P3(argv[4]);
	}

	return 0;

}
