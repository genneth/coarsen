#include <boost/random.hpp>
#include <numeric>
#include <iostream>

boost::mt19937 rng;

template <size_t N>
struct grid_lattice {
	explicit grid_lattice(): g({{0}}) {g[N/2][N/2] = 1;}

	unsigned int operator()(int i, int j) const {return g[i][j];}
	unsigned int & operator()(int i, int j) {return g[i][j];}

	size_t size() const {
		size_t s = 0;
		for(size_t i = 0; i < N; ++i) {
			for(size_t j = 0; j < N; ++j)
				s += g[i][j];
		}
		return s;
	}
private:
	unsigned int g[N][N];
};

template <size_t N>
std::ostream & operator<<(std::ostream & os, const grid_lattice<N> & g)
{
	for(unsigned int i = 0; i < N; ++i) {
		for(unsigned int j = 0; j < N; ++j)
			os << g(i,j) << " ";
		os << std::endl;
	}
	return os;
}

int main()
{
	using namespace std;
	using namespace boost;

	const size_t N = 21;

	exponential_distribution<> time_distribution(N*N);
	variate_generator<mt19937 &, exponential_distribution<> >
		next_event(rng, time_distribution);
	uniform_int<> x_distribution(0, N-1), y_distribution(0, N-1);
	variate_generator<mt19937 &, uniform_int<> >
		next_x(rng, x_distribution), next_y(rng, y_distribution);
	uniform_int<> direction(0,3);
	variate_generator<mt19937 &, uniform_int<> >
		neighbour(rng, direction);

	for(int count = 0; count < 10000;) {
		grid_lattice<N> grid;
		double t = 0.0;
		do {
			double dt = next_event();
			unsigned int x = next_x(), y = next_y();

			int dir = neighbour();
			int dx =    (dir & 1)  * ((dir & 2) - 1);
			int dy = (1-(dir & 1)) * ((dir & 2) - 1);

			grid(x, y) = grid((x+dx)%N, (y+dy)%N);

			t += dt;
		} while(t < 10.0);

		if(grid.size() > 0) {
			cout << grid.size() << endl;
			cout << grid << endl;
			count++;
		}
	}

	return 0;
}
 
