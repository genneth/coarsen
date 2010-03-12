#include <boost/random.hpp>
#include <numeric>
#include <iostream>

boost::mt19937 rng;

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
		unsigned int grid[N][N] = {{0}};
		grid[N/2][N/2] = 1;
		double t = 0.0;
		do {
			double dt = next_event();
			unsigned int x = next_x(), y = next_y();

			int dir = neighbour();
			int dx = (dir & 1) * 2 - 1;
			int dy = (dir & 2) - 1;

			grid[x][y] = grid[(x+dx)%N][(y+dy)%N];

			t += dt;
			// cout << t << " " << x << " " << y << " " << dx << " " << dy << endl;
		} while(t < 10.0);

		size_t size = 0;
		for(size_t i = 0; i < N; ++i) {
			for(size_t j = 0; j < N; ++j) {
				size += grid[i][j];
			}
		}
		if(size > 0) {
			cout << size << endl;
			count++;
		}
	}

	return 0;
}
 
