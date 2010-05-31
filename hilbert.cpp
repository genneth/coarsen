#include "hilbert.hpp"
#include <boost/assert.hpp>
#include <cmath>

using namespace std;

static uint32_t graycode(uint32_t x)
{
	return x^(x>>1U);
}

static uint32_t igraycode(uint32_t x)
{
	if(x == 0) return x;
	
	uint32_t m = ceil(log((double)x) / log(2.0))+1;
	uint32_t i = x, j = 1;
	while(j < m) {
		i = i ^ (x>>j);
		j += 1;
	}
	return i;
}

static uint32_t rrot(uint32_t x, uint32_t i, uint32_t width)
{
	BOOST_ASSERT(x < (1U << width));
	i = i%width;
	x = (x>>i) | (x<<(width-i));
	return x&((1U<<width)-1U);
}

static uint32_t lrot(uint32_t x, uint32_t i, uint32_t width)
{
	BOOST_ASSERT(x < (1U << width));
	i = i%width;
	x = (x<<i) | (x>>(width-i));
	return x&((1U<<width)-1U);
}

static uint32_t tsb(uint32_t x, uint32_t width)
{
	BOOST_ASSERT(x < (1U << width));
	uint32_t i = 0;
	while(x&1 && i <= width) {
		x >>= 1;
		i++;
	}
	return i;
}

static uint32_t setbit(uint32_t x, uint32_t w, uint32_t i, uint32_t b)
{
	BOOST_ASSERT(b == 0 || b == 1);
	BOOST_ASSERT(i < w);
	if(b) return x | (1U << (w-i-1));
	else  return x & (~(1U << (w-i-1)));
}

static uint32_t bitrange(uint32_t x, uint32_t width, uint32_t start, uint32_t end)
{
	return x >> (width - end) & ((1 << (end - start)) - 1);
}

static uint32_t transform(uint32_t entry, uint32_t direction, uint32_t width, uint32_t x)
{
	BOOST_ASSERT(x < (1U << width));
	BOOST_ASSERT(entry < (1U << width));
	return rrot((x^entry), direction + 1, width);
}

static uint32_t itransform(uint32_t entry, uint32_t direction, uint32_t width, uint32_t x)
{
	BOOST_ASSERT(x < (1U << width));
	BOOST_ASSERT(entry < (1U << width));
	return lrot(x, direction + 1, width)^entry;
}

static uint32_t direction(uint32_t x, uint32_t n)
{
	BOOST_ASSERT(x < (1U << n));
	if(x == 0) 
		return 0;
	else if(x % 2 == 0)
		return tsb(x-1, n)%n;
	else
		return tsb(x, n)%n;
}

static uint32_t entry(uint32_t x)
{
	if(x == 0) return 0;
	else
		return graycode(2*((x-1)/2));
}

uint64_t hilbert_index(uint32_t dimension, uint32_t order, vector<uint32_t> p)
{
	BOOST_ASSERT(dimension == p.size());
	
	uint64_t h = 0;
	uint32_t e = 0, d = 1;
	
	for(uint32_t i = 0; i < order; ++i) {
		uint32_t l = 0;
		for(uint32_t x = 0; x < dimension; ++x) {
			uint32_t b = bitrange(p[x], order, i, i+1);
			l |= b<<x;
		}
		l = transform(e, d, dimension, l);
		uint32_t w = igraycode(l);
		e = e ^ lrot(entry(w), d+1, dimension);
		d = (d + direction(w, dimension) + 1)%dimension;
		h = (h<<dimension)|w;
	}
	return h;
}

vector<uint32_t> hilbert_point(uint32_t dimension, uint32_t order, uint64_t h)
{
	uint32_t hwidth = order*dimension;
	uint32_t e = 0, d = 0;
	vector<uint32_t> p(dimension, 0);
	for(uint32_t i = 0; i < order; ++i) {
		uint32_t w = bitrange(h, hwidth, i*dimension, i*dimension+dimension);
		uint32_t l = graycode(w);
		l = itransform(e, d, dimension, l);
		for(uint32_t j = 0; j < dimension; ++j) {
			uint32_t b = bitrange(l, dimension, j, j+1);
			p[j] = setbit(p[j], order, i, b);
		}
		e = e ^ lrot(entry(w), d+1, dimension);
		d = (d + direction(w, dimension) + 1)%dimension;
	}
	return p;
}

