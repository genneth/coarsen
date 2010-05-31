#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <cstdint>
#include <vector>

uint64_t hilbert_index(uint32_t dim, uint32_t order, std::vector<uint32_t> p);
std::vector<uint32_t> hilbert_point(uint32_t dim, uint32_t order, uint64_t h);

#endif // HILBERT_HPP
