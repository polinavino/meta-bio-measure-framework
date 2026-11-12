// Random-to-random insertion Gibbs sampler for uniform linear extensions.
// Conditional on the other objects' order, every legal slot for a selected
// object has equal probability. Choosing the object uniformly gives detailed
// balance with the uniform distribution. Self transitions make it aperiodic.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 6) return 2;
    const int n = std::atoi(argv[1]), burn = std::atoi(argv[2]);
    const int draws = std::atoi(argv[3]), thin = std::atoi(argv[4]);
    std::mt19937_64 rng(std::strtoull(argv[5], nullptr, 10));
    std::vector<uint8_t> before(static_cast<size_t>(n) * n);
    std::vector<int32_t> order(n);
    std::cin.read(reinterpret_cast<char*>(before.data()), before.size());
    std::cin.read(reinterpret_cast<char*>(order.data()), n * sizeof(int32_t));
    if (!std::cin) return 3;
    std::uniform_int_distribution<int> select(0, n - 1);
    std::vector<double> ranks(n);
    for (int sweep = 0; sweep < burn + draws * thin; ++sweep) {
        for (int move = 0; move < n; ++move) {
            int old = select(rng), item = order[old], lo = 0, hi = n - 1;
            for (int j = old - 1; j >= 0; --j) {
                if (before[static_cast<size_t>(order[j]) * n + item]) {
                    lo = j + 1; break;
                }
            }
            for (int j = old + 1; j < n; ++j) {
                if (before[static_cast<size_t>(item) * n + order[j]]) {
                    hi = j - 1; break;
                }
            }
            int dest = std::uniform_int_distribution<int>(lo, hi)(rng);
            if (dest < old) std::move_backward(order.begin() + dest, order.begin() + old, order.begin() + old + 1);
            if (dest > old) std::move(order.begin() + old + 1, order.begin() + dest + 1, order.begin() + old);
            order[dest] = item;
        }
        if (sweep >= burn && (sweep - burn) % thin == 0) {
            for (int j = 0; j < n; ++j) ranks[order[j]] = j + 1;
            std::cout.write(reinterpret_cast<char*>(ranks.data()), n * sizeof(double));
        }
    }
    return std::cout ? 0 : 4;
}
