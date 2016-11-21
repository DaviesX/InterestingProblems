#include <cmath>
#include <ctime>
#include <cstring>
#include <random>
#include <limits>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

//////////////////////////////////////////// Algorithms /////////////////////////////////////////////
struct DiffRecord
{
        DiffRecord(int diff, const std::vector<int>& b_k):
                diff(diff), b_k(b_k)
        {
        }

        DiffRecord(int diff):
                diff(diff)
        {
        }

        int                     diff;
        std::vector<int>        b_k;
};

DiffRecord closest_sequence_dc(const int* a, const int* b, int i, int j, int len_a, int len_b)
{
        // C(n, k) = 0, no pairing is possible.
        if (i == len_a)
                return DiffRecord(0);
        // The number of elements left in "a" is bigger than that in "b."
        if (len_a - i > len_b - j)
                return DiffRecord(std::numeric_limits<int>::max()/2);

        // Running on recurrence relation.
        // C(i, j) = min{|a[i] - b[j]| + C(i + 1, j + 1), C(i, j + 1)}
        DiffRecord min(std::numeric_limits<int>::max());
        for (int k = j; k < len_b; k ++) {
                const DiffRecord& record = closest_sequence_dc(a, b, i + 1, k + 1, len_a, len_b);
                int diff_ik = abs(a[i] - b[k]) + record.diff;

                if (diff_ik < min.diff) {
                        min.diff = diff_ik;
                        min.b_k = record.b_k;
                        min.b_k.push_back(b[k]);
                }
        }
        return min;
}

int closest_sequence_dc(const int* a, const int* b, int len_a, int len_b, int* b_prime)
{
        const DiffRecord& record = closest_sequence_dc(a, b, 0, 0, len_a, len_b);
        for (int i = 0; i < len_a; i ++) {
                b_prime[i] = record.b_k[len_a - 1 - i];
        }
        return record.diff;
}

int closest_sequence_dp(const int* a, const int* b, int len_a, int len_b, int* b_prime)
{
        int* min_diff = new int [len_a*len_b];

        // Boundary conditions.
        // C(n - i, m - i) = sum(k from i to n - 1) {|a[n - k]  - b[m - k]|}
        int partial_sum = 0;
        for (int i = len_a - 1; i >= 0; i --) {
                int j = len_b - len_a + i;
                min_diff[i + j*len_a] = partial_sum + abs(a[i] - b[j]);
                partial_sum += abs(a[i] - b[j]);
        }
        // C(n - 1, i) = min(k from i to m) {|a[n - 1] - b[k]|}
        int min = std::numeric_limits<int>::max();
        for (int j = len_b - 1; j >= len_a - 1; j --) {
                int diff = abs(a[(len_a - 1)] - b[j]);
                if (diff < min)
                        min = diff;
                min_diff[(len_a - 1) + j*len_a] = min;
        }
        // Build up recurrence.
        // C(i, j) = min{|a[i] - b[j]| + C(i + 1, j + 1), C(i, j + 1)}
        int d = 3;
        for (int i = len_a - 2; i >= 0; i --) {
                for (int j = len_b - d; j >= i; j --) {
                        min_diff[i + j*len_a] = std::min(min_diff[(i + 1) + (j + 1)*len_a] + abs(a[i] - b[j]),
                                        min_diff[i + (j + 1)*len_a]);
                }
                d ++;
        }

        // Trace the path.
        int shortest_j = 0;
        for (int i = 0; i < len_a; i ++) {
                int min_i = min_diff[i + shortest_j*len_a];
                for (int j = shortest_j; j <= len_b - len_a + i; j ++) {
                        int diff = min_diff[i + j*len_a];
                        if (diff > min_i) {
                                b_prime[i] = b[j - 1];
                                shortest_j = j;
                                break;
                        }
                }
        }

        int abs_min = min_diff[0];
        delete [] min_diff;
        return abs_min;
}

void random_start(int* b_k, int len_a, int len_b)
{
        std::vector<int> candidates(len_b);
        for (int i = 0; i < len_b; i ++)
                candidates[i] = i;
        std::random_shuffle(candidates.begin(), candidates.end());
        for (int i = 0; i < len_a; i ++)
                b_k[i] = candidates[i];
        std::sort(b_k, &b_k[len_a]);
}

int compute_abs_diff(const int* a, const int* b, int len_a, const int* b_k)
{
        int diff = 0;
        for (int i = 0; i < len_a; i ++)
                diff += abs(a[i] - b[b_k[i]]);
        return diff;
}

void copy_bk_to_bprime(int* b_prime, const int* b, int len_a, const int* b_k)
{
        for (int j = 0; j < len_a; j ++)
                b_prime[j] = b[b_k[j]];
}

int neighbor_of(const int* b_k, int len_a, int len_b, int& selected)
{
retry:
        selected = rand()%len_a;
        int min_k = (selected == 0 ? 0 : b_k[selected - 1]) + 1;
        int max_k = (selected == len_a - 1 ? len_b - 1 : b_k[selected + 1]) - 1;
        if (min_k >= max_k)
                goto retry;
        else
                return min_k + rand()%(max_k - min_k + 1);
}

int closest_sequence_ls(const int* a, const int* b, int len_a, int len_b, int* b_prime)
{
        int* b_k = new int [len_a];

        int glb_min = std::numeric_limits<int>::max();
        for (int i = 0; i < 500; i ++) {
                random_start(b_k, len_a, len_b);

                int min = compute_abs_diff(a, b, len_a, b_k);
                bool has_improvment;
                do {
                        has_improvment = false;
                        for (int j = 0; j < len_b; j ++) {
                                int selected;
                                int new_neighbor = neighbor_of(b_k, len_a, len_b, selected);
                                if (new_neighbor == -1)
                                        continue;

                                int old = b_k[selected];

                                // Hill climbing.
                                int improvement = abs(a[selected] - b[old]) - abs(a[selected] - b[new_neighbor]);
                                if (improvement > 0) {
                                        // Confirm neighbor.
                                        has_improvment = true;
                                        min = min - improvement;
                                        b_k[selected] = new_neighbor;
                                }
                        }
                } while (has_improvment);

                if (min < glb_min) {
                        copy_bk_to_bprime(b_prime, b, len_a, b_k);
                        glb_min = min;
                }
        }
        delete [] b_k;
        return glb_min;
}

bool bernoulli_0(float p_1, float& p)
{
        p = (float) (rand()%RAND_MAX)/(float) RAND_MAX;
        return p <= p_1;
}

int annealing(const int* a, const int* b, int len_a, int len_b, int* b_k,
              uint64_t& k, const uint64_t interval, const uint64_t MAX_STEPS, const int MAX_TOLERANCE)
{
        int min = compute_abs_diff(a, b, len_a, b_k);

        float temperature;
        float t0 = (float) min/len_a;
        const float alpha = 0.95f;
        const float gamma = 1.0f/MAX_STEPS;

        int tolerance = 0;

        int max_impro = 0;
        int max_regress = 0;
        uint64_t e = k + interval;
        while (k < e || tolerance < MAX_TOLERANCE) {
                uint64_t l;
                max_regress = 0;
                max_impro = 0;
                for (l = k; l < k + len_b; l ++) {
                        temperature = t0*pow(alpha, l*gamma);

                        int selected = rand()%len_a;
                        int new_neighbor = neighbor_of(b_k, len_a, len_b, selected);
                        if (new_neighbor == -1)
                                continue;

                        int old = b_k[selected];

                        // Hill climbing.
                        int improvement = abs(a[selected] - b[old]) - abs(a[selected] - b[new_neighbor]);
                        if (improvement > 0) {
                                // Confirm neighbor.
                                min = min - improvement;
                                b_k[selected] = new_neighbor;

                                max_impro = std::max(improvement, max_impro);
                        } else {
                                // Accept bad move with the probability exp(-|ds|/T).
                                float prob = exp((float) improvement/temperature);
                                float p;
                                bool is_accept = bernoulli_0(prob, p);
                                if (is_accept) {
                                        min = min - improvement;
                                        b_k[selected] = new_neighbor;

                                        max_regress = std::min(improvement, max_regress);
                                }
                        }
                }
                //std::cout << "temp " << temperature << ", "
                //          << "max impro " << max_impro << ", "
                //          << "max regress " << max_regress << std::endl;
                if (max_impro == 0 && max_regress == 0)
                        tolerance ++;
                else
                        tolerance = 0;
                k = l;
        }
        return min;
}

int closest_sequence_sa(const int* a, const int* b, int len_a, int len_b, int* b_prime)
{
        const int64_t MAX_STEPS = 10000000;

        int* b_k = new int [len_a];
        random_start(b_k, len_a, len_b);

        uint64_t s = 0;
        int min = annealing(a, b, len_a, len_b, b_k, s, MAX_STEPS, MAX_STEPS, 100);

        copy_bk_to_bprime(b_prime, b, len_a, b_k);
        delete [] b_k;
        return min;
}

void reproduce(int* population, int pop_size, float mutation_rate,
               const int* a, const int* b, int len_a, int len_b,
               std::pair<int, int>* fitness,
               bool* is_taken,
               const std::map<int, int>& survivors)
{
        for (int i = 0; i < pop_size; i ++)
                is_taken[i] = false;
        int num_mutations = ceil(mutation_rate*(float) len_a);
        for (std::pair<int, int> survivor: survivors) {
                int n_reproduced = 0;
                for (int i = 0; i < pop_size && n_reproduced < survivor.second; i ++) {
                        if (survivors.find(i) == survivors.end() && !is_taken[i]) {
                                is_taken[i] = true;

                                // Copy DNA and mutate.
                                memcpy(&population[i*len_a], &population[survivor.first*len_a], sizeof(int)*len_a);
                                fitness[i].second = fitness[survivor.first].second;
                                for (int m = 0; m < num_mutations; m ++) {
                                        int selected;
                                        int neighbor = neighbor_of(&population[i*len_a], len_a, len_b, selected);
                                        int old = population[i*len_a + selected];
                                        fitness[i].first = i;
                                        fitness[i].second = fitness[i].second -
                                                        abs(a[selected] - b[old]) + abs(a[selected] - b[neighbor]);
                                        population[i*len_a + selected] = neighbor;
                                }

                                n_reproduced ++;
                        }
                }
        }
}

int closest_sequence_ga(const int* a, const int* b, int len_a, int len_b, int* b_prime)
{
        int NUM_RESTARTS = 10;
        int MAX_ITERATIONS = 2000;
        int MAX_POPULATION = 1000;
        int MAX_SURVIVAL = 10;
        const float MUTATION_RATE = 0.01f;

        int* population = new int [len_a*MAX_POPULATION];
        std::pair<int, int>* fitness = new std::pair<int, int> [MAX_POPULATION];
        std::pair<int, int>* sorted_fitness = new std::pair<int, int> [MAX_POPULATION];
        bool* is_taken = new bool [MAX_POPULATION];
        std::map<int, int> survivors;

        int glb_min = std::numeric_limits<int>::max();

        for (int r = 0; r < NUM_RESTARTS; r ++) {
                // Initial population.
                for (int i = 0; i < MAX_POPULATION; i ++) {
                        random_start(&population[i*len_a], len_a, len_b);
                        fitness[i] = std::pair<int, int>(i, compute_abs_diff(a, b, len_a, &population[i*len_a]));
                        is_taken[i] = true;
                }

                int i = 0;
                while (true) {
                        // Selection.
                        memcpy(sorted_fitness, fitness, sizeof(std::pair<int, int>)*MAX_POPULATION);
                        std::stable_sort(sorted_fitness, sorted_fitness + MAX_POPULATION,
                                         [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) -> bool
                        {
                                return p1.second < p2.second;
                        });

                        if (i == MAX_ITERATIONS - 1)
                                break;

                        survivors.clear();
                        int total_score = 0;
                        for (int j = 0; j < MAX_SURVIVAL; j ++) {
                                total_score += sorted_fitness[j].second;
                        }

                        for (int j = 0; j < MAX_SURVIVAL; j ++) {
                                int num_succ = MAX_POPULATION*
                                                (float) (1.0f - (float) sorted_fitness[j].second/(float) total_score);
                                survivors[sorted_fitness[j].first] = num_succ;
                        }

                        // Reproduce.
                        reproduce(population, MAX_POPULATION, MUTATION_RATE,
                                  a, b, len_a, len_b, fitness, is_taken, survivors);
                        i ++;
                }

                int min = sorted_fitness[0].second;
                if (min < glb_min) {
                        copy_bk_to_bprime(b_prime, b, len_a, &population[sorted_fitness[0].first*len_a]);
                        glb_min = min;
                }
        }

        delete [] population;
        delete [] fitness;
        delete [] sorted_fitness;
        delete [] is_taken;

        return glb_min;
}

///////////////////////////////////////////// Tests //////////////////////////////////////////////

template<typename T>
void print_array(T* a, int len, std::ostream& os)
{
        for (int i = 0; i < len; i ++)
                os << a[i] << " ";
}

template<typename T>
bool is_equal_sequence(T* a, T* b, int len)
{
        for (int i = 0; i < len; i ++) {
                if (a[i] != b[i])
                        return false;
        }
        return true;
}

int main()
{
        int seq_a[] = {1, -2, 4, 6, 8};
        int seq_b[] = {100, 64, 8, 12, 54, 29, -12, 54, 86, 12, 120};
        int closest[sizeof(seq_a)/sizeof(int)];
        int abs_diff = closest_sequence_dc(seq_a, seq_b,
                                           sizeof(seq_a)/sizeof(int), sizeof(seq_b)/sizeof(int),
                                           closest);

        std::cout << "a = ";
        print_array(seq_a, sizeof(seq_a)/sizeof(int), std::cout);
        std::cout << std::endl;

        std::cout << "b = ";
        print_array(seq_b, sizeof(seq_b)/sizeof(int), std::cout);
        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << "Absolute difference of (a, b) using DC: " << abs_diff << std::endl;
        std::cout << "b' = ";
        print_array(closest, sizeof(seq_a)/sizeof(int), std::cout);
        std::cout << std::endl;

        abs_diff = closest_sequence_dp(seq_a, seq_b,
                                       sizeof(seq_a)/sizeof(int), sizeof(seq_b)/sizeof(int),
                                       closest);
        std::cout << "Absolute difference of (a, b) using DP: " << abs_diff << std::endl;
        std::cout << "b' = ";
        print_array(closest, sizeof(seq_a)/sizeof(int), std::cout);
        std::cout << std::endl;

        std::cout << std::endl << "Speed comparison: " << std::endl;
        const float ratio = 3;
        const int SIZE_LARGE_A = 1000;
        const int SIZE_LARGE_B = SIZE_LARGE_A*ratio;
        int* large_a = new int [SIZE_LARGE_A];
        int* large_b = new int [SIZE_LARGE_B];
        int* large_closest = new int [SIZE_LARGE_A];
        for (unsigned i = 0; i < SIZE_LARGE_A; i ++)
                large_a[i] = std::rand()%2000 - 1000;
        for (unsigned i = 0; i < SIZE_LARGE_B; i ++)
                large_b[i] = std::rand()%2000 - 1000;

        std::cout << "|a| = " << SIZE_LARGE_A << std::endl;
        std::cout << "|b| = " << SIZE_LARGE_B << std::endl;

        clock_t t0, t1;
        if (SIZE_LARGE_A*SIZE_LARGE_B <= 300) {
                t0 = std::clock();
                abs_diff = closest_sequence_dc(large_a, large_b,
                                               SIZE_LARGE_A, SIZE_LARGE_B,
                                               large_closest);
                t1 = std::clock();
                std::cout << "Absolute difference of (a, b) using DC: " << abs_diff << std::endl;
                std::cout << "Using time " << (float) (t1 - t0)/CLOCKS_PER_SEC << "s" << std::endl;
        } else {
                std::cout << "DC is not possible for such problem size." << std::endl;
        }

        t0 = std::clock();
        abs_diff = closest_sequence_dp(large_a, large_b,
                                       SIZE_LARGE_A, SIZE_LARGE_B,
                                       large_closest);
        t1 = std::clock();

        std::cout << "Absolute difference of (a, b) using DP: " << abs_diff << std::endl;
        std::cout << "Using time " << (float) (t1 - t0)/CLOCKS_PER_SEC << "s" << std::endl;

        t0 = std::clock();
        abs_diff = closest_sequence_ls(large_a, large_b,
                                       SIZE_LARGE_A, SIZE_LARGE_B,
                                       large_closest);
        t1 = std::clock();
        std::cout << "Absolute difference of (a, b) using LS: " << abs_diff << std::endl;
        std::cout << "Using time " << (float) (t1 - t0)/CLOCKS_PER_SEC << "s" << std::endl;

        t0 = std::clock();
        abs_diff = closest_sequence_sa(large_a, large_b,
                                       SIZE_LARGE_A, SIZE_LARGE_B,
                                       large_closest);
        t1 = std::clock();
        std::cout << "Absolute difference of (a, b) using SA: " << abs_diff << std::endl;
        std::cout << "Using time " << (float) (t1 - t0)/CLOCKS_PER_SEC << "s" << std::endl;

        t0 = std::clock();
        abs_diff = closest_sequence_ga(large_a, large_b,
                                       SIZE_LARGE_A, SIZE_LARGE_B,
                                       large_closest);
        t1 = std::clock();

        std::cout << "Absolute difference of (a, b) using GA: " << abs_diff << std::endl;
        std::cout << "Using time " << (float) (t1 - t0)/CLOCKS_PER_SEC << "s" << std::endl;

        delete [] large_a;
        delete [] large_b;
        delete [] large_closest;
        return 0;
}
