# Closest Sequence
The difference between two sequences of the same length a1, a2, a3,..., an and b1, b2, b3,..., bn can be defined as the sum of absolute differences between their respective elements:

diff(a, b) = |a1 - b1| + |a2 - b2| + ... + |an - bn|.

For the given sequences a and b (not necessarily having the same lengths) find a subsequence b' of b such that diff(a, b') is minimal. Return this difference.

Example

For a = [1, 2, 6] and b = [0, 1, 3, 4, 5], the output should be
closestSequence2(a, b) = 2.

The best subsequence will be b' = [1, 3, 5] which has a difference of 2 with a.

The problem was solved using recursion and dynamic programming. Dynamic programming to this problem can be optimized to use only O(n) space (where n = ||b||), however, to actually obtain the best subsequence, it would require O(mn) space where m = ||a||, n = ||b||. It becomes impractical for large sequences. Therefore, three types of local search algorithms are implemented to find the approximate solution. In particular, local search with random restart is the fastest but only finds a solution that is around 100% worse than the minimum. Simulated annealing runs the slowest but is managed to obtain a sequence that is only around 7% bigger than the minimum absolute differences. Genetics algorithm doesn't perform well when the ratio ||b||/||a|| is large (large search space) but is sometimes able to find the exact solution when the search space is small.
