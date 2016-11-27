# InterestingProblems
A collection of solutions to some problems I am interested in.

# Closest Sequence
The difference between two sequences of the same length a1, a2, a3,..., an and b1, b2, b3,..., bn can be defined as the sum of absolute differences between their respective elements:

diff(a, b) = |a1 - b1| + |a2 - b2| + ... + |an - bn|.

For the given sequences a and b (not necessarily having the same lengths) find a subsequence b' of b such that diff(a, b') is minimal. Return this difference.

Example

For a = [1, 2, 6] and b = [0, 1, 3, 4, 5], the output should be
closestSequence2(a, b) = 2.

The best subsequence will be b' = [1, 3, 5] which has a difference of 2 with a.

The problem was solved using recursion and dynamic programming. Dynamic programming to this problem can be optimized to use only O(n) space (where n = ||b||), however, to actually obtain the best subsequence, it would require O(mn) space where m = ||a||, n = ||b||. It becomes impractical for large sequences. Therefore, three types of local search algorithms are implemented to find the approximate solution. In particular, local search with random restart is the fastest but only finds a solution that is around 100% worse than the minimum. Simulated annealing runs the slowest but is managed to obtain a sequence that is only around 7% bigger than the minimum absolute differences. Genetics algorithm doesn't perform well when the ratio ||b||/||a|| is large (large search space) but is sometimes able to find the exact solution when the search space is small.

# Neat Printer
Consider the problem of neatly printing a paragraph with a monospaced font (all characters having the same width) on a printer. The input text is a sequence of n words of lengths l1, l2, ..., ln, measured in characters. We want to print this paragraph neatly on a number of lines that hold a maximum of M characters each. Our criteron of  “neatness” is as follows. If a given line contains words i through  j, where i ≤ j, and we leave exactly one space between words, the number of extra space characters at the end of the line is M − (j - i) − ∑k=i to j (lk) (that is, "character limit" - "the number of spaces taken between words" - "the number of characters taken by words from i to j"), which must be nonnegative so that the words fit on the line. We wish to minimize the sum, over all lines except the last, of the cubes of the numbers of extra space characters at the ends of lines. Give an algorithm to print a paragraph of n words neatly on a printer. Analyze the running time and space requirements of your algorithm.

# Count of Range Sum
Given an integer array nums, return the number of range sums that lie in [lower, upper] inclusive.
Range sum S(i, j) is defined as the sum of the elements in nums between indices i and j (i ≤ j), inclusive.

Note:
A naive algorithm of O(n2) is trivial. You MUST do better than that.

Example:
Given nums = [-2, 5, -1], lower = -2, upper = 2,
Return 3.
The three ranges are : [0, 0], [2, 2], [0, 2] and their respective sums are: -2, -1, 2.

# Binary Tree Maximum Path Sum
Given a binary tree, find the maximum path sum.

For this problem, a path is defined as any sequence of nodes from some starting node to any node in the tree along the parent-child connections. The path must contain at least one node and does not need to go through the root.

For example:
Given the below binary tree,
```
       1
      / \
     2   3
```
Return 6.
