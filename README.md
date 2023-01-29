# Sequence-Alignment-Problem
This project aims to find tme and memory efficient ways to solve sequence alignment problem by implementing the [Hirschberg's Algorithm](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm#:~:text=Hirschberg's%20algorithm%20is%20simply%20described,of%20DNA%20and%20protein%20sequences.).

# Introduction

A sequence alignment problem focuses on finding the similarity between two sequences n and m
with the idea of minimizing the cost of alignment. The most common use case of the Sequence Alignment
problem is to find how closely are two DNA sequences related. In the following project, our input consists
of sequences formed by characters ‘A’, ‘C’, ‘T’, ‘G’ representing 4 bases of amino acids. The cost of
alignment is found by adding a gap penalty and mismatch cost. The gap penalty is augmented to the total cost
when there is a character in one sequence and a gap (missing character) in another sequence. For the
following project gap penalty is fixed to a value of 30. Mismatch cost is when a character in one sequence is
different from the corresponding character in another sequence. The mismatch cost is 0 for matching characters.
For the following project, a fixed mismatch cost is pre-defined for all 4 base pairs. We have implemented
the solution using two different approaches using Python Language: -

1) Dynamic Programming (basic approach)
2) Divide and Conquer plus Dynamic Programming (efficient approach)

# CPU Time Plot

![timPlot](https://user-images.githubusercontent.com/22619455/215299415-767fa5a8-5987-4cee-ba07-39ccb2c90e89.jpg)

# Analysis of Results

There is a tradeoff between CPU Time and memory consumed by the program. As the String length
increases the memory consumed by the Dynamic Programming (basic solution) becomes significantly
high as compared to D&C + Dynamic Programming(efficient) solution implementation though CPU run
Time continues to be double, hence for genomic sequences where lengths of sequences range into
billions the approach of divide and conquer is extremely efficient.

# Observations

**CPU Time: -** The Tme complexity of both the approaches is O(nm) asymptotically, where m is the length
of the first generated string and n is the length of the second generated string. There is a difference in CPU run
Time, where the memory efficient approach is taking approximately double the time to execute
compared to the basic approach.

**Memory Consump 4 on: -** For input sequences of small length, we see the memory consumed by efficient
solution is equal to or a little more than that of the basic approach. Once the string size increases a particular
threshold, the memory consumed by the basic program increases exponentially and the memory consumed
by the efficient solution is significantly less than the basic solution.

# Memory Consumption Plot

![memPlot](https://user-images.githubusercontent.com/22619455/215299419-990a592e-5b2f-451d-9439-cab5d82e81d0.jpg)

# Insights

**CPU Time: -** When we implement the Efficient (D&C plus Dynamic Programming) solution with every
proceeding divide step, the size of the subproblem becomes half as we keep finding the point of split for
the second sequence. Only during the first divide step, we would be doing the same amount of work as
done by the basic approach (dynamic programming solution). So, the progression is as follows

Cmn + ½ Cmn + ¼ Cmn...... = 2Cmn

That is CPU Tme is double for efficient solution (2Cmn) though asymptotic Tme complexity is O(mn)

**Memory Consump 4 on: -** Memory for basic solution requires O(nm) space to find the optimal solution.
Even though we only require the previously filled column to find the minimum penalty, we require enTre
an array of size nm to backtrack to find the actual string alignment and not just the penalty. The idea of
saving the alignment as we move forward won’t be helpful as when we are using dynamic programming,
we must check the entire solution space before arriving at the optimum solution. When we are finding
the solution using the divide and conquer strategy, we use the efficient version of the Dynamic programming
approach (using the last 2 columns) for finding the optimal split and constructing the optimal solution using the
basic DP algorithm when we have sequences of very small size. Hence, we only allocate memory for only
the current and last column of the 2d array, therefore, making the space complexity O(m).

For input sequences of small length, we see the memory consumed by the efficient solution is equal to or
a little more than that of the basic approach. This is because of function overhead for the recursive calls within
the divide and conquers algorithm. For longer sequences, the memory consumption by the dynamic
programming matrix quickly outgrows this function overhead in the divide and conquer algorithm.
Hence, asymptotically, the divide and conquer algorithm has a beBer space complexity than the dynamic
programming algorithm at the expense of Time complexity.
