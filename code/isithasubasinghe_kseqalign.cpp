/*
**
** TESTED ON: GNU 10.3.0 & GNU 11.1.0
** COMPILING COMMAND: g++ isithasubasinghe_kseqalign.cpp -I. -fopenmp -static-libstdc++ -O3 -fno-exceptions -fno-rtti -march=native -funroll-loops -pedantic -Wall -o kseqalign
**
** OMP_NUM_THREADS should ideally be chosen such that entire NUMA nodes are dedicated to the program, given NUMA nodes: 0,1,2 with CPUs {8,8,1}, OMP_NUM_THREADS should be 16 not 17 in this instance, some effort is placed on the user to determine
** the cost/benefit ratio here. Refer to the report on our NUMA mitigation strategy to understand why. It isn't a major limitation however, roughly 1/17ths of the memory would be present on memory closer to NUMA node 2, it simply is not likely to provide much benefit
**
** OMP_PROC_BIND=close
** OMP_PLACES=cores
**
** SLURM DETAILS:
** PARTITION: physical
** CONSTRAINT: physg5
** CPU: Problem size dependant, left to the user
** MEMORY: 32GB
**
 */

#include "sha512.hh"
#include <cstring>
#include <iostream>
#include <string>
#include <sys/time.h>

using namespace std;

static std::string getMinimumPenalties(std::string *genes, int k, int pxy,
                                       int pgap, int *penalties);
static int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                             int *xans, int *yans);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current wallclock time, for performance measurement
uint64_t GetTimeStamp() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

int main(int argc, char **argv) {
  int misMatchPenalty;
  int gapPenalty;
  int k;
  std::cin >> misMatchPenalty;
  std::cin >> gapPenalty;
  std::cin >> k;
  std::string *genes = new std::string[k];
  for (int i = 0; i < k; i++)
    std::cin >> genes[i];

  int numPairs = k * (k - 1) / 2;

  int *penalties = new int[numPairs];

  uint64_t start = GetTimeStamp();

  // return all the penalties and the hash of all allignments
  std::string alignmentHash =
      getMinimumPenalties(genes, k, misMatchPenalty, gapPenalty, penalties);

  // print the time taken to do the computation
  printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

  // print the alginment hash
  std::cout << alignmentHash << std::endl;

  for (int i = 0; i < numPairs; i++) {
    std::cout << penalties[i] << " ";
  }
  std::cout << std::endl;
  delete[] genes;
  delete[] penalties;
  return 0;
}

[[gnu::always_inline]] [[gnu::hot]] inline int min3(int a, int b, int c) {
  // Return the minimum of 3 variables
  return a < b ? (a < c ? a : c) : (b < c ? b : c);
}

[[gnu::cold]] static int **new2dSingle(int width, int height) {
  int **dp = new int *[width];
  size_t size = width;
  size *= height;
  // not catching malloc error
  // as it is not reliable due to overcommitting etc.
  int *dp0 = new int[size];
  dp[0] = dp0;

  for (int i = 1; i < width; i++) {
    dp[i] = i * height + dp0;
  }
  return dp;
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
[[gnu::always_inline]] [[gnu::hot]] inline int **new2d(int width, int height) {
  int **dp = new int *[width];

  // #pragma omp parallel for default(none) shared(dp, width)
  // for(int i=0; i < width; i++) {
  //     dp[0] = 0;
  // }
  size_t size = width;
  size *= height;
  int *dp0 = new int[size];

  #pragma omp parallel default(none) shared(dp, dp0, size, width)
  {
    #pragma omp for nowait
    for (int i = 0; i < width; i++) {
      dp[0] = 0;
    }
    #pragma omp for
    for (size_t i = 0; i < size; i++) {
      dp0[i] = 0;
    }
  }

  dp[0] = dp0;
  // #pragma omp parallel for default(none) shared(dp0, size)
  // for(size_t i=0; i < size; i++) {
  //     dp0[i] = 0;
  // }

  for (int i = 1; i < width; i++)
    dp[i] = dp[i - 1] + height;

  return dp;
}

static std::string getMinimumPenalties(std::string *genes, int k, int pxy,
                                       int pgap, int *penalties) {
  int probNum = 0;
  std::string alignmentHash = "";
  for (int i = 1; i < k; i++) {
    for (int j = 0; j < i; j++) {
      std::string gene1 = genes[i];
      std::string gene2 = genes[j];
      int m = gene1.length(); // length of gene1
      int n = gene2.length(); // length of gene2
      int l = m + n;
      int *xans = new int[l + 1];
      int *yans = new int[l + 1];
      penalties[probNum] =
          getMinimumPenalty(gene1, gene2, pxy, pgap, xans, yans);
      // Since we have assumed the answer to be n+m long,
      // we need to remove the extra gaps in the starting
      // id represents the index from which the arrays
      // xans, yans are useful
      int id = 1;
      int a;
      for (a = l; a >= 1; a--) {
        if ((char)yans[a] == '_' && (char)xans[a] == '_') {
          id = a + 1;
          break;
        }
      }
      std::string align1 = "";
      std::string align2 = "";
      for (a = id; a <= l; a++) {
        align1.append(1, (char)xans[a]);
      }
      for (a = id; a <= l; a++) {
        align2.append(1, (char)yans[a]);
      }
      std::string align1hash = sw::sha512::calculate(align1);
      std::string align2hash = sw::sha512::calculate(align2);
      std::string problemhash =
          sw::sha512::calculate(align1hash.append(align2hash));
      alignmentHash = sw::sha512::calculate(alignmentHash.append(problemhash));

      // Uncomment for testing purposes
      // std::cout << penalties[probNum] << std::endl;
      // std::cout << align1 << std::endl;
      // std::cout << align2 << std::endl;
      // std::cout << std::endl;

      probNum++;
      delete[] xans;
      delete[] yans;
    }
  }
  return alignmentHash;
}

[[gnu::always_inline]] [[gnu::hot]] inline void
diagonalise(int **dp, int width, int height, int di, int dj,
            const std::string &x, const std::string &y, int pxy, int pgap) {
  int new_width = width / di;
  int new_height = height / dj;

  int diagonals =
      new_width + new_height + ((width % di + height % dj) > 0 ? 1 : 0);

  int outer_i = 0;
  int outer_j = 0;

  int diff_i = 0;
  int diff_j = 0;
  int diag_i = 0;
  int diag_j = 0;
  int length = 0;

  {
    for (int d = 0; d < diagonals; d++) {
      diff_i = outer_i;
      diff_j = height - outer_j - 1;

      diag_i = 1 + (diff_i / di);
      diag_j = 1 + (diff_j / dj);
      length = std::min(diag_i, diag_j);

      #pragma omp parallel for schedule(dynamic) default(none) shared(x, y, dp, di, dj, pgap, pxy, length, width, height, outer_i, outer_j)
      for (int tile = 0; tile < length; tile++) {
        int inner_i = std::max(1, outer_i - (tile * di));
        int inner_j = std::max(1, outer_j + (tile * dj));

        int imax = std::min(inner_i + di, width);
        int jmax = std::min(inner_j + dj, height);
        for (int j = inner_j; j < jmax; j++) {
          for (int i = inner_i; i < imax; i++) {
            if (x[j - 1] == y[i - 1]) {
              dp[j][i] = dp[j - 1][i - 1];
            } else {
              dp[j][i] = min3(dp[j - 1][i - 1] + pxy, dp[j - 1][i] + pgap,
                              dp[j][i - 1] + pgap);
            }
          }
        }
      }

      if (outer_i + di < width) {
        outer_i += di;
      } else {
        outer_j += dj;
      }
    }
  }
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
static int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                             int *__restrict__ xans, int *__restrict__ yans) {

  int i, j; // intialising variables

  int m = x.length(); // length of gene1
  int n = y.length(); // length of gene2

  // table for storing optimal substructure answers
  int **dp;

  if (std::min(m, n) > 1000) {
    dp = new2d(m + 1, n + 1);
    #pragma omp parallel default(none) shared(dp, pgap, m, n)
    {

      #pragma omp for nowait
      for (int i = 0; i <= m; i++) {
        dp[i][0] = i * pgap;
      }
      #pragma omp for
      for (int i = 0; i <= n; i++) {
        dp[0][i] = i * pgap;
      }
    }

    // intialising the table

    diagonalise(dp, n + 1, m + 1, (n + 1) / (8 + 5) + 1, (m + 1) / (8 + 5) + 1,
                x, y, pxy, pgap);
  } else {
    dp = new2dSingle(m + 1, n + 1);

    for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
        if (x[i - 1] == y[j - 1]) {
          dp[i][j] = dp[i - 1][j - 1];
        } else {
          dp[i][j] = min3(dp[i - 1][j - 1] + pxy, dp[i - 1][j] + pgap,
                          dp[i][j - 1] + pgap);
        }
      }
    }
  }

  // calcuting the minimum penalty

  // Reconstructing the solution
  int l = n + m; // maximum possible length

  i = m;
  j = n;

  int xpos = l;
  int ypos = l;

  while (!(i == 0 || j == 0)) {
    if (x[i - 1] == y[j - 1]) {
      xans[xpos--] = (int)x[i - 1];
      yans[ypos--] = (int)y[j - 1];
      i--;
      j--;
    } else if (dp[i - 1][j - 1] + pxy == dp[i][j]) {
      xans[xpos--] = (int)x[i - 1];
      yans[ypos--] = (int)y[j - 1];
      i--;
      j--;
    } else if (dp[i - 1][j] + pgap == dp[i][j]) {
      xans[xpos--] = (int)x[i - 1];
      yans[ypos--] = (int)'_';
      i--;
    } else if (dp[i][j - 1] + pgap == dp[i][j]) {
      xans[xpos--] = (int)'_';
      yans[ypos--] = (int)y[j - 1];
      j--;
    }
  }
  while (xpos > 0) {
    if (i > 0)
      xans[xpos--] = (int)x[--i];
    else
      xans[xpos--] = (int)'_';
  }
  while (ypos > 0) {
    if (j > 0)
      yans[ypos--] = (int)y[--j];
    else
      yans[ypos--] = (int)'_';
  }

  int ret = dp[m][n];

  delete[] dp[0];
  delete[] dp;

  return ret;
}
