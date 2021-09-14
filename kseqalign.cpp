// CPP program to solve the sequence alignment
// problem. Adapted from
// https://www.geeksforgeeks.org/sequence-alignment-problem/
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <set>
#include <immintrin.h>
#include <numa.h>
#include "sha512.hh"

using namespace std;

static std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
                                int *penalties);

static int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                      int *xans, int *yans);



/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current wallclock time, for performance measurement
[[gnu::cold]]
static uint64_t GetTimeStamp() {
  return chrono::duration_cast<chrono::microseconds>(
             chrono::system_clock::now().time_since_epoch())
      .count();
}

[[gnu::cold]]
int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
  int misMatchPenalty;
  int gapPenalty;
  int k;
  FILE *fp = freopen("mseq.dat", "r", stdin);
  // #ifdef DEBUG_BUILD
  // freopen("mseq-simple.dat", "r", stdin);
  // #endif
  std::cin >> misMatchPenalty;
  std::cin >> gapPenalty;
  std::cin >> k;

  auto genes = new std::string[k];

  for (int i = 0; i < k; i++) {
    std::cin >> genes[i];
  }

  int numPairs = k * (k - 1) / 2;

  int *penalties = new int[numPairs];

  uint64_t start = GetTimeStamp();

  std::string alignmentHash = getMinimumPenalties(genes, k, misMatchPenalty, gapPenalty, penalties);
  std::cout << GetTimeStamp() - start << " us\n";
  std::cout << "ALIGNMENT: " << alignmentHash << std::endl;

  // return all the penalties and the hash of all allignments

  for (int i = 0; i < numPairs; i++) {
    std::cout << penalties[i] << " ";
  }
  std::cout << std::endl;

  fclose(fp);
  delete[] genes;
  delete[] penalties;
  return 0;
}

[[gnu::always_inline]] [[gnu::hot]]
inline int min3(int a, int b, int c) {
  // Return the minimum of 3 variables
  return a < b ? (a < c ? a : c) : (b < c ? b : c);
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
static int **new2d(int width, int height) {
  int **dp = new int *[width];
  size_t size = width;
  size *= height;
  // not catching malloc error
  // as it is not reliable due to overcommitting etc.
  int *dp0 = new int[size];
  dp[0] = dp0;

  #pragma omp parallel default(none) shared(dp0, size, dp, width, height)
  {

    #pragma omp for nowait
    for(size_t i=0; i < size; i++) {
      dp0[i] = 0;
    }

    #pragma omp for
    for (int i = 1; i < width; i++) {
      dp[i] = i*height + dp0;
    }
  }
  return dp;
}

[[gnu::cold]]
static std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
                                int *penalties) {
  std::string alignmentHash;

  int probNum = 0;

  {
    for (int i = 1; i < k; i++) {
      for (int j = 0; j < i; j++) {
        auto gene1 = genes[i];
        auto gene2 = genes[j];
        auto m = gene1.length(); // length of gene1
        auto n = gene2.length(); // length of gene2
        auto l = m + n;
        int *xans = new int[l + 1];
        int *yans = new int[l + 1];

        penalties[probNum] =
            getMinimumPenalty(gene1, gene2, pxy, pgap, xans, yans);
        // Since we have assumed the answer to be n+m long,
        // we need to remove the extra gaps in the starting
        // id represents the index from which the arrays
        // xans, yans are useful
        unsigned long id = 1;
        unsigned long a;
        for (a = l; a >= 1; a--) {
          if ((char)yans[a] == '_' && (char)xans[a] == '_') {
            id = a + 1;
            break;
          }
        }
        std::string align1;
        std::string align2;
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
        {

          alignmentHash =
              sw::sha512::calculate(alignmentHash.append(problemhash));

          // Uncomment for testing purposes
          //					std::cout << penalties[probNum] <<
          //std::endl; 					std::cout << align1 << std::endl; 					std::cout << align2 <<
          //std::endl; 					std::cout << std::endl;

          probNum++;
        }

        delete[] xans;
        delete[] yans;
      }
    }
  }

  return alignmentHash;
}


[[gnu::always_inline]] [[gnu::hot]]
inline void diagonalise(int **dp, int width, int height, int di, int dj, const std::string &x, const std::string &y, int pxy, int pgap) {
  int new_width = width / di;
  int new_height = height / dj;

  int diagonals = new_width + new_height + ((width%di  + height%dj) > 0 ? 1 : 0);

  int outer_i = 0;
  int outer_j = 0;

  int diff_i = 0;
  int diff_j = 0;
  int diag_i = 0;
  int diag_j = 0;
  int length = 0;

  {
    for(int d = 0; d < diagonals; d++) {
      diff_i = outer_i;
      diff_j = height - outer_j -1;

      diag_i = 1 + (diff_i / di);
      diag_j = 1 + (diff_j / dj);
      length = std::min(diag_i, diag_j);

      #pragma omp parallel for schedule(dynamic) default(none) shared(x, y, dp, di, dj, pgap, pxy, length, width, height, outer_i, outer_j)
      for(int tile=0; tile < length; tile++) {
        int inner_i = std::max(1, outer_i - (tile * di));
        int inner_j = std::max(1, outer_j + (tile * dj));

        int imax = std::min(inner_i + di, width);
        int jmax = std::min(inner_j + dj, height);
        for(int i = inner_i; i < imax; i++) {
          for(int j = inner_j; j < jmax; j++) {
            if(x[j-1] == y[i-1]) {
              dp[j][i] = dp[j-1][i-1];
            }else {
              dp[j][i] = min3(dp[j-1][i-1] + pxy, dp[j-1][i] + pgap, dp[j][i-1] + pgap);
            }
          }
        }
      }

      if (outer_i + di < width) {
        outer_i += di;
      }else {
        outer_j += dj;
      }
    }
  }
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
[[gnu::hot]]
static int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                      int *__restrict__ xans, int *__restrict__ yans) {


  int m = static_cast<int>(x.length()); // length of gene1
  int n = static_cast<int>(y.length()); // length of gene2

  // table for storing optimal substructure answers
  int **dp = new2d(m + 1, n + 1);
  size_t size = m + 1;
  size *= n + 1;
  memset(dp[0], 0, size);

  #pragma omp parallel default(none) shared(dp, m, n, pgap)
  {
    #pragma omp for simd nowait
    for (int i = 0; i <= m; i++) {
      dp[i][0] = i * pgap;
    }
    #pragma omp for
    for (int i = 0; i <= n; i++) {
      dp[0][i] = i * pgap;
    }
  }
  const int N = m + 1;
  const int M = n + 1;

  diagonalise(dp, M, N, M/(8*2) , N/(8*2) , x, y, pxy, pgap);
	// Reconstructing the solution
	int l = n + m; // maximum possible length

	int i = m; int j = n;

	int xpos = l;
	int ypos = l;

	while ( !(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)x[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)y[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;
  return ret;
}
