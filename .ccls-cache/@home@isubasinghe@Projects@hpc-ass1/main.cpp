#include <iostream>
#include <chrono>
#include <cstdio>
#include <new>
#include <memory>
using namespace std;

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d(int width, int height) {
  int **dp = new int *[width];
  size_t size = width;
  size *= height;
  // not catching malloc error
  // as it is not reliable due to overcommitting etc.
  int *dp0 = new int[size];
  dp[0] = dp0;
  for (int i = 1; i < width; i++) {
    dp[i] = dp[i - 1] + height;
  }
  return dp;
}


void setMatrix(int **dp, int width, int height, int value) {
    for(int i =0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            dp[i][j] = value;
        }
    }
}


void diagonalise(int **dp, int width, int height, int di, int dj) {
  int new_width = width / di;
  int new_height = height / dj;

  int diagonals = new_width + new_height + (width%di  + height%dj > 0 ? 1 : 0);

  int outer_i = 0;
  int outer_j = 0;

  for(int d = 0; d < diagonals; d++) {
    std::cout << "DIAG_START\n";
    int diff_i = outer_i;
    int diff_j = height - outer_j -1;

    int diag_i = 1 + (diff_i / di);
    int diag_j = 1 + (diff_j / dj);
    int length = std::min(diag_i, diag_j);

    for(int tile=0; tile < length; tile++) {
      std::cout << "TILE_START" << std::endl;
      int inner_i = std::max(1, outer_i - (tile * di));
      int inner_j = std::max(1, outer_j + (tile * dj));

      int imax = std::min(inner_i + di, width);
      int jmax = std::min(inner_j + dj, height);
      for(int i = inner_i; i < imax; i++) {
        for(int j = inner_j; j < jmax; j++) {
          std::cout << "(" << i << ", " << j << ")" << std::endl;
          dp[i][j] = tile;
        }
      }
      std::cout << "TILE_END" << std::endl;
    }

    outer_i += di;
    if(outer_i >= width) {
      outer_i = width - di;
      outer_j += dj;
    }
  }
}

uint64_t GetTimeStamp() {
  return chrono::duration_cast<chrono::microseconds>(
             chrono::system_clock::now().time_since_epoch())
      .count();
}
int main(int argc, char *argv[]) {
  if(argc != 3) {
    return 0;
  }
  int pgap = atoi(argv[1]);
  int sz = atoi(argv[2]);

  int *arr = static_cast<int *>(aligned_alloc(32, sz));
  arr[0] = 1;

  auto start = GetTimeStamp();

  // #pragma omp parallel for simd
  for(int i = 1; i < sz; i++) {
    arr[i] = i * pgap + (arr[i-1]);
  }
  std::cout << GetTimeStamp() - start << std::endl;
  std::cout << arr[sz-1] << std::endl;
  return 0;
}
