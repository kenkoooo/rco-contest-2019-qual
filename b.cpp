#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <queue>
#include <tuple>
#include <vector>

using namespace std;

#define REP(i, n) for (int(i) = 0; (i) < (n); (i)++)

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  if (!v.empty()) {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& p) {
  out << "[" << p.first << ", " << p.second << "]";
  return out;
}

uint64_t xor_shift64(uint64_t* state) {
  uint64_t x = *state;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  *state = x;
  return x;
}

// --- template ---

class UnionFind {
 private:
  struct nodeinfo {
    int par;
    int rank;
    int size;
    nodeinfo(int p) : par(p), rank(0), size(1) {}
  };
  std::vector<nodeinfo> node;

 public:
  UnionFind(int n) : node() {
    node.reserve(n);
    for (int i = 0; i < n; ++i) {
      node.push_back(nodeinfo(i));
    }
  }
  int root(int x) {
    if (node[x].par == x) {
      return x;
    }
    return node[x].par = root(node[x].par);
  }
  void unite(int x, int y) {
    x = root(x);
    y = root(y);
    if (x == y) {
      return;
    }
    if (node[x].rank < node[y].rank) {
      node[x].par = y;
      node[y].size += node[x].size;
      node[x].size = 0;
    } else {
      node[y].par = x;
      node[x].size += node[y].size;
      node[y].size = 0;
      if (node[x].rank == node[y].rank) {
        ++node[x].rank;
      }
    }
  }

  bool is_same_set(int x, int y) { return root(x) == root(y); }

  int size(int x) {
    x = root(x);
    return node[x].size;
  }
};

using State = tuple<double, int, int>;
const int DX[] = {0, 0, 1, -1};
const int DY[] = {1, -1, 0, 0};
const int INF = 1e9;

int N, M;
int A[50][50];

auto backward(vector<vector<double>>& dist, int i, int j) {
  vector<tuple<int, int>> ans;
  while (dist[i][j] > 0) {
    auto cur_min = dist[i][j];
    assert(cur_min < INF);
    tuple<int, int> next = make_tuple(i, j);

    REP(d, 4) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;

      if (dist[ni][nj] < cur_min) {
        next = make_tuple(ni, nj);
        cur_min = dist[ni][nj];
      }
    }

    int ni, nj;
    tie(ni, nj) = next;

    assert(ni != i || nj != j);
    dist[i][j] = INF;
    ans.emplace_back(ni, nj);
    i = ni;
    j = nj;
  }

  return ans;
}

vector<tuple<int, int>> dijkstra(int si, int sj, UnionFind& uf) {
  const int start_size = uf.size(si * N + sj);

  priority_queue<State, vector<State>, greater<State>> q;
  vector<vector<double>> dist(N, vector<double>(N, INF));
  dist[si][sj] = 0;
  q.emplace(0, si, sj);
  while (q.size() != 0) {
    int cost, i, j;
    tie(cost, i, j) = q.top();
    q.pop();

    if (A[i][j] == 9 && uf.size(i * N + j) + start_size >= 9) {
      return backward(dist, i, j);
    }

    REP(d, 4) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;
      double w = 9.00001 - A[ni][nj];
      if (dist[ni][nj] > dist[i][j] + w) {
        dist[ni][nj] = dist[i][j] + w;
        q.emplace(dist[ni][nj], ni, nj);
      }
    }
  }

  return {};
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  vector<tuple<int, int, int>> ans;

  cin >> N >> M;
  REP(i, N) REP(j, N) {
    cin >> A[i][j];
    auto cost = 9 - A[i][j];
    if (cost <= 3 && cost + ans.size() <= M) {
      A[i][j] += cost;
      REP(z, cost) ans.emplace_back(1, i, j);
    }
  }

  UnionFind uf(N * N);
  REP(i, N) REP(j, N) REP(d, 4) {
    int ni = i + DX[d];
    int nj = j + DY[d];
    if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;
    if (A[i][j] != A[ni][nj]) continue;
    int from = i * N + j;
    int to = ni * N + nj;
    uf.unite(from, to);
  }

  vector<bool> started(N * N, false);
  vector<vector<tuple<int, int>>> paths;
  vector<tuple<double, int, int>> costs;
  REP(i, N) REP(j, N) {
    if (A[i][j] != 9) continue;
    int x = i * N + j;
    x = uf.root(x);

    if (started[x]) continue;
    started[x] = true;

    int size = uf.size(i * N + j);
    if (size >= 9) continue;

    auto path = dijkstra(i, j, uf);
    int cost = 0;
    for (const auto& t : path) {
      int pi, pj;
      tie(pi, pj) = t;
      cost += 9 - A[pi][pj];
    }

    double cost_per_size = cost;
    cost_per_size /= size;
    costs.emplace_back(cost_per_size, cost, paths.size());
    paths.emplace_back(path);
  }

  sort(costs.begin(), costs.end());
  for (const auto& t : costs) {
    double cps;
    int cost;
    int index;
    tie(cps, cost, index) = t;
    const auto& path = paths[index];
    cost = 0;
    for (const auto& s : path) {
      int i, j;
      tie(i, j) = s;
      int turn = 9 - A[i][j];
      cost += turn;
    }

    if (cost + ans.size() >= M - 200) continue;
    for (const auto& s : path) {
      int i, j;
      tie(i, j) = s;
      int turn = 9 - A[i][j];
      REP(z, turn) ans.emplace_back(1, i, j);
      A[i][j] = 9;
    }
  }

  REP(i, N) REP(j, N) REP(d, 4) {
    int ni = i + DX[d];
    int nj = j + DY[d];
    if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;
    if (A[i][j] != A[ni][nj]) continue;
    int from = i * N + j;
    int to = ni * N + nj;
    uf.unite(from, to);
  }

  vector<tuple<int, int, int>> candidates;
  vector<bool> done(N * N, false);
  REP(i, N) REP(j, N) {
    int x = i * N + j;
    x = uf.root(x);
    if (done[x]) continue;
    done[x] = true;
    if (uf.size(x) < A[i][j]) continue;
    int value = uf.size(x) * A[i][j];
    candidates.emplace_back(value, i, j);
  }

  sort(candidates.rbegin(), candidates.rend());
  for (const auto& t : candidates) {
    int v, i, j;
    tie(v, i, j) = t;
    ans.emplace_back(2, i, j);
    if (ans.size() == M) break;
  }

  for (const auto& t : ans) {
    int o, i, j;
    tie(o, i, j) = t;
    cout << o << " " << i << " " << j << "\n";
  }
}
