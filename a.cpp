#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <tuple>
#include <vector>

using namespace std;

#define REP(i, n) for (size_t(i) = 0; (i) < (n); (i)++)

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

const double INF = 1e18;

size_t N;
double X[200], Y[200];
double dist[200][200];

double calc_distance(double x1, double y1, double x2, double y2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  return sqrt(dx * dx + dy * dy);
}

double improve_2opt(vector<size_t>& path) {
  double cur_improvement = 0.0;
  size_t index_from = 0;
  size_t index_to = 0;

  REP(i, N) {
    auto from1 = path[i];
    auto to1 = path[i + 1];
    auto cost1 = dist[from1][to1];
    for (size_t j = i + 1; j < N; j++) {
      auto from2 = path[j];
      auto to2 = path[j + 1];
      auto cost2 = dist[from2][to2];

      auto prev_cost = cost1 + cost2;
      auto next_cost = dist[from1][from2] + dist[to1][to2];
      if (prev_cost > next_cost && prev_cost - next_cost > cur_improvement) {
        cur_improvement = prev_cost - next_cost;
        index_from = i + 1;
        index_to = j + 1;
      }
    }
  }
  if (index_from == index_to) {
    return 0.0;
  }

  reverse(path.begin() + index_from, path.begin() + index_to);
  return cur_improvement;
}

tuple<vector<size_t>, double> greedy() {
  double total_cost = 0.0;
  vector<size_t> res;
  vector<bool> vis(N, false);

  size_t cur = 0;
  vis[cur] = true;
  res.push_back(cur);

  REP(i, (N - 1)) {
    double cur_min = INF;
    size_t cur_next = N;
    REP(next, N) {
      if (vis[next]) {
        continue;
      }
      if (dist[cur][next] < cur_min) {
        cur_next = next;
        cur_min = dist[cur][next];
      }
    }

    total_cost += dist[cur][cur_next];
    res.push_back(cur_next);
    vis[cur_next] = true;
    cur = cur_next;
  }

  res.push_back(0);
  total_cost += dist[cur][0];

  return make_tuple(res, total_cost);
}

tuple<vector<size_t>, double> solve(double average) {
  REP(i, N) REP(j, N) {
    auto d = calc_distance(X[i], Y[i], X[j], Y[j]);
    dist[i][j] = (d - average) * (d - average);
    dist[j][i] = dist[i][j];
  }

  vector<size_t> path;
  double total_cost;
  tie(path, total_cost) = greedy();

  REP(i, 300) {
    double improvement = improve_2opt(path);
    if (improvement < 1e-9) {
      break;
    }
    assert(total_cost >= improvement);
    total_cost -= improvement;
  }
  return make_tuple(path, total_cost);
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  auto start = chrono::system_clock::now();

  cin >> N;
  REP(i, N) { cin >> X[i] >> Y[i]; }

  uint64_t random_state = 71;

  vector<vector<size_t>> paths;
  vector<tuple<double, size_t>> costs;
  while (true) {
    double r = xor_shift64(&random_state);
    r /= numeric_limits<uint64_t>::max();

    double average = 200.0 + r * 200.0;
    vector<size_t> path;
    double cost;
    tie(path, cost) = solve(average);
    costs.emplace_back(cost, paths.size());
    paths.emplace_back(path);

    auto ms = chrono::duration_cast<chrono::milliseconds>(
                  chrono::system_clock::now() - start)
                  .count();
    if (ms > 1900) {
      break;
    }
  }
  sort(costs.begin(), costs.end());

  size_t min_cost_index = get<1>(costs[0]);

  const auto& path = paths[min_cost_index];
  REP(i, N) { cout << path[i] << endl; }
}
