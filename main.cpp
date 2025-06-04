#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;

constexpr double pow2(double x) {
    return x * x;
};
constexpr int roundint(double x) {
    return static_cast<int>(round(x));
}

template<typename F, typename S>
inline ostream& operator<<(ostream& os, const pair<F, S>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}
template<typename T>
inline ostream& operator<<(ostream& os, const vector<T>& vec) {
    os << "[";
    if (!vec.empty()) {
        os << vec[0];
        for (size_t i = 1; i < vec.size(); ++i) {
            os << ", " << vec[i];
        }
    }
    os << "]";
    return os;
}

struct Timer {
    const chrono::system_clock::time_point start;
    Timer() : start(chrono::system_clock::now()) {}
    double get() const {
        chrono::system_clock::time_point end = chrono::system_clock::now();
        chrono::duration<double> sec = end - start;
        return sec.count();
    }
};

struct Color {
    double C, M, Y;
    Color(double C_ = 0, double M_ = 0, double Y_ = 0) : C(C_), M(M_), Y(Y_) {}

    double distance(const Color& other) const {
        return sqrt(pow2(C - other.C) + pow2(M - other.M) + pow2(Y - other.Y));
    }

    friend ostream& operator<<(ostream& os, const Color& color) {
        return os << "C: " << color.C << ", M: " << color.M << ", Y: " << color.Y;
    }
    friend istream& operator>>(istream& is, Color& color) {
        return is >> color.C >> color.M >> color.Y;
    }
};

inline Color mixed_color(const Color& c1, double r1, const Color& c2, double r2) {
    double total_ratio = r1 + r2;
    return Color(
        (c1.C * r1 + c2.C * r2) / total_ratio,
        (c1.M * r1 + c2.M * r2) / total_ratio,
        (c1.Y * r1 + c2.Y * r2) / total_ratio
    );
}
inline Color mixed_color(const Color& c1, int r1, const Color& c2, int r2) {
    int total_ratio = r1 + r2;
    return Color(
        (c1.C * r1 + c2.C * r2) / total_ratio,
        (c1.M * r1 + c2.M * r2) / total_ratio,
        (c1.Y * r1 + c2.Y * r2) / total_ratio
    );
}

struct Point3D {
    double x, y, z;
    vector<pair<int, int>> ratio;

    Point3D(Color color = Color(), vector<pair<int, int>> ratio_ = vector<pair<int, int>>())
        : x(color.C), y(color.M), z(color.Y), ratio(move(ratio_)) {}

    double getCoord(int axis) const {
        if (axis == 0) return x;
        if (axis == 1) return y;
        return z;
    }
};

constexpr double distanceSquared(const Point3D& p1, const Point3D& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}

struct KdNode {
    Point3D point;
    unique_ptr<KdNode> left;
    unique_ptr<KdNode> right;
    int axis; // 分割軸 (0 for x, 1 for y, 2 for z)

    KdNode(const Point3D& pt, int ax) : point(pt), axis(ax) {}
};

class KdTree {
public:
    KdTree() {
        random_device rd;
        gen = mt19937(rd());
    }

    void build(vector<Point3D>& points) {
        root = buildRecursive(points, 0, points.size() - 1, 0);
    }

    Point3D findClosestPoint(const Point3D& target) const {
        Point3D bestPoint = root->point;
        double bestDistSq = distanceSquared(target, root->point);
        findClosestRecursive(root.get(), target, bestPoint, bestDistSq);
        return bestPoint;
    }

private:
    unique_ptr<KdNode> root;
    mt19937 gen;

    unique_ptr<KdNode> buildRecursive(vector<Point3D>& points, int start, int end, int depth) {
        if (start > end) {
            return nullptr;
        }

        int axis = depth % 3;

        if (start == end) { // 葉ノード
            return make_unique<KdNode>(points[start], axis);
        }
        
        uniform_int_distribution<int> distrib(start, end);
        int pivot_initial_idx = distrib(gen);

        // ピボットを範囲の最後に移動
        swap(points[pivot_initial_idx], points[end]);
        const Point3D& pivot_element = points[end];
        double pivot_coord_val = pivot_element.getCoord(axis);

        // パーティション処理: points[start...end-1] をピボット基準で分割
        // store_idx は「ピボットより小さい要素グループ」の末尾の次を指す
        int store_idx = start; 
        for (int i = start; i < end; ++i) {
            if (points[i].getCoord(axis) < pivot_coord_val) {
                swap(points[i], points[store_idx]);
                store_idx++;
            }
        }
        
        // ピボットを最終的な位置（store_idx）に移動
        swap(points[store_idx], points[end]);
        int final_pivot_idx = store_idx;

        auto node = make_unique<KdNode>(points[final_pivot_idx], axis);
        node->left = buildRecursive(points, start, final_pivot_idx - 1, depth + 1);
        node->right = buildRecursive(points, final_pivot_idx + 1, end, depth + 1);

        return node;
    }

    void findClosestRecursive(const KdNode* node, const Point3D& target, Point3D& bestPoint, double& bestDistSq) const {
        if (!node) {
            return;
        }

        double currentDistSq = distanceSquared(target, node->point);
        if (currentDistSq < bestDistSq) {
            bestDistSq = currentDistSq;
            bestPoint = node->point;
        }

        int axis = node->axis;
        double targetCoord = target.getCoord(axis);
        double nodeCoord = node->point.getCoord(axis);

        const unique_ptr<KdNode>* firstChild = (targetCoord < nodeCoord) ? &(node->left) : &(node->right);
        const unique_ptr<KdNode>* secondChild = (targetCoord < nodeCoord) ? &(node->right) : &(node->left);

        if (firstChild->get()) {
            findClosestRecursive(firstChild->get(), target, bestPoint, bestDistSq);
        }

        // もう一方のサブツリーも探索する必要があるかチェック
        // (分割平面までの距離の2乗が現在の最小距離の2乗より小さい場合)
        if (secondChild->get()) {
            double distToPlaneSq = (targetCoord - nodeCoord) * (targetCoord - nodeCoord);
            if (distToPlaneSq < bestDistSq) {
                findClosestRecursive(secondChild->get(), target, bestPoint, bestDistSq);
            }
        }
    }
};

struct Input {
    int N, K, H, T, D;
    vector<Color> own;
    vector<Color> target;

    friend istream& operator>>(istream& is, Input& input) {
        is >> input.N >> input.K >> input.H >> input.T >> input.D;
        input.own.resize(input.K);
        for (Color& color : input.own) {
            is >> color;
        }
        input.target.resize(input.H);
        for (Color& color : input.target) {
            is >> color;
        }
        return is;
    }
};

struct Output {
    vector<vector<bool>> v, h;
    vector<vector<int>> operation;

    Output(int N) : v(N, vector<bool>(N - 1)), h(N - 1, vector<bool>(N)) {}

    void add(int i, int j, int k) {
        operation.push_back({1, i, j, k});
    }

    void submission(int i, int j) {
        operation.push_back({2, i, j});
    }

    void remove(int i, int j) {
        operation.push_back({3, i, j});
    }

    void partition(int i1, int j1, int i2, int j2) {
        operation.push_back({4, i1, j1, i2, j2});
    }

    friend ostream& operator<<(ostream& os, const Output& output) {
        for (const auto& row : output.v) {
            for (int i = 0; i < (int)row.size(); ++i) {
                os << row[i] << (i < (int)row.size() - 1 ? ' ' : '\n');
            }
        }
        for (const auto& row : output.h) {
            for (int i = 0; i < (int)row.size(); ++i) {
                os << row[i] << (i < (int)row.size() - 1 ? ' ' : '\n');
            }
        }
        for (const auto& op : output.operation) {
            for (int i = 0; i < (int)op.size(); ++i) {
                os << op[i] << (i < (int)op.size() - 1 ? ' ' : '\n');
            }
        }
        return os;
    }
};

struct ContainerStatus {
    static constexpr int LIMIT = 10;

    int turn;
    int score;
    vector<int> result;
    vector<int> capacity;

    ContainerStatus() : turn(0), score(0) {}
    ContainerStatus(const Input &input) : turn(0), score(0), capacity(input.K, 0) {}

    void add(const Input &input, const vector<vector<vector<pair<int, int>>>> &closestRatio, int index) {
        result.push_back(index);
        int targetIndex = result.size() - 1;
        int ratioSize = closestRatio[targetIndex][index].size();
        if (ratioSize == 1) {
            turn += 2;
            score += roundint(input.target[targetIndex].distance(input.own[closestRatio[targetIndex][index][0].first]) * 10000);
            return;
        }
        for (int i = 0; i < ratioSize; ++i) {
            int colorIndex = closestRatio[targetIndex][index][i].first;
            int ratio = closestRatio[targetIndex][index][i].second;
            if (capacity[colorIndex] < ratio) {
                if (capacity[colorIndex] == 0) {
                    turn += 3;
                } else {
                    turn += 5;
                }
                capacity[colorIndex] += LIMIT - ratio;
            } else {
                capacity[colorIndex] -= ratio;
                turn += 2;
            }
        }
        turn += 1;

        Color mixedColor;
        int totalRatio = 0;
        for (int i = 0; i < ratioSize; ++i) {
            int colorIndex = closestRatio[targetIndex][index][i].first;
            int ratio = closestRatio[targetIndex][index][i].second;
            mixedColor = mixed_color(mixedColor, totalRatio, input.own[colorIndex], ratio);
            totalRatio += ratio;
        }
        score += roundint(input.target[result.size() - 1].distance(mixedColor) * 10000);
    }

    int getScore(const Input &input, bool is_final) const {
        if (!is_final) {
            return score;
        }
        else {
            return score + accumulate(capacity.begin(), capacity.end(), 0) / 10 * input.D;
        }
    }

    Output toOutput(const Input &input, const vector<vector<vector<pair<int, int>>>> &closestRatio) const {
        Output output(input.N);
        for (int i = 0; i < input.N - 1; ++i) {
            for (int j = 0; j < min(input.N - 1, input.K); ++j) {
                output.v[i][j] = true;
            }
        }
        for (int i = 0; i < input.K; ++i) {
            output.h[LIMIT - 1][i] = true;
        }

        vector<int> now_capacity(input.K, 0);
        for (int ti = 0; ti < (int)result.size(); ++ti) {
            int index = result[ti];
            const auto& ratio = closestRatio[ti][index];
            if ((int)ratio.size() == 1) {
                output.add(input.N - 1, input.N - 1, ratio[0].first);
                output.submission(input.N - 1, input.N - 1);
                continue;
            }

            for (int i = 0; i < (int)ratio.size(); ++i) {
                int color_index = ratio[i].first;
                int ratio_value = ratio[i].second;
                if (now_capacity[color_index] < ratio_value) {
                    if (now_capacity[color_index] == 0) {
                        output.add(LIMIT - 1, color_index, color_index);
                    } else {
                        output.partition(now_capacity[color_index] + LIMIT - 1, color_index, now_capacity[color_index] + LIMIT, color_index);
                        output.add(now_capacity[color_index] + LIMIT - 1, color_index, color_index);
                        output.partition(now_capacity[color_index] - 1, color_index, now_capacity[color_index], color_index);
                    }
                    output.partition(now_capacity[color_index] + LIMIT - ratio_value - 1, color_index, now_capacity[color_index] + LIMIT - ratio_value, color_index);
                } else {
                    if (now_capacity[color_index] - ratio_value != 0) {
                        output.partition(now_capacity[color_index] - ratio_value - 1, color_index, now_capacity[color_index] - ratio_value, color_index);
                    }
                }
            }

            for (int i = 0; i < (int)ratio.size(); ++i) {
                int color_index = ratio[i].first;
                int ratio_value = ratio[i].second;
                if (now_capacity[color_index] < ratio_value) {
                    output.partition(now_capacity[color_index] + LIMIT - 1, color_index, now_capacity[color_index] + LIMIT, color_index);
                } else {
                    output.partition(now_capacity[color_index] - 1, color_index, now_capacity[color_index], color_index);
                }
            }
            output.submission(input.N - 1, input.N - 1);
            for (int i = 0; i < (int)ratio.size(); ++i) {
                int color_index = ratio[i].first;
                int ratio_value = ratio[i].second;
                if (now_capacity[color_index] < ratio_value) {
                    now_capacity[color_index] += LIMIT - ratio_value;
                } else {
                    if (now_capacity[color_index] - ratio_value == 0) {
                        output.partition(LIMIT - 1, color_index, LIMIT, color_index);
                    }
                    now_capacity[color_index] -= ratio_value;
                }
            }
        }

        return output;
    }
};

struct Solution {
    static constexpr int DEPTHS[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 7, 6, 6, 5, 5, 5};
    static constexpr int MAX_WIDTH = 400;

    Input input;
    vector<vector<vector<pair<int, int>>>> closestRatio;

    Solution(const Input& input_) : input(input_), closestRatio(input_.H) {}

    Output solve() {
        get_closestRatio();

        ContainerStatus containerStatus = beam_search(MAX_WIDTH);

        return containerStatus.toOutput(input, closestRatio);
    }

    void get_closestRatio() {
        Timer timer;
        int max_depth = DEPTHS[input.K - 1];
        vector<vector<Point3D>> points(max_depth);
        vector<pair<int, int>> ratio;
        auto dfs = [&](auto&& self, int pt, int sum, Color color, int depth) -> void {
            if (sum == 10) {
                points[depth - 1].emplace_back(color, ratio);
                return;
            }
            if (depth >= max_depth) {
                return;
            }
            if (pt >= input.K) {
                return;
            }
            for (int npt = pt + 1; npt < input.K; ++npt) {
                for (int k = 1; k <= 10 - sum; ++k) {
                    ratio.emplace_back(npt, k);
                    Color mixed = mixed_color(input.own[npt], k, color, sum);
                    self(self, npt, sum + k, mixed, depth + 1);
                    ratio.pop_back();
                }
            }
        };
        dfs(dfs, -1, 0, Color(), 0);

        vector<KdTree> kdTrees(max_depth);
        for (int i = 0; i < max_depth; ++i) {
            kdTrees[i].build(points[i]);
        }

        for (int i = 0; i < max_depth; ++i) {
            for (int j = 0; j < input.H; ++j) {
                Point3D targetPoint(input.target[j]);
                closestRatio[j].push_back(kdTrees[i].findClosestPoint(targetPoint).ratio);
            }
        }
    }

    ContainerStatus beam_search(int beam_width) {
        vector<ContainerStatus> nodes;
        nodes.push_back(ContainerStatus(input));
        for (int i = 0; i < input.H; ++i) {
            bool is_final = (i >= input.H - input.K);
            unordered_map<int, int> next_nodes_map;
            vector<ContainerStatus> next_nodes;

            Point3D targetPoint(input.target[i]);
            for (const auto& node : nodes) {
                for (int j = 0; j < (int)closestRatio[i].size(); ++j) {
                    ContainerStatus new_node = node;
                    new_node.add(input, closestRatio, j);
                    if (next_nodes_map.count(new_node.turn) == 0 
                            || next_nodes[next_nodes_map[new_node.turn]].getScore(input, is_final) > new_node.getScore(input, is_final)) {
                        next_nodes.push_back(new_node);
                        next_nodes_map[new_node.turn] = next_nodes.size() - 1;
                    }
                }
            }
            vector<pair<int, int>> sorted_nodes;
            for (const auto& [turn, next_node] : next_nodes_map) {
                if (turn < min(input.T, input.T * (i + 1) / input.H + beam_width / 2)) {
                    sorted_nodes.emplace_back(next_nodes[next_node].getScore(input, is_final), next_node);
                }
            }
            sort(sorted_nodes.begin(), sorted_nodes.end());
            nodes.clear();
            for (int j = 0; j < min(beam_width, (int)sorted_nodes.size()); ++j) {
                nodes.push_back(next_nodes[sorted_nodes[j].second]);
            }
        }

        int best_score = numeric_limits<int>::max();
        ContainerStatus best_node(input);
        for (const auto& node : nodes) {
            if (best_score > node.getScore(input, true)) {
                best_node = node;
                best_score = node.getScore(input, true);
            }
        }

        return best_node;
    }
};

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    
    Input input;
    cin >> input;

    Solution solution(input);
    Output output = solution.solve();
    cout << output << flush;

    return 0;
}