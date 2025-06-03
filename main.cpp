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

// 3次元座標を表す構造体
struct Point3D {
    double x, y, z;
    vector<int> ratio;

    Point3D(Color color = Color(), vector<int> ratio_ = vector<int>(20))
        : x(color.C), y(color.M), z(color.Y), ratio(std::move(ratio_)) {}

    // 指定された軸の座標値を返す
    double getCoord(int axis) const {
        if (axis == 0) return x;
        if (axis == 1) return y;
        return z;
    }
};

// 2点間のユークリッド距離の2乗を計算する関数
constexpr double distanceSquared(const Point3D& p1, const Point3D& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}

// k-d木のノード構造体
struct KdNode {
    Point3D point;
    std::unique_ptr<KdNode> left;
    std::unique_ptr<KdNode> right;
    int axis; // 分割軸 (0 for x, 1 for y, 2 for z)

    KdNode(const Point3D& pt, int ax) : point(pt), axis(ax) {}
};

// k-d木を構築するクラス
class KdTree {
public:
    KdTree() = default;

    // 点のリストからk-d木を構築
    void build(std::vector<Point3D>& points) {
        if (points.empty()) {
            root = nullptr;
            return;
        }
        root = buildRecursive(points, 0, points.size() - 1, 0);
    }

    // 最も近い点を見つける
    Point3D findClosestPoint(const Point3D& target) const {
        if (!root) {
            throw std::runtime_error("KdTree is empty.");
        }
        Point3D bestPoint = root->point; // 初期値をルートノードの点とする
        double bestDistSq = distanceSquared(target, root->point);
        findClosestRecursive(root.get(), target, bestPoint, bestDistSq);
        return bestPoint;
    }

private:
    std::unique_ptr<KdNode> root;

    std::unique_ptr<KdNode> buildRecursive(std::vector<Point3D>& points, int start, int end, int depth) {
        if (start > end) {
            return nullptr;
        }

        int axis = depth % 3; // 0: x, 1: y, 2: z
        int mid = start + (end - start) / 2;

        // std::nth_elementを使って中央値を見つけ、その位置に配置する
        // これにより、ソート全体を行うよりも効率的に中央値を得られる
        std::nth_element(points.begin() + start, points.begin() + mid, points.begin() + end + 1,
            [axis](const Point3D& a, const Point3D& b) {
                return a.getCoord(axis) < b.getCoord(axis);
            });

        auto node = std::make_unique<KdNode>(points[mid], axis);
        node->left = buildRecursive(points, start, mid - 1, depth + 1);
        node->right = buildRecursive(points, mid + 1, end, depth + 1);

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

        // ターゲットが分割平面のどちら側にあるかで、探索する子ノードの順番を決める
        const std::unique_ptr<KdNode>* firstChild = (targetCoord < nodeCoord) ? &(node->left) : &(node->right);
        const std::unique_ptr<KdNode>* secondChild = (targetCoord < nodeCoord) ? &(node->right) : &(node->left);

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
    static constexpr int SIZE = 20;
    static constexpr int LIMIT = 10;

    int turn;
    int score;
    vector<vector<int>> result;
    vector<int> capacity;

    ContainerStatus() : turn(0), score(0), capacity(SIZE, 0) {}

    void add(const Input &input, vector<int> ratio) {
        result.push_back(ratio);
        for (int i = 0; i < SIZE; ++i) {
            if (ratio[i] == LIMIT) {
                turn += 2;
                score += roundint(input.target[result.size() - 1].distance(input.own[i]) * 10000);
                return;
            }
        }
        for (int i = 0; i < SIZE; ++i) {
            if (ratio[i] == 0) continue;
            if (capacity[i] < ratio[i]) {
                capacity[i] += LIMIT - ratio[i];
                if (capacity[i] == 0) {
                    turn += 3;
                } else {
                    turn += 5;
                }
            } else {
                capacity[i] -= ratio[i];
                turn += 2;
            }
        }
        turn += 1;

        Color mixedColor;
        int totalRatio = 0;
        for (int i = 0; i < SIZE; ++i) {
            if (ratio[i] == 0) continue;
            mixedColor = mixed_color(mixedColor, totalRatio, input.own[i], ratio[i]);
            totalRatio += ratio[i];
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

    Output toOutput(const Input &input) const {
        Output output(SIZE);
        for (int i = 0; i < SIZE - 1; ++i) {
            for (int j = 0; j < min(SIZE - 1, input.K); ++j) {
                output.v[i][j] = true;
            }
        }
        for (int i = 0; i < input.K; ++i) {
            output.h[LIMIT - 1][i] = true;
        }

        vector<int> now_capacity(SIZE, 0);
        for (const auto& ratio : result) {
            bool is_full = false;
            for (int i = 0; i < SIZE; ++i) {
                if (ratio[i] == LIMIT) {
                    output.add(SIZE - 1, SIZE - 1, i);
                    output.submission(SIZE - 1, SIZE - 1);
                    is_full = true;
                    continue;
                }
            }
            if (is_full) continue;

            for (int i = 0; i < SIZE; ++i) {
                if (ratio[i] == 0) continue;
                if (now_capacity[i] < ratio[i]) {
                    if (now_capacity[i] == 0) {
                        output.add(LIMIT - 1, i, i);
                    } else {
                        output.partition(now_capacity[i] + LIMIT - 1, i, now_capacity[i] + LIMIT, i);
                        output.add(now_capacity[i] + LIMIT - 1, i, i);
                        output.partition(now_capacity[i] - 1, i, now_capacity[i], i);
                    }
                    output.partition(now_capacity[i] + LIMIT - ratio[i] - 1, i, now_capacity[i] + LIMIT - ratio[i], i);
                } else {
                    if (now_capacity[i] - ratio[i] != 0) {
                        output.partition(now_capacity[i] - ratio[i] - 1, i, now_capacity[i] - ratio[i], i);
                    }
                }
            }
            
            for (int i = 0; i < SIZE; ++i) {
                if (ratio[i] == 0) continue;
                if (now_capacity[i] < ratio[i]) {
                    output.partition(now_capacity[i] + LIMIT - 1, i, now_capacity[i] + LIMIT, i);
                } else {
                    output.partition(now_capacity[i] - 1, i, now_capacity[i], i);
                }
            }
            output.submission(SIZE - 1, SIZE - 1);
            for (int i = 0; i < SIZE; ++i) {
                if (ratio[i] == 0) continue;
                if (now_capacity[i] < ratio[i]) {
                    now_capacity[i] += LIMIT - ratio[i];
                } else {
                    if (now_capacity[i] - ratio[i] == 0) {
                        output.partition(LIMIT - 1, i, LIMIT, i);
                    }
                    now_capacity[i] -= ratio[i];
                }
            }
        }

        return output;
    }
};

struct Solution {
    Input input;
    vector<KdTree> kdTrees;
    int DEPTH;

    Solution(const Input& input_) : input(input_), DEPTH(5 <= input.K && input.K <= 18 ? (6 <= input.K && input.K <= 15 ? 6 : 5) : 4) {}

    Output solve() {
        init_kdtree();

        int beam_width = 4;
        ContainerStatus containerStatus = beam_search(beam_width);

        return containerStatus.toOutput(input);
    }

    void init_kdtree() {
        vector<vector<Point3D>> points(DEPTH);
        vector<int> ratio(20, 0);
        auto dfs = [&](auto&& self, int pt, int sum, Color color, int depth) -> void {
            if (sum == 10) {
                points[depth - 1].emplace_back(color, ratio);
                return;
            }
            if (depth >= DEPTH) {
                return;
            }
            if (pt >= input.K) {
                return;
            }
            for (int npt = pt + 1; npt < input.K; ++npt) {
                for (int k = 1; k <= 10 - sum; ++k) {
                    ratio[npt] = k;
                    Color mixed = mixed_color(input.own[npt], k, color, sum);
                    self(self, npt, sum + k, mixed, depth + 1);
                    ratio[npt] = 0;
                }
            }
        };
        dfs(dfs, -1, 0, Color(), 0);
        for (int i = 0; i < DEPTH; ++i) {
            kdTrees.emplace_back();
            kdTrees.back().build(points[i]);
        }
    }

    ContainerStatus beam_search(int beam_width) {
        vector<ContainerStatus> nodes;
        nodes.push_back(ContainerStatus());
        for (int i = 0; i < input.H; ++i) {
            bool is_final = (i >= input.H - input.K);
            map<int, ContainerStatus> next_nodes_map;
            Point3D targetPoint(input.target[i]);
            for (const auto& node : nodes) {
                for (int j = 0; j < DEPTH; ++j) {
                    Point3D closestPoint = kdTrees[j].findClosestPoint(targetPoint);
                    ContainerStatus new_node = node;
                    new_node.add(input, closestPoint.ratio);
                    if (next_nodes_map.count(new_node.turn) == 0 || next_nodes_map[new_node.turn].getScore(input, is_final) > new_node.getScore(input, is_final)) {
                        next_nodes_map[new_node.turn] = new_node;
                    }
                }
            }

            vector<ContainerStatus> next_nodes;
            vector<pair<int, int>> sorted_nodes;
            for (const auto& [turn, next_node] : next_nodes_map) {
                if (turn < min(input.T, input.T * (i + 1) / input.H + beam_width / 2)) {
                    next_nodes.push_back(next_node);
                    sorted_nodes.emplace_back(next_node.getScore(input, is_final), next_nodes.size() - 1);
                }
            }
            sort(sorted_nodes.begin(), sorted_nodes.end());
            nodes.clear();
            for (int j = 0; j < min(beam_width, (int)sorted_nodes.size()); ++j) {
                nodes.push_back(next_nodes[sorted_nodes[j].second]);
            }
        }

        double best_score = numeric_limits<double>::max();
        ContainerStatus best_node;
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