#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;

constexpr double pow2(double x) {
    return x * x;
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
    vector<vector<int>> result;
    vector<int> capacity;
    int turn;

    ContainerStatus() : capacity(SIZE, 0), turn(0) {}

    void add(vector<int> ratio) {
        result.push_back(ratio);
        for (int i = 0; i < SIZE; ++i) {
            if (ratio[i] == LIMIT) {
                turn += 2;
                return;
            }
        }
        for (int i = 0; i < SIZE; ++i) {
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
            turn += 1;
        }
    }

    Output toOutput() const {
        Output output(SIZE);
        for (int i = 0; i < SIZE - 1; ++i) {
            for (int j = 0; j < SIZE - 1; ++j) {
                output.v[i][j] = true;
            }
        }
        for (int i = 0; i < SIZE; ++i) {
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

    const int SIZE = 20;
    const int LIMIT = 10;
};

struct Solution {
    Input input;

    Solution(const Input& input_) : input(input_) {}

    Output solve() {
        vector<Point3D> points1;
        for (int i = 0; i < input.K; ++i) {
            vector<int> ratio(20, 0);
            ratio[i] = 10;
            points1.emplace_back(input.own[i], ratio);
        }
        KdTree kdTree1;
        kdTree1.build(points1);

        vector<Point3D> points2;
        for (int i = 0; i < input.K - 1; ++i) {
            for (int j = i + 1; j < input.K; ++j) {
                for (int k = 1; k < 10; ++k) {
                    Color mixed = mixed_color(input.own[i], k, input.own[j], 10 - k);
                    vector<int> ratio(20, 0);
                    ratio[i] = k;
                    ratio[j] = 10 - k;
                    Point3D point(mixed, ratio);
                    points2.push_back(point);
                }
            }
        }
        KdTree kdTree2;
        kdTree2.build(points2);

        vector<Point3D> points3;
        for (int i = 0; i < input.K - 2; ++i) {
            for (int j = i + 1; j < input.K - 1; ++j) {
                for (int k = j + 1; k < input.K; ++k) {
                    for (int l = 1; l < 9; ++l) {
                        for (int m = 1; m < 10 - l; ++m) {
                            vector<int> ratio(20, 0);
                            ratio[i] = l;
                            ratio[j] = m;
                            ratio[k] = 10 - l - m;
                            Color mixed = mixed_color(input.own[i], l, input.own[j], m);
                            mixed = mixed_color(mixed, l + m, input.own[k], 10 - l - m);
                            Point3D point(mixed, ratio);
                            points3.push_back(point);
                        }
                    }
                }
            }
        }
        KdTree kdTree3;
        kdTree3.build(points3);

        ContainerStatus containerStatus;
        for (int i = 0; i < input.H; ++i) {
            Point3D targetPoint(input.target[i]);
            Point3D closestPoint1 = kdTree1.findClosestPoint(targetPoint);
            Point3D closestPoint2 = kdTree2.findClosestPoint(targetPoint);
            Point3D closestPoint3 = kdTree3.findClosestPoint(targetPoint);

            if (input.T - containerStatus.turn - 20 > (input.H - i) * 2) {
                containerStatus.add(closestPoint3.ratio);
            } else {
                containerStatus.add(closestPoint1.ratio);
            }
        }

        return containerStatus.toOutput();
    }
};

int main() {
    Input input;
    cin >> input;

    Solution solution(input);
    Output output = solution.solve();
    cout << output << flush;

    return 0;
}