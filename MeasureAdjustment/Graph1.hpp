#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <algorithm>
#include <string>
#include <set>

using namespace std;

class Graph {
private:
    // 数据结构定义
    unordered_map<string, int> node_map;  // 点名字符串到整型的映射表
    vector<string> node_list;             // 整型到点名字符串的映射
    unordered_map<int, unordered_set<int>> adj;  // 邻接表（整型存储）
    vector<vector<int>> cycles;           // 存储所有最小独立环
    int max_search_depth;                 // 最大搜索深度限制

public:
    Graph() : max_search_depth(50) {}  // 默认最大搜索深度为50

    // 添加边（自动处理点名字符串到整型的转换）
    void addEdge(const string& u, const string& v) {
        // 将点名字符串映射为整型
        if (!node_map.count(u)) {
            node_map[u] = node_list.size();
            node_list.push_back(u);
        }
        if (!node_map.count(v)) {
            node_map[v] = node_list.size();
            node_list.push_back(v);
        }
        int u_id = node_map[u], v_id = node_map[v];

        // 无向图添加边

        adj[u_id].insert(v_id);
        adj[v_id].insert(u_id);
    }

    // 主入口：搜索所有最小独立环
    vector<vector<string>> findMinimalCycles() {
        vector<bool> visited(node_list.size(), false);
        vector<int> parent(node_list.size(), -1);

        // 遍历所有节点作为起点
        for (const auto& node : node_list) {
            int u = node_map[node];
            if (!visited[u]) {
                stack<pair<int, int>> dfs_stack;  // <当前节点, 父节点>
                dfs_stack.push({ u, -1 });
                vector<int> path;

                // 迭代DFS代替递归
                while (!dfs_stack.empty()) {
                    auto [current, prev] = dfs_stack.top();
                    dfs_stack.pop();

                    // 回溯机制：如果当前节点已访问过，需要回退路径
                    if (visited[current]) {
                        while (!path.empty() && path.back() != prev) {
                            visited[path.back()] = false;
                            path.pop_back();
                        }
                        continue;
                    }

                    // 标记访问并记录路径
                    visited[current] = true;
                    path.push_back(current);

                    // 遍历邻接节点
                    for (int neighbor : adj[current]) {
                        if (neighbor == prev) continue;  // 跳过父节点

                        if (!visited[neighbor]) {
                            dfs_stack.push({ neighbor, current });
                        }
                        else {
                            // 发现回边，提取环
                            auto pos = find(path.begin(), path.end(), neighbor);
                            if (pos != path.end()) {
                                vector<int> cycle(pos, path.end());
                                cycle.push_back(neighbor);  // 闭合环
                                normalizeCycle(cycle);      // 标准化环
                                if (isMinimalCycle(cycle)) { // 验证是否为最小环
                                    cycles.push_back(cycle);
                                }
                            }
                        }
                    }// end for neighbor
                }// end while stack
            }// end if !visited
        }// end for all nodes

        // 转换为字符串形式返回
        return convertCyclesToString();
    }

private:
    // 标准化环：最小节点作为起点，字典序排列
    void normalizeCycle(vector<int>& cycle) {
        if (cycle.size() < 3) {
            cycle.clear();
            return;
        }

        // 找到最小节点的位置
        auto min_it = min_element(cycle.begin(), cycle.end() - 1);
        cycle.back() = *min_it;
        rotate(cycle.begin(), min_it, cycle.end() - 1);

        // 确保顺时针或逆时针方向唯一
        if (cycle[1] > cycle[cycle.size() - 2]) {
            reverse(cycle.begin() + 1, cycle.end() - 1);
        }
    }

    // 验证是否为最小独立环（未被其他环包含）
    bool isMinimalCycle(const vector<int>& cycle) {
        for (int node : cycle) {
            std::cout << "-" << node_list[node];
        }
        std::cout << std::endl;
        //return true;
        /*unordered_set<int> nodes(cycle.begin(), cycle.end());
        for (const auto& exist_cycle : cycles) {
            unordered_set<int> exist_nodes(exist_cycle.begin(), exist_cycle.end());
            if (includes(exist_nodes.begin(), exist_nodes.end(),
                nodes.begin(), nodes.end())) {
                return false;
            }
        }
        return true;*/

        std::set<int> nodes(cycle.begin(), cycle.end());
        for (const auto& exist_cycle : cycles) {
            std::set<int> exist_nodes(exist_cycle.begin(), exist_cycle.end());

            // 检查当前节点集合是否包含在已存在的节点集合中
            if (std::includes(exist_nodes.begin(), exist_nodes.end(), nodes.begin(), nodes.end())) {
                return false;
            }
        }
        return true;
    }

    // 将整型环转换为点名字符串环
    vector<vector<string>> convertCyclesToString() {
        vector<vector<string>> result;
        for (const auto& cycle : cycles) {
            vector<string> str_cycle;
            for (int node : cycle) {
                str_cycle.push_back(node_list[node]);
            }
            result.push_back(str_cycle);
        }
        return result;
    }
};

/*int main() {
    // 示例测试（论文中CPIII网示例）
    Graph g;
    g.addEdge("A", "B");
    g.addEdge("B", "C");
    g.addEdge("C", "D");
    g.addEdge("D", "A");  // 环1: A-B-C-D
    g.addEdge("B", "D");
    g.addEdge("D", "E");
    g.addEdge("E", "B");  // 环2: B-D-E

    auto cycles = g.findMinimalCycles();

    cout << "Found " << cycles.size() << " minimal independent cycles:" << endl;
    for (const auto& cycle : cycles) {
        for (const string& node : cycle) {
            cout << node << " ";
        }
        cout << endl;
    }

    return 0;
}*/