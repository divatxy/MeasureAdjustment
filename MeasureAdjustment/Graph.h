#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <stack>

using namespace std;

class Graph {
private:
    int V;                      // 顶点数
    vector<vector<int>> adj;    // 邻接表（无向图）
    vector<vector<int>> cycles; // 存储所有最小独立环
    vector<int> parent;         // 生成树的父节点
    vector<bool> visited;       // 标记顶点是否被访问过

public:
    Graph(int vertices) : V(vertices), adj(vertices), parent(vertices, -1), visited(vertices, false) {}

    // 添加边（无向图）
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // 查找所有最小独立环
    vector<vector<int>> findMinimalCycles() {
        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                dfs(i);
            }
        }
        return cycles;
    }

private:
    // DFS遍历生成树，并检测回边形成的最小环
    void dfs(int u) {
        stack<int> stack;
        stack.push(u);
        visited[u] = true;

        while (!stack.empty()) {
            int current = stack.top();
            bool hasUnvisited = false;

            for (int v : adj[current]) {
                if (v == parent[current]) continue; // 忽略父节点

                if (!visited[v]) {
                    parent[v] = current;
                    visited[v] = true;
                    stack.push(v);
                    hasUnvisited = true;
                    break;
                }
                else {
                    // 发现回边，提取路径形成环
                    vector<int> cycle = findCycle(current, v);
                    if (!cycle.empty()) {
                        normalizeCycle(cycle);
                        cycles.push_back(cycle);
                    }
                }
            }

            if (!hasUnvisited) {
                stack.pop();
            }
        }
    }

    // 从回边中提取环的路径
    vector<int> findCycle(int u, int v) {
        vector<int> cycle;

        // 从u向上回溯到公共祖先，找到环的路径
        int ancestor = u;
        while (ancestor != -1 && ancestor != v) {
            cycle.push_back(ancestor);
            ancestor = parent[ancestor];
        }

        if (ancestor == v) {
            cycle.push_back(v);
            reverse(cycle.begin(), cycle.end());
            int temp = u;
            while (temp != v) {
                cycle.push_back(temp);
                temp = parent[temp];
            }
        }
        else {
            cycle.clear(); // 非环路径（如横向边）
        }

        return cycle;
    }

    // 标准化环：按最小顶点排序，确保唯一性
    void normalizeCycle(vector<int>& cycle) {
        if (cycle.size() < 3) {
            cycle.clear();
            return;
        }

        // 找到最小顶点的位置
        auto minIt = min_element(cycle.begin(), cycle.end());
        rotate(cycle.begin(), minIt, cycle.end());

        // 确保环的方向唯一（例如，选择相邻节点较小的方向）
        if (cycle[1] > cycle.back()) {
            reverse(cycle.begin() + 1, cycle.end());
        }
    }
};

/*
int main() {
    // 示例：创建一个包含多个环的图
    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 0);    // 环1：0-1-2-3
    g.addEdge(2, 4);
    g.addEdge(4, 1);    // 环2：1-2-4
    g.addEdge(4, 5);
    g.addEdge(5, 2);    // 环3：2-4-5

    // 查找所有最小独立环
    vector<vector<int>> cycles = g.findMinimalCycles();

    // 输出结果
    cout << "Found " << cycles.size() << " minimal independent cycles:\n";
    for (const auto& cycle : cycles) {
        for (int v : cycle) {
            cout << v << " ";
        }
        cout << endl;
    }

    return 0;
};
*/