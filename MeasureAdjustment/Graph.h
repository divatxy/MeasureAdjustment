#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <stack>

using namespace std;

class Graph {
private:
    int V;                      // ������
    vector<vector<int>> adj;    // �ڽӱ�����ͼ��
    vector<vector<int>> cycles; // �洢������С������
    vector<int> parent;         // �������ĸ��ڵ�
    vector<bool> visited;       // ��Ƕ����Ƿ񱻷��ʹ�

public:
    Graph(int vertices) : V(vertices), adj(vertices), parent(vertices, -1), visited(vertices, false) {}

    // ��ӱߣ�����ͼ��
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // ����������С������
    vector<vector<int>> findMinimalCycles() {
        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                dfs(i);
            }
        }
        return cycles;
    }

private:
    // DFS�����������������ر��γɵ���С��
    void dfs(int u) {
        stack<int> stack;
        stack.push(u);
        visited[u] = true;

        while (!stack.empty()) {
            int current = stack.top();
            bool hasUnvisited = false;

            for (int v : adj[current]) {
                if (v == parent[current]) continue; // ���Ը��ڵ�

                if (!visited[v]) {
                    parent[v] = current;
                    visited[v] = true;
                    stack.push(v);
                    hasUnvisited = true;
                    break;
                }
                else {
                    // ���ֻرߣ���ȡ·���γɻ�
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

    // �ӻر�����ȡ����·��
    vector<int> findCycle(int u, int v) {
        vector<int> cycle;

        // ��u���ϻ��ݵ��������ȣ��ҵ�����·��
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
            cycle.clear(); // �ǻ�·���������ߣ�
        }

        return cycle;
    }

    // ��׼����������С��������ȷ��Ψһ��
    void normalizeCycle(vector<int>& cycle) {
        if (cycle.size() < 3) {
            cycle.clear();
            return;
        }

        // �ҵ���С�����λ��
        auto minIt = min_element(cycle.begin(), cycle.end());
        rotate(cycle.begin(), minIt, cycle.end());

        // ȷ�����ķ���Ψһ�����磬ѡ�����ڽڵ��С�ķ���
        if (cycle[1] > cycle.back()) {
            reverse(cycle.begin() + 1, cycle.end());
        }
    }
};

/*
int main() {
    // ʾ��������һ�������������ͼ
    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 0);    // ��1��0-1-2-3
    g.addEdge(2, 4);
    g.addEdge(4, 1);    // ��2��1-2-4
    g.addEdge(4, 5);
    g.addEdge(5, 2);    // ��3��2-4-5

    // ����������С������
    vector<vector<int>> cycles = g.findMinimalCycles();

    // ������
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