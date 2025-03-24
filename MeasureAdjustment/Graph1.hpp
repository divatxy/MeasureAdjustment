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
    // ���ݽṹ����
    unordered_map<string, int> node_map;  // �����ַ��������͵�ӳ���
    vector<string> node_list;             // ���͵������ַ�����ӳ��
    unordered_map<int, unordered_set<int>> adj;  // �ڽӱ����ʹ洢��
    vector<vector<int>> cycles;           // �洢������С������
    int max_search_depth;                 // ��������������

public:
    Graph() : max_search_depth(50) {}  // Ĭ������������Ϊ50

    // ��ӱߣ��Զ���������ַ��������͵�ת����
    void addEdge(const string& u, const string& v) {
        // �������ַ���ӳ��Ϊ����
        if (!node_map.count(u)) {
            node_map[u] = node_list.size();
            node_list.push_back(u);
        }
        if (!node_map.count(v)) {
            node_map[v] = node_list.size();
            node_list.push_back(v);
        }
        int u_id = node_map[u], v_id = node_map[v];

        // ����ͼ��ӱ�

        adj[u_id].insert(v_id);
        adj[v_id].insert(u_id);
    }

    // ����ڣ�����������С������
    vector<vector<string>> findMinimalCycles() {
        vector<bool> visited(node_list.size(), false);
        vector<int> parent(node_list.size(), -1);

        // �������нڵ���Ϊ���
        for (const auto& node : node_list) {
            int u = node_map[node];
            if (!visited[u]) {
                stack<pair<int, int>> dfs_stack;  // <��ǰ�ڵ�, ���ڵ�>
                dfs_stack.push({ u, -1 });
                vector<int> path;

                // ����DFS����ݹ�
                while (!dfs_stack.empty()) {
                    auto [current, prev] = dfs_stack.top();
                    dfs_stack.pop();

                    // ���ݻ��ƣ������ǰ�ڵ��ѷ��ʹ�����Ҫ����·��
                    if (visited[current]) {
                        while (!path.empty() && path.back() != prev) {
                            visited[path.back()] = false;
                            path.pop_back();
                        }
                        continue;
                    }

                    // ��Ƿ��ʲ���¼·��
                    visited[current] = true;
                    path.push_back(current);

                    // �����ڽӽڵ�
                    for (int neighbor : adj[current]) {
                        if (neighbor == prev) continue;  // �������ڵ�

                        if (!visited[neighbor]) {
                            dfs_stack.push({ neighbor, current });
                        }
                        else {
                            // ���ֻرߣ���ȡ��
                            auto pos = find(path.begin(), path.end(), neighbor);
                            if (pos != path.end()) {
                                vector<int> cycle(pos, path.end());
                                cycle.push_back(neighbor);  // �պϻ�
                                normalizeCycle(cycle);      // ��׼����
                                if (isMinimalCycle(cycle)) { // ��֤�Ƿ�Ϊ��С��
                                    cycles.push_back(cycle);
                                }
                            }
                        }
                    }// end for neighbor
                }// end while stack
            }// end if !visited
        }// end for all nodes

        // ת��Ϊ�ַ�����ʽ����
        return convertCyclesToString();
    }

private:
    // ��׼��������С�ڵ���Ϊ��㣬�ֵ�������
    void normalizeCycle(vector<int>& cycle) {
        if (cycle.size() < 3) {
            cycle.clear();
            return;
        }

        // �ҵ���С�ڵ��λ��
        auto min_it = min_element(cycle.begin(), cycle.end() - 1);
        cycle.back() = *min_it;
        rotate(cycle.begin(), min_it, cycle.end() - 1);

        // ȷ��˳ʱ�����ʱ�뷽��Ψһ
        if (cycle[1] > cycle[cycle.size() - 2]) {
            reverse(cycle.begin() + 1, cycle.end() - 1);
        }
    }

    // ��֤�Ƿ�Ϊ��С��������δ��������������
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

            // ��鵱ǰ�ڵ㼯���Ƿ�������Ѵ��ڵĽڵ㼯����
            if (std::includes(exist_nodes.begin(), exist_nodes.end(), nodes.begin(), nodes.end())) {
                return false;
            }
        }
        return true;
    }

    // �����ͻ�ת��Ϊ�����ַ�����
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
    // ʾ�����ԣ�������CPIII��ʾ����
    Graph g;
    g.addEdge("A", "B");
    g.addEdge("B", "C");
    g.addEdge("C", "D");
    g.addEdge("D", "A");  // ��1: A-B-C-D
    g.addEdge("B", "D");
    g.addEdge("D", "E");
    g.addEdge("E", "B");  // ��2: B-D-E

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