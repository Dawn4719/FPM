#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <malloc.h>
#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"

#define qw cout<<endl;
#define dl endl
#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define NC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))
#define gti std::chrono::high_resolution_clock::now()
#define cti(y, x) std::chrono::duration_cast<std::chrono::nanoseconds>(y - x).count()
#define fi first
#define se second
string QUE = "5";
string PATH = "ER100";
int QUERY_NUMS = 3;
int D = 1;
int SPECIAL = 2;
int Meth = 1;
double tim1, tim2, tim3, tim4;
bool test = false;
bool in_file = false;
bool cot = false;
int path_nums;
int debug_cnt = 0;
bool RE;
bool topk = false;
int top = 0x3f3f3f3f;

std::string input_query_graph_file = "../../test/querys/" + QUE;
std::string input_data_graph_file = "../../test/graphs/Dataset/" + PATH + "K/" + PATH + "K.txt";
std::string input_filter_type = "CECI";
std::string input_order_type = "CECI";
std::string input_engine_type = "CECI";
std::string input_max_embedding_num = "MAX";
std::string input_time_limit = "60";
std::string input_order_num = "100";
std::string input_distribution_file_path = "temp.distribution";
std::string input_csr_file_path;

double sum_load_graphs_time_in_ns = 0;
size_t sum_memory_cost_in_bytes = 0;
size_t init_memory_cost_in_bytes = 0;
size_t sum_embedding_count = 0;
size_t sum_call_count = 0;

size_t try_count;
size_t dfs_all_count;
vector<pair<int, vector<int>>> match_result;
vector<int> ORD;
bool* query_is_matched;
size_t sp_match_res;
size_t cur_sp_match_res;
TreeNode *ceci_tree;
ui *ceci_order;
ui *candidates_count;
bool* is_show;
ui *idx ;
ui *idx_count ;
ui *embedding ;
vector<vector<ui>> valid_candidates ;
bool *visited_vertices ;

size_t all_memory;
size_t all_query_vertex_nums;
size_t time_cnt;
const int B = (1 << QUERY_NUMS) - 2;
string getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
    return tmp;
}

double time1;
double time2;
double time3;
double time4;

vector<long long> memory;
ui **candidates = nullptr;
ui* matching_order;
ui* pivots;
map<string, map<string, map<string, double>>> time_limit;
int mx = -1;
vector<bool> dfs_st;
vector<int> dfs_dist;
auto total_start = std::chrono::high_resolution_clock::now();
map<string, size_t> get_index_mem() {
    FILE* fp = fopen("/proc/self/status", "r");
    char line[128];
    map<string, size_t> res;
    while (fgets(line, 128, fp) != NULL)
    {
//        if (strncmp(line, "VmPeak", 2) == 0)
//        {
//            cout << line << endl;
////            printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
//        }
        if (strncmp(line, "VmRSS:", 6) == 0) {
            string p = line;
            res["now"] = size_t(stoull(p.substr(6)));
//            cout << line;
        }
        if (strncmp(line, "VmPeak:", 7) == 0) {
            string p = line;
            res["pk"] = size_t(stoull(p.substr(7)));
//            cout << line;
        }
    }
    fclose(fp);
    return res;
}
vector<int> pst;
vector<int> path;
vector<pair<int, vector<int>>> res_graph;
auto LB(vector<pair<int, vector<int>>>& V, uint val) {
    size_t L = 0, R = V.size();
    while (L < R) {
        size_t MID = (L + R) >> 1;
        if (V[MID].first >= val) R = MID;
        else L = MID + 1;
    }
    return L;
}
void add2graph(Graph* data_graph, vector<int>& path) {
    for (int i = 0; i < path.size() - 1; ++i) {
        auto v = path[i];
        auto u = path[i + 1];
        ui nbrs_cnt;
        auto nbrs = data_graph->getVertexDNeighbors(v, nbrs_cnt); // v_f's neighbors
        auto idx_ = std::lower_bound(nbrs, nbrs + nbrs_cnt, u);
        if (nbrs_cnt == 0) continue;
        if (*idx_ != u) continue;

        auto lb = LB(res_graph, v);
        if (lb >= res_graph.size() || res_graph[lb].first != v) res_graph.insert(res_graph.begin() + lb, {v, vector<int>{}});
        auto lb2 = lower_bound(res_graph[lb].second.begin(), res_graph[lb].second.end(), u);
        if (lb2 == res_graph[lb].second.end() || *lb2 != u) res_graph[lb].second.insert(lb2, u);
    }
    for (int i = path.size() - 1; i > 0; --i) {
        auto v = path[i];
        auto u = path[i - 1];
        ui nbrs_cnt;
        auto nbrs = data_graph->getVertexDNeighbors(v, nbrs_cnt); // v_f's neighbors
        auto idx_ = std::lower_bound(nbrs, nbrs + nbrs_cnt, u);
        if (nbrs_cnt == 0 || *idx_ != u) continue;
        auto lb = LB(res_graph, v);
        if (lb >= res_graph.size() || res_graph[lb].first != v) res_graph.insert(res_graph.begin() + lb, {v, vector<int>{}});
        auto lb2 = lower_bound(res_graph[lb].second.begin(), res_graph[lb].second.end(), u);
        if (lb2 == res_graph[lb].second.end() || *lb2 != u) res_graph[lb].second.insert(lb2, u);
    }
}
bool get_path3D(Graph* data_graph, int b, int u, int fa, int dd, vector<pair<int, int>>& path) {
    if (dd >= D) {
        return false;
    }
    ui nbrs_cnt;
    auto nbrs = data_graph->getVertexNeighbors(u, nbrs_cnt); // v_f's neighbors
    bool this_edge = false;
    for (int i = 0; i < nbrs_cnt; ++i) {
        auto v = nbrs[i];

        if (v == fa) continue;
        if (pst[v] < 0) {
            if (-pst[v] <= QUERY_NUMS) {
                if (-pst[v] != b) {
                    // cout << u << " " << v << endl;
                    path.emplace_back(u, (int)v);
                    this_edge = true;
                }
            }
            continue;
        }
        if (pst[v] > 0) {
            if (pst[v] > dd + 1) {
                pst[v] = dd + 1;
                if (get_path3D(data_graph, b, v, u, dd + 1, path)) {
                    path.emplace_back(u, (int)v);
                    this_edge = true;
                }
                else
                    pst[v] = 0;
            }
        }
        if (pst[v] == 0) {
            pst[v] = dd + 1;
            if (get_path3D(data_graph, b, v, u, dd + 1, path)) {
                path.emplace_back(u, (int)v);
                this_edge = true;
            }
            else {
                pst[v] = 0;
            }
        }
    }
    return this_edge;
}
map<int, vector<vector<vector<int>>>> his;
bool check(Graph* data_graph, Graph** query_graph, int cnt, vector<pair<int, vector<vector<int>>>>& all_match_res, int which_i) {
    dfs_all_count++;
    if (topk) return true;

    if (cnt == QUERY_NUMS) {
        int sum = 0;
        int cnt_ = 0;
        vector<int> all_node;
        auto tmp = match_result;
        for (auto& i : tmp) {
            sort(i.second.begin(), i.second.end());
            for (auto j : i.second)
                sum += j;
        }
        sort(tmp.begin(), tmp.end());
        vector<vector<int>> nw_res;
        for (auto i : tmp) nw_res.emplace_back(i.second);

        if (his.find(sum) != his.end()) {
            auto it = his[sum];
            bool pa = true;
            int pa_cnt = 0;
            for (int k = 0; k < it.size(); ++k) {
                pa = true;
                for (int i = 0; i < QUERY_NUMS; ++i) {
                    if (nw_res[i].size() != it[k][i].size()) {
                        pa = false;
                        break;
                    }
                    for (int j = 0; j < nw_res[i].size(); ++j) {
                        if (nw_res[i][j] != it[k][i][j]) {
                            pa = false;
                            break;
                        }
                    }
                    if (!pa) break;
                }
                if (pa)
                    return true;
            }
        }
        his[sum].emplace_back(nw_res);
    //    fstream fs = fstream("../result1.txt", ios::app);
        // auto bk = match_result;
        // for (auto &i : bk) sort(i.second.begin(), i.second.end());
        // map<int, bool> neiMp;
        // for (const auto& i : match_result) {
        //     cout << "< " << i.first << ": ";
        //     for (auto v : i.second) {
        //         cout << v << " ";
        //         ui nbrs_cnt;
        //         const VertexID *nbrs = data_graph->getVertexNeighbors(v, nbrs_cnt);
        //         for (int ii = 0; ii < nbrs_cnt; ii++) {
        //             VertexID v = nbrs[ii];
        //             neiMp[v] = true;
        //         }
        //     }
        //     cout << "> ";
        // }
        // cout << endl;
        // for (auto vv : neiMp) {
        //     cout << vv.first << " (" << data_graph->getVertexDegree(vv.first) << ") ";
        // }
        // cout << endl;
        // int x; cin >> x;

        auto sstart = std::chrono::high_resolution_clock::now();

        fill(pst.begin(), pst.end(), 0);
        for (const auto& i : match_result) {
            for (auto v : i.second) {
                if (pst[v] >= 0) pst[v] = -i.first - 1;
            }
        }
        vector<vector<pair<int, int>>> Path;
        for (const auto& i : match_result) {
            vector<pair<int, int>> pa;
            for (auto v : i.second) {
                get_path3D(data_graph, i.first + 1, v, -1, 0, pa);
            }
            Path.emplace_back(pa);
            for (const auto& v : pa) {
                if (pst[v.first] > 0) pst[v.first] = 0;
                if (pst[v.second] > 0) pst[v.second] = 0;
            }
        }

        auto eend = std::chrono::high_resolution_clock::now();
        time4 += std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
        if (NC(total_time) >= 3600) {
            topk = true;
            cout << "time exced" << endl;
        }

        cur_sp_match_res++;
        if (cur_sp_match_res % 10000 == 0)
            cout <<  "asd " << cur_sp_match_res << " " << NC(total_time) << "ms " << get_index_mem()["pk"] << "kb" << endl;
        //cout << cur_sp_match_res << endl;
        if (cur_sp_match_res >= top) {
            cout<< "Sdasdasda" << " " << NC(total_time) << endl;
            topk = true;
        }
        if (topk) return true;
        return true;
    }

    bool* st = new bool[data_graph->vertices_count_];
    int* dist = new int[data_graph->vertices_count_];
    bool* in_queue = new bool[data_graph->vertices_count_];

    memset(st, 0, sizeof(bool) * data_graph->vertices_count_);
    memset(dist, 0, sizeof(int) * data_graph->vertices_count_);
    memset(in_queue, 0, sizeof(bool) * data_graph->vertices_count_);

    queue<int> q;
    for (const auto& pp : match_result) {
        for (int p : pp.second) {
            q.emplace(p);
            in_queue[p] = true;
        }
    }
    bool this_turn = false;
    while (!q.empty()) {
        int vertex = q.front(); q.pop();
        // if (which_i == 1)
        //     cout << vertex << endl;
        int d = dist[vertex];
        if (st[vertex])
            continue;
        // cout << vertex << " " << d << endl;
        for (int i = 1; i < QUERY_NUMS; ++i) {
            if (!query_is_matched[i]) {
                for (int j = 0; j < all_match_res[i].second.size(); ++j) {
                    try_count++;
                    if (find(all_match_res[i].second[j].begin(), all_match_res[i].second[j].end(), vertex) != all_match_res[i].second[j].end()) {
                        match_result.emplace_back(i, all_match_res[i].second[j]);
                        query_is_matched[i] = true;
                        if (check(data_graph, query_graph, cnt + 1, all_match_res, which_i)) {
                            this_turn = true;
                        }
                        if (topk) return true;
                        match_result.pop_back();
                        query_is_matched[i] = false;
                    }
                }
            }
        }
        st[vertex] = true;
        if (d + 1 > D) continue;
        // 邻居加入队列
        ui nbrs_cnt;
        const VertexID *nbrs = data_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (int ii = 0; ii < nbrs_cnt; ii++) {
            VertexID v = nbrs[ii];
            if (!st[v]) // csr寸的是点值
            {
                if (!in_queue[v]) // 不在队列中
                {
                    q.emplace(v);
                    in_queue[v] = true;
                    dist[v] = d + 1;
                }
            }
        }
    }
    delete []st;
    delete []dist;
    delete []in_queue;
    return this_turn;
}
map<ui, map<ui, bool>> dont_need_match;
void base_solve(double load_graphs_time_in_ns, Graph** query_graph, Graph* data_graph) {
//    std::cout << "Start queries..." << std::endl;
//    std::cout << "-----" << std::endl;
//    std::cout << "Filter candidates..." << std::endl;
    auto sstart = std::chrono::high_resolution_clock::now();

    ui max_vetex_nums = 0;
    for (int i = 0; i < QUERY_NUMS; ++i) {
        if (query_graph[i]->vertices_count_ > max_vetex_nums) {
            max_vetex_nums = query_graph[i]->vertices_count_;
        }
    }
    size_t MEM_SIZE = 0;

    vector<pair<int, vector<vector<int>>>> all_match_res;
    all_match_res.reserve(QUERY_NUMS);
    auto cecitime1 = std::chrono::high_resolution_clock::now();

    for (int current = 0; current < QUERY_NUMS; ++current) {
        auto each_begin = std::chrono::high_resolution_clock::now();
        candidates_count = new VertexID[mx];
        ceci_order = new VertexID[mx];
        ceci_tree = new TreeNode[mx];
        for (int tree_i = 0; tree_i < mx; ++tree_i) {
            ceci_tree[tree_i].initialize(mx);
            ceci_tree[tree_i].clear();
        }
        matching_order = new ui[mx];
        pivots = new ui[mx];
        candidates = new ui*[mx];
        for (int i = 0; i < mx; ++i) {
            candidates[i] = new ui[data_graph->vertices_count_];
        }
        idx = new ui[mx];
        idx_count = new ui[mx];
        embedding = new ui[mx];
        valid_candidates.resize(mx);
        visited_vertices = new bool[data_graph->vertices_count_];
        std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
        std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
        bool filter_res;

        filter_res = FilterVertices::CECIFilter(data_graph, query_graph[current], candidates, candidates_count, ceci_order,
                                                ceci_tree, TE_Candidates, NTE_Candidates, -1, -1, nullptr, SPECIAL);
        if (!filter_res) {
            cout << "filter not match" << endl;
            return;
        }
        // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
        std::vector<ui> optimal_candidates_count;
        double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                              candidates_count, optimal_candidates_count);
        FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
        if (!test) {
            std::cout << "-----" << std::endl;
            std::cout << "Build indices..." << std::endl;
        }

        Edges ***edge_matrix = NULL;

        size_t memory_cost_in_bytes = 0;

        if (TE_Candidates.empty()) {
            cout << "no result" << endl;
            return;
        }
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph[current], candidates_count, ceci_order,
                                                                    ceci_tree,
                                                                    TE_Candidates, NTE_Candidates);
//        BuildTable::printTableCardinality(query_graph[current], ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);
        if (!test) {
            std::cout << "-----" << std::endl;
            std::cout << "Generate a matching order..." << std::endl;
        }

//        ui *matching_order = NULL;
//        ui *pivots = NULL;
//        ui **weight_array = NULL;

        size_t order_num = 0;

        sscanf(input_order_num.c_str(), "%zu", &order_num);

        std::vector<std::vector<ui>> spectrum;

        GenerateQueryPlan::generateCECIQueryPlan(query_graph[current], ceci_tree, ceci_order, matching_order, pivots);

        if (input_order_type != "Spectrum") {
//            GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[current], matching_order, pivots);
//            GenerateQueryPlan::printSimplifiedQueryPlan(query_graph[current], matching_order);
        } else {
//            std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
        }
        if (!test) {
            std::cout << "-----" << std::endl;
            std::cout << "Enumerate..." << std::endl;
        }
        size_t output_limit = 0;
        size_t embedding_count = 0;
        if (input_max_embedding_num == "MAX") {
            output_limit = std::numeric_limits<size_t>::max();
        } else {
            sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
        }


#if ENABLE_QFLITER == 1
        EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

        size_t call_count = 0;
        size_t time_limit = 0;
        sscanf(input_time_limit.c_str(), "%zu", &time_limit);

        vector<vector<int>> cur_match_res;

        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph[current], ceci_tree, candidates,
                                                          candidates_count, TE_Candidates,
                                                          NTE_Candidates, ceci_order, output_limit, call_count, cur_match_res,idx,idx_count,embedding,valid_candidates,visited_vertices);
        query_graph[current]->emb_count = embedding_count;
//        if (50 <= embedding_count) {
//            vector<int> lb;
//            for (auto iii: ppa)
//                lb.emplace_back(data_graph->getVertexLabel(iii));
//            sort(lb.begin(), lb.end());
//            for (auto iii : lb)
//                cout << iii << " ";
//            cout << "[] ";
//            for (auto i : ppa)
//                cout << data_graph->getVertexLabel(i) << " ";
//            cout << "[] ";
//            for (auto i : ppa)
//                cout << i << " ";
//            cout << "-embedding_count: " << embedding_count << endl;
//        }
        if (embedding_count == 0) {
            cout << current << " no result" << endl;
            return;
        }
//        return;
        all_match_res.emplace_back(current, cur_match_res);
#ifdef DISTRIBUTION
        std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
        outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
        delete[] EvaluateQuery::distribution_count_;
#endif
        if (!test) {
            std::cout << "--------------------------------------------------------------------" << std::endl;
            std::cout << "Release memories..." << std::endl;
        }
        /**
         * Release the allocated memories.
         */
        delete[] candidates_count;
        delete[] ceci_tree;
        delete[] ceci_order;
        delete[] matching_order;
        delete[] pivots;
        for (ui i = 0; i < query_graph[current]->getVerticesCount(); ++i) {
            delete[] candidates[i];
        }
        delete[] candidates;
        delete[] idx;
        delete[] idx_count;
        delete[] embedding;
        // delete[] valid_candidates;
        delete[] visited_vertices;
//        std::cout << "--------------------------------------------------------------------" << std::endl;
        printf("#Embeddings: %zu\n", embedding_count);
        auto each_end = std::chrono::high_resolution_clock::now();
        double each_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - each_begin).count();
        cout << NC(each_time) << endl;
    }

    auto cecitime2 = std::chrono::high_resolution_clock::now();
    double cecitime = std::chrono::duration_cast<std::chrono::nanoseconds>(cecitime2 - cecitime1).count();
    if (test) return;
    std::cout << "query_graph match result:" << std::endl;
    MEM_SIZE += all_match_res.size();
    for (auto& i : all_match_res) {
        mx = std::max(mx, (int)i.second.size());
        std::cout << i.second.size() << std::endl;
        for (auto& j : i.second) {
            MEM_SIZE += j.size();
        }
    }

    // sort(query_graph, query_graph + QUERY_NUMS, [](const Graph* a, const Graph* b) {
    //     return a->emb_count < b->emb_count;
    // });

    cout << "Query order: " << endl;
    for (int j = 0; j < QUERY_NUMS; ++j) {
        cout << query_graph[j]->id_ << " " << query_graph[j]->emb_count << endl;
        all_match_res[query_graph[j]->id_].first = query_graph[j]->emb_count;
    }
    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
//    for (auto& i : all_match_res) sort(i.second.begin(), i.second.end());
//    sort(all_match_res.begin(), all_match_res.end(), [](const pair<int, vector<vector<int>>>& a, const pair<int, vector<vector<int>>>& b) {
//        return a.first < b.first;
//    });
//    for (int i = 0; i < QUERY_NUMS; ++i)
//        sort(all_match_res[i].second.begin(), all_match_res[i].second.end());
    if (!test) {
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Begin SP Match!" << std::endl;
    }
    cur_sp_match_res = 0;
    topk = false;
//    for (int i = 0; i < all_match_res[0].second.size(); ++i) {
//        for (auto j :all_match_res[0].second[i])
//            cout << j << " ";
//        cout << endl;
//    }
    for (int i = 0; i < all_match_res[0].second.size(); ++i) {
        for (int j = 0; j < QUERY_NUMS; ++j) query_is_matched[j] = false;
        topk = false;
        his.clear();
        auto each_begin = std::chrono::high_resolution_clock::now();
        match_result.clear();
        query_is_matched[0] = true;
        match_result.emplace_back(0, all_match_res[0].second[i]);

//        cout << i << " " << topk << " " << cur_sp_match_res << endl;
        check(data_graph, query_graph, 1, all_match_res, i);
        auto each_end = std::chrono::high_resolution_clock::now();
        double each_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - each_begin).count();
//        sp_match_res += cur_sp_match_res;
        cout << i << " " << cur_sp_match_res << " " << cur_sp_match_res << " " << NANOSECTOSEC(each_time) << endl;
//        cur_sp_match_res = 0;
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - total_start).count();
        if (topk) break;
    }
    auto eend = std::chrono::high_resolution_clock::now();
    double MatchTime = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
    std::cout << "--------------------------------------------------------------------" << std::endl;

    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Memory cost (MB): %.4lf\n", BYTESTOMB(all_memory + sum_memory_cost_in_bytes));
    printf("#Embeddings: %zu\n", cur_sp_match_res);
    printf("DFS call count: %zu\n", dfs_all_count);
    printf("try count: %zu\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(MatchTime));
//    printf("Avg queue size: %llu\n", sizeof(int) * (avg_res + avg_res2) / calculate_count);

    std::cout << "End." << std::endl;
}

void CECI(Graph* data_graph, Graph** query_graph, vector<vector<int>>& cur_match_res, int i, int u, int vertex) {
    if (u != -1)
    if (data_graph->getVertexDegree(vertex) < query_graph[i]->getVertexDegree(u))
        return;
    const std::unordered_map<LabelID, ui> *u_nlf = query_graph[i]->getVertexNLF(u);
    const std::unordered_map<LabelID, ui> *v_nlf = data_graph->getVertexNLF(vertex);
    if (u != -1) {
        if (v_nlf->size() >= u_nlf->size()) {
            bool is_valid = true;
            for (auto element: *u_nlf) {
                auto iter = v_nlf->find(element.first);
                if (iter == v_nlf->end() || iter->second < element.second) {
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid) {
                return;
            }
        } else {
            return;
        }
    }

    for (int tree_i = 0; tree_i < mx; ++tree_i) { ceci_tree[tree_i].clear(); }

    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    bool filter_res;

    filter_res = FilterVertices::CECIFilter(data_graph, query_graph[i],
                                            candidates,
                                            candidates_count, ceci_order, ceci_tree,
                                            TE_Candidates,
                                            NTE_Candidates, u, vertex,
                                            nullptr, SPECIAL);
    if (!filter_res) {
        return;
    }
#ifdef OPTIMAL_CANDIDATES
    std::vector<ui> optimal_candidates_count;
     double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                               candidates_count, optimal_candidates_count);
     FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
    if (TE_Candidates.empty()) {
        return;
    }
    GenerateQueryPlan::generateCECIQueryPlan(query_graph[i], ceci_tree,
                                             ceci_order,
                                             matching_order,
                                             pivots);
    GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[i],matching_order,pivots);

    size_t output_limit = 0;
    size_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<size_t>::max();
    } else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }
#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif
    size_t call_count = 0;
    size_t time_limit = 0;
    //                            sscanf(input_time_limit.c_str(), "%zu", &time_limit);
    //                            auto t3e = gti;
    //                            time3 += cti(t3e, t3b);

    cur_match_res.reserve(10);
//    cout << "begin extd" << endl;
    embedding_count = EvaluateQuery::exploreCECIStyle(data_graph,
                                                      query_graph[i],
                                                      ceci_tree,
                                                      candidates,
                                                      candidates_count,
                                                      TE_Candidates,
                                                      NTE_Candidates, ceci_order,
                                                      output_limit,
                                                      call_count, cur_match_res,
                                                      idx, idx_count, embedding,
                                                      valid_candidates,
                                                      visited_vertices);
//    cout << "end extd " << embedding_count << endl;
#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
     outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
     delete[] EvaluateQuery::distribution_count_;
#endif
}

bool BFS(Graph* data_graph, Graph** query_graph, int cnt)
{
    if (topk) return true;
    dfs_all_count++;
    if (cnt == QUERY_NUMS) {
        int sum = 0;
        int cnt_ = 0;
        vector<int> all_node;
        auto tmp = match_result;
        for (auto& i : tmp) {
            sort(i.second.begin(), i.second.end());
            for (auto j : i.second)
                sum += j;
        }
        sort(tmp.begin(), tmp.end());
        vector<vector<int>> nw_res;
        for (const auto& i : tmp) nw_res.emplace_back(i.second);

        if (his.find(sum) != his.end()) {
            auto it = his[sum];
            bool pa = true;
            int pa_cnt = 0;
            for (int k = 0; k < it.size(); ++k) {
                pa = true;
                for (int i = 0; i < QUERY_NUMS; ++i) {
                    if (nw_res[i].size() != it[k][i].size()) {
                        pa = false;
                        break;
                    }
                    for (int j = 0; j < nw_res[i].size(); ++j) {
                        if (nw_res[i][j] != it[k][i][j]) {
                            pa = false;
                            break;
                        }
                    }
                    if (!pa) break;
                }
                if (pa)
                    return true;
            }
        }
        his[sum].emplace_back(nw_res);

        auto sstart = std::chrono::high_resolution_clock::now();
        fill(pst.begin(), pst.end(), 0);
        for (const auto& i : match_result) {
            for (auto v : i.second) {
                if (pst[v] >= 0) pst[v] = -i.first - 1;
            }
        }
        // for (const auto& i : match_result) {
        //     cout << "< " << i.first << ": ";
        //     for (auto v : i.second) {
        //         cout << v << " ";
        //     }
        //     cout << "> ";
        // }
        // cout << endl;
        // int x; cin >> x;
        vector<vector<pair<int, int>>> Path;
        for (const auto& i : match_result) {
            vector<pair<int, int>> pa;
            for (auto v : i.second) {
                get_path3D(data_graph, i.first + 1, v, -1, 0, pa);
            }
            Path.emplace_back(pa);
            for (const auto& v : pa) {
                if (pst[v.first] > 0) pst[v.first] = 0;
                if (pst[v.second] > 0) pst[v.second] = 0;
            }
        }
        auto eend = std::chrono::high_resolution_clock::now();
        time4 += std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
        if (NC(total_time) >= 3600) {
            topk = true;
            cout << "time exced" << endl;
        }
        cur_sp_match_res++;
        if (cur_sp_match_res % 10000 == 0)
            cout <<  "asd " << cur_sp_match_res << " " << NC(total_time) << "ms " << get_index_mem()["pk"] << "kb" << endl;
        if (cur_sp_match_res >= top)
            topk = true;
        if (topk) return true;
        return true;
    }

    bool this_turn = false;

    vector<bool> st(data_graph->vertices_count_, false);
    vector<bool> in_queue(data_graph->vertices_count_, false);
    vector<int> dist(data_graph->vertices_count_, 0);

    queue<int> q;
    for (const auto& pp : match_result) {
        for (int p : pp.second)
        {q.emplace(p); in_queue[p] = true;}
    }

    while (!q.empty()) {
        int vertex = q.front(); q.pop();
        int d = dist[vertex];
        // if (vertex == 1907)
        //     cout << endl;
        if (st[vertex])
            continue;
        // 拓展开始
        for (int i = 1; i < QUERY_NUMS; ++i)
        {
            if (!query_is_matched[i]) { // 查询图i未匹配

                ui count = 0;

                auto itr = query_graph[i]->getVerticesByLabel(data_graph->getVertexLabel(vertex), count);
                try_count += count;
                if (count > 0) {
                    // vertex是否被匹配过

                    // 查询图i中这个点的候选集，可以通过度数排序之类
                    // 在数据图上进行该查询图的匹配
                    for (int j = 0; j < count; ++j) { // 遍历查询图i与vertex同类别的点的下标
                        vector<vector<int>> cur_match_res;

                        CECI(data_graph, query_graph, cur_match_res, i, itr[j], vertex);
                        if (cnt == 1 && i == 2 && !cur_match_res.empty())
                            cout << endl;
                        for (const auto& cur_match: cur_match_res) {
                            match_result.emplace_back(i, cur_match);
                            query_is_matched[i] = true;

                            if (BFS(data_graph, query_graph, cnt + 1)) {
                                this_turn = true;
                            }
                            if (topk) return true;
                            match_result.pop_back();
                            query_is_matched[i] = false;
                        }
                    }
                }
            }
        }
        st[vertex] = true;
        if (d + 1 > D) continue;
        // 邻居加入队列

        ui nbrs_cnt;
        const VertexID *nbrs = data_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (int ii = 0; ii < nbrs_cnt; ii++) {
            VertexID v = nbrs[ii];
            if (!st[v]) // csr寸的是点值
            {
                if (!in_queue[v]) // 不在队列中
                {
                    q.emplace(v);
                    in_queue[v] = true;
                    dist[v] = d + 1;
                }
            }
        }
    }
    return this_turn;
}
vector<bool> dfs_in_queue;
bool DFS(Graph* data_graph, Graph** query_graph, ui vertex, int cnt, bool from_vertex, bool need_match)
{
    if (topk) return true;
    dfs_all_count++;
    if (cnt == QUERY_NUMS)
    {
        int sum = 0;
        int cnt_ = 0;
        vector<int> all_node;
        auto tmp = match_result;
        for (auto& i : tmp) {
            sort(i.second.begin(), i.second.end());
            for (auto j : i.second)
                sum += j;
        }
        sort(tmp.begin(), tmp.end());
        vector<vector<int>> nw_res;
        for (auto i : tmp) nw_res.emplace_back(i.second);

        if (his.find(sum) != his.end()) {
            auto it = his[sum];
            bool pa = true;
            int pa_cnt = 0;
            for (int k = 0; k < it.size(); ++k) {
                pa = true;
                for (int i = 0; i < QUERY_NUMS; ++i) {
                    if (nw_res[i].size() != it[k][i].size()) {
                        pa = false;
                        break;
                    }
                    for (int j = 0; j < nw_res[i].size(); ++j) {
                        if (nw_res[i][j] != it[k][i][j]) {
                            pa = false;
                            break;
                        }
                    }
                    if (!pa) break;
                }
                if (pa)
                    return true;
            }
        }
        his[sum].emplace_back(nw_res);

        auto sstart = std::chrono::high_resolution_clock::now();
        fill(pst.begin(), pst.end(), 0);
        for (const auto& i : match_result) {
            for (auto v : i.second) {
                if (pst[v] >= 0) pst[v] = -i.first - 1;
            }
        }
        vector<vector<pair<int, int>>> Path;
        for (const auto& i : match_result) {
            vector<pair<int, int>> pa;
            for (auto v : i.second) {
                get_path3D(data_graph, i.first + 1, v, -1, 0, pa);
            }
            Path.emplace_back(pa);
            for (const auto& v : pa) {
                if (pst[v.first] > 0) pst[v.first] = 0;
                if (pst[v.second] > 0) pst[v.second] = 0;
            }
        }
        auto eend = std::chrono::high_resolution_clock::now();
        time4 += std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
        if (NC(total_time) >= 3600) {
            topk = true;
            cout << "time exced" << endl;
        }
        cur_sp_match_res++;
        if (cur_sp_match_res % 10000 == 0)
            cout <<  "asd " << cur_sp_match_res << " " << NC(total_time) << "ms " << get_index_mem()["pk"] << "kb" << endl;
        if (cur_sp_match_res >= top)
            topk = true;
        if (topk) return true;
        return true;
    }

    bool this_turn = false;
    if (!from_vertex) {
        for (const auto& cur : match_result) {
            for (auto ve : cur.second) {
                if (dfs_st[ve]) continue;
                for (int i = 0; i < QUERY_NUMS; ++i) {
                    if (query_is_matched[i]) continue;
                    ui count;
                    auto itr = query_graph[i]->getVerticesByLabel(data_graph->getVertexLabel(ve), count);
                    try_count += count;
                    if (count > 0) {
                        for (int j = 0; j < count; ++j) { // 遍历查询图i与vertex同类别的点的下标
                            if (data_graph->getVertexDegree(vertex) < query_graph[i]->getVertexDegree(itr[j])) continue;
                            const std::unordered_map<LabelID, ui>* u_nlf = query_graph[i]->getVertexNLF(itr[j]);
                            const std::unordered_map<LabelID, ui>* v_nlf = data_graph->getVertexNLF(ve);

                            if (v_nlf->size() >= u_nlf->size()) {
                                bool is_valid = true;
                                for (auto element: *u_nlf) {
                                    auto iter = v_nlf->find(element.first);
                                    if (iter == v_nlf->end() || iter->second < element.second) {
                                        is_valid = false;
                                        break;
                                    }
                                }
                                if (!is_valid) continue;
                            } else {
                                continue;
                            }
                            for (int tree_i = 0; tree_i < mx; ++tree_i) {
                                ceci_tree[tree_i].clear();
                            }
                            std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
                            bool filter_res;

                            filter_res = FilterVertices::CECIFilter(data_graph, query_graph[i], candidates,
                                                                    candidates_count, ceci_order,ceci_tree,
                                                                    TE_Candidates,
                                                                    NTE_Candidates,itr[j], ve,
                                                                    nullptr, SPECIAL);
                            if (!filter_res) {
                                continue;
                            }


                            // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
                            std::vector<ui> optimal_candidates_count;
        double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                              candidates_count, optimal_candidates_count);
        FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
//                        std::cout << "-----" << std::endl;
//                        std::cout << "Build indices..." << std::endl;

//                            size_t memory_cost_in_bytes = 0;

                            if (TE_Candidates.empty()) {
                                continue;
                            }
//                            memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph[i], candidates_count,ceci_order,ceci_tree,TE_Candidates, NTE_Candidates);
//        BuildTable::printTableCardinality(query_graph[current], ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);

//                        std::cout << "-----" << std::endl;
//                        std::cout << "Generate a matching order..." << std::endl;


//                            ui *matching_order = NULL;
//                            ui *pivots = NULL;

                            size_t order_num = 0;

                            sscanf(input_order_num.c_str(), "%zu", &order_num);

                            std::vector<std::vector<ui>> spectrum;

                            GenerateQueryPlan::generateCECIQueryPlan(query_graph[i], ceci_tree, ceci_order, matching_order,
                                                                     pivots);

                            if (input_order_type != "Spectrum") {
                                GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[i], matching_order, pivots);
//                            GenerateQueryPlan::printSimplifiedQueryPlan(query_graph[0], matching_order);
                            } else {
                                std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
                            }

//                        std::cout << "-----" << std::endl;
//                        std::cout << "Enumerate..." << std::endl;
                            size_t output_limit = 0;
                            size_t embedding_count = 0;
                            if (input_max_embedding_num == "MAX") {
                                output_limit = std::numeric_limits<size_t>::max();
                            } else {
                                sscanf(input_max_embedding_num.c_str(),"%zu", &output_limit);
                            }

#if ENABLE_QFLITER == 1
                            EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

                            size_t call_count = 0;
                            size_t time_limit = 0;
                            sscanf(input_time_limit.c_str(), "%zu", &time_limit);

                            vector<vector<int>> cur_match_res;
                            cur_match_res.reserve(10);
                            embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph[i], ceci_tree,
                                                                              candidates,
                                                                              candidates_count, TE_Candidates,
                                                                              NTE_Candidates, ceci_order, output_limit,
                                                                              call_count, cur_match_res,idx,idx_count,embedding,valid_candidates,visited_vertices);
                            //                            memory_cost_in_bytes += sizeof(int) * cur_match_res.size();


#ifdef DISTRIBUTION
                            std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
        outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
        delete[] EvaluateQuery::distribution_count_;
#endif

                            for (const auto& cur_match: cur_match_res) {
                                match_result.emplace_back(i, cur_match);
                                query_is_matched[i] = true;
                                auto dfs_st_bk = dfs_st;
                                auto in_queue_bk = dfs_in_queue;
                                auto dist_bk = dfs_dist;
                                fill(dfs_dist.begin(), dfs_dist.end(), 0);
                                for (const auto& iii : match_result) {for (auto jjj : iii.second) {dfs_in_queue[jjj] = true;}}
                                if (DFS(data_graph, query_graph, ve, cnt + 1, false, true)) {
                                    this_turn = true;
                                }
                                if (topk) return true;
                                dfs_st = dfs_st_bk;
                                dfs_dist = dist_bk;
                                dfs_in_queue = in_queue_bk;
                                match_result.pop_back();
                                query_is_matched[i] = false;
                            }
                        }
                    }
                }
            }
        }

        for (const auto& i : match_result) {
            for (auto ve : i.second) {
                if (dfs_dist[ve] + 1 > D) continue;
                ui nbrs_cnt;
                const VertexID *nbrs = data_graph->getVertexNeighbors(ve, nbrs_cnt);
                for (int ii = 0; ii < nbrs_cnt; ii += SPECIAL) {
                    VertexID v = nbrs[ii];
                    if (!dfs_st[v]) // csr寸的是点值
                    {
                        assert(dfs_dist[v] == 0);
                        if (!dfs_in_queue[v]) // 不在队列中
                        {
                            dfs_dist[v] = dfs_dist[ve] + 1;
                            if (DFS(data_graph, query_graph, v, cnt, true, true))
                                this_turn = true;
                            if (topk) return true;

                        }
                    }
                    else {
                        if (!dfs_in_queue[v]) {
                            if (dfs_dist[ve] + 1 < dfs_dist[v] || dfs_dist[v] == 0) {
                                dfs_dist[v] = dfs_dist[ve] + 1;
                                if (DFS(data_graph, query_graph, v, cnt, true, false)) {
                                    this_turn = true;
                                }
                                if (topk) return true;

                            }
                        }
                    }
                }
            }
        }
        dfs_st[vertex] = true;
        return this_turn;
    }
    this_turn = false;
    if (from_vertex) {
        if (need_match) {
            for (int i = 1; i < QUERY_NUMS; ++i) {
                if (query_is_matched[i]) continue;
                ui count;
                auto itr = query_graph[i]->getVerticesByLabel(data_graph->getVertexLabel(vertex), count);
                bool this_vertex = false;
                try_count += count;
                if (count > 0) {
                    for (int j = 0; j < count; ++j) { // 遍历查询图i与vertex同类别的点的下标

//                    if (nei_st.find(itr[j]) != nei_st.end()) continue;
//                    nei_st[itr[j]] = true;
                        if (data_graph->getVertexDegree(vertex) < query_graph[i]->getVertexDegree(itr[j])) continue;
                        const std::unordered_map<LabelID, ui>* u_nlf = query_graph[i]->getVertexNLF(itr[j]);
                        const std::unordered_map<LabelID, ui>* v_nlf = data_graph->getVertexNLF(vertex);

                        if (v_nlf->size() >= u_nlf->size()) {
                            bool is_valid = true;
                            for (auto element: *u_nlf) {
                                auto iter = v_nlf->find(element.first);
                                if (iter == v_nlf->end() || iter->second < element.second) {
                                    is_valid = false;
                                    break;
                                }
                            }
                            if (!is_valid) continue;
                        } else {
                            continue;
                        }
//                    auto start = std::chrono::high_resolution_clock::now();

//                    ui **candidates = nullptr; // candidates[i][j] the i-th node's j-th candidate node
//                    ui *candidates_count = nullptr;
//
//                    TreeNode *ceci_tree = nullptr;
//                    ui *ceci_order = nullptr;
                        for (int tree_i = 0; tree_i < mx; ++tree_i) {
                            ceci_tree[tree_i].clear();
                        }
                        std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
                        std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
                        bool filter_res;

                        filter_res = FilterVertices::CECIFilter(data_graph, query_graph[i], candidates,
                                                                candidates_count, ceci_order,
                                                                ceci_tree, TE_Candidates, NTE_Candidates, itr[j], vertex,
                                                                nullptr, SPECIAL);
                        if (!filter_res) {
                            continue;
                        }

//                    auto end = std::chrono::high_resolution_clock::now();
//                    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

                        // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
                        std::vector<ui> optimal_candidates_count;
        double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                              candidates_count, optimal_candidates_count);
        FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
//                        std::cout << "-----" << std::endl;
//                        std::cout << "Build indices..." << std::endl;

//                    size_t memory_cost_in_bytes = 0;

                        if (TE_Candidates.empty()) {
                            continue;
                        }

                        size_t order_num = 0;

                        sscanf(input_order_num.c_str(), "%zu", &order_num);

                        std::vector<std::vector<ui>> spectrum;

                        GenerateQueryPlan::generateCECIQueryPlan(query_graph[i], ceci_tree, ceci_order, matching_order,
                                                                 pivots);

                        if (input_order_type != "Spectrum") {
                            GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[i], matching_order, pivots);
//                            GenerateQueryPlan::printSimplifiedQueryPlan(query_graph[0], matching_order);
                        } else {
                            std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
                        }

//                        std::cout << "-----" << std::endl;
//                        std::cout << "Enumerate..." << std::endl;
                        size_t output_limit = 0;
                        size_t embedding_count = 0;
                        if (input_max_embedding_num == "MAX") {
                            output_limit = std::numeric_limits<size_t>::max();
                        } else {
                            sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
                        }

#if ENABLE_QFLITER == 1
                        EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

                        size_t call_count = 0;
                        size_t time_limit = 0;
                        sscanf(input_time_limit.c_str(), "%zu", &time_limit);

//                    start = std::chrono::high_resolution_clock::now();
                        vector<vector<int>> cur_match_res;
                        cur_match_res.reserve(100);
                        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph[i], ceci_tree,
                                                                          candidates,
                                                                          candidates_count, TE_Candidates,
                                                                          NTE_Candidates, ceci_order, output_limit,
                                                                          call_count, cur_match_res,idx,idx_count,embedding,valid_candidates,visited_vertices);

#ifdef DISTRIBUTION
                        std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
        outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
        delete[] EvaluateQuery::distribution_count_;
#endif

                        for (const auto& cur_match: cur_match_res) {
                            match_result.emplace_back(i, cur_match);
                            query_is_matched[i] = true;

                            auto dfs_st_bk = dfs_st;
                            auto in_queue_bk = dfs_in_queue;
                            auto dist_bk = dfs_dist;
                            fill(dfs_dist.begin(), dfs_dist.end(), 0);

                            for (const auto& iii : match_result) {for (auto jjj : iii.second) {dfs_in_queue[jjj] = true;}}
                            if (DFS(data_graph, query_graph, vertex, cnt + 1, false, true)) {
                                this_turn = true;
                            }
                            if (topk) return true;

                            dfs_st = dfs_st_bk;
                            dfs_dist = dist_bk;
                            dfs_in_queue = in_queue_bk;
                            match_result.pop_back();
                            query_is_matched[i] = false;
                        }
                    }
                }
            }
        }

        dfs_st[vertex] = true;
        if (dfs_dist[vertex] + 1 > D) return this_turn;
        ui nbrs_cnt;
        const VertexID *nbrs = data_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (int ii = 0; ii < nbrs_cnt; ii++) {
            VertexID v = nbrs[ii];
            if (!dfs_st[v]) // csr寸的是点值
            {
                assert(dfs_dist[v] == 0);
                if (!dfs_in_queue[v]) // 不在队列中
                {
                    dfs_dist[v] = dfs_dist[vertex] + 1;
                    if (DFS(data_graph, query_graph, v, cnt, true, true))
                        this_turn = true;
                    if (topk) return true;
                }
            }
            else {
                if (!dfs_in_queue[v]) {
                    if (dfs_dist[vertex] + 1 < dfs_dist[v] || dfs_dist[v] == 0) {
                        dfs_dist[v] = dfs_dist[vertex] + 1;
                        if (DFS(data_graph, query_graph, v, cnt, true, false)) {
                            this_turn = true;
                        }
                        if (topk) return true;
                    }
                }
            }
        }
        return this_turn;
    }
}
void bfs_solve(double total_time_in_ns, vector<vector<int>>& all_match_res, Graph** query_graph, Graph* data_graph){
    auto sstart = std::chrono::high_resolution_clock::now();

    std::cout << "Begin SP Match!" << std::endl;
    cout << "PATH: " << PATH << " Query_Num: " << QUERY_NUMS << " D: " << D << " Meth: bfs" << dl;

    for (int i = 0; i < all_match_res.size(); ++i) {
        his.clear();
        auto each_begin = std::chrono::high_resolution_clock::now();
        query_is_matched[0] = true;
        match_result.emplace_back(0, all_match_res[i]);
        BFS(data_graph, query_graph, 1);
        auto each_end = std::chrono::high_resolution_clock::now();
        double each_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - each_begin).count();
        cout << i << " " << cur_sp_match_res << " " << NC(each_time) << " " << NC(time1) << " " << NC(time2) << " " << NC(time3) <<endl;
        match_result.clear();
//        sp_match_res += cur_sp_match_res;
//        cur_sp_match_res = 0;
        if (topk) break;
    }
    auto eend = std::chrono::high_resolution_clock::now();
    double MatchTime = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
    double totime = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
    printf("#Embeddings: %zu\n", cur_sp_match_res);
    printf("DFS call count: %zu\n", dfs_all_count);
    printf("try count: %zu\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(MatchTime));
//    printf("Avg queue size: %llu\n", sizeof(int) * (avg_res + avg_res2) / calculate_count);

    std::cout << "End." << std::endl;
}
void dfs_solve(double total_time_in_ns, vector<vector<int>>& all_match_res, Graph** query_graph, Graph* data_graph){
    std::cout << "Begin SP Match!" << std::endl;
    auto sstart = std::chrono::high_resolution_clock::now();

    dfs_st.resize(data_graph->vertices_count_);
    dfs_dist.resize(data_graph->vertices_count_);
    dfs_in_queue.resize(data_graph->vertices_count_);
    int i = 0;
    for (auto & all_match_re : all_match_res) {
        topk = false;
        his.clear();
        fill(dfs_st.begin(), dfs_st.end(), 0);
        fill(dfs_dist.begin(), dfs_dist.end(), 0);
        fill(dfs_in_queue.begin(), dfs_in_queue.end(), 0);
        match_result.clear();
        query_is_matched[0] = true;
        match_result.emplace_back(0, all_match_re);
        for (auto j : all_match_re) dfs_in_queue[j] = true;
        auto each_begin = std::chrono::high_resolution_clock::now();
        DFS(data_graph, query_graph, match_result[0].second[0], 1, false, true);
        auto each_end = std::chrono::high_resolution_clock::now();
        double each_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - each_begin).count();
//        sp_match_res += cur_sp_match_res;
//        cur_sp_match_res = 0;
        cout << i++ << " " << NC(each_time) << " " << get_index_mem()["pk"] - memory[0] << endl;
        if (topk) break;
    }
    auto eend = std::chrono::high_resolution_clock::now();
    double MatchTime = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
    double totime = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
    printf("#Embeddings: %zu\n", cur_sp_match_res);
    printf("DFS call count: %zu\n", dfs_all_count);
    printf("try count: %zu\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(MatchTime));

    std::cout << "End." << std::endl;
}

int khop;

map<string, vector<int>> biasMp;
vector<pair<int, int>> queryR;

bool begin_extend(Graph* data_graph, Graph** query_graph, int cnt) {
    if (topk) return true;
    dfs_all_count++;
    if (cnt == QUERY_NUMS)
    {
        int sum = 0;
        int cnt_ = 0;
        vector<int> all_node;
        auto tmp = match_result;
        for (auto& i : tmp) {
            sort(i.second.begin(), i.second.end());
            for (auto j : i.second)
                sum += j;
        }
        sort(tmp.begin(), tmp.end());
        vector<vector<int>> nw_res;
        for (const auto& i : tmp) nw_res.emplace_back(i.second);

        if (his.find(sum) != his.end()) {
            auto it = his[sum];
            bool pa = true;
            int pa_cnt = 0;
            for (int k = 0; k < it.size(); ++k) {
                pa = true;
                for (int i = 0; i < QUERY_NUMS; ++i) {
                    if (nw_res[i].size() != it[k][i].size()) {
                        pa = false;
                        break;
                    }
                    for (int j = 0; j < nw_res[i].size(); ++j) {
                        if (nw_res[i][j] != it[k][i][j]) {
                            pa = false;
                            break;
                        }
                    }
                    if (!pa) break;
                }
                if (pa)
                    return true;
            }
        }
        his[sum].emplace_back(nw_res);
// for (auto i : match_result) {
//            auto ma = i.second;
//            sort(ma.begin(), ma.end());
//            for (auto j : ma)
//                cout << j << " ";
//            cout << " ";
//        }
//        cout << endl;
        auto sstart = std::chrono::high_resolution_clock::now();
        fill(pst.begin(), pst.end(), 0);
        for (const auto& i : match_result) {
            for (auto v : i.second) {
                if (pst[v] >= 0) pst[v] = -i.first - 1;
            }
        }
        vector<vector<pair<int, int>>> Path;
        for (const auto& i : match_result) {
            vector<pair<int, int>> pa;
            for (auto v : i.second) {
                get_path3D(data_graph, i.first + 1, v, -1, 0, pa);
            }
            Path.emplace_back(pa);
            for (const auto& v : pa) {
                if (pst[v.first] > 0) pst[v.first] = 0;
                if (pst[v.second] > 0) pst[v.second] = 0;
            }
        }

        // for (const auto& i : match_result) {
        //     cout << "< " << i.first << ": ";
        //     for (auto v : i.second) {
        //         cout << v << " ";
        //     }
        //     cout << "> ";
        // }
        // cout << endl;
        auto eend = std::chrono::high_resolution_clock::now();
        time4 += std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - total_start).count();
        if (NC(total_time) >= 3600) {
            topk = true;
            cout << "time exced" << endl;
        }
        cur_sp_match_res++;
        if (cur_sp_match_res % 10000 == 0)
            cout <<  "asd " << cur_sp_match_res << " " << NC(total_time) << "ms " << get_index_mem()["pk"] << "kb" << endl;
        if (cur_sp_match_res >= top)
            topk = true;
        if (topk) return true;
        return true;
    }

    bool this_turn = false;

    vector<bool> st(data_graph->vertices_count_, false);
    vector<bool> in_queue(data_graph->vertices_count_, false);
    vector<bool> canNextVertex(data_graph->vertices_count_, false);
    vector<int> dist(data_graph->vertices_count_, 0);
    vector<bool> canMatch(QUERY_NUMS, false);
    // queue<int> q1;
    queue<int> q;
    for (const auto& pp : match_result) {
        for (int p : pp.second)
        {q.emplace(p); in_queue[p] = true; canNextVertex[p] = true;}
    }
    Graph* subgraph = new Graph(true, false);
    subgraph->move(data_graph);
    // Graph subgraph(*data_graph);

    memset(subgraph->offsets_, 0, (subgraph->vertices_count_ + 1) * sizeof(int));

    while (!q.empty()) {
        int vertex = q.front(); q.pop();
        int d = dist[vertex];

        if (st[vertex]) {
            // q.pop();
            continue;
        }
        st[vertex] = true;
        if (d > D + queryR.back().first) {
            continue;
        }
        // q.pop();
        // 邻居加入队列
        assert(vertex < subgraph->vertices_count_);
        // cout << vertex << endl;
        subgraph->offsets_[vertex] = data_graph->offsets_[vertex];
        subgraph->offsets_[vertex + 1] = data_graph->offsets_[vertex + 1];

        ui nbrs_cnt;
        const VertexID *nbrs = data_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (int ii = 0; ii < nbrs_cnt; ii++) {
            VertexID v = nbrs[ii];
            if (!st[v]) // csr寸的是点值 不在队列中
            {
                if ( !in_queue[v]) {
                    in_queue[v] = true;
                    dist[v] = d + 1;
                    q.emplace(v);
                    if (dist[v] <= D) {
                        for (int jj = 1; jj < QUERY_NUMS; ++jj) {
                            if (!query_is_matched[jj]) {
                                ui count = 0;
                                assert(v < data_graph->vertices_count_);
                                auto itr = query_graph[jj]->getVerticesByLabel(data_graph->getVertexLabel(v), count);
                                if (count > 0) {
                                    canMatch[jj] = true;
                                    canNextVertex[v] = true;
                                    // cout << jj << " " << v << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 1; i <= subgraph->vertices_count_; i++) {
        if (subgraph->offsets_[i] == 0)
            subgraph->offsets_[i] = subgraph->offsets_[i - 1];
    }

    for (int i = 1; i < QUERY_NUMS; ++i) {
        if (query_is_matched[i]) continue;
        if (!canMatch[i]) continue;

        vector<vector<int>> cur_match_res;
        CECI(subgraph, query_graph, cur_match_res, i, -1, -1);
        // tmp_match_res[i] = cur_match_res;
        query_is_matched[i] = true;

        for (const auto eachMatch : cur_match_res) {
            bool canNext = false;
            for (auto ve : eachMatch) {
                if (canNextVertex[ve]) {
                    canNext = true;
                    break;
                }
            }
            if (!canNext) continue;
            match_result.emplace_back(i, eachMatch);
            begin_extend(data_graph, query_graph, cnt + 1);
            match_result.pop_back();
            if (topk) return true;
        }
        query_is_matched[i] = false;
    }
    if (topk) return true;
    delete subgraph;

    /*for (int i = 0; i < queryR.size(); i++) {
        if (query_is_matched[queryR[i].second]) continue;
        if (!canMatch[queryR[i].second]) continue;
        while (!q.empty()) {
            int vertex = q.front();
            int d = dist[vertex];
            if (st[vertex]) {
                q.pop();
                continue;
            }

            /#1#/ 拓展开始
            for (int i = 1; i < QUERY_NUMS; ++i)
            {
                if (!query_is_matched[i]) { // 查询图i未匹配
                    ui count = 0;

                    auto itr = query_graph[i]->getVerticesByLabel(data_graph->getVertexLabel(vertex), count);
                    try_count += count;
                    if (count > 0) {
                        // vertex是否被匹配过

                        // 查询图i中这个点的候选集，可以通过度数排序之类
                        // 在数据图上进行该查询图的匹配
                        for (int j = 0; j < count; ++j) { // 遍历查询图i与vertex同类别的点的下标
                            vector<vector<int>> cur_match_res;

                            CECI(data_graph, query_graph, cur_match_res, i, itr[j], vertex);

                            for (const auto& cur_match: cur_match_res) {
                                match_result.emplace_back(i, cur_match);
                                query_is_matched[i] = true;
                                if (begin_extend(data_graph, query_graph, cnt + 1)) {
                                    this_turn = true;
                                }
                                if (topk) return true;
                                match_result.pop_back();
                                query_is_matched[i] = false;
                            }
                        }
                    }
                }
            }
            st[vertex] = true;
            if (d > D + queryR[i].first) {
                break;
                // continue;
            }
            q.pop();
            // 邻居加入队列
            assert(vertex < subgraph->vertices_count_);
            // cout << vertex << endl;
            subgraph->offsets_[vertex] = data_graph->offsets_[vertex];
            subgraph->offsets_[vertex + 1] = data_graph->offsets_[vertex + 1];

            ui nbrs_cnt;
            const VertexID *nbrs = data_graph->getVertexNeighbors(vertex, nbrs_cnt);
            for (int ii = 0; ii < nbrs_cnt; ii++) {
                VertexID v = nbrs[ii];
                if (!st[v]) // csr寸的是点值
                {
                    if (!in_queue[v]) // 不在队列中
                    {
                        q.emplace(v);
                        in_queue[v] = true;
                        dist[v] = d + 1;
                    }
                }
            }
        }
        vector<vector<int>> cur_match_res;
        CECI(subgraph, query_graph, cur_match_res, queryR[i].second, -1, -1);
        // tmp_match_res[i] = cur_match_res;
        query_is_matched[queryR[i].second] = true;

        for (const auto eachMatch : cur_match_res) {
            bool canNext = false;
            for (auto ve : eachMatch) {
                if (canNextVertex[ve]) {
                    canNext = true;
                    break;
                }
            }
            // if (!canNext) continue;
            match_result.emplace_back(queryR[i].second, eachMatch);
            begin_extend(data_graph, query_graph, cnt + 1);
            match_result.pop_back();
        }
        query_is_matched[queryR[i].second] = false;
    }*/

    // vector<vector<vector<int>>> tmp_match_res(QUERY_NUMS);
    //
    // for (int i = 1; i < QUERY_NUMS; ++i) {
    //     if (!query_is_matched[i]) {
    //         // 查询图i未匹配
    //         vector<vector<int>> cur_match_res;
    //         CECI(data_graph, query_graph, cur_match_res, i, -1, -1);
    //         tmp_match_res[i] = cur_match_res;
    //     }
    // }

    return this_turn;
}

void subgraph_solve(double total_time_in_ns, vector<vector<int>>& all_match_res, Graph** query_graph, Graph* data_graph) {
    auto sstart = std::chrono::high_resolution_clock::now();
    std::cout << "Begin SP Match!" << std::endl;
    // top = 100;
    cout << "PATH: " << PATH << " Query_Num: " << QUERY_NUMS << " D: " << D << " Meth: subgraph" << " top: " << top << dl;
    // pst = new int[data_graph->vertices_count_];
    // pst2 = new int[data_graph->vertices_count_];
    int epoch = 0;
//    get_nei_label_vertex.resize(data_graph->vertices_count_);
    cout << "max_size " << (get_index_mem()["pk"] - memory[0]) << "KB" << endl;
    khop=D+1;
    // for (auto & all_match_re : all_match_res) {
    for (int i = 0; i < all_match_res.size(); i++){
        his.clear();
        auto & all_match_re = all_match_res[i];
        match_result.clear();
        match_result.emplace_back(0, all_match_re);
        auto each_begin = std::chrono::high_resolution_clock::now();
        // topk = false;

        begin_extend(data_graph, query_graph, 1);

        auto each_end = std::chrono::high_resolution_clock::now();
        double each_time = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end - each_begin).count();

        cout << "each " << i << " " << cur_sp_match_res << " " << NC(each_time) << "ms " << endl;
        epoch++;
    }
    sp_match_res = cur_sp_match_res;
    auto eend = std::chrono::high_resolution_clock::now();
    double ALL_TIME = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
    sum_memory_cost_in_bytes += sizeof (bool)* data_graph->vertices_count_;
    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    printf("Memory cost (MB): %.4lf, %.4lf\n", BYTESTOMB(all_memory + sum_memory_cost_in_bytes + init_memory_cost_in_bytes), BYTESTOMB(init_memory_cost_in_bytes));
    printf("#Embeddings: %d\n", sp_match_res);
    printf("DFS call count: %d\n", dfs_all_count);
    printf("try count: %d\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(ALL_TIME + total_time_in_ns));
    cout << NC(ALL_TIME) << " " << NC(total_time_in_ns) << " " << NC(tim1) << " " << NC(tim2) << endl;
    std::cout << "End." << std::endl;
    cout << "max_size " << (get_index_mem()["pk"] - memory[0]) << "KB" << endl;
    // delete [] pst;
    // delete [] pst2;
}

int main(int argc, char** argv) {
    cout << getTime() << endl;
    test = false;
    auto mp = get_index_mem();
    memory.emplace_back(mp["pk"]);
//  CL-10M-1d8-L5

    QUE = "1";
    PATH = "MSRC-21C";
    D = 1;
    QUERY_NUMS = 6;
    Meth = 0;
//    input_query_graph_file = "../../../test/querys/" + QUE;
//    input_data_graph_file = "../../../test/graphs/Dataset/" + PATH + "/" + PATH + ".txt";

    MatchingCommand command(argc, argv);
    QUE = command.getQueryGraphFilePath();
    PATH = command.getDataGraphFilePath();
    D = stoi(command.getDist());
    QUERY_NUMS = stoi(command.getQueryNums());
    Meth = stoi(command.getMeth());
    top = 100;
    input_query_graph_file = "/home/dbia-graph/qsl/FPM/q/" + QUE;
    input_data_graph_file = "/home/dbia-graph/qsl/Dataset/FPMData/" + PATH + "/" + PATH + ".txt";
    double value = 3600;
    time_limit["CL"]["1"]["2"] = value;
    time_limit["CL"]["1"]["3"] = value;
    time_limit["CL"]["1"]["4"] = value;
    time_limit["CL"]["2"]["2"] = value;
    time_limit["CL"]["2"]["3"] = value;
    time_limit["CL"]["2"]["4"] = value;

    SPECIAL = 1;
    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph CSR: " << input_csr_file_path << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tEngine Type: " << input_engine_type << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tOrder Num: " << input_order_num << std::endl;
    std::cout << "\tDistribution File Path: " << input_distribution_file_path << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    query_is_matched = new bool[QUERY_NUMS];
    for (int i = 0; i < QUERY_NUMS; ++i) query_is_matched[i] = false;
    is_show = new bool[QUERY_NUMS];

    /**
     * Load input graphs.
     */

    std::cout << "Load graphs..." << std::endl;

//    Graph* query_graph = new Graph(true);
//    query_graph->loadGraphFromFile(input_query_graph_file);
//    query_graph->buildCoreTable();
    bool has_dir_ = false;
    if (PATH == "CL") has_dir_ = true;
    Graph** query_graph = new Graph*[QUERY_NUMS];

    for (int i = 0; i < QUERY_NUMS; ++i) {
        query_graph[i] = new Graph(true, has_dir_);
        query_graph[i]->id_ = i;
        query_graph[i]->loadGraphFromFile(input_query_graph_file + "/q" + to_string(i) + ".txt");
//        query_graph[i]->loadGraphFromFile("../../test/querys/4/q0.txt");
        // query_graph[i]->buildCoreTable();
        mx = max(mx, (int)query_graph[i]->vertices_count_);
        all_query_vertex_nums += query_graph[i]->vertices_count_;
    }

    Graph* data_graph = new Graph(true, has_dir_);

    if (input_csr_file_path.empty()) {
        cout << "Load csr file" << endl;
        data_graph->loadGraphFromFile(input_data_graph_file);
    }
    else {
        cout << "NOT Load csr file" << endl;
        std::string degree_file_path = input_csr_file_path + "_deg.bin";
        std::string edge_file_path = input_csr_file_path + "_adj.bin";
        std::string label_file_path = input_csr_file_path + "_label.bin";
        data_graph->loadGraphFromFileCompressed(degree_file_path, edge_file_path, label_file_path);
    }
    mp = get_index_mem();
    memory.emplace_back(mp["pk"]);

    biasMp["MSRC-21C"] = {3, 1, 2, 5, 2, 3};
    biasMp["MSRC-21"] = {2, 5, 1, 6, 2, 2};
    biasMp["DD"] = {2, 1, 3, 4, 2, 2};
    for (int i = 1; i < QUERY_NUMS; i++) {
        queryR.emplace_back(biasMp[PATH][i], i);
    }
    sort(queryR.begin(), queryR.end());

//    auto end1 = std::chrono::high_resolution_clock::now();
//    double load_graphs_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start).count();

//    if (test) {
////    cst2 = new bool[data_graph->vertices_count_];
////    cst3 = new bool[data_graph->vertices_count_];
//
////        for (const auto& i : all_match_res) {
////            for (auto j : i) {
////                memset(cst, 0, sizeof(bool) * (data_graph->vertices_count_));
////                memset(cst2, 0, sizeof(bool) * (data_graph->vertices_count_));
////                memset(cst3, 0, sizeof(bool) * (data_graph->vertices_count_));
//////            for (auto j: all_match_re) cst[j] = true;
//////            auto v = all_match_re[0];
//////            cout << i << endl;
////                ppa.clear();
////                ppa2.clear();
////                dfs_get(data_graph,  j, j, -1, 0);
////                EDGE.clear();
////            }
////        }
//        cst = new bool[data_graph->vertices_count_];
////        cst2 = new bool[data_graph->vertices_count_];
////        cst3 = new bool[data_graph->vertices_count_];
//        for (int i = 0; i < data_graph->vertices_count_; i ++) {
//            memset(cst, 0, sizeof(bool) * (data_graph->vertices_count_));
////            memset(cst2, 0, sizeof(bool) * (data_graph->vertices_count_));
////            memset(cst3, 0, sizeof(bool) * (data_graph->vertices_count_));
////            for (auto j: all_match_re) cst[j] = true;
////            auto v = all_match_re[0];
////            cout << i << endl;
//            if (i % 5000 == 0)
//                cout << i << endl;
//            ppa.clear();
//            ppa2.clear();
//            dfs(data_graph, query_graph, i, i, -1, 0, 0);
//            EDGE.clear();
//        }
//        return 0;
//    }
    auto start = std::chrono::high_resolution_clock::now();
    total_start = std::chrono::high_resolution_clock::now();
//    cst = new bool[data_graph->vertices_count_];
    idx = new ui[mx];
    idx_count = new ui[mx];
    embedding = new ui[mx];
    valid_candidates.resize(mx);
    visited_vertices = new bool[data_graph->vertices_count_];

//    cst = new bool[data_graph->vertices_count_];
//    cst3 = new bool[data_graph->vertices_count_];
//    bst = new int[data_graph->vertices_count_];
    pst.resize(data_graph->vertices_count_);
    if (Meth == 0) {
        base_solve(0, query_graph, data_graph);
        delete[] is_show;
        delete[] query_is_matched;
//        delete[] cst;
        for (int i = 0; i < QueryNums; ++i)
            delete query_graph[i];
        delete[] query_graph;
        delete data_graph;
        return 0;
    }
//    dont_need_match = new bool*[data_graph->vertices_count_ + 1];
//    for (int i = 0; i <= data_graph->vertices_count_; ++i)
//        dont_need_match[i] = new bool [QUERY_NUMS];

    candidates_count = new VertexID[mx];
    ceci_order = new VertexID[mx];
    ceci_tree = new TreeNode[mx];
    for (ui i = 0; i < mx; ++i) {
        ceci_tree[i].initialize(mx);
    }
    matching_order = new ui[mx];
    pivots = new ui[mx];
    candidates = new ui*[mx];
    for (int i = 0; i < mx; ++i) {
        candidates[i] = new ui[data_graph->vertices_count_];
    }

//    sort(query_graph, query_graph + QUERY_NUMS, [](const Graph* a, const Graph* b) {
//        if (a->vertices_count_ == b->vertices_count_) {
//            if (a->edges_count_ == b->edges_count_)
//                return a->l_count > b->l_count;
//            else
//                return a->edges_count_ > b->edges_count_;
//        }
//        else
//            return a->vertices_count_ > b->vertices_count_;
//    });

    cout << "Query order: ";
    for (int i = 0; i < QUERY_NUMS; ++i)
        cout << query_graph[i]->id_ << " ";
    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Metha Information" << std::endl;
    for (int i = 0; i < QUERY_NUMS; ++i)
        query_graph[i]->printGraphMetaData();
//    query_graph->printGraphMethaData();
    std::cout << "-----" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;
    match_result.reserve(QUERY_NUMS);

    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;
    ui max_vetex_nums = 0;
    for (int i = 0; i < QUERY_NUMS; ++i) { if (query_graph[i]->vertices_count_ > max_vetex_nums) {max_vetex_nums = query_graph[i]->vertices_count_;}}

    vector<vector<int>> all_match_res;
    all_match_res.reserve(1000);

//    ui **candidates = NULL; // candidates[i][j] the i-th node's j-th candidate node
//    ui *candidates_count = NULL;

//    TreeNode *ceci_tree = NULL;
//    ui *ceci_order = NULL;

    GenerateFilteringPlan::generateCECIFilterPlan(data_graph, query_graph[0], ceci_tree, ceci_order, -1);

    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    bool filter_res;
    double tmp;
    filter_res = FilterVertices::CECIFilter(data_graph, query_graph[0], candidates, candidates_count, ceci_order,
                                            ceci_tree, TE_Candidates, NTE_Candidates, -1, -1, nullptr, SPECIAL);

    if (!filter_res) {
        cout << "filter not match" << endl;
        return 0;
    }

    // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
    std::vector<ui> optimal_candidates_count;
        double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                              candidates_count, optimal_candidates_count);
        FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
    std::cout << "-----" << std::endl;
    std::cout << "Build indices..." << std::endl;

    Edges ***edge_matrix = NULL;

    size_t memory_cost_in_bytes = 0;

    if (TE_Candidates.empty()) {
        cout << "no result" << endl;
        return 0;
    }
    memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph[0], candidates_count, ceci_order,
                                                                ceci_tree,
                                                                TE_Candidates, NTE_Candidates);
//        BuildTable::printTableCardinality(query_graph[current], ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);

    std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order..." << std::endl;

//    ui *matching_order = NULL;
//    ui *pivots = NULL;
    ui **weight_array = NULL;

    size_t order_num = 0;

    sscanf(input_order_num.c_str(), "%zu", &order_num);

    std::vector<std::vector<ui>> spectrum;

    GenerateQueryPlan::generateCECIQueryPlan(query_graph[0], ceci_tree, ceci_order, matching_order, pivots);

    if (input_order_type != "Spectrum") {
        GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[0], matching_order, pivots);
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph[0], matching_order);
    } else {
        std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
    }

    std::cout << "-----" << std::endl;
    std::cout << "Enumerate..." << std::endl;
    size_t output_limit = 0;
    size_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<size_t>::max();
    } else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }

#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

    size_t call_count = 0;
    size_t time_limit = 0;
    sscanf(input_time_limit.c_str(), "%zu", &time_limit);

    embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph[0], ceci_tree, candidates,
                                                      candidates_count, TE_Candidates,
                                                      NTE_Candidates, ceci_order, output_limit, call_count,
                                                      all_match_res,idx,idx_count,embedding,valid_candidates,visited_vertices);

    memory_cost_in_bytes += sizeof(int) * all_match_res.size() * all_match_res[0].size();

#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
        outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
        delete[] EvaluateQuery::distribution_count_;
#endif

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;

    // sort(all_match_res.begin(), all_match_res.end(), [](vector<int>& a, vector<int>& b) {
    //     for (int i = 0; i < a.size(); ++i) {
    //         if (a[i] != b[i])
    //             return a[i] < b[i];
    //     }
    //     return true;
    // });

    cout << "finish" << endl;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    double total_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    sum_memory_cost_in_bytes += memory_cost_in_bytes;
    printf("First match Memory cost (MB): %.4lf\n", BYTESTOMB(sum_memory_cost_in_bytes));
    sum_embedding_count += embedding_count;
    sum_call_count += call_count;
    printf("#Embeddings: %zu\n", embedding_count);

    std::cout << "query_graph match result:" << std::endl;
    cout << all_match_res.size() << endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
//    for (auto i : all_match_res) {
//        map<int, int> nbr_lab_cnt;
//        for (auto j : i) {
//            cout << j << "(" << data_graph->getVertexLabel(j) << ") ";
//            ui nbrs_cnt;
//            auto nbrs = data_graph->getVertexDNeighbors(j, nbrs_cnt); // v_f's neighbors
//            for (int k = 0; k < nbrs_cnt; k++) {
//                nbr_lab_cnt[data_graph->getVertexLabel(nbrs[k])]++;
//            }
//        }
//        for (auto k : nbr_lab_cnt) {
//            cout << "[" << k.first << " " << k.second << "]" << endl;
//        }
//        cout << endl;
//    }

    if (Meth == 1)
        bfs_solve(total_time_in_ns, all_match_res, query_graph, data_graph);
    if (Meth == 2)
        subgraph_solve(total_time_in_ns, all_match_res, query_graph, data_graph);
    if (Meth == 3)
        dfs_solve(total_time_in_ns, all_match_res, query_graph, data_graph);

//    delete[] bst;
    delete[] is_show;
    delete[] query_is_matched;
    // delete[] valid_candidates;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;

    for (int i = 0; i < QueryNums; ++i)
        delete query_graph[i];
    delete[] query_graph;
    delete data_graph;
    return 0;
}
// auto tt1 = std::chrono::high_resolution_clock::now();
// auto tt2 = std::chrono::high_resolution_clock::now();
// auto tt = std::chrono::duration_cast<std::chrono::nanoseconds>(tt2 - tt1).count();
// cout << p << endl;