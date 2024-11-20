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
int SPECIAL = 1;
int Meth = 1;
int LL, RR;
double tim1, tim2, tim3, tim4;
bool test = false;

fstream fin = fstream("/home/dbia/qsl/FINISHED/PSM/nm/result.csv", ios::app);
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
vector<long long> memory;
string getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
    return tmp;
}

double sum_load_graphs_time_in_ns = 0;
size_t sum_memory_cost_in_bytes = 0;
size_t init_memory_cost_in_bytes = 0;
size_t sum_embedding_count = 0;
size_t sum_call_count = 0;

map<int,bool> all_lab;
auto sstart = std::chrono::high_resolution_clock::now();
size_t sp_match_res;
int thread_num;
int mx = -1;
vector<vector<int>> all_match_res;
mutex mtx;
mutex mtx_in_match;
mutex mtx2;
vector<bool> thread_state;
int all_match_idx;
vector<vector<int>> thread2idx;

int top = 3e6;
vector<pair<int, int>> new_idx;
size_t la;
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
size_t max_size{}, max_cal_size{};
class Worker {
public:
    size_t try_count{};
    size_t dfs_all_count{};
    bool* query_is_matched{};
    TreeNode *ceci_tree{};
    ui *ceci_order{};
    ui *candidates_count{};
    ui *idx{} ;
    ui *idx_count{} ;
    ui *embedding{} ;
    vector<vector<ui>> valid_candidates ;
    bool *visited_vertices{} ;
    size_t int_size{};
    size_t all_memory{};
    size_t all_query_vertex_nums{};
    ui **candidates = nullptr;
    ui* matching_order{};
    ui* pivots{};
    int cur_sp_match_res;

    vector<vector<int>> match_result;
    vector<pair<int, vector<pair<int, vector<int>>>>> get_nei_label_vertex;

    struct Node {
        vector<int> match;
        int step;
        bool flag;
        pair<int, int> fa;
        Node(){
            step = QUERY_NUMS;
            flag = false;
            fa = {-1, -1};
        }
    };
    struct SN {
        vector<Node> n;
        map<int, vector<int>> sum2idx;
        map<int, vector<int>> has_extend;
    };
    vector<SN> N;
    int space_cnt{};
    vector<bool> spanning_st;
    map<int, int> lab;
    map<int, map<int, bool>> have_caled;
    int* pst{};
    int* pst2{};
    vector<int> path;
    vector<pair<int, vector<int>>> res_graph;
    vector<pair<int, int>> now;
    map<int, vector<map<int, vector<int>>>> tree;
    size_t nei_all_size = 0;
    vector<bool> qim;
    vector<int> match_num;
    map<int, vector<map<int, map<int, bool>>>> Nst;
    map<int, map<int, bool>> n_st;
    map<int, vector<vector<pair<int, int>>>> unq;

    bool topk = false;
    map<int, map<int, int>> mp_cnt;
    vector<vector<pair<int, int>>> Path;
    int Worker_idx;
    double this_time{};
    Worker(Graph* data_graph, Graph** query_graph, int Worker_idx_) {
        cur_sp_match_res = 0;
        Worker_idx = Worker_idx_;
        idx = new ui[mx];
        idx_count = new ui[mx];
        embedding = new ui[mx];
        valid_candidates.resize(mx);
        visited_vertices = new bool[data_graph->vertices_count_];

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
        query_is_matched = new bool[QUERY_NUMS];
        pst = new int[data_graph->vertices_count_];
        pst2 = new int[data_graph->vertices_count_];
    }
    size_t LB(vector<pair<int, vector<int>>>& V, uint val) {
        size_t L = 0, R = V.size();
        while (L < R) {
            size_t MID = (L + R) >> 1;
            if (V[MID].first >= val) R = MID;
            else L = MID + 1;
        }
        return L;
    }
    size_t LB2(vector<pair<int, vector<pair<int, vector<int>>>>>& V, uint val) {
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
            if (idx_ >= nbrs + nbrs_cnt || *idx_ != u) {
                swap(u, v);
                auto lb = LB(res_graph, v);
                if (lb >= res_graph.size() || res_graph[lb].first != v) res_graph.insert(res_graph.begin() + lb, {v, vector<int>{}});
                auto lb2 = lower_bound(res_graph[lb].second.begin(), res_graph[lb].second.end(), u);
                if (lb2 == res_graph[lb].second.end() || *lb2 != u) res_graph[lb].second.insert(lb2, u);
            }
            else {
                auto lb = LB(res_graph, v);
                if (lb >= res_graph.size() || res_graph[lb].first != v) res_graph.insert(res_graph.begin() + lb, {v, vector<int>{}});
                auto lb2 = lower_bound(res_graph[lb].second.begin(), res_graph[lb].second.end(), u);
                if (lb2 == res_graph[lb].second.end() || *lb2 != u) res_graph[lb].second.insert(lb2, u);
            }
        }
        // for (int i = path.size() - 1; i > 0; --i) {
        //     auto v = path[i];
        //     auto u = path[i - 1];
        //     ui nbrs_cnt;
        //     auto nbrs = data_graph->getVertexDNeighbors(v, nbrs_cnt); // v_f's neighbors
        //     auto idx_ = std::lower_bound(nbrs, nbrs + nbrs_cnt, u);
        //     if (nbrs_cnt == 0 || *idx_ != u) continue;
        //     auto lb = LB(res_graph, v);
        //     if (lb >= res_graph.size() || res_graph[lb].first != v) res_graph.insert(res_graph.begin() + lb, {v, vector<int>{}});
        //     auto lb2 = lower_bound(res_graph[lb].second.begin(), res_graph[lb].second.end(), u);
        //     if (lb2 == res_graph[lb].second.end() || *lb2 != u) res_graph[lb].second.insert(lb2, u);
        // }
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
                        if (pst[u] >= 0) pst2[u] = 1;
                        this_edge = true;
                    }
                }
                continue;
            }
            if (pst[v] > 0) {
                if (pst[v] == dd + 1) {
                    // cout << u << " " << v << endl;
                    path.emplace_back(u, (int)v);
                    if (u == 1917)
                        if (pst2[u] == 0) pst2[u] = pst2[v] + 1;
                        else pst2[u] = min(pst2[u], pst2[v] + 1);
                    this_edge = true;
                    continue;
                }
                if (pst[v] > dd + 1) {
                    pst[v] = dd + 1;
                    if (get_path3D(data_graph, b, v, u, dd + 1, path)) {
                        // cout << u <<  " " << v << endl;
                        path.emplace_back(u, (int)v);
                        if (pst2[u] == 0) pst2[u] = pst2[v] + 1;
                        else pst2[u] = min(pst2[u], pst2[v] + 1);
                        this_edge = true;
                    }
                    else
                        pst[v] = 0;
                }
                if (pst[v] < dd + 1) {
                    if (pst[u] + pst2[v] < D) {
                        path.emplace_back(u, (int)v);
                        this_edge = true;
                    }
                }
            }
            if (pst[v] == 0) {
                pst[v] = dd + 1;
                if (get_path3D(data_graph, b, v, u, dd + 1, path)) {
                    // cout << u << " " << v << endl;
                    path.emplace_back(u, (int)v);
                    this_edge = true;
                    if (pst2[u] == 0) pst2[u] = pst2[v] + 1;
                    else pst2[u] = min(pst2[u], pst2[v] + 1);
                }
                else {
                    pst[v] = 0;
                }
            }
        }
        return this_edge;
    }
    bool add(int vertex, size_t lower_) {
        fstream fss = fstream("/home/dbia/qsl/FINISHED/PSM/" + PATH +"D" + to_string(D) + "/" + to_string(vertex) + ".txt");
        pair<int, vector<pair<int, vector<int>>>> nodes;
        nodes.first = vertex;
        vertex = -1;
        string line;
        while (fss >> vertex) {
            int idx_ = -1, lab_ = -1, num_ = -1, x;
            fss >> idx_ >> lab_ >> num_;

            if (all_lab.find(lab_) == all_lab.end()) {
                while (fss >> vertex)
                    if (vertex == -1)
                        break;
                continue;
            }
            nodes.second.resize(nodes.second.size() + 1);
            nodes.second.back().first = lab_;
            nodes.second.back().second.reserve(num_);
            for (int i = 0; i < num_; ++i) {
                fss >> x;
                nodes.second.back().second.emplace_back(x);
            }
            fss >> vertex;
        }
        fss.close();
        if (!nodes.second.empty()) {
            get_nei_label_vertex.insert(get_nei_label_vertex.begin() + lower_, nodes);
            return true;
        }
        return false;
    }
    void extd(Graph* data_graph) {
        if (topk) return;
        if (now.size() == QUERY_NUMS)   {
            int sum = 0;
            auto tmp = now;
            for (int i = 1; i < now.size(); ++i) {
                sum += now[i].first + now[i].second;
            }
            // mp_cnt[now[2].first][now[2].second]++;
            // cout << now[2].first << " " << now[2].second << " " << now[3].first << " " << now[3].second << endl;
            sort(tmp.begin(), tmp.end());
            if (unq.find(sum) != unq.end()) {
                bool same = true;
                for (auto ve : unq[sum]) {
                    same = true;
                    for (int i = 0; i < ve.size(); ++i) {
                        if (ve[i] != tmp[i]) {
                            same = false;
                            break;
                        }
                    }
                    if (same)
                        return;
                }
            }
            unq[sum].emplace_back(tmp);
//        mtx.lock();
            cur_sp_match_res++;
//             cout << Worker_idx << " " << cur_sp_match_res << endl;
//        mtx.unlock();

            if (cur_sp_match_res >= top) {
                // cout << "topk break" << endl;
                // mtx.lock();
                // cout << Worker_idx << " " << cur_sp_match_res << endl;
                // mtx.unlock();
                topk = true;
            }

//            fstream fs = fstream("../../result2.txt", ios::app);
//            for (auto i : now) {
//                auto ma = N[i.first].n[i.second].match;
//                sort(ma.begin(), ma.end());
//                for (auto j : ma)
//                    fs << j << " ";
//                fs << " ";
//            }
//            fs << endl;
//            for (auto i : Path) {
//                for (auto j : i) {
//                    fs << "[" << j.first << " " << j.second << "] ";
//                }
//            }
//            fs << endl;
//            fs.close();

            res_graph.clear();
            // memset(pst, 0, 4 * data_graph->vertices_count_);
            // memset(pst2, 0, 4 * data_graph->vertices_count_);
            map<int, bool> has_used;

            // for (auto i : now) {
            //     for (auto v : N[i.first].n[i.second].match) {
            //         if (has_used[v]) return;
            //         has_used[v] = true;
            //     }
            // }
            // for (auto i : now) {
//            auto i = now.back();
//            for (auto v : N[i.first].n[i.second].match) {
//                pst[v] = -i.first - 1;
//            }
            // }

//            vector<pair<int, int>> pa;
//            for (auto v : N[i.first].n[i.second].match) get_path3D(data_graph, i.first + 1, v, -1, 0, pa);
//            Path.emplace_back(pa);
//            //         assert(Path.size() + 1 == now.size());
//            for (auto v : N[i.first].n[i.second].match) pst[v] = 0;
//            for (const auto& v : Path.back()) {pst[v.first] = 0, pst[v.second] = 0, pst2[v.first] = 0, pst2[v.second] = 0;}
//            Path.pop_back();
            // }
            auto eeend = std::chrono::high_resolution_clock::now();
            auto timmm = std::chrono::duration_cast<std::chrono::nanoseconds>(eeend - sstart).count();
            if (cur_sp_match_res % 1000 == 0) {
                cout << cur_sp_match_res << " " << NC(timmm) << " " << endl;
            }

            return;
        }
        for (auto i : now) {
            if (i.first < tree.size())
                for (const auto& j : tree[i.first][i.second]) {
                    if (qim[j.first]) continue;
                    for (auto k : j.second) {
                        if (!Nst[i.first][i.second][j.first][k]) {
                            // dbg_cnt++;
                            // cout << dbg_cnt << endl;
                            // if (dbg_cnt == 26395)
                            //     cout << endl;
                            // assert(Nst[1][39][2].size() == 8272);
                            qim[j.first] = true;
                            now.emplace_back(j.first, k);
                            if (QUERY_NUMS >= 4 && now.size() == QUERY_NUMS - 1)
                                Nst[i.first][i.second][j.first][k] = true;

                            auto ii = now.back();
                            for (auto v : N[ii.first].n[ii.second].match) if (pst[v] >= 0) pst[v] = -ii.first - 1;
                            vector<pair<int, int>> pa;
                            for (auto v : N[ii.first].n[ii.second].match)
                                get_path3D(data_graph, ii.first + 1, v, -1, 0, pa);

                            Path.emplace_back(pa);

                            // cout << "1 " << now.size() << " " << Path.size() << endl;
                            // assert(Path.size() + 1 == now.size());
                            extd(data_graph);
                            Nst[i.first][i.second][j.first][k] = false;
                            if (topk) return;
                            qim[j.first] = false;
                            // if (j_ == 2)

                            for (auto v : N[ii.first].n[ii.second].match) if (pst[v] == -ii.first - 1) pst[v] = 0;
                            for (const auto& v : Path.back()) {
                                if (pst[v.first] > 0) pst[v.first] = 0;
                                if (pst[v.second] > 0) pst[v.second] = 0;
                                if (pst2[v.first] > 0) pst2[v.first] = 0;
                                if (pst2[v.second] > 0) pst2[v.second] = 0;
                            }
                            now.pop_back();
                            Path.pop_back();
                            // cout << "2 " << now.size() << " " << Path.size() << endl;
                        }
                    }
                }
        }
    }
    bool dfs(Graph* data_graph, pair<int, int> u, int cnt) {
        if (topk) return true;
        if (cnt >= QUERY_NUMS) return true;
        for (auto i : tree[u.first][u.second]) {
            if (qim[i.first]) continue;
            for (auto j : i.second) {
                // dbg_cnt++;
                // if (dbg_cnt == 169023)
                //     cout << endl;
                // cout << dbg_cnt << endl;
                // assert(Nst[1][39][2].size() == 8272);
                now.emplace_back(i.first, j);

                auto ii = now.back();
                for (auto v : N[ii.first].n[ii.second].match) if (pst[v] >= 0) pst[v] = -ii.first - 1;
                vector<pair<int, int>> pa;
                for (auto v : N[ii.first].n[ii.second].match)
                    get_path3D(data_graph, ii.first + 1, v, -1, 0, pa);
                Path.emplace_back(pa);

                // cout << "3 " << now.size() << " " << Path.size() << endl;
                qim[i.first] = true;
                if (i.first < tree.size() && !dfs(data_graph, {i.first, j}, cnt+1))
                    return false;
                if (topk) return true;
                // assert(Nst[1][39][2].size() > 0);
                Nst[u.first][u.second][i.first][j] = true;

                n_st.clear();
                extd(data_graph);
                if (topk) return true;
                qim[i.first] = false;
                now.pop_back();

                for (auto v : N[ii.first].n[ii.second].match) if (pst[v] == -ii.first - 1) pst[v] = 0;
                for (const auto& v : Path.back()) {
                    if (pst[v.first] > 0) pst[v.first] = 0;
                    if (pst[v.second] > 0) pst[v.second] = 0;
                    if (pst2[v.first] > 0) pst2[v.first] = 0;
                    if (pst2[v.second] > 0) pst2[v.second] = 0;
                }
                Path.pop_back();

                // cout << "4 " << now.size() << " " << Path.size() << endl;
                match_num[i.first]--;
                for (auto each_match_num : match_num) if (each_match_num == 0) return false;
            }
        }
        return true;
    }
    void CECI(Graph* data_graph, Graph** query_graph, vector<vector<int>>& cur_match_res, int i, int u, int vertex) {
        if (data_graph->getVertexDegree(vertex) < query_graph[i]->getVertexDegree(u))
            return;
        const std::unordered_map<LabelID, ui> *u_nlf = query_graph[i]->getVertexNLF(u);
        const std::unordered_map<LabelID, ui> *v_nlf = data_graph->getVertexNLF(vertex);

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
    map<int, map<int, size_t>> epoch_cnt;
    void begin_extend(Graph* data_graph, Graph** query_graph, int match_idx) {
        auto each_begin1 = std::chrono::high_resolution_clock::now();
//        cout<< 'bb' << endl;
        N.clear();
        unq.clear();
        // eg2id.clear();
        //    g.clear();
        have_caled.clear();
        lab.clear();
        spanning_st.clear();
        now.clear();
        match_num.clear();
        qim.clear();
        Nst.clear();
        tree.clear();
        epoch_cnt.clear();
        match_num.resize(QUERY_NUMS);
        match_num[0] = 1;
        N.resize(QUERY_NUMS);
        space_cnt = 0;

        size_t sum = 0;
        // eg2id.resize(QUERY_NUMS);
        for (auto i : all_match_res[match_idx]) sum += i;
        N[0].n.resize(1);
        N[0].n[0].match = all_match_res[match_idx];
        N[0].n[0].step = 0;
        N[0].sum2idx[sum].emplace_back(0);
        bool enough_break = false;
        size_t bound = 1;
        for (int epoch = 0; epoch < QUERY_NUMS - 1; ++epoch) {
//            cout << "N size" << endl;
//            for (int j = 0; j < QUERY_NUMS; ++j) {
//                cout << "j: " << N[j].n.size() << endl;
//            }
            for (int j = 0; j < QUERY_NUMS; ++j) {
                // cout << epoch << " " << j << " N size" << endl;
                // if (epoch == 1) {
                //     cout << "thread: " << Worker_idx << endl;
                //     for (int jj = 0; jj < QUERY_NUMS; ++jj) {
                //         cout << "i: " << N[jj].n.size() << endl;
                //     }
                // }
                tree[j].resize(N[j].n.size());
                for (int n_j = 0; n_j < N[j].n.size(); ++n_j) {
                    if (top != 0x3f3f3f3f && epoch == 0) {
                        size_t ubt = 1;
                        for (int j_ = 0; j_ < QUERY_NUMS; ++j_) {
                            ubt *= N[j_].n.size();
                        }
                        if (ubt >= top) {
                            enough_break = true;
                            // cout << "break bound: " << ubt << endl;
                            break;
                        }
                    }
                    auto& n = N[j].n[n_j];
                    if (n.step == epoch) {
                        if (!n.flag && n.step < QUERY_NUMS - 1) {
                            for (auto match_vertex : n.match) {
                                auto lower_ = LB2(get_nei_label_vertex, match_vertex);
                                if (lower_ >= get_nei_label_vertex.size() || get_nei_label_vertex[lower_].first != match_vertex) {
                                    if (!add(match_vertex, lower_))
                                        continue;
                                }
                                auto lab_ver = get_nei_label_vertex[lower_].second;
                                for (int i = 1; i < QUERY_NUMS; ++i) {
                                    if (i == j) continue;
                                    for (auto LB_ : query_graph[i]->all_labs) {
                                        auto lb = LB_.first;

                                        ui u_cnt;
                                        auto lb2v = query_graph[i]->getVerticesByLabel(lb, u_cnt);
                                        auto lower = LB(lab_ver, lb);
                                        if (lower < lab_ver.size() && lab_ver[lower].first == lb) {
                                            for (auto vertex : lab_ver[lower].second) {
                                                if (bound >= top) {
                                                    enough_break = true;
                                                    break;
                                                }
                                                if (epoch_cnt[1][i] * epoch_cnt[1][j] > top) continue;

                                                size_t sm = 1;
                                                for (int ii_ = 1; ii_ < QUERY_NUMS; ++ii_) {
                                                    if (epoch_cnt[1][ii_]) {
                                                        sm *= epoch_cnt[1][ii_];
                                                    }
                                                }
                                                if (sm >= top && epoch_cnt[1][i]) {
                                                    break;
                                                }
                                                // cout << sm << " " << epoch_cnt[1][i] << endl;
                                                // if (!epoch) for (int p_p = 1; p_p < QUERY_NUMS; p_p++) cout << N[p_p].n.size() << " \n"[p_p + 1 == QUERY_NUMS];
                                                if (top != 0x3f3f3f3f && epoch == 0) {
                                                    size_t ubt = 1;
                                                    for (int j_ = 0; j_ < QUERY_NUMS; ++j_) {
                                                        ubt *= N[j_].n.size();
                                                    }
                                                    if (ubt >= top) {
                                                        enough_break = true;
                                                        // cout << "break bound: " << ubt << endl;
                                                        break;
                                                    }
                                                }
                                                if (N[i].has_extend.find(vertex) != N[i].has_extend.end()) {
                                                    if (!N[i].has_extend[vertex].empty()) {
                                                        // add_edge
                                                        for (auto edges : N[i].has_extend[vertex]) {
                                                            int n_i = edges;
                                                            if (n.step < N[i].n[n_i].step) {
//                                                                if (i > tree[j][n_j].size()) tree[j][n_j].resize(QUERY_NUMS);
                                                                auto lb_ = lower_bound(tree[j][n_j][i].begin(), tree[j][n_j][i].end(), n_i);
                                                                if (lb_ == tree[j][n_j][i].end() || *lb_ != n_i) {
//                                                                    tree[j][n_j][i].first = i;
                                                                    tree[j][n_j][i].insert(lb_, n_i);
                                                                    if (top != 0x3f3f3f3f) {
                                                                        if (epoch == 0) {
                                                                            size_t ubt = 1;
                                                                            for (int j = 0; j < QUERY_NUMS; ++j) {
                                                                                ubt *= N[j].n.size();
                                                                                // cout << j << " " << N[j].n.size() << endl;
                                                                            }
                                                                            if (ubt >= top) {
                                                                                enough_break = true;
                                                                                cout << "break bound: " << ubt << endl;
                                                                                break;
                                                                            }
                                                                        }
                                                                        else {
                                                                            map<int, bool> st_;
                                                                            bound++;
                                                                            st_[j] = true;
                                                                            st_[i] = true;
                                                                            for (int qi_ = 1; qi_ < QUERY_NUMS; ++qi_) {
                                                                                if (!st_[qi_]) {
                                                                                    bound += tree[j][n_j][qi_].size();
                                                                                }
                                                                                if (bound >= top) {
                                                                                    // cout << "break 0 epoch: " << epoch << endl;
                                                                                    enough_break = true;
                                                                                    break;
                                                                                }if (enough_break) break;
                                                                            }
                                                                        }
                                                                        if (enough_break) break;
                                                                    }if (enough_break) break;
                                                                    match_num[i]++;
                                                                }if (enough_break) break;
                                                            }if (enough_break) break;
                                                        }if (enough_break) break;
                                                    }if (enough_break) break;
                                                    continue;
                                                }if (enough_break) break;

                                                N[i].has_extend[vertex] = vector<int> {};
                                                N[i].has_extend[vertex].reserve(10);
                                                for (int u_ = 0; u_ < u_cnt; ++u_) {
                                                    auto u = lb2v[u_];
                                                    try_count++;
                                                    vector<vector<int>> cur_match_res;
                                                    dfs_all_count++;
                                                    CECI(data_graph, query_graph, cur_match_res, i, u, vertex);
                                                    //                                                 cout << i << " " << u << " " << vertex << endl;
                                                    if (cur_match_res.empty()) continue;
                                                    for (auto each_res: cur_match_res) {
                                                        sum = 0;
                                                        for (auto each_res_i: each_res) sum += each_res_i;
                                                        if (N[i].sum2idx.find(sum) == N[i].sum2idx.end()) {
                                                            //                                                        N[i].sum2idx[sum].reserve(3);
                                                            N[i].sum2idx[sum].emplace_back(N[i].n.size());
                                                            sort(each_res.begin(), each_res.end());
                                                            Node a;
                                                            a.match = each_res;
                                                            a.step = N[j].n[n_j].step + 1;
                                                            N[i].n.emplace_back(a);
                                                            epoch_cnt[a.step][i]++;
                                                            match_num[i]++;
                                                            map<int, vector<int>> tmp;
                                                            tmp[j].emplace_back(n_j);
//                                                            if (i > tree[j][n_j].size()) {
//                                                                tree[j][n_j].resize(QUERY_NUMS);
//                                                                for (int tree_resize_idx = 0; tree_resize_idx < QUERY_NUMS; ++tree_resize_idx)
//                                                                    tree[j][n_j][tree_resize_idx].first = tree_resize_idx;
//                                                            }
                                                            tree[j][n_j][i].emplace_back(N[i].n.size() - 1);
                                                            N[i].has_extend[vertex].emplace_back(N[i].n.size() - 1);

                                                            // if (n_j >= eg2id[j].size()) {
                                                            //     eg2id[j].resize(n_j + 1);
                                                            // }
                                                            // if (N[i].n.size() - 1 >= eg2id[i].size()) eg2id[i].resize(N[i].n.size() - 1 + 1);
                                                            // if (eg2id[j][n_j] == 0) {
                                                            //     lab[space_cnt] = j;
                                                            //     eg2id[j][n_j] = space_cnt++;
                                                            // }
                                                            // if (eg2id[i][N[i].n.size() - 1] == 0) {
                                                            //     lab[space_cnt] = i;
                                                            //     eg2id[i][N[i].n.size() - 1] = space_cnt++;
                                                            // }
                                                            // if (eg2id[j].find(n_j) == eg2id[j].end()) {
                                                            //     lab[space_cnt] = j;
                                                            //     eg2id[j][n_j] = space_cnt++;
                                                            // }
                                                            // if (eg2id[i].find(N[i].n.size() - 1) == eg2id[i].end()) {
                                                            //     lab[space_cnt] = i;
                                                            //     eg2id[i][N[i].n.size() - 1] = space_cnt++;
                                                            // }
                                                            map<int, bool> st_;
                                                            st_[i] = true;
                                                            st_[j] = true;
                                                            bound++;
                                                            if (top != 0x3f3f3f3f) {
                                                                if(epoch == 0) {
                                                                    size_t ubt = 1;
                                                                    for (int j = 0; j < QUERY_NUMS; ++j) {
                                                                        ubt *= N[j].n.size();
                                                                        // cout << j << " " << N[j].n.size() << endl;
                                                                    }
                                                                    if (ubt >= top) {
                                                                        enough_break = true;
                                                                        cout << "break bound: " << ubt << endl;
                                                                        break;
                                                                    }
                                                                }
                                                                else{
                                                                    for (int qi_ = 1; qi_ < QUERY_NUMS; ++qi_) {
                                                                        if (!st_[qi_]) {
                                                                            bound += epoch_cnt[1][qi_];
                                                                            bound += tree[j][n_j][qi_].size();
                                                                        }
                                                                        if (bound >= top) {
                                                                            // cout << "break 1 epoch: " << epoch << endl;
                                                                            enough_break = true;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                                if (enough_break) break;
                                                            }
                                                            if (enough_break) break;
                                                        }
                                                        else {
                                                            bool f = true;
                                                            int n_i = -1;
                                                            for (auto sum_i_idx : N[i].sum2idx[sum]) {
                                                                auto& history_match = N[i].n[sum_i_idx].match;
                                                                if (each_res.size() == history_match.size()) {
                                                                    f = true;
                                                                    sort(each_res.begin(), each_res.end());
                                                                    sort(history_match.begin(), history_match.end());
                                                                    for (int idx_ = 0; idx_ < each_res.size(); ++idx_) {
                                                                        if (each_res[idx_] != history_match[idx_]) {
                                                                            f = false;
                                                                            break;
                                                                        }
                                                                    }
                                                                    if (!f)
                                                                        continue;
                                                                    n_i = sum_i_idx;
                                                                    break;
                                                                }
                                                                else f = false;
                                                            }
                                                            if (!f) {
                                                                N[i].sum2idx[sum].emplace_back(N[i].n.size());

                                                                Node a;
                                                                sort(each_res.begin(), each_res.end());
                                                                a.match = each_res;
                                                                a.step = N[j].n[n_j].step + 1;
                                                                N[i].n.emplace_back(a);
                                                                match_num[i]++;
                                                                map<int, vector<int>> tmp;
                                                                tmp[j].emplace_back(n_j);
//                                                                if (i > tree[j][n_j].size()) {
//                                                                    tree[j][n_j].resize(QUERY_NUMS);
//                                                                    for (int tree_resize_idx = 0; tree_resize_idx < QUERY_NUMS; ++tree_resize_idx)
//                                                                        tree[j][n_j][tree_resize_idx].first = tree_resize_idx;
//                                                                }
                                                                tree[j][n_j][i].emplace_back(N[i].n.size() - 1);
                                                                epoch_cnt[a.step][i]++;
                                                                N[i].has_extend[vertex].emplace_back(N[i].n.size() - 1);
                                                                // if (n_j >= eg2id[j].size()) eg2id[j].resize(n_j + 1);
                                                                // if (N[i].n.size() - 1 >= eg2id[i].size()) eg2id[i].resize(N[i].n.size() - 1 + 1);
                                                                // if (eg2id[j][n_j] == 0) {
                                                                //     lab[space_cnt] = j;
                                                                //     eg2id[j][n_j] = space_cnt++;
                                                                // }
                                                                // if (eg2id[i][N[i].n.size() - 1] == 0) {
                                                                //     lab[space_cnt] = i;
                                                                //     eg2id[i][N[i].n.size() - 1] = space_cnt++;
                                                                // }
                                                                // if (eg2id[j].find(n_j) == eg2id[j].end()) {
                                                                //     lab[space_cnt] = j;
                                                                //     eg2id[j][n_j] = space_cnt++;
                                                                // }
                                                                // if (eg2id[i].find(N[i].n.size() - 1) == eg2id[i].end()) {
                                                                //     lab[space_cnt] = i;
                                                                //     eg2id[i][N[i].n.size() - 1] = space_cnt++;
                                                                // }
                                                                bound++;
                                                                map<int, bool> st_;
                                                                st_[i] = true;
                                                                st_[j] = true;
                                                                if (top != 0x3f3f3f3f) {
                                                                    if(epoch == 0) {
                                                                        size_t ubt = 1;
                                                                        for (int j = 0; j < QUERY_NUMS; ++j) {
                                                                            ubt *= N[j].n.size();
                                                                            // cout << j << " " << N[j].n.size() << endl;
                                                                        }
                                                                        if (ubt >= top) {
                                                                            enough_break = true;
                                                                            cout << "break bound: " << ubt << endl;
                                                                            break;
                                                                        }
                                                                    }
                                                                    else {
                                                                        for (int qi_ = 1; qi_ < QUERY_NUMS; ++qi_) {
                                                                            if (!st_[qi_]) {
                                                                                bound += epoch_cnt[1][qi_];
                                                                                bound += tree[j][n_j][qi_].size();
                                                                            }
                                                                            if (bound >= top) {
                                                                                // cout << "break 2 epoch: " << epoch << endl;
                                                                                enough_break = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                    if (enough_break) break;
                                                                }
                                                            }
                                                            else {
                                                                if (n.step < N[i].n[n_i].step) {
//                                                                    if (i > tree[j][n_j].size()) {
//                                                                        tree[j][n_j].resize(QUERY_NUMS);
//                                                                        for (int tree_resize_idx = 0; tree_resize_idx < QUERY_NUMS; ++tree_resize_idx)
//                                                                            tree[j][n_j][tree_resize_idx].first = tree_resize_idx;
//                                                                    }
                                                                    auto lb_ = lower_bound(tree[j][n_j][i].begin(), tree[j][n_j][i].end(), n_i);
                                                                    if (lb_ == tree[j][n_j][i].end() || *lb_ != n_i) {
                                                                        tree[j][n_j][i].insert(lb_, n_i);
                                                                        match_num[i]++;
                                                                        //                                                                    g.emplace_back(j, n_j, i, n_i);
                                                                        //                                                                     if (!N[i].has_extend[vertex].empty() && N[i].has_extend[vertex][0] == n_j) continue;
                                                                        N[i].has_extend[vertex].emplace_back(n_i);
                                                                        map<int, bool> st_;
                                                                        bound++;
                                                                        st_[j] = true;
                                                                        st_[i] = true;
                                                                        for (int qi_ = 1; qi_ < QUERY_NUMS; ++qi_) {
                                                                            if (!st_[qi_]){
                                                                                bound += tree[j][n_j][qi_].size();
                                                                            }
                                                                            if (bound >= top) {
                                                                                // cout << "break 3 epoch: " << epoch << endl;
                                                                                enough_break = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                        if (enough_break) break;
                                                                    }
                                                                    if (enough_break) break;
                                                                    //                                                                 if (N[i].has_extend[vertex].empty())
                                                                    //                                                                 N[i].has_extend[vertex].emplace_back(n_j);
                                                                }
                                                                if (enough_break) break;
                                                            }
                                                            if (enough_break) break;
                                                        }
                                                        if (enough_break) break;
                                                    }
                                                    if (enough_break) break;
                                                }
                                                if (enough_break) break;
                                                //                                             if (N[i].has_extend[vertex].empty())
                                                //                                                 N[i].has_extend.erase(vertex);
                                            }
                                            if (enough_break) break;
                                        }
                                        if (enough_break) break;
                                    }
                                    if (enough_break) break;
                                }
                                if (enough_break) break;
                            }
                            if (enough_break) break;
                        }
                        if (enough_break) break;
                        n.flag = true;
                    }
                    if (enough_break) break;
                }
                if (enough_break) break;
            }
            if (top != 0x3f3f3f3f && epoch == 0) {
                size_t ubt = 1;
                for (int j = 0; j < QUERY_NUMS; ++j) {
                    ubt *= N[j].n.size();
                    // cout << j << " " << N[j].n.size() << endl;
                }
                if (ubt >= top) {
                    enough_break = true;
                    cout << "break bound: " << ubt << endl;
                    break;
                }
            }
            if (enough_break) break;
            if (epoch == 0) {
                bound = 1;
                for (int j = 1; j < QUERY_NUMS; ++j) {
                    bound *= N[j].n.size();
                }
            }
        }

        // cout << "enough_break: " << enough_break << endl;
        auto each_end1 = std::chrono::high_resolution_clock::now();
        tim1 += std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count();
//    this_time += std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count();
        // cout << "build time: " << NC(std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count()) << endl;
        //    g.clear();
        auto ti = std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - sstart).count();
        //    cout << NC(ti) << "s" << endl;
        for (auto & i : N) {
            if (i.n.empty())
                return;
        }

//         cout << "N size" << endl;
//         for (int j = 0; j < QUERY_NUMS; ++j) {
//             cout << "j: " << N[j].n.size() << endl;
//         }
        if (thread_num == 1) {
            auto new_la = get_index_mem()["pk"];
            if (new_la > la) {
                size_t all_cnt = tree.size();
                for (auto i : tree) {
                    all_cnt += i.second.size();
                    for (auto j : i.second) {
                        all_cnt += j.size();
                        for (auto l : j) {
                            all_cnt += l.second.size();
                        }
                    }
                }
                all_cnt *= 4;
                for (auto & i : N) {
                    for (const auto& j : i.n) {
                        all_cnt += 4 * j.match.size();
                        all_cnt += 8;
                    }
                    for (const auto& j : i.has_extend) {
                        all_cnt += sizeof(int);
                        all_cnt += sizeof(int) * j.second.size();
                        all_cnt += 3 * sizeof(void*) + sizeof(char);
                    }
                    for (const auto& j : i.sum2idx) {
                        all_cnt += sizeof(int);
                        all_cnt += sizeof(int) * j.second.size();
                        all_cnt += 3 * sizeof(void*) + sizeof(char);
                    }
                }
                max_cal_size = all_cnt;
            }
        }

//    for (auto& i : tree) {
//        for (auto& j : i.second) {
//            sort(j.begin(), j.end(), [](pair<int, vector<int>>& a, pair<int, vector<int>>& b) {
//                return a.second.size() < b.second.size();
//            });
//        }
//    }

        //    sort(N.begin(), N.end(), [](SN& a, SN& b) {
        //        return a.n.size() < b.n.size();
        //    });

        // Nst.resize(tree.size());

        for (const auto& i : tree) {
            Nst[i.first].resize(i.second.size() + 1);
        }

        // cout << endl;
        // for (int j = 0; j < QUERY_NUMS; ++j) {
        //     cout << j << ": " << match_num[j] << endl;
        // }
        dfs_all_count = 0;

        for (auto & i : N) {
            i.sum2idx.clear();
            i.has_extend.clear();
        }

        auto each_begin2 = std::chrono::high_resolution_clock::now();

        qim.resize(QUERY_NUMS);
        //    qim[0] = true;

        now.emplace_back(0, 0);
        for (auto v : N[0].n[0].match) pst[v] = -1;
        // cout << "------------ build over ------------" << endl;
        // cout << "edge_cnt: " << space_cnt << endl;
        dfs(data_graph, {0, 0}, 1);
        // cout << "------------ finish extd ------------" << endl;
        auto each_end2 = std::chrono::high_resolution_clock::now();
        tim2 += std::chrono::duration_cast<std::chrono::nanoseconds>(each_end2 - each_begin2).count();
//        mtx2.lock();
        // top -= cur_sp_match_res;
//        mtx2.unlock();
        Path.clear();
        if (topk) return;
    }
    void Excute(double total_time_in_ns, Graph** query_graph, Graph* data_graph) {
//        for (; all_match_idx < all_match_res.size(); ) {
        for (int i = Worker_idx; i < all_match_res.size(); i += thread_num) {
            nei_all_size = 0;
            get_nei_label_vertex.clear();

//            mtx_in_match.lock();
////            match_result.clear();
////            match_result.emplace_back(all_match_res[all_match_idx]);
//            all_match_idx++;
//            if (all_match_idx >= all_match_res.size()) break;
//            mtx_in_match.unlock();

//            thread2idx[Worker_idx].emplace_back(all_match_idx);
//             cout << "Worker_idx: " << Worker_idx << " " << i << endl;

            topk = false;
//            cout << all_match_idx << endl;
            begin_extend(data_graph, query_graph, i);
            if (topk) break;
        }
//        delete []pst;
//        delete []pst2;
//        delete []query_is_matched;
    }
};
vector<double> times;

void exc(Graph* data_graph, Graph** query_graph, int Worker_idx) {
    Worker worker(data_graph, query_graph, Worker_idx);
    worker.Excute(0, query_graph, data_graph);
//    times.emplace_back(worker.this_time);
    if (thread_num == 1) {
        auto new_la = get_index_mem()["pk"];
        size_t nei_all_size;
        if (new_la > la) {
            nei_all_size = worker.get_nei_label_vertex.size();
            for (const auto &iii: worker.get_nei_label_vertex) {
                nei_all_size += iii.second.size();
                for (const auto &jii: iii.second) {
                    nei_all_size += jii.second.size();
                }
            }
            max_size = nei_all_size;
            la = new_la;
        }
    }
}
int add(int vertex, size_t lower_) {
    fstream fss = fstream("/home/dbia/qsl/FINISHED/PSM/" + PATH +"D" + to_string(D) + "/" + to_string(vertex) + ".txt");
//    pair<int, vector<pair<int, vector<int>>>> nodes;
//    nodes.first = vertex;
    vertex = -1;
    string line;
    int cnt_ = 0;
    while (fss >> vertex) {
        int idx_ = -1, lab_ = -1, num_ = -1, x;
        fss >> idx_ >> lab_ >> num_;

        if (all_lab.find(lab_) == all_lab.end()) {
            while (fss >> vertex)
                if (vertex == -1)
                    break;
            continue;
        }
//        nodes.second.resize(nodes.second.size() + 1);
//        nodes.second.back().first = lab_;
//        nodes.second.back().second.reserve(num_);
        for (int i = 0; i < num_; ++i) {
            fss >> x;
            cnt_++;
//            nodes.second.back().second.emplace_back(x);
        }
        fss >> vertex;
    }
    fss.close();

    return cnt_;
}
void bfs_init_solve(double total_time_in_ns, vector<vector<int>>& all_match_res, Graph** query_graph, Graph* data_graph) {
    auto sstart = std::chrono::high_resolution_clock::now();
    std::cout << "Begin SP Match!" << std::endl;
    top /= thread_num;
    cout << "PATH: " << PATH << " Query_Num: " << QUERY_NUMS << " D: " << D << " Meth: bfs_init " << thread_num << " top " << top << dl;

    thread_num = min(thread_num, (int)all_match_res.size());
    // thread_res.resize((int)all_match_res.size());
//    new_idx.resize(all_match_res.size());
//    for (int i = 0; i < all_match_res.size(); ++i) {
//        int sum = 0;
//        for (auto v : all_match_res[i]) sum += add(v, 0);
////        cout << sum << endl;
//        new_idx[i] = {sum, i};
//    }
//    sort(new_idx.begin(), new_idx.end(), [](pair<int, int>& a, pair<int, int>& b) {
//        return a.first > b.first;
//    });
    std::vector<std::thread> threads;
    thread_state.reserve(thread_num);
    thread2idx.resize(thread_num);
    all_match_idx = 0;
    for (int i = 0; i < thread_num; ++i) {
        threads.emplace_back(exc, data_graph, query_graph, i);
    }
    for (std::thread& t : threads) {
        t.join();
    }
//    for (int i = 0; i < thread_num;++i) {
//        cout << i << " " << thread2idx[i].size() << endl;
//    }
    auto eend = std::chrono::high_resolution_clock::now();
    double ALL_TIME = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart).count();
    sum_memory_cost_in_bytes += sizeof (bool)* data_graph->vertices_count_;
    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    // printf("Memory cost (MB): %.4lf, %.4lf\n", BYTESTOMB(all_memory + sum_memory_cost_in_bytes + init_memory_cost_in_bytes), BYTESTOMB(init_memory_cost_in_bytes));
    printf("#Embeddings: %d\n", sp_match_res);
    // printf("DFS call count: %d\n", dfs_all_count);
    // printf("try count: %d\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(ALL_TIME + total_time_in_ns));
    cout << NC(ALL_TIME) << " " << NC(total_time_in_ns) << " " << NC(tim1) << " " << NC(tim2) << endl;
    std::cout << "End." << std::endl;
    // cout << "max_size " << (get_index_mem()["pk"] - memory[0]) << "KB" << endl;
    // cout << "max_cal_size " << (double)max_cal_size / 1024 << "kB" << endl;
    // cout << "nei_all_size " << (double)nei_all_size * 4 / 1024 << "kB" << endl;
    fin << "new" << "," << PATH << "," << QUE << ","<< QUERY_NUMS << "," << D << ","
        << sp_match_res << ","
        << NC(ALL_TIME) << ","
        << NC(total_time_in_ns) << ","
        << NC(ALL_TIME + total_time_in_ns) << ","
        << 0 << ","
        << NC(tim1) << ","
        << NC(tim2) << ","
        << (double)max_size * 4 / 1024 << ","
        << (double)max_cal_size / 1024 << endl;
    fin.close();
}

// void Init(Graph* data_graph, string pat, int l, int r) {
//     vector<bool> cst5(data_graph->vertices_count_, false);
//     vector<int> gdist(data_graph->vertices_count_, 0);
//     assert(l < data_graph->vertices_count_);
//     assert(r <= data_graph->vertices_count_);
//     for (int u = l; u < r; ++u) {
//         fill(cst5.begin(), cst5.end(), 0);
//         fill(gdist.begin(), gdist.end(), 0);
//
//         get_nei_label_vertex[u][data_graph->getVertexLabel(u)].emplace_back(u);
//
//         queue<int> q;
//         q.push(u);
//         while (!q.empty()) {
//             auto ver = q.front(); q.pop();
//             if (cst5[ver] || gdist[ver] >= D) {
//                 continue;
//             }
//             cst5[ver] = true;
//             ui nbrs_cnt;
//             const VertexID *nbrs = data_graph->getVertexNeighbors(ver, nbrs_cnt); // 获取点vertex的邻居
//             for (int j = 0; j < nbrs_cnt; ++j) { // 遍历点vertex的邻居
//                 int v = nbrs[j];
//                 int lb = data_graph->getVertexLabel(v);
//                 if (!cst5[v]) {
//                     if (gdist[ver] + 1 <= D)
//                     q.push(v);
//                     if (gdist[v] == 0) {
//                         get_nei_label_vertex[u][lb].emplace_back(v);
// //                        cout << v << endl;
//                         gdist[v] = gdist[ver] + 1;
//                     }
//                 }
//             }
//         }
//         fstream fin2 = fstream(pat + to_string(u)+".txt", ios::out);
//         int j_ = 0;
//         for (const auto& j : get_nei_label_vertex[u]) {
//             fin2 << u << " " << j_++ << " ";
//             fin2 << j.first << " " << j.second.size() << " ";
//             for (auto k: j.second) {
//                 fin2 << k << " ";
//             }
//             fin2 << -1 << endl;
//         }
//         fin2.close();
//         get_nei_label_vertex[u].clear();
//     }
//     char path_[50] = "ls /home/qsl/CL";
//     strcat(path_, to_string(D).c_str());
//     char path_end[10] = "/ | wc -w";
//     strcat(path_, path_end);
//     system(path_);
// //    cout << l << " " << r << " finish" << endl;
// }
// void WR(Graph* data_graph, string pat, bool f) {
//     get_nei_label_vertex.resize(data_graph->vertices_count_);
//     for (int i = 0; i < data_graph->vertices_count_; i += data_graph->vertices_count_ / 10) {
//         vector<thread> thrd;
//         int thrd_num = 32;
//         int l = 0, r = data_graph->vertices_count_ / 10 / 32;
//         for (int j = 0; j < thrd_num; ++j) {
//             thrd.emplace_back(Init, data_graph, pat, i + l, i + l + r);
//             l += r;
//         }
//         for (std::thread& t : thrd) {
//             t.join();
//         }
//         cout << i << " " << i + l << " finish " << getTime() << endl;
//     }
// //    Init(data_graph, pat, 0, data_graph->vertices_count_);
//
//     // if (f) {
//     //     for (int i = 0; i < data_graph->vertices_count_; ++i) {
//     //         fstream fin2 = fstream(pat + to_string(i)+".txt", ios::out);
//     //         for (int j_ = 0; j_ < get_nei_label_vertex[i].second.size(); ++j_) {
//     //             auto j = get_nei_label_vertex[i].second[j_];
//     //             fin2 << i << " " << j_ << " ";
//     //             fin2 << j.first << " " << j.second.size() << " ";
//     //             for (auto k: j.second) {
//     //                 fin2 << k << " ";
//     //             }
//     //             fin2 << -1 << endl;
//     //         }
//     //         fin2.close();
//     //     }
//     // }
// }
// void RD(Graph* data_graph, string pat) {
//     cout << "RD" << endl;
//     fstream file = fstream(pat, ios::in);
//     std::string line;
//     int la = -1;
//     char type;
//     bool ps = false;
//     get_nei_label_vertex.resize(1);
//     while (std::getline(file, line)) { // 逐行读取文件内容
// //        std::cout << line << std::endl; // 输出读取的行
//         int id = -1, lb = -1, ix = -1, ct = -1, nm = -1;
//         int x = 0;
//         for (int i = 0; i < line.size(); ++i) {
//             if (line[i] == ' ' || line[i] == '\n') {
//                 if (id == -1) {
//                     if (la != x && la != -1) {
//                         if (get_nei_label_vertex.back().second.size() != 0)
//                             get_nei_label_vertex.resize(get_nei_label_vertex.size() + 1);
//                     }
//                     id = x;
//                     la = x;
//                     x = 0;
//                 }
//                 else if (ix == -1) {
//                     ix = x;
//                     x = 0;
//                 }
//                 else if (nm == -1 && ix == 0) {
//                     nm = x;
//                     get_nei_label_vertex.back().first = id;
//                     x = 0;
//                 }
//                 else if (lb == -1) {
//                     lb = x;
//                     if (all_lab.find(lb) == all_lab.end()) break;
//                     x = 0;
//                     get_nei_label_vertex.back().second.resize(get_nei_label_vertex.back().second.size() + 1);
//                     get_nei_label_vertex.back().second.back().first = lb;
//                 }
//                 else if (ct == -1) {
//                     ct = x;
//                     get_nei_label_vertex.back().second.back().second.reserve(ct);
//                     x = 0;
//                 }
//                 else {
//                     get_nei_label_vertex.back().second.back().second.emplace_back(x);
//                     x = 0;
//                 }
//             }
//             else {
//                 x = x * 10 + line[i] - '0';
//             }
//         }
//     }
//     file.close();
// }

int main(int argc, char** argv) {
    cout << getTime() << endl;
    cout << "begin" << endl;
    // auto mp = get_index_mem();
    // memory.emplace_back(mp["pk"]);
    test = false;
    // la = get_index_mem()["pk"];
//    QUE = "1";
//    PATH = "MSRC-21C";
//    D = 3;
//    QUERY_NUMS = 3;
//    Meth = 2;
//    input_query_graph_file = "../../../test/querys/" + QUE;
//    input_data_graph_file = "../../../test/graphs/Dataset/" + PATH + "/" + PATH + ".txt";
//    thread_num = 1;
    MatchingCommand command(argc, argv);
    QUE = command.getQueryGraphFilePath();
    PATH = command.getDataGraphFilePath();
    D = stoi(command.getDist());
    QUERY_NUMS = stoi(command.getQueryNums());
    thread_num = stoi(command.getMeth());
    Meth = 2;
    input_query_graph_file = "../../test/querys/" + QUE;
    input_data_graph_file = "../../test/graphs/Dataset/" + PATH + "/" + PATH + ".txt";
    bool has_dir_ = false;
    if (PATH == "CL") has_dir_ = true;
//    srand(time(0));
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
    std::cout << "\tL: " << LL << " R: " << RR << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;



    /**
     * Load input graphs.
     */

    std::cout << "Load graphs..." << std::endl;
//    auto start = std::chrono::high_resolution_clock::now();

//    Graph* query_graph = new Graph(true);
//    query_graph->loadGraphFromFile(input_query_graph_file);
//    query_graph->buildCoreTable();

    Graph** query_graph = new Graph*[QUERY_NUMS];
    for (int i = 0; i < QUERY_NUMS; ++i) {
        query_graph[i] = new Graph(true, has_dir_);
        query_graph[i]->id_ = i;
        query_graph[i]->loadGraphFromFile(input_query_graph_file + "/q" + to_string(i) + ".txt");
//        query_graph[i]->loadGraphFromFile("../../test/querys/4/q0.txt");
//        query_graph[i]->buildCoreTable();
        mx = max(mx, (int)query_graph[i]->vertices_count_);
        // all_query_vertex_nums += query_graph[i]->vertices_count_;
        for (auto j : query_graph[i]->all_labs)
            all_lab[j.first] = true;
        // cout << "Q" << i << endl;
        // for (int j = 0; j < query_graph[i]->vertices_count_; ++j) {
        //     ui ncnt;
        //     auto nei = query_graph[i]->getVertexNeighbors(j, ncnt);
        //     for (int k = 0; k < ncnt; k++) {
        //         cout << j << " " << nei[k] << endl;
        //     }
        // }
        // cout << endl;
        // for (int j = 0; j < query_graph[i]->vertices_count_; ++j) {
        //     ui ncnt;
        //     auto nei = query_graph[i]->getVertexDNeighbors(j, ncnt);
        //     for (int k = 0; k < ncnt; k++) {
        //         cout << j << " " << nei[k] << endl;
        //     }
        // }
        // cout << endl;
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
//    cout << "After data graph " << endl;
    // mp = get_index_mem();
    // memory.emplace_back(mp["pk"]);
//    auto end1 = std::chrono::high_resolution_clock::now();
//    double load_graphs_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start).count();

    TreeNode *ceci_tree;
    ui *ceci_order;
    ui *candidates_count;

    ui *idx ;
    ui *idx_count ;
    ui *embedding ;
    vector<vector<ui>> valid_candidates ;
    bool *visited_vertices ;

    size_t int_size;
    size_t all_memory;
    size_t all_query_vertex_nums;

    ui **candidates = nullptr;
    ui* matching_order;
    ui* pivots;

    auto start = std::chrono::high_resolution_clock::now();
    idx = new ui[mx];
    idx_count = new ui[mx];
    embedding = new ui[mx];
    valid_candidates.resize(mx);
    visited_vertices = new bool[data_graph->vertices_count_];

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

    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;

    ui max_vetex_nums = 0;
    for (int i = 0; i < QUERY_NUMS; ++i) { if (query_graph[i]->vertices_count_ > max_vetex_nums) {max_vetex_nums = query_graph[i]->vertices_count_;}}


    all_match_res.reserve(1000);
    size_t memory_cost_in_bytes = 0;
    size_t embedding_count = 0;
    size_t call_count = 0;
    {
//        GenerateFilteringPlan::generateCECIFilterPlan(data_graph, query_graph[0], ceci_tree, ceci_order, -1);
//        cout << "After generateCECIFilterPlan" << endl;
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

        if (TE_Candidates.empty()) {
            cout << "no result" << endl;
            return 0;
        }

        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph[0], candidates_count, ceci_order,
                                                                    ceci_tree,
                                                                    TE_Candidates, NTE_Candidates);

        std::cout << "-----" << std::endl;
        std::cout << "Generate a matching order..." << std::endl;

        size_t order_num = 0;

        sscanf(input_order_num.c_str(), "%zu", &order_num);

        GenerateQueryPlan::generateCECIQueryPlan(query_graph[0], ceci_tree, ceci_order, matching_order, pivots);

//        GenerateQueryPlan::checkQueryPlanCorrectness(query_graph[0], matching_order, pivots);
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph[0], matching_order);

        std::cout << "-----" << std::endl;
        std::cout << "Enumerate..." << std::endl;
        size_t output_limit = 0;
        if (input_max_embedding_num == "MAX") {
            output_limit = std::numeric_limits<size_t>::max();
        } else {
            sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
        }
//        output_limit = 10; // opt
#if ENABLE_QFLITER == 1
        EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

        size_t time_limit = 0;
        sscanf(input_time_limit.c_str(), "%zu", &time_limit);

        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph[0], ceci_tree, candidates,
                                                          candidates_count, TE_Candidates,
                                                          NTE_Candidates, ceci_order, output_limit, call_count,
                                                          all_match_res,idx,idx_count,embedding, valid_candidates,visited_vertices);
        memory_cost_in_bytes += sizeof(int) * all_match_res.size() * all_match_res[0].size();

#ifdef DISTRIBUTION
        std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
        outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
        delete[] EvaluateQuery::distribution_count_;
#endif

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Release memories..." << std::endl;
    }

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
//    cout << "Matchs: " << endl;
//    for (auto i : all_match_res) {
//        for (auto j : i) {
//            cout << j << " ";
//        }
//        cout << endl;
//    }
//    std::cout << "--------------------------------------------------------------------" << std::endl;

    // mp = get_index_mem();
    // memory.emplace_back(mp["pk"]);
//     if (Meth != 1) {
//         string pat = "/home/qsl/" + PATH +"D" + to_string(D) + "/";
//         WR(data_graph, pat, true);
// //        RD(data_graph, pat);
//         cout << getTime() << endl;
//         return 0;
//     }

//    for (auto i : get_nei_label_vertex){
//        for (auto j : i.second) {
//            nei_all_size ++;
//            nei_all_size ++;
//            nei_all_size += j.second.size();
//        }
//    }
//    nei_all_size *= 4;
    // mp = get_index_mem();
    // memory.emplace_back(mp["pk"]);
    sstart = std::chrono::high_resolution_clock::now();
    if (Meth == 2)
        bfs_init_solve(total_time_in_ns, all_match_res, query_graph, data_graph);

//    delete[] valid_candidates;
//    delete[] idx;
//    delete[] idx_count;
//    delete[] embedding;
//    delete[] visited_vertices;
//
//    for (int i = 0; i < QueryNums; ++i)
//        delete query_graph[i];
//    delete[] query_graph;
//    delete data_graph;
//    fss.close();
    return 0;
}
