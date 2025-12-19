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

fstream fin = fstream("/home/qsl/exp/FPM/result.csv", ios::app);
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
// mutex mtx;
// mutex mtx_in_match;
// mutex mtx2;
// mutex mtx3;
double pathtime;

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
vector<pair<int, int>> to_match_idx;
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
    vector<int> chooedNode;
    map<int, vector<pair<int, vector<int>>>> get_nei_label_vertex;

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
    map<int, vector<vector<vector<int>>>> his;


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
//                    if (u == 1917)
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
        fstream fss = fstream("/home/qsl/exp/FPM/Dataset/" + PATH +"D" + to_string(D) + "/" + to_string(vertex) + ".txt");
        vector<pair<int, vector<int>>> nodes;
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
            nodes.resize(nodes.size() + 1);
            nodes.back().first = lab_;
            nodes.back().second.reserve(num_);
            for (int i = 0; i < num_; ++i) {
                fss >> x;
                nodes.back().second.emplace_back(x);
            }
            fss >> vertex;
        }
        fss.close();
        if (!nodes.empty()) {
            get_nei_label_vertex[lower_] = std::move(nodes);
            return true;
        }
        return false;
    }
    bool add(int vertex, size_t lower_, Graph* data_graph, vector<pair<int, vector<int>>>& ve) {
        // cout << "add" << endl;
        queue<int> q;
        q.emplace(vertex);
        map<int, bool> st;
        map<int, int> dist;

        dist[vertex] = 0;
        map<int, vector<int>> lab;
        while (!q.empty()) {
            auto p = q.front(); q.pop();
            auto d = dist[p];
            if (st[p]) continue;
            if (d >= D) continue;
            st[p] = true;

            ui nbrs_cnt;
            const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbrs_cnt); // 获取点vertex的邻居
            for (int j = 0; j < nbrs_cnt; ++j) { // 遍历点vertex的邻居
                int v = nbrs[j];
                if (!st[v]) {
                    if (d + 1 <= D)
                        q.emplace(v);
                    if (dist.find(v) == dist.end()) {
                        auto lable = data_graph->getVertexLabel(v);
                        lab[lable].emplace_back(v);
                        dist[v] = d + 1;
                    }
                }
            }
        }
        // get_nei_label_vertex[lower_].reserve(lab.size());
        for (auto i : lab) {
            ve.emplace_back(i);
        }
        return true;
    }
    void extd(Graph* data_graph) {

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

            /*vector<vector<int>> matchres;
            for (auto i : now) {
                matchres.emplace_back(N[i.first].n[i.second].match);
                for (auto v : N[i.first].n[i.second].match)
                    sum += v;
            }
            sort(matchres.begin(), matchres.end());

            if (his.find(sum) != his.end()) {
                auto it = his[sum];
                bool pa = true;
                int pa_cnt = 0;
                for (int k = 0; k < it.size(); ++k) {
                    pa = true;
                    for (int i = 0; i < QUERY_NUMS; ++i) {
                        if (matchres[i].size() != it[k][i].size()) {
                            pa = false;
                            break;
                        }
                        for (int j = 0; j < matchres[i].size(); ++j) {
                            if (matchres[i][j] != it[k][i][j]) {
                                pa = false;
                                break;
                            }
                        }
                        if (!pa) break;
                    }
                    if (pa)
                        return;
                }
            }
            his[sum].emplace_back(matchres);
            */


            // mtx.lock();
            sp_match_res++;
//            cout << sp_match_res << endl;
            // mtx.unlock();
            cur_sp_match_res++;
            return;
        }
        for (auto i : now) {
            if (i.first < tree.size())
                for (const auto& j : tree[i.first][i.second]) {
                    if (qim[j.first]) continue;
                    for (auto k : j.second) {
                        if (chooedNode[j.first] < k) {
//                            if (sp_match_res >= 3000000) {
//                                topk = true;
//                                return;
//                            }
                            qim[j.first] = true;
                            now.emplace_back(j.first, k);
                            // if (QUERY_NUMS >= 4 && now.size() == QUERY_NUMS - 1)
                            //     Nst[i.first][i.second][j.first][k] = true;
                            auto ti = gti;
                            auto ii = now.back();
                            for (auto v : N[ii.first].n[ii.second].match) if (pst[v] >= 0) pst[v] = -ii.first - 1;
                            vector<pair<int, int>> pa;
                            for (auto v : N[ii.first].n[ii.second].match) {
//                                if (sp_match_res >= 3000000) {
//                                    topk = true;
//                                    return;
//                                }
                                get_path3D(data_graph, ii.first + 1, v, -1, 0, pa);
                            }

                            Path.emplace_back(pa);
                            auto tie = gti;
                            pathtime += cti(tie, ti);
                            // cout << "1 " << now.size() << " " << Path.size() << endl;
                            // assert(Path.size() + 1 == now.size());
                            auto oldNode = chooedNode[j.first];
                            chooedNode[j.first] = k;
                            extd(data_graph);
                            chooedNode[j.first] = oldNode;
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

        if (cnt >= QUERY_NUMS) return true;
        for (const auto& i : tree[u.first][u.second]) {
//            cout << cnt << endl;
            if (qim[i.first]) continue;
            for (auto j : i.second) {


                if (j < chooedNode[i.first]) continue;
                auto old = chooedNode[i.first];
                chooedNode[i.first] = j;
                now.emplace_back(i.first, j);

                auto ii = now.back();
                auto ti = gti;
                for (auto v : N[ii.first].n[ii.second].match) if (pst[v] >= 0) pst[v] = -ii.first - 1;
                vector<pair<int, int>> pa;
                for (auto v : N[ii.first].n[ii.second].match) {
//                    if (sp_match_res >= 3000000) {
//                        topk = true;
//                        return true;
//                    }
                    get_path3D(data_graph, ii.first + 1, v, -1, 0, pa);
                }
                Path.emplace_back(pa);
                auto tie = gti;
                pathtime += cti(tie, ti);
                // cout << "3 " << now.size() << " " << Path.size() << endl;
                qim[i.first] = true;
                if (i.first < tree.size() && !dfs(data_graph, {i.first, j}, cnt+1)) {
                    chooedNode[i.first] = old;
                    return false;
                }



                n_st.clear();
                extd(data_graph);
                chooedNode[i.first] = old;


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
                // for (auto each_match_num : match_num) if (each_match_num == 0) return false;
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
//    map<int, map<int, size_t>> epoch_cnt;
    void begin_extend(Graph* data_graph, Graph** query_graph, int match_idx) {
        auto each_begin1 = std::chrono::high_resolution_clock::now();
//        cout<< 'bb' << endl;
        N.clear();
        unq.clear();
        his.clear();
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
//        epoch_cnt.clear();
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
        vector<pair<int, vector<int>>> lab_ver;
        for (int epoch = 0; epoch < QUERY_NUMS - 1; ++epoch) {
            for (int j = 0; j < QUERY_NUMS; ++j) {
                tree[j].resize(N[j].n.size());
                for (int n_j = 0; n_j < N[j].n.size(); ++n_j) {
                    auto& n = N[j].n[n_j];
                    if (n.step == epoch)
                        {
                        if (!n.flag && n.step < QUERY_NUMS - 1) {
                            for (auto match_vertex : n.match) {
//                            auto lower_ = LB2(get_nei_label_vertex, match_vertex);
                                if (get_nei_label_vertex[match_vertex].empty()) {
                                    if (!add(match_vertex, match_vertex))
                                        continue;
                                }
                                lab_ver = get_nei_label_vertex[match_vertex];

                                // // get_nei_label_vertex.clear();
                                // lab_ver.clear();
                                // add(match_vertex, match_vertex, data_graph, lab_ver);

                                for (int i = 1; i < QUERY_NUMS; ++i) {
                                    if (i == j) continue;
                                    for (auto LB_ : query_graph[i]->all_labs) {
                                        auto lb = LB_.first;

                                        ui u_cnt;
                                        auto lb2v = query_graph[i]->getVerticesByLabel(lb, u_cnt);
                                        auto lower = LB(lab_ver, lb);
                                        if (lower < lab_ver.size() && lab_ver[lower].first == lb) {
                                            for (auto vertex : lab_ver[lower].second) {
                                                if (N[i].has_extend.find(vertex) != N[i].has_extend.end()) {
                                                    if (!N[i].has_extend[vertex].empty()) {
                                                        // add_edge
                                                        for (auto edges : N[i].has_extend[vertex]) {
                                                            int n_i = edges;
                                                            if (n.step < N[i].n[n_i].step) {
//                                                                if (i > tree[j][n_j].size()) tree[j][n_j].resize(QUERY_NUMS);
                                                                auto lb_ = lower_bound(tree[j][n_j][i].begin(), tree[j][n_j][i].end(), n_i);
                                                                if (lb_ == tree[j][n_j][i].end() || *lb_ != n_i) {
                                                                    tree[j][n_j][i].insert(lb_, n_i);
                                                                    match_num[i]++;
                                                                }
                                                            }
                                                        }
                                                    }
                                                    continue;
                                                }

                                                // N[i].has_extend[vertex] = vector<int> {};
                                                // N[i].has_extend[vertex].reserve(10);
                                                for (int u_ = 0; u_ < u_cnt; ++u_) {
                                                    auto u = lb2v[u_];
                                                    try_count++;
                                                    vector<vector<int>> cur_match_res;
                                                    dfs_all_count++;
                                                    CECI(data_graph, query_graph, cur_match_res, i, u, vertex);
                                                    //                                                 cout << i << " " << u << " " << vertex << endl;
                                                    if (cur_match_res.empty()) continue;
                                                    for (auto& each_res: cur_match_res) {
                                                        sum = 0;
                                                        sort(each_res.begin(), each_res.end());
                                                        for (auto each_res_i: each_res) sum += each_res_i;
                                                        if (N[i].sum2idx.find(sum) == N[i].sum2idx.end())
                                                            {
                                                            N[i].sum2idx[sum].emplace_back(N[i].n.size());
                                                            sort(each_res.begin(), each_res.end());
                                                            Node a;
                                                            a.match = each_res;
                                                            a.step = N[j].n[n_j].step + 1;
                                                            N[i].n.emplace_back(a);
//                                                            epoch_cnt[a.step][i]++;
                                                            match_num[i]++;
                                                            map<int, vector<int>> tmp;
                                                            tmp[j].emplace_back(n_j);

                                                            tree[j][n_j][i].emplace_back(N[i].n.size() - 1);
                                                            N[i].has_extend[vertex].emplace_back(N[i].n.size() - 1);

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
                                                                tree[j][n_j][i].emplace_back(N[i].n.size() - 1);
                                                                N[i].has_extend[vertex].emplace_back(N[i].n.size() - 1);
                                                            }
                                                            else {
                                                                if (n.step <= N[i].n[n_i].step) {
                                                                    auto lb_ = lower_bound(tree[j][n_j][i].begin(), tree[j][n_j][i].end(), n_i);
                                                                    if (lb_ == tree[j][n_j][i].end() || *lb_ != n_i) {
                                                                        tree[j][n_j][i].insert(lb_, n_i);
                                                                        match_num[i]++;
                                                                   if (!N[i].has_extend[vertex].empty() && N[i].has_extend[vertex][0] == n_j) continue;
                                                                        N[i].has_extend[vertex].emplace_back(n_i);
                                                                    }

                                                                }
                                                                else {
                                                                    if (N[i].has_extend.find(vertex) == N[i].has_extend.end()) {
                                                                        auto p = N[i].has_extend[vertex];
                                                                        sort(p.begin(), p.end());
                                                                        auto lb2_ = lower_bound(p.begin(), p.end(), n_i); //
                                                                        if (lb2_ == p.end() || *lb2_ != n_i) {
                                                                            // cout << endl;
                                                                            N[i].has_extend[vertex].insert(lb2_, n_i);
                                                                        }
                                                                        // assert(lb__ < p.end() && *lb__ == n_i);
                                                                    }
                                                                    else {
                                                                        N[i].has_extend[vertex].emplace_back(n_i);
                                                                    }
                                                                }

                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                        n.flag = true;
                    }

                }

            }
        }

//        cout << "build over: " << Worker_idx << endl;
        auto each_end1 = std::chrono::high_resolution_clock::now();
        // mtx2.lock();

        tim1 += cti(each_end1, each_begin1);
//        cout << "Worker_idx: " << Worker_idx << " idx_time: " << NC(std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count()) << endl;
        // mtx2.unlock();
//    this_time += std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count();
        // cout << "build time: " << NC(std::chrono::duration_cast<std::chrono::nanoseconds>(each_end1 - each_begin1).count()) << endl;
        //    g.clear();
        //    cout << NC(ti) << "s" << endl;
        for (auto & i : N) {
            if (i.n.empty())
                return;
        }

        get_nei_label_vertex.clear();
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
        // mtx3.lock();
        // tim2 = max(tim2, (double)std::chrono::duration_cast<std::chrono::nanoseconds>(each_end2 - each_begin2).count());
        tim2 += cti(each_end2, each_begin2);
//        cout << "Work_idx: " << Worker_idx << " graph time: " << NC(std::chrono::duration_cast<std::chrono::nanoseconds>(each_end2 - each_begin2).count()) << "ms " << cur_sp_match_res << endl;
        // mtx3.unlock();

        Path.clear();

    }
    void Excute(double total_time_in_ns, Graph** query_graph, Graph* data_graph) {
        chooedNode.resize(QUERY_NUMS, -1);
        for (int i = Worker_idx; i < all_match_res.size(); i += thread_num) {
//            get_nei_label_vertex.clear();
            cout << "Worker_idx: " << Worker_idx << " " << i << endl;

            begin_extend(data_graph, query_graph, to_match_idx[i].second);
        }
    }
};
vector<double> times;

void exc(Graph* data_graph, Graph** query_graph, int Worker_idx) {
    Worker worker(data_graph, query_graph, Worker_idx);
    worker.Excute(0, query_graph, data_graph);
}

void bfs_init_solve(double total_time_in_ns, vector<vector<int>>& all_match_res, Graph** query_graph, Graph* data_graph) {
    auto sstart2 = std::chrono::high_resolution_clock::now();
    std::cout << "Begin SP Match!" << std::endl;
//    top /= thread_num;

    thread_num = min(thread_num, (int)all_match_res.size());
    cout << "PATH: " << PATH << " Query_Num: " << QUERY_NUMS << " D: " << D << " Meth: bfs_init " << thread_num << " thread_num " << thread_num << dl;

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
    // exc(data_graph, query_graph, 0);
    std::vector<std::thread> threads;

    for (int i = 0; i < thread_num; ++i) {
        threads.emplace_back(exc, data_graph, query_graph, i);
    }
    for (std::thread& t : threads) {
        t.join();
    }

    auto eend = std::chrono::high_resolution_clock::now();
    double ALL_TIME = std::chrono::duration_cast<std::chrono::nanoseconds>(eend - sstart2).count();
    sum_memory_cost_in_bytes += sizeof (bool)* data_graph->vertices_count_;
    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    // printf("Memory cost (MB): %.4lf, %.4lf\n", BYTESTOMB(all_memory + sum_memory_cost_in_bytes + init_memory_cost_in_bytes), BYTESTOMB(init_memory_cost_in_bytes));
    printf("#Embeddings: %d\n", sp_match_res);
    // printf("DFS call count: %d\n", dfs_all_count);
    // printf("try count: %d\n", try_count);
    printf("ALL Time: %.4lf\n", NANOSECTOSEC(ALL_TIME + total_time_in_ns));
    cout << NC(ALL_TIME) << " " << NC(total_time_in_ns) << " " << NC(tim1) << " " << NC(tim2) << endl;
    std::cout << "End." << std::endl;

    fin.close();
}


int main(int argc, char** argv) {
    cout << getTime() << endl;
    cout << "begin" << endl;
    auto mp = get_index_mem();
    memory.emplace_back(mp["pk"]);
    test = false;
    // la = get_index_mem()["pk"];
    QUE = "6";
    PATH = "DD";
    D = 1;
    QUERY_NUMS = 3;
    Meth = 2;
//    input_query_graph_file = "../../../test/querys/" + QUE;
//    input_data_graph_file = "../../../test/graphs/Dataset/" + PATH + "/" + PATH + ".txt";
    thread_num = 1;
    MatchingCommand command(argc, argv);
    QUE = command.getQueryGraphFilePath();
    PATH = command.getDataGraphFilePath();
    D = stoi(command.getDist());
    QUERY_NUMS = stoi(command.getQueryNums());
    thread_num = stoi(command.getThread());
    Meth = 2;
    input_query_graph_file = "/home/qsl/exp/FPM/q/" + QUE;
    input_data_graph_file = "/home/qsl/exp/FPM/Dataset/" + PATH + ".txt";
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
    sort(all_match_res.begin(), all_match_res.end());
    std::cout << "query_graph match result:" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    to_match_idx.resize(all_match_res.size());
    for (int i = 0; i < all_match_res.size(); ++i) {
        to_match_idx[i].first = 0;
        to_match_idx[i].second = i;
        for (auto v : all_match_res[i]) {
//            auto p = data_graph->getVertexNLF(v);
//            if (p.f)
            to_match_idx[i].first += data_graph->getVertexDegree(v);
        }
    }
    // if (thread_num == 1)
    //     sort(to_match_idx.begin(), to_match_idx.end(), [](pair<int, int>& a, pair<int, int>& b){
    //         return a.first > b.first;
    //     });
    // else
    //     sort(to_match_idx.begin(), to_match_idx.end(), [](pair<int, int>& a, pair<int, int>& b){
    //         return a.first < b.first;
    //     });
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

    // delete[] valid_candidates;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;

    for (int i = 0; i < QueryNums; ++i)
        delete query_graph[i];
    delete[] query_graph;
    delete data_graph;
//    fss.close();
    return 0;
}
