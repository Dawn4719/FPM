//
// Created by ssunah on 6/22/18.
//

#ifndef SUBGRAPHMATCHING_GRAPH_H
#define SUBGRAPHMATCHING_GRAPH_H

#include <unordered_map>
#include <iostream>
#include <cmath>
#include "configuration/types.h"
#include "configuration/config.h"
#include <cstring>
/**
 * A graph is stored as the CSR format.
 */

class Graph {
public:
    bool enable_label_offset_;
    int id_;
    ui l_count;
    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;
    int emb_count;
    ui* offsets_;
    VertexID * neighbors_;
    LabelID* labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;
    bool has_dir_;
    bool is_q_;
    // direct
    ui* doffsets_;  // degree sum
    VertexID * dneighbors_;

    int* core_table_;
    ui core_length_;
    std::unordered_map<int, bool> all_labs;
    std::unordered_map<LabelID, ui> labels_frequency_;

#if OPTIMIZED_LABELED_GRAPH == 1
    ui* labels_offsets_;
    std::unordered_map<LabelID, ui>* nlf_;
#endif

private:
    void BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    void BuildNLF();
    void BuildLabelOffset();
#endif

public:
    Graph(const bool enable_label_offset, bool has_dir) {
        enable_label_offset_ = enable_label_offset;
        id_ = -1;
        l_count = 0;
        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        core_length_ = 0;
        has_dir_ = has_dir;
        offsets_ = NULL;
        neighbors_ = NULL;
        doffsets_ = NULL;  // degree sum
        dneighbors_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();

#if OPTIMIZED_LABELED_GRAPH == 1
        labels_offsets_ = NULL;
        nlf_ = NULL;
#endif
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        if (has_dir_ || is_q_) {
            delete[] doffsets_;
            delete[] dneighbors_;
        }
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
#if OPTIMIZED_LABELED_GRAPH == 1
        delete[] labels_offsets_;
        delete[] nlf_;
#endif
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromFileCompressed(const std::string& degree_path, const std::string& edge_path,
                                     const std::string& label_path);
    void storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
                               const std::string& label_path);
    void printGraphMetaData();
public:
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        if (offsets_[id + 1] < offsets_[id]) return 0;
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getVertexDDegree(const VertexID id) const {
        if (has_dir_) {
            if (doffsets_[id + 1] < doffsets_[id]) return 0;
            return doffsets_[id + 1] - doffsets_[id];
        }
        else {
            if (offsets_[id + 1] < offsets_[id]) return 0;
            return offsets_[id + 1] - offsets_[id];
        }
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const ui getCoreValue(const VertexID id) const {
        return core_table_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        if (offsets_[id + 1] < offsets_[id]) count = 0;
        return neighbors_ + offsets_[id];
    }
    // d
    const ui * getVertexDNeighbors(const VertexID id, ui& count) const {
        if (has_dir_) {
            count = doffsets_[id + 1] - doffsets_[id];
            if (doffsets_[id + 1] < doffsets_[id]) count = 0;
            return dneighbors_ + doffsets_[id];
        }
        else {
            count = offsets_[id + 1] - offsets_[id];
            if (offsets_[id + 1] < offsets_[id]) count = 0;
            return neighbors_ + offsets_[id];
        }
    }

    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        if (id + 1 > labels_count_) {
            count = 0;
            return reverse_index_;
        }
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

#if OPTIMIZED_LABELED_GRAPH == 1
    const ui * getNeighborsByLabel(const VertexID id, const LabelID label, ui& count) const {
        ui offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return nlf_ + id;
    }

    bool checkEdgeExistence(const VertexID u, const VertexID v, const LabelID u_label) const {
        ui count = 0;
        const VertexID* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
#endif

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors = getVertexNeighbors(v, count);
        std::cout << u << " " << v << std::endl;
        for (int i = 0; i < count; ++i) {
            std::cout << neighbors[i] << " ";
        }
        std::cout << std::endl;
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void buildCoreTable();

    void move(const Graph* other) {
        enable_label_offset_ = other->enable_label_offset_;
        id_ = other->id_;
        l_count = other->l_count;
        vertices_count_ = other->vertices_count_;
        edges_count_ = other->edges_count_;
        labels_count_ = other->labels_count_;
        max_degree_ = other->max_degree_;
        max_label_frequency_ = other->max_label_frequency_;
        emb_count = other->emb_count;
        has_dir_ = other->has_dir_;
        is_q_ = other->is_q_;
        core_length_ = other->core_length_;
        all_labs = other->all_labs;
        labels_frequency_ = other->labels_frequency_;


        offsets_ = new ui[vertices_count_ +  1]; // csr
        memcpy(offsets_, other->offsets_, (vertices_count_ + 1) * sizeof(ui));
        // std::fill_n(offsets_, vertices_count_ + 1, 0);
        if (has_dir_) {
            neighbors_ = new VertexID[edges_count_ * 2];
            memcpy(neighbors_, other->neighbors_, (edges_count_ * 2) * sizeof(VertexID));
            dneighbors_ = new VertexID[edges_count_];
            memcpy(dneighbors_, other->dneighbors_, (edges_count_) * sizeof(VertexID));
            doffsets_ = new ui[vertices_count_ +  1];
            doffsets_[0] = 0;
            std::fill_n(doffsets_, vertices_count_ + 1, 0);
        }
        else {
            neighbors_ = new VertexID[edges_count_ + 10];
            memcpy(neighbors_, other->neighbors_, (edges_count_ + 10) * sizeof(VertexID));
        }
        labels_ = new LabelID[vertices_count_];
        memcpy(labels_, other->labels_, (vertices_count_) * sizeof(VertexID));
        reverse_index_ = new ui[vertices_count_];
        memcpy(reverse_index_, other->reverse_index_, (vertices_count_) * sizeof(VertexID));
        reverse_index_offsets_= new ui[labels_count_ + 1];
        memcpy(reverse_index_offsets_, other->reverse_index_offsets_, (labels_count_ + 1) * sizeof(VertexID));

#if OPTIMIZED_LABELED_GRAPH == 1
        nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];
        for (ui i = 0; i < vertices_count_; ++i) {
            nlf_[i] = other->nlf_[i];
        }
        // labels_offsets_ = new ui[(size_t)vertices_count_ * labels_count_ + 1];
        // memcpy(labels_offsets_, other->labels_offsets_, ((size_t)vertices_count_ * labels_count_ + 1) * sizeof(ui));
#endif
    }
};


#endif //SUBGRAPHMATCHING_GRAPH_H
