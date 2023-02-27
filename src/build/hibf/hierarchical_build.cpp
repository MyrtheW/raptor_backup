// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/build/hibf/construct_ibf.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/initialise_max_bin_kmers.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/build/hibf/loop_over_children.hpp>
#include <raptor/build/hibf/update_user_bins.hpp>

namespace raptor::hibf
{

/*!\brief Recursively builds the HIBF from the precomputed layout.
* \details Recursive function to build the Hierarchical Interleaved Bloom Filter from the layout produced by chopper,
 * read into the `node_map`.
 * This version supports the Dynamic HIBF, by taking into account empty bins
 * and by updating the supporting tables (e.g. the FPR and occupancy table).
* \param[in] empty_bin_kmers The variable is recursively passed on, such that merged bins are allocated suffiencient
 * space to accomodate all empty bins in its subtree.
* \param[out] data The data object sores the HIBF and supporting tables.
* \author Adapted by Myrthe Willemsen
*/
template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build_data<data_layout_mode> & data,
                          build_arguments const & arguments,
                          bool is_root,
                          size_t & empty_bin_kmers)
{
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<size_t> kmers{};

    // initialize lower level IBF
    size_t const max_bin_tbs =
    initialise_max_bin_kmers(kmers, ibf_positions, filename_indices, current_node, data, arguments, empty_bin_kmers);
    auto lower_ibf_idx = ibf_positions[data.node_map[current_node].max_bin_index];
    auto && ibf = construct_ibf(kmers, max_bin_tbs, current_node, data, arguments, empty_bin_kmers);
    data.hibf.occupancy_table[ibf_pos].resize(ibf.bin_count());
    data.hibf.fpr_table[ibf_pos].resize(ibf.bin_count()); // Update the FPR and occupancy table for the dynamic HIBF.
    data.hibf.ibf_vector[ibf_pos] = ibf;// This is required for updating the fpr table during insert_into_ibf.
    //insert_into_ibf(parent_kmers, kmers, max_bin_tbs, data.node_map[current_node].max_bin_index, ibf, is_root); // previously the insertion was part of the construct_ibf, however it makes more sense to seperate this.
    insert_into_ibf(parent_kmers, kmers, std::make_tuple((uint64_t) ibf_pos, data.node_map[current_node].max_bin_index,
                            (uint64_t) max_bin_tbs), data.hibf, ibf, is_root); // previously the insertion was part of the construct_ibf, however it makes more sense to seperate this.
    kmers.clear(); // reduce memory peak

    // Parse all other children (merged bins) of the current ibf
    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, arguments, is_root, empty_bin_kmers);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];
        if (std::filesystem::path(record.filenames[0]).extension() !=".empty_bin"){ // Only the first entry of `filenames` stores the actual filename, hence it has to be indexed with [0]
            compute_kmers(kmers, arguments, record);
            insert_into_ibf(parent_kmers, kmers, std::make_tuple((uint64_t) ibf_pos, (uint64_t) record.bin_indices.back(),
                            (uint64_t) record.number_of_bins.back()), data.hibf, ibf, is_root);
            kmers.clear();
        }else{ // if we are dealing with an empty bin, their size will be extracted and and added to `empty_bin_kmers`, which is to be passed on to parent merged bins.
            std::string kmer_count = (std::string) std::filesystem::path(record.filenames[0]).stem(); // The empty bin's intended size will be extracted from the filename
            double kmer_count_double = ::atof(kmer_count.c_str());
            empty_bin_kmers += static_cast<size_t>(kmer_count_double); // The empty bin's size is added to `empty_bin_kmers`, which is to be passed on to the parent merged bins of the empty bin.
        }
        update_user_bins(data, filename_indices, record);
    }

    data.hibf.ibf_vector[ibf_pos] = std::move(ibf); // assign the ibf again to the ibf_vector to make sure that the updates are incorporated.
    data.hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
    data.hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}


template size_t hierarchical_build<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                      lemon::ListDigraph::Node const &,
                                                                      build_data<seqan3::data_layout::uncompressed> &,
                                                                      build_arguments const &,
                                                                      bool,
                                                                      size_t &);

template size_t hierarchical_build<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                    lemon::ListDigraph::Node const &,
                                                                    build_data<seqan3::data_layout::compressed> &,
                                                                    build_arguments const &,
                                                                    bool,
                                                                    size_t &);

} // namespace raptor::hibf

