// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/initialise_max_bin_kmers.hpp>
#include <raptor/build/hibf/update_user_bins.hpp>

namespace raptor::hibf
{
/*!\brief Computes the k-mers of the maximum bin of the current IBF to the kmers object.
* \details If the max we are dealing with an empty bin, no k-mers will be computed.
 * Instead the k-mer-estimate of the empty bin will be used within the `construct_ibf` function.
* \param[out] kmers The set object to which the k-mers of the max bin will be inserted.
* \param[in] ibf_positions
* \param[in] filename_indices
* \param[in] node
* \param[out] data The data object sores the HIBF and supporting tables.
* \param[in] arguments
* \param[in] empty_bin_kmers The variable is recursively passed on, such that merged bins are allocated suffiencient
 * space to accomodate all empty bins in its subtree.
 * \return The number of bins needed to store the max bin in this IBF.
* \author Adapted by Myrthe Willemsen //TODO
*/
template <seqan3::data_layout data_layout_mode>
size_t initialise_max_bin_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                                std::vector<int64_t> & ibf_positions,
                                std::vector<int64_t> & filename_indices,
                                lemon::ListDigraph::Node const & node,
                                build_data<data_layout_mode> & data,
                                build_arguments const & arguments,
                                size_t & empty_bin_kmers)
{
    auto & node_data = data.node_map[node];

    if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
    {
        // recursively initialize favourite child first
        ibf_positions[node_data.max_bin_index] =
            hierarchical_build(kmers, node_data.favourite_child, data, arguments, false, empty_bin_kmers);
        return 1;
    }
    else // max bin is not a merged bin
    {
        // we assume that the max record is at the beginning of the list of remaining records.
        auto const & record = node_data.remaining_records[0];
        if (std::filesystem::path(record.filenames[0]).extension() !=".empty_bin"){ // If instead we are dealing with an empty bin, no kmers will be computed. Instead the kmer-estimate of the empty bin will be used within `construct_ibf`.
            compute_kmers(kmers, arguments, record);
            update_user_bins(data, filename_indices, record);
        }
        return record.number_of_bins.back();
    }
}

template size_t initialise_max_bin_kmers<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                            std::vector<int64_t> &,
                                                                            std::vector<int64_t> &,
                                                                            lemon::ListDigraph::Node const &,
                                                                            build_data<seqan3::data_layout::uncompressed> &,
                                                                            build_arguments const &,
                                                                            size_t &);

template size_t initialise_max_bin_kmers<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                          std::vector<int64_t> &,
                                                                          std::vector<int64_t> &,
                                                                          lemon::ListDigraph::Node const &,
                                                                          build_data<seqan3::data_layout::compressed> &,
                                                                          build_arguments const &,
                                                                          size_t &);

} // namespace raptor::hibf
