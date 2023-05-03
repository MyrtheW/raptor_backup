// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
robin_hood::unordered_flat_set<size_t> create_ibfs_from_chopper_pack(build_data<data_layout_mode> & data, build_arguments const & arguments, bool is_root)
{
    read_chopper_pack_file(data, arguments.bin_file);
    lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
    robin_hood::unordered_flat_set<size_t> root_kmers{};

    size_t const t_max{data.node_map[root].number_of_technical_bins};
    data.compute_fp_correction(t_max, arguments.hash, arguments.fpr);

    size_t empty_bin_kmers=0;
    hierarchical_build(root_kmers, root, data, arguments, is_root, empty_bin_kmers);

    data.hibf.fpr_max = arguments.fpr; // This parameter is added here as member to the HIBF datastructure, such that they will be stored and can be used during search or updating. Note: `window_size`, `compressed`, `parts`,  `bin_path` and `shape` stored as part of the `index` datastructure. `kmer_size` could be obtained from `shape`,
    data.hibf.t_max = t_max;
    // prepare index with data structures to allow for updates
    data.hibf.user_bins.initialize_filename_position_to_ibf_bin();
    data.hibf.initialize_previous_ibf_id();
    data.hibf.initialize_ibf_sizes();
    return root_kmers; // this is used when rebuilding, to insert into the parent merged bin.
}

template robin_hood::unordered_flat_set<size_t>
create_ibfs_from_chopper_pack<seqan3::data_layout::uncompressed>(build_data<seqan3::data_layout::uncompressed> &,
                                                                 build_arguments const &, bool is_root);
template robin_hood::unordered_flat_set<size_t>
create_ibfs_from_chopper_pack<seqan3::data_layout::compressed>(build_data<seqan3::data_layout::compressed> &,
                                                               build_arguments const & arguments, bool is_root);

} // namespace raptor::hibf
