// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/bin_size_in_bits.hpp>
#include <raptor/build/hibf/construct_ibf.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data<data_layout_mode> & data,
                                                 build_arguments const & arguments,
                                                 bool is_root)
{
    auto & node_data = data.node_map[node];
    auto const & record = node_data.remaining_records[0]; //myrthe, check if it is always the case that initialized max bin kmers is called before.
    unsigned long kmers_per_bin{};

//    assert (record.filenames[idx] != "empty_bin"){ // or kmers.size() == 0. myrthe
//       //myrthe varify with Svenja that this is correct //size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(record.estimated_sizes[0]) / number_of_bins))}; //assuming estimated_sizes[0] belongs to the first bin?
    if (node_data.favourite_child == lemon::INVALID //not a merged bin
    and std::filesystem::path(record.filenames[0]).extension() ==".empty_bin"){ // empty bin
        std::string kmer_count = (std::string) std::filesystem::path(record.filenames[0]).stem();
        double kmer_count_double = ::atof(kmer_count.c_str());
        kmers_per_bin = static_cast<size_t>(std::ceil(static_cast<double>(kmer_count_double/ number_of_bins)));
    }else{
        double kmer_size = static_cast<double> (kmers.size());
        if (node_data.favourite_child == lemon::INVALID){ //merged bin
            // how many empty bin children? use index.ibf().occupancy_table(kmers.size(), ibf_idx, start_bin_idx, number_of_bins);

        }
        kmers_per_bin=static_cast<size_t>(std::ceil(static_cast<double>(kmer_size) / number_of_bins));
    }

    double const bin_bits{static_cast<double>(bin_size_in_bits(arguments, kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{arguments.hash}};

    insert_into_ibf(parent_kmers, kmers, number_of_bins, node_data.max_bin_index, ibf, is_root);

    return ibf;
}

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                 robin_hood::unordered_flat_set<size_t> &,
                                                 size_t const,
                                                 lemon::ListDigraph::Node const &,
                                                 build_data<seqan3::data_layout::uncompressed> &,
                                                 build_arguments const &,
                                                 bool);

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                               robin_hood::unordered_flat_set<size_t> &,
                                               size_t const,
                                               lemon::ListDigraph::Node const &,
                                               build_data<seqan3::data_layout::compressed> &,
                                               build_arguments const &,
                                               bool);

template <seqan3::data_layout data_layout_mode>
seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data<data_layout_mode> & data,
                                                 build_arguments const & arguments,
                                                 bool is_root,
                                                 size_t lower_ibf_idx,
                                                 size_t & empty_bin_kmers)
{
    auto & node_data = data.node_map[node];
    auto const & record = node_data.remaining_records[0]; //myrthe, check if it is always the case that initialized max bin kmers is called before.
    unsigned long kmers_per_bin{};

//    assert (record.filenames[idx] != "empty_bin"){ // or kmers.size() == 0. myrthe
//       //myrthe varify with Svenja that this is correct //size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(record.estimated_sizes[0]) / number_of_bins))}; //assuming estimated_sizes[0] belongs to the first bin?
    if (node_data.favourite_child == lemon::INVALID //not a merged bin
    and std::filesystem::path(record.filenames[0]).extension() ==".empty_bin"){ // empty bin
        std::string kmer_count = (std::string) std::filesystem::path(record.filenames[0]).stem();
        double kmer_count_double = ::atof(kmer_count.c_str());
        kmers_per_bin = static_cast<size_t>(std::ceil(static_cast<double>(kmer_count_double/ number_of_bins)));
        empty_bin_kmers += static_cast<size_t>(kmer_count_double);
    }else{
        double kmer_size = static_cast<double> (kmers.size());
        if (node_data.favourite_child != lemon::INVALID){ //merged bin
                    // if it is a merged bin, check if it has empty bin childeren. if so, increase kmers_per_bins appropiately.

//            lower_ibf_idx = data.hibf.ibf().next_ibf_id(ibf_pos)[0]; // is this already filled out??
//             //how many empty bin children? use
//            int number_empty_bins =0;
//            auto& lower_ibf_occupancy = data.hibf.occupancy_table[lower_ibf_idx];
//            for (auto bin_idx=0; bin_idx< (int) lower_ibf_occupancy.size(); bin_idx++){
//                if (lower_ibf_occupancy[bin_idx] == 0){number_empty_bins++;};
//            }
//            kmer_size += (1+number_empty_bins/(lower_ibf_occupancy.size() - number_empty_bins));
//            // use ibf_pos as ibf_idx and start_bin_idx = 0 and next_ibf_id. data.

            kmer_size += empty_bin_kmers;
        }
        kmers_per_bin=static_cast<size_t>(std::ceil(static_cast<double>(kmer_size) / number_of_bins));
    }

    double const bin_bits{static_cast<double>(bin_size_in_bits(arguments, kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{arguments.hash}};

    //insert_into_ibf(parent_kmers, kmers, number_of_bins, node_data.max_bin_index, ibf, is_root);

    return ibf;
}

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                 robin_hood::unordered_flat_set<size_t> &,
                                                 size_t const,
                                                 lemon::ListDigraph::Node const &,
                                                 build_data<seqan3::data_layout::uncompressed> &,
                                                 build_arguments const &,
                                                 bool,
                                                 size_t,
                                                 size_t &);

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                               robin_hood::unordered_flat_set<size_t> &,
                                               size_t const,
                                               lemon::ListDigraph::Node const &,
                                               build_data<seqan3::data_layout::compressed> &,
                                               build_arguments const &,
                                               bool,
                                               size_t,
                                               size_t &);

} // namespace raptor::hibf
