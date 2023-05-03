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

/*!\brief Creates a configuration object which is passed to chopper's layout algorithm.
* \details Calculates the number of kmers that have to be stored per bloom filter (i.e. bin) and creates an IBF with this bin size.
 * Previously the insertion of the maximal bin was also part of the construct_ibf function, however as this is another operation, it makes more sense to seperate this.
* \param[in] kmers The kmers
* \param[in] number_of_bins
* \param[in] node
* \param[out] data The data object sores the HIBF and supporting tables.
* \param[in] arguments
* \param[in] empty_bin_kmers The variable is recursively passed on, such that merged bins are allocated suffiencient
 * space to accomodate all empty bins in its subtree.
* \return ibf; the IBF that has been created.
* \author Adapted by Myrthe Willemsen
*/
template <seqan3::data_layout data_layout_mode>
seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data<data_layout_mode> & data,
                                                 build_arguments const & arguments,
                                                 size_t & empty_bin_kmers)
{
    auto & node_data = data.node_map[node];
    auto const & record = node_data.remaining_records[0];

    // Calculate the number of kmers that have to be stored per bloom filter (i.e. bin).
    unsigned long kmers_per_bin{};
    if (node_data.favourite_child == lemon::INVALID // not a merged bin
    and std::filesystem::path(record.filenames[0]).extension() ==".empty_bin"){ // we are dealing with an empty bin
        std::string kmer_count = (std::string) std::filesystem::path(record.filenames[0]).stem(); // the empty bin's intended size will be extracted from the filename
        double kmer_count_double = ::atof(kmer_count.c_str()); // TODO, perhaps it would be better to take the size from the next UB, to ensure that it is not larger than the HLL size estimate, which could cause ver unnesessary rebuilding.
        empty_bin_kmers += static_cast<size_t>(kmer_count_double); // the empty bin's size is added to `empty_bin_kmers`, which is to be passed on to the parent merged bins of the empty bin.
        kmers_per_bin = static_cast<size_t>(std::ceil(static_cast<double>(kmer_count_double / number_of_bins))); // For the construction of the IBF we need to calculate the length of a single bloom filter (bin size), by dividing the number of kmers of the UB by the number of bins among which it will be stored.
    }else{
        double kmer_size = static_cast<double> (kmers.size());
        if (node_data.favourite_child != lemon::INVALID){ // merged bin
            kmer_size += empty_bin_kmers; // in case we are initializing an IBF where the largest bin is a merged bin, we must take into account the `empty_bin_kmers` that it should be able to store after future updates, next to all its children `kmers`.
        }
        kmers_per_bin=static_cast<size_t>(std::ceil(static_cast<double>(kmer_size) / number_of_bins));
    }

    double const bin_bits{static_cast<double>(bin_size_in_bits(arguments, kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{arguments.hash}};

    return ibf;
}

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                 size_t const,
                                                 lemon::ListDigraph::Node const &,
                                                 build_data<seqan3::data_layout::uncompressed> &,
                                                 build_arguments const &,
                                                 size_t &);

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                               size_t const,
                                               lemon::ListDigraph::Node const &,
                                               build_data<seqan3::data_layout::compressed> &,
                                               build_arguments const &,
                                               size_t &);

} // namespace raptor::hibf

