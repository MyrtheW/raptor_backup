// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>
#include <raptor/adjust_seed.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/dna4_traits.hpp>
#include "raptor/index.hpp"

namespace raptor::hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index, // Bin_index: xth bin in an IBF
                     seqan3::interleaved_bloom_filter<> & ibf,
                     bool is_root)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            ibf.emplace(value, bin_idx);
            if (!is_root)
                parent_kmers.insert(value);
        }
    }
}

template <typename arguments_t> //Myrthe 14.10
void insert_into_ibf(arguments_t const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf)
{
    auto const bin_index = seqan3::bin_index{static_cast<size_t>(record.bin_indices.back())};

    if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : record.filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                ibf.emplace(minimiser_value, bin_index);
        }
    }
    else
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                       seqan3::window_size{arguments.window_size},
                                                       seqan3::seed{adjust_seed(arguments.shape.count())});

        for (auto const & filename : record.filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash : seq | hash_view)
                    ibf.emplace(hash, bin_index);
    }
}
template void insert_into_ibf<build_arguments>(build_arguments const & arguments, //Myrthe 14.10
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf);
template void insert_into_ibf<upgrade_arguments>(upgrade_arguments const & arguments, //Myrthe 14.10
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf);




} // namespace raptor::hibf

namespace raptor{
    // only if IBF is uncompressed.
     //Myrthe
//     using index_structure_t = std::conditional_t<seqan3::uncompressed, index_structure::hibf_compressed, index_structure::hibf>;

void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers, // kmers or minimizers
                    std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                    raptor_index<index_structure::hibf> & index)
{ // change this, so that if you have multiple technical bins of different sizes, you can create the chunks in accordance to these sizes.

    size_t const ibf_idx = std::get<0>(index_triple);
    size_t const start_bin_idx = std::get<1>(index_triple);
    size_t const number_of_bins = std::get<2>(index_triple);
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF , or data.hibf.ibf_vector[]
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{start_bin_idx + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            auto const bin_index = seqan3::bin_index{static_cast<size_t>(bin_idx)}; //  seqan3::bin_index const bin_idx{bin
            ibf.emplace(value, bin_index);
        }
    }
}

}