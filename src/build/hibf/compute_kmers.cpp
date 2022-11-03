// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/dna4_traits.hpp>

namespace raptor::hibf
{
template <typename arguments_t> //Myrthe 14.10
void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   arguments_t const & arguments,
                   std::vector<std::string> const & filenames)//std::basic_string<char>
{
    if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                kmers.insert(minimiser_value);
                // track a kmer_count here, Myrthe 14.10
        }
    }
    else
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        for (auto const & filename : filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape.count())}))
                    kmers.insert(hash);
    }
}
template <typename arguments_t> //Myrthe 14.10
void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   arguments_t const & arguments,
                   chopper_pack_record const & record){
    compute_kmers(kmers, arguments, record.filenames);
};

template void compute_kmers<build_arguments>(robin_hood::unordered_flat_set<size_t> & kmers, //Myrthe 14.10
                   build_arguments const & arguments,
                   chopper_pack_record const & record);
template void compute_kmers<upgrade_arguments>(robin_hood::unordered_flat_set<size_t> & kmers,
                   upgrade_arguments const & arguments,
                   chopper_pack_record const & record);
} // namespace raptor::hibf
