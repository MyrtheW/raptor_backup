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
/*!\brief Insertions of an UB into an IBF //TODO DOC
 * \details The algorithm inserts the k-mers of a UB in one or more TBs
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \param[in] ibf_idx, bin_idx, number_of_bins A index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \param[in] ibf a  reference to the IBF to insert the
 * \author Myrthe Willemsen
 */

// alternative function that uses index as input,
template <seqan3::data_layout data_layout_mode>
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers, // kmers or minimizers
                     std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                     raptor::hierarchical_interleaved_bloom_filter<data_layout_mode> & index,
                     seqan3::interleaved_bloom_filter<> & ibf,
                     bool is_root)
{
    size_t const ibf_idx = std::get<0>(index_triple);
    size_t const start_bin_idx = std::get<1>(index_triple);
    size_t const number_of_bins = std::get<2>(index_triple);
    //auto& ibf = index.ibf_vector[ibf_idx]; //  select the IBF , or data.hibf.ibf_vector[] or auto&&
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};


    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{start_bin_idx + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            ibf.emplace(value, bin_idx);
            if (!is_root)
                parent_kmers.insert(value);
        }
    }
    index.update_occupancy_table((size_t) kmers.size(), ibf_idx, start_bin_idx, number_of_bins);
    auto fpr = index.update_fpr(ibf_idx, start_bin_idx, number_of_bins); // this should be done after updating the occupancy table.

}

template void insert_into_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> & ,
                     robin_hood::unordered_flat_set<size_t> const & , // kmers or minimizers
                     std::tuple <uint64_t, uint64_t, uint16_t> ,
                     raptor::hierarchical_interleaved_bloom_filter<seqan3::data_layout::uncompressed> & ,
                     seqan3::interleaved_bloom_filter<> & ,
                     bool );

template void insert_into_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> & ,
                     robin_hood::unordered_flat_set<size_t> const & , // kmers or minimizers
                     std::tuple <uint64_t, uint64_t, uint16_t> ,
                     raptor::hierarchical_interleaved_bloom_filter<seqan3::data_layout::compressed> &,
                              seqan3::interleaved_bloom_filter<> & ,
                     bool);

template <typename arguments_t> //Myrthe 14.10
void insert_into_ibf(arguments_t const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf) //  edge case when there is no splitting at the root and the function can be simplified (and is probably more efficient than the first, more generic one)
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



/*!\brief This version of insert_into_ibf is only used for updating.
 * \param[in] kmers The kmers to be inserted.
 * \param[in] index_triple consisting of the ibf_idx, start_bin_idx and number_of_bins.
 * \param[in] index The HIBF.
 * \param[in] rebuild_index_tuple A tuple of the index of IBF and TB that need to be rebuild.
 * If the current TB reaches the FPR_max, it will be equal to the input value (when using the layout-rebuild method)
 * \attention This function is only available for **uncompressed** Interleaved Bloom Filters.
 * \author Myrthe Willemsen
 */
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers, // kmers or minimizers
                    std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                    raptor_index<index_structure::hibf> & index,
                    std::tuple <uint64_t, uint64_t> rebuild_index_tuple) // change to std::set rebuild_indexes if using the splitting method.
                    // only if IBF is uncompressed.
{ // change this, so that if you have multiple technical bins of different sizes, you can create the chunks in accordance to these sizes.
    int union_count = 0; //You should somehow check if the bits set to 1 were already set to 1.
    size_t const ibf_idx = std::get<0>(index_triple);
    size_t const start_bin_idx = std::get<1>(index_triple);
    size_t const number_of_bins = std::get<2>(index_triple);
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF , or data.hibf.ibf_vector[]
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    for (auto chunk : kmers | seqan3::views::chunk(chunk_size)) //QUESTION: Why cut it into chunks?
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{start_bin_idx + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            auto const bin_index = seqan3::bin_index{static_cast<size_t>(bin_idx)}; //  seqan3::bin_index const bin_idx{bin
            union_count += ibf.emplace_exists(value, bin_index); // Union count value will remain 0 if inserted in an empty bin.
        }
    }

    // to improve the implementation, Perhaps do the FPR calculations for all bins to which kmers will be inserted before actually inserting.
    index.ibf().update_occupancy_table(kmers.size()-union_count, ibf_idx, start_bin_idx, number_of_bins);
    auto fpr = index.ibf().update_fpr(ibf_idx, start_bin_idx, number_of_bins); // this should be done after updating the occupancy table.
    if (fpr > index.ibf().fpr_max){
        assert(index.ibf().next_ibf_id[ibf_idx][start_bin_idx] != ibf_idx); //assert that fpr should not reach fpr max for the leaf bin ibf, because we search a location such that it is 'feasible', i.e. binsize should be sufficient to accomodate new UB.
        rebuild_index_tuple = std::make_tuple(ibf_idx, start_bin_idx);
    }
}

} // namespace raptor