// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <raptor/argument_parsing/upgrade_arguments.hpp> //Myrthe 14.10
#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/chopper_pack_record.hpp>
#include "raptor/index.hpp"

namespace raptor::hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf,
                     bool is_root);

template <typename arguments_t> //Myrthe 14.10
void insert_into_ibf(arguments_t const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf);

void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers, // kmers or minimizers
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf); //Myrthe 16.10
} // namespace raptor::hibf
namespace raptor
{
using index_structure_t = std::conditional_t<seqan3::uncompressed, index_structure::hibf_compressed, index_structure::hibf>;

void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers, // kmers or minimizers
                    std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                    raptor_index<index_structure_t> & index);  //Myrthe 20.10
} // namespace raptor