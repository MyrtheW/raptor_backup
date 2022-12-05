// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/chopper_pack_record.hpp>
#include <raptor/argument_parsing/update_arguments.hpp> //Myrthe 14.10

namespace raptor::hibf
{
template <typename arguments_t> //Myrthe 14.10
void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   arguments_t const & arguments,
                   chopper_pack_record const & record);

template <typename arguments_t> //Myrthe 24.10
void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   arguments_t const & arguments,
                   std::vector<std::string> const & filenames); //remove const?
} // namespace raptor::hibf
