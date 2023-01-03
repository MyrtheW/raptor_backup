// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor
{
//void split_ibf(size_t ibf_idx,
//                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments);
//void recall_dp(size_t ibf_idx,
//                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments);
void recall_layout_2(size_t ibf_idx,
                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments);

} // namespace raptor
