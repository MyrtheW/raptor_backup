// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/build_data.hpp>
#include <raptor/index.hpp>
#include <chopper/configuration.hpp>

namespace raptor
{
//void split_ibf(size_t ibf_idx,
//                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments);
chopper::configuration  layout_config(std::string subtree_bin_paths, update_arguments const & arguments);
void recall_layout_1(size_t ibf_idx,
                  raptor_index<index_structure::hibf> & index, chopper::configuration & arguments);

void recall_layout_2(size_t ibf_idx,
                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments);
void write_filenames(std::string bin_path, std::set<std::string> user_bin_filenames);
build_arguments build_config(std::string subtree_bin_paths, update_arguments const & arguments);
template <seqan3::data_layout data_layout_mode>
void recall_build_1(build_arguments & arguments);
template <typename T> void remove_indices(std::vector<T> & arg) ;
void merge_indexes(raptor_index<index_structure::hibf> & index, raptor_index<index_structure::hibf> & subindex, size_t ibf_idx);

} // namespace raptor
