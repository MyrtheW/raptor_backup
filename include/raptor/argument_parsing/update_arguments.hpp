// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor
{

struct update_arguments
{  //TODO documentation
    // TODO which arguments are needed for updating, and which can be stored as part of the IBF
    //!\brief
    uint32_t window_size{};
    //!\brief
    uint8_t kmer_size{};
    //!\brief
    seqan3::shape shape{};
    //!\brief
    uint8_t parts{1u};

    //!\brief The index is compressed
    bool compressed{false};

    //!\brief The index is an HIBF
    bool is_hibf{false};

    //!\brief The percentage of empty bins sampled during layout computation.
    double empty_bin_percentage{0.1};

    //!\brief Should updates account for sequence similarities.
    bool similarity{false};

    // Following arguments are mutually exclusive and indicate what update operation should be performed
    bool delete_ubs{false};
    bool insert_ubs{false};
    bool insert_sequences{false};
    bool delete_sequences{false};

    // Filenames
    std::filesystem::path bin_file{};
    std::filesystem::path in_file{};
    std::filesystem::path out_file{};
    std::filesystem::path sketch_directory{"sketches"};
    std::string insert_sequence_appendix{"_insertsequences"};

    //!\brief
    std::vector<std::vector<std::string>> bin_path{};

    //!\brief
    bool is_minimiser{false}; //Myrthe 14.10
};

} // namespace raptor
