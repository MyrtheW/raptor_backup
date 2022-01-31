// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/search/detail/destroyed_indirectly_by_error.hpp>

TEST(issue, 142)
{
    std::vector<double> const destroyed_by_error =
        raptor::detail::destroyed_indirectly_by_error(15u, 6u, seqan3::ungapped{4u});
}
