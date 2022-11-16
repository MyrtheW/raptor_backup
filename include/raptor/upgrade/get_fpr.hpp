//
//  myrthe on 17/10/2022.

#pragma once

#include <robin_hood.h>

#include <raptor/argument_parsing/upgrade_arguments.hpp>
#include "raptor/index.hpp"

namespace raptor
{
    std::tuple <uint64_t, uint64_t, uint16_t> get_location(std::vector<std::string> const & filename,
                                                       size_t kmer_count,
                                                       raptor_index<index_structure::hibf> & index);
    void upgrade_hibf(upgrade_arguments const & arguments, raptor_index<index_structure::hibf> & index);
    void delete_ub(std::vector<std::string> const & filename, raptor_index<index_structure::hibf> & index);
    size_t find_ibf_idx_traverse_by_fpr(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx);
    std::tuple <uint64_t, uint64_t>  find_empty_bin_idx(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx);

}