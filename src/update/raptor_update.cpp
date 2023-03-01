// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/update/update.hpp>
#include <raptor/update/load_hibf.hpp>
#include "raptor/index.hpp"
#include "raptor/update/rebuild.hpp"
#include <raptor/build/store_index.hpp>
#include <raptor/update/insertions.hpp>

namespace raptor
{
/*!\brief Prunes subtree from the original HIBF
 * \details One should remove the IBFs in the original index which were part of the subtree that had to be rebuild.
 * If using some sort of splitting, then removing only needs to happen once since both new subindexes share the same original ibfs.
 * \param[in|out] index the original HIBF
 * \param[in] ibf_idx the index of the IBF where the subtree needs to be removed, including the ibf_idx itself.
 * \author Myrthe Willemsen
 */
void raptor_update(update_arguments const & arguments)
{    if (arguments.is_hibf) // and not arguments.compressed remove arguments.compressed later, when decompressing is effective.
    {
        auto index = raptor_index<index_structure::hibf>{}; // Does not do anything with arguments? Strangely seems only done in store_index.
        load_hibf(index, arguments);

        std::tuple<size_t,size_t> index_tuple = std::make_tuple(0, 14);
        partial_rebuild(index_tuple, index, arguments);
        // TODO after rebuilding, the '            assert(current_filename_index < 0);' becomes false in hibf.hpp
        if //constexpr
        (not arguments.compressed){ // should be constexpr, otherwise it will try for all vlaues of compressed
            if (arguments.insert_ubs==true){
                insert_ubs(arguments, index); // currently requires uncompressed type?  requires (compressed == false)
            }else if(arguments.insert_sequences==true){
                insert_sequences(arguments, index);
            }else if(arguments.delete_ubs==true){
                delete_ubs(arguments, index);
            }else if(arguments.delete_sequences==true){

            }else{
                std::cout << "please select one of the options {delete-ubs, insert-ubs, insert-sequences, delete-sequences}";
            }

        }
    // if arguments.compressed, it should be compressed again
    store_index(arguments.in_file, index, arguments);
    //store_index(arguments.in_file, std::move(index), arguments); // store index
    //todo: there is a bug AFTER storing the index. Ask Svenja. free(): invalid pointer, as part of ~user_bins. In chopper build, the given arguments are similar.
    // problem is specifically in destructing filename_position_to_ibf_bin of the user bin class, die ik heb gemaakt .
    }
       return;
}

} // namespace raptor
