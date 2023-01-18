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

void raptor_update(update_arguments const & arguments)
{    if (arguments.is_hibf) // and not arguments.compressed remove arguments.compressed later, when decompressing is effective.
    {
        auto index = raptor_index<index_structure::hibf>{}; // Does not do anything with arguments? Strangely seems only done in store_index.
        load_hibf(index, arguments);

        partial_rebuild(0, index, arguments);

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
    //todo: there is a bug in storing the index. Ask Svenja.
    }
       return;
}

} // namespace raptor
