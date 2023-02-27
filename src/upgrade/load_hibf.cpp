#include <raptor/search/search_single.hpp> // to make sure index_structure::hibf_compresse exists here.
#include <raptor/upgrade/load_hibf.hpp>
#include <raptor/upgrade/get_fpr.hpp>
#include "raptor/build/store_index.hpp"

// TODO this file is old  and can be deleted.
//#include upgrade file
//namespace raptor::hibf{
//    void some_function(std::string filename){
//        bool x = exists_filename(filename); //index_structure::};
//}
namespace raptor
{
template <bool compressed>
void load_hibf(upgrade_arguments const & arguments)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::hibf_compressed, index_structure::hibf>;
    auto index = raptor_index<index_structure_t>{}; // Does not do anything with arguments? Strangely seems only done in store_index.
    double index_io_time{0.0};
    load_index(index, arguments, index_io_time); // add: arguments.parts - 1, if HIBF consists of multiple
    // prepare index to allow for updates
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin();
    index.ibf().initialize_previous_ibf_id(); //
    // create datastructure for perfect ns.

    if constexpr (not compressed){ // should be constexpr, otherwise it will try for all vlaues of compressed
        upgrade_hibf(arguments, index); //
    } // requires uncompressed type?  requires (compressed == false)
    else{ // todo try uncompressing a IBF
        auto some_compressed_ibf = index.ibf().ibf_vector[0];
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_uncompressed{some_compressed_ibf};
    }

    // if arguments.compressed, it should be compressed again (?)
    store_index(arguments.in_file, std::move(index), arguments); // store index
}

template void load_hibf<false>(upgrade_arguments const & arguments);
template void load_hibf<true>(upgrade_arguments const & arguments);

}


// for each of the bins in bin_paths: (maybe sort? or do this smartly somehow?


// request_ibf_idx
// request_user_bin_idx

// when using such an existing function, make sure the arguments have the right arguments type, by using a template (see store_index.hpp how to) , typename arguments_t