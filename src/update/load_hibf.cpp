#include <raptor/search/search_single.hpp>
#include <raptor/update/load_hibf.hpp>

#include <raptor/build/store_index.hpp>

//


namespace raptor
{
void load_hibf(raptor_index<index_structure::hibf> & index, update_arguments const & arguments) // perhaps better to have index as in and output of this function, because now it is calling the update function within.
{   if (not arguments.compressed){
    double index_io_time{0.0};
    load_index(index, arguments, index_io_time); // add: arguments.parts - 1, if HIBF consists of multiple
    // prepare index to allow for updates
    index.ibf().fpr_max = 0.05; //todo store as part of IBF or make part of argument.
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin();
    index.ibf().initialize_previous_ibf_id(); //
    index.ibf().ibf_sizes.resize(index.ibf().ibf_vector.size()); //todo make part of initialize_ibf_sizes
    index.ibf().initialize_ibf_sizes();
    // occupancy table, dummy
//    index.ibf().occupancy_table.resize(index.ibf().ibf_vector.size());
//    index.ibf().fpr_table.resize(index.ibf().ibf_vector.size());
//    for (size_t ibf_idx=0; ibf_idx< index.ibf().ibf_vector.size(); ibf_idx++){
//            index.ibf().occupancy_table[ibf_idx].resize(index.ibf().ibf_vector[ibf_idx].bin_count());
//            index.ibf().fpr_table[ibf_idx].resize(index.ibf().ibf_vector[ibf_idx].bin_count());
//                for (size_t bin_idx=0; bin_idx< index.ibf().ibf_vector[ibf_idx].bin_count(); bin_idx++){
//                    index.ibf().fpr_table[ibf_idx][bin_idx] = 0.5;
//                     index.ibf().occupancy_table[ibf_idx][bin_idx] = 10;
//                }
//    }
                //index.ibf().occupancy_table[0].resize(index.ibf().ibf_vector[ibf_idx].bin_count());

    // create datastructure for perfect ns.

    // make sure that the necessary data structures are loaded/created.
    }else{ // todo try decompressing a IBF
        auto some_compressed_ibf = index.ibf().ibf_vector[0];
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_uncompressed{some_compressed_ibf};
    }
    return;
}
//
//template void load_hibf<false>(raptor_index<index_structure::hibf> & index, update_arguments const & arguments);
//template void load_hibf<true>(raptor_index<index_structure::hibf_compressed> & index, update_arguments const & arguments);

}


// for each of the bins in bin_paths: (maybe sort? or do this smartly somehow?


// request_ibf_idx
// request_user_bin_idx

// when using such an existing function, make sure the arguments have the right arguments type, by using a template (see store_index.hpp how to) , typename arguments_t