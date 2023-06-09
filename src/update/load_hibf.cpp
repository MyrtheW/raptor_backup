#include <raptor/search/search_single.hpp>
#include <raptor/update/load_hibf.hpp>
#include <raptor/update/insertions.hpp>
#include <raptor/build/store_index.hpp>

//


namespace raptor
{
template <bool compressed>
void load_hibf(update_arguments const & arguments) // perhaps better to have index as in and output of this function, because now it is calling the update function within.
{
    using index_structure_t = std::conditional_t<compressed, index_structure::hibf_compressed, index_structure::hibf>;
    auto index = raptor_index<index_structure_t>{}; // Does not do anything with arguments? Strangely seems only done in store_index.
    double index_io_time{0.0};
    load_index(index, arguments, index_io_time); // add: arguments.parts - 1, if HIBF consists of multiple
    // prepare index to allow for updates
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin();
    index.ibf().initialize_previous_ibf_id(); //
    // occupancy table, dummy
    index.ibf().occupancy_table.resize(index.ibf().ibf_vector.size());
    index.ibf().fpr_table.resize(index.ibf().ibf_vector.size());
    for (size_t ibf_idx=0; ibf_idx< index.ibf().ibf_vector.size(); ibf_idx++){
            index.ibf().occupancy_table[ibf_idx].resize(index.ibf().ibf_vector[ibf_idx].bin_count());
            index.ibf().fpr_table[ibf_idx].resize(index.ibf().ibf_vector[ibf_idx].bin_count());
                for (size_t bin_idx=0; bin_idx< index.ibf().ibf_vector[ibf_idx].bin_count(); bin_idx++){
                    index.ibf().fpr_table[ibf_idx][bin_idx] = 0.5;
                     index.ibf().occupancy_table[ibf_idx][bin_idx] = 10;
                }

    }
                //index.ibf().occupancy_table[0].resize(index.ibf().ibf_vector[ibf_idx].bin_count());

    // create datastructure for perfect ns.

    // make sure that the necessary data structures are loaded/created.

    if constexpr (not compressed){ // should be constexpr, otherwise it will try for all vlaues of compressed
        if (arguments.insert_ubs==true){
            update_hibf(arguments, index); // currently requires uncompressed type?  requires (compressed == false)
        }else if(arguments.insert_sequences==true){
            insert_sequences(arguments, index);
        }else if(arguments.delete_ubs==true){
            delete_ubs(arguments, index);
        }else if(arguments.delete_sequences==true){

        }else{
            std::cout << "please select one of the options {delete-ubs, insert-ubs, insert-sequences, delete-sequences}";
        }
    }
    else{ // todo try uncompressing a IBF
        auto some_compressed_ibf = index.ibf().ibf_vector[0];
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_uncompressed{some_compressed_ibf};
    }

    // if arguments.compressed, it should be compressed again
    store_index(arguments.in_file, std::move(index), arguments); // store index
}

template void load_hibf<false>(update_arguments const & arguments);
template void load_hibf<true>(update_arguments const & arguments);

}


// for each of the bins in bin_paths: (maybe sort? or do this smartly somehow?


// request_ibf_idx
// request_user_bin_idx

// when using such an existing function, make sure the arguments have the right arguments type, by using a template (see store_index.hpp how to) , typename arguments_t