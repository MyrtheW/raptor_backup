#include <raptor/search/search_single.hpp>
#include <raptor/update/load_hibf.hpp>

#include <raptor/build/store_index.hpp>

//


namespace raptor
{

/*!\brief Calls the layout algorithm from the chopper library
 * \param[in]
 * \param[in] update_arguments configuration object with parameters required for calling an update operation
 * \author Myrthe
 */
void load_hibf(raptor_index<index_structure::hibf> & index, update_arguments const & arguments) // perhaps better to have index as in and output of this function, because now it is calling the update function within.
{   if (not arguments.compressed){
    double index_io_time{0.0};
    load_index(index, arguments, index_io_time); // add: arguments.parts - 1, if HIBF consists of multiple
    // prepare index to allow for updates
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin();
    index.ibf().initialize_previous_ibf_id(); //
    index.ibf().initialize_ibf_sizes();
    }else{ // try decompressing a IBF
        auto some_compressed_ibf = index.ibf().ibf_vector[0];
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_uncompressed{some_compressed_ibf};
    }
    return;
}

} // end namespace

