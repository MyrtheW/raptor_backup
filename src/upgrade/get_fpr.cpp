#include <raptor/upgrade/get_fpr.hpp>
//#include <raptor/search/search_single.hpp>  // test: to get some functionalities
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>

namespace raptor
{
   // using index_structure_t = std::conditional_t<seqan3::compressed, index_structure::hibf_compressed, index_structure::hibf>;
//template<bool compressed>
std::tuple <uint64_t, uint64_t, uint16_t> get_location(std::vector<std::string> const & filename,
                                                       size_t kmer_count,
                                                       raptor_index<index_structure::hibf> & index){
    // test:  only return {0,0,0} see if the issue still occurs without calling all functions.
//bool x = index.ibf().user_bins.exists_filename(filename);
//using user_bins = index.ibf().user_bins;
if (index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin, inserts it if it does not exist yet.
    return (index.ibf().user_bins.find_filename(filename[0]));
    //  or use request_ibf_idx and request_user_bin_idx
    // check if the UB still has enough space for the extra inserts.
}else{
    // FIND LOCATION FOR NEW UB
    // one possible outcome should be that a new layout is made.

    // no resizing of (H)IBF
    std::tuple <uint64_t, uint64_t, uint16_t> index_triple = {0,0,0} ; //or index_pair_list
//    ... find best location
//    // resizing of (H)IBF
//    ... rebuild, create new empty bins and directly insert new bins,
//    // update
//    update_filename_indices(filename, bin_idx, ibf_idx); //,number_of_bins
//    update_empty_bin_index(); // an index pointing to the empty bins.

    return index_triple;
    }
}
//template std::tuple <uint64_t, uint64_t, uint16_t> get_location<true>(std::vector<std::string> const & filename,
//                                                       size_t kmer_count,
//                                                       raptor_index<hierarchical_interleaved_bloom_filter<seqan3::data_layout::compressed>> index);
//template std::tuple <uint64_t, uint64_t, uint16_t> get_location<false>(std::vector<std::string> const & filename,
//                                                       size_t kmer_count,
//                                                       raptor_index<hierarchical_interleaved_bloom_filter<seqan3::data_layout::uncompressed>> index);
    void upgrade_hibf(upgrade_arguments const & arguments, raptor_index<index_structure::hibf> index)
//void upgrade_hibf(upgrade_arguments const & arguments, raptor_index<index_structure_t> index)
{ //how to parse pointer to HIBF in function?
    //std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    // Add new bins
    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
    //size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for  (auto &filename: arguments.bin_path){// loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
            raptor::hibf::compute_kmers(kmers, arguments, filename); // or  std::vector<std::string> filename2 = {{filename}};
//            index_structure::hibf::compute_kmers(kmers, arguments, test); // computes kmers that are stored in the yet empty "kmers" object, based on the record.file
//            index_structure::hibf::compute_kmers(kmers, arguments, std::vector<std::string>{{filename}}); // computes kmers that are stored in the yet empty "kmers" object, based on the record.file
            size_t kmer_count = kmers.size();
            std::tuple <uint64_t, uint64_t, uint16_t> index_triple = get_location(filename, kmer_count, index); //TODO: std::move is not correct to use here. https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used
            // also return number of bins or list of indices, for if bin is split. index_triple; bin_idx, ibf_idx, number_of_bins
            while (std::get<0>(index_triple) != index.ibf().ibf_vector.size()){             // loop over parents to also insert the kmers in merged bins.
                insert_into_ibf(kmers, index_triple, index); // add to HIBF.
                auto index_tuple = index.ibf().previous_ibf_id[std::get<0>(index_triple)]; // update index triple:
                index_triple = std::make_tuple(std::get<0>(index_tuple), std::get<1>(index_tuple), 1);
            }
            kmers.clear();
    }
}



//DELETE UB
void delete_ub(std::vector<std::string> const & filename,
                    raptor_index<index_structure::hibf> & index){

    if (not index.ibf().user_bins.exists_filename(filename[0])) // first find index
    {
        std::cout << "Warning: the user bin you want to delete is not present in the HIBF: "; // + filename; //--> if not: return error , make sure by doing this, we dont accidently add the filename..
    }else{
        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);

        // empty UB
        size_t const ibf_idx = std::get<0>(index_triple);
        size_t const start_bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF

        for (size_t chunk_number=0; chunk_number<number_of_bins; ++chunk_number)
        {
            seqan3::bin_index const bin_idx{start_bin_idx + chunk_number};
            auto const bin_index = seqan3::bin_index{static_cast<size_t>(bin_idx)}; //  seqan3::bin_index const bin_idx{bin
            ibf.clear(bin_index);
        }
    }

    index.ibf().user_bins.delete_filename(filename[0]);  // update filename tables. even if the UB did not exist, it might have been added through the STL .find() function.

     // update FPR table and #kmer table.
    //...

}

} // end namespace


//NOTES
// compute_fp_correction and bin_size_in_bits from ibf_fpr can help with this. Note that ibf_fpr seems not to be added as a target yet.

//UPDATING
                    //Myrthe Question 14.10: what does this update exactely?
            //        data.hibf.ibf_vector[ibf_pos] = std::move(ibf);
            //        data.hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
            //        data.hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);
            //        update_user_bins(data, filename_indices, record); //Myrthe Question 14.10: what does this update exactely? I think it does too much, like remaking fiename_indices.
                        // counting :
                            // get bin size: size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored);
                            // cardinality; //!< The size/weight of the bin (either a kmer count or hll sketch estimation) , e.g. using inline void count_kmers


//TRAVERSE HIBF: The following should go into hpp?
//    template <std::ranges::forward_range value_range_t>
//    void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
//    {
//        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>(); // search the ibf. How does this work? Can I write an equivalent to obtain FPR/occupancy?
//        auto & result = agent.bulk_count(values);
//
//        uint16_t sum{};
//
//        for (size_t bin{}; bin < result.size(); ++bin)
//        {
//            sum += result[bin];
//
//            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);
//
//            if (current_filename_index < 0) // merged bin
//            {
//                if (sum >= threshold) // if ibf contains the query.
//                    bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold); //recursive
//                sum = 0u;
//            }
//            else if (bin + 1u == result.size() ||                                                    // last bin
//                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
//            {
//                if (sum >= threshold)  // if ibf contains the query.
//                    result_buffer.emplace_back(current_filename_index);
//                sum = 0u;
//            }
//        }
//    }



    // function that recalls layout
    // function that recalls hibf construction
        // new params:
        // seq_inserts
        // ub_inserts

    // first function in store_index, calling index, with last argument bin_path, actually updates the IBF.
    // this we might use for filling an empty UB.
//    void insert_sequences(index, arguments){
//        // bin paths and index are already loaded in the raptor_upgrad.cpp and upgrade_parsing.cpp.
//        // processed arguments: .bin_paths.
//        // arguments: .upgrade_ubs vs .upgrade_seqs (or decide yourself by checking if exists.
//    }// arguments: HIBF, path to sequences
//    void insert_sequences_tb(){ // arguments: TB, one file path to sequences
//    // for empty user bins: this file should contain all sequences
//    // for existing user bins: (does the provided file also contain sequences that were added? Then we should not add them again)
//
//    // use: insert_into_ibf
//
//    }
