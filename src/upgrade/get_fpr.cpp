//#include <raptor/upgrade/get_fpr.hpp>
//#include <raptor/build/hibf/compute_kmers.hpp>
//#include <raptor/build/hibf/insert_into_ibf.hpp>
////#include <algorithm> // for std::max
//// THIS FILE IS OLD!!!
//// CAN BE DELETED
//
//
//namespace raptor
//{
//   // using index_structure_t = std::conditional_t<seqan3::compressed, index_structure::hibf_compressed, index_structure::hibf>;
////template<bool compressed>
//std::tuple <uint64_t, uint64_t, uint16_t> get_location(std::vector<std::string> const & filename,
//                                                       size_t kmer_count,
//                                                       raptor_index<index_structure::hibf> & index){
//    // test:  only return {0,0,0} see if the issue still occurs without calling all functions.
////bool x = index.ibf().user_bins.exists_filename(filename);
////using user_bins = index.ibf().user_bins;
//if (index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin, inserts it if it does not exist yet.
//    return (index.ibf().user_bins.find_filename(filename[0]));
//    //  or use request_ibf_idx and request_user_bin_idx
//    // check if the UB still has enough space for the extra inserts.
//}else{
//    // FIND LOCATION FOR NEW UB
//    // one possible outcome might be that a new layout is made.this now happens in insert_into_ibf.
//    //METHOD 1
//    size_t root_idx = 0; // would assume so based on bulk_contains_impl
//    size_t ibf_idx = find_ibf_idx_traverse_by_fpr(kmer_count, index, root_idx);
//    std::tuple <uint64_t, uint64_t> index_tuple = find_empty_bin_idx(kmer_count, index, ibf_idx);
//    std::tuple <uint64_t, uint64_t, uint16_t> index_triple = {ibf_idx, std::get<0>(index_tuple), std::get<1>(index_tuple)};             // dummy: std::tuple <uint64_t, uint64_t, uint16_t> index_triple = {0,0,0} ; //or index_pair_list
//
//    return index_triple;
//    }
//}
//
//void upgrade_hibf(upgrade_arguments const & arguments,
//                  raptor_index<index_structure::hibf> & index)   //std::move is not correct to use here. https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used
//{
//    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
//    for  (auto &filename: arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
//            raptor::hibf::compute_kmers(kmers, arguments, filename); // or  std::vector<std::string> filename2 = {{filename}};
//            size_t kmer_count = kmers.size();
//            std::tuple <uint64_t, uint64_t, uint16_t> index_triple = get_location(filename, kmer_count, index);       //  index_triple; bin_idx, ibf_idx, number_of_bins
//            while (std::get<0>(index_triple) != index.ibf().ibf_vector.size()){             // loop over parents to also insert the kmers in merged bins.
//                insert_into_ibf(kmers, index_triple, index); // add to HIBF. If at some point the IBF will be rebuild, then make sure that kmers stay the same. . The next step will be fine, as you only need the index of the IBF, and the IBF that is rebuild only grows. However, this approach is inefficient if it would trigger again a rebuild on the next level.
//                auto index_tuple = index.ibf().previous_ibf_id[std::get<0>(index_triple)]; // update index triple:
//                index_triple = std::make_tuple(std::get<0>(index_tuple), std::get<1>(index_tuple), 1); // number of bins will be 1 for the merged bins (i assume)
//            }
//            kmers.clear();
//    }
//}
//
//
//
////DELETE UB
//void delete_ub(std::vector<std::string> const & filename,
//                    raptor_index<index_structure::hibf> & index){
//
//    if (not index.ibf().user_bins.exists_filename(filename[0])) // first find index
//    {
//        std::cout << "Warning: the user bin you want to delete is not present in the HIBF: "; // + filename; //--> if not: return error , make sure by doing this, we dont accidently add the filename..
//    }else{
//        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
//
//        // empty UB
//        size_t const ibf_idx = std::get<0>(index_triple);
//        size_t const start_bin_idx = std::get<1>(index_triple);
//        size_t const number_of_bins = std::get<2>(index_triple);
//        auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF
//
//        for (size_t chunk_number=0; chunk_number<number_of_bins; ++chunk_number)
//        {
//            seqan3::bin_index const bin_idx{start_bin_idx + chunk_number}; //auto const bin_index = seqan3::bin_index{static_cast<size_t>(bin_idx)}; //  seqan3::bin_index const bin_idx{bin
//            ibf.clear(bin_idx);
//        }
//
//         for (size_t offset=0; offset < number_of_bins; offset++){ // update FPR table and occupancy=#kmer table.
//            index.ibf().fpr_table[ibf_idx][start_bin_idx+offset] = 0;
//            index.ibf().occupancy_table[ibf_idx][start_bin_idx+offset] = 0;
//         }
//    }
//
//    index.ibf().user_bins.delete_filename(filename[0]);  // update filename tables. even if the UB did not exist, it might have been added through the STL .find() function.
//}
//
//
////TRAVERSE HIBF
//size_t find_ibf_idx_traverse_by_fpr(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx){
//    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF , or data.hibf.ibf_vector[]
//    if (index.ibf().ibf_max_kmers(ibf_idx) > kmer_count){ // kmer-capacity of IBF > bin size new UB, go down if possible. Instead of maximal capcity, you can calculate the optimal kmer_ size.
//        size_t best_mb_idx = ibf.bin_count(); double best_fpr = 1; // initialize the best idx outside of the ibf, such that we can use this after the loop to check if a MB was found.
//         for (size_t bin_idx; bin_idx < ibf.bin_count(); ++bin_idx){ //loop over bins to find the bext merged bin
//            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
//                auto fpr = index.ibf().get_fpr(ibf_idx, bin_idx);
//                if (fpr < best_fpr){
//                    best_fpr = fpr;
//                    best_mb_idx = bin_idx;
//                }
//            }
//         }
//         if (best_mb_idx > ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
//             return ibf_idx;
//         }else{
//             auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //next_ibf_id[ibf_id_high].size()
//             return (find_ibf_idx_traverse_by_fpr(kmer_count, index, next_ibf_idx));
//         }
//    }else{ // kmer-capacity of IBF < bin size new UB, go up if possible
//        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
//    }
//    }
//
//std::tuple <uint64_t, uint64_t>  find_empty_bin_idx(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx){
//    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
//
//    size_t number_of_bins = 1; // calculate number of user bins needed.
//    if (ibf_idx==0 and index.ibf().ibf_max_kmers(ibf_idx) < kmer_count){ // if we are at the root we might have to split bins.
//              number_of_bins = index.ibf().number_of_bins(ibf_idx, kmer_count);     // calculate among how many bins we should split
//    }
//    // insert in the first EB encountered. improve using empty bin datastructure and rank operation.
//    size_t bin_idx=0;
//    for (; bin_idx < ibf_bin_count; bin_idx++){
//        if (index.ibf().occupancy_table[ibf_idx][bin_idx] == 0){
//            if (index.ibf().occupancy_table[ibf_idx][bin_idx + number_of_bins - 1] == 0){ // you can speed this up for number_of_bins ==1
//                break;
//            }else{
//               bin_idx += number_of_bins - 1;
//            }
//        }
//    } // if the whole loop is run through, no appropiate empty bin is found and the bin idx will be the size of the IBF.
//    if (bin_idx == ibf_bin_count){// then there is no empty bin. Resize IBF .  or ibf.bin_count()
//        double EB_percentage = 0.1;
//        size_t new_ibf_bin_count = std::max((size_t) std::round(EB_percentage*ibf_bin_count), ibf_bin_count + number_of_bins);
//        index.ibf().ibf_vector[ibf_idx].increase_bin_number_to((seqan3::bin_count) new_ibf_bin_count);
//    }
//    return std::make_tuple(bin_idx, number_of_bins); //index_tuple
//}
//
//} // end namespace
//
//
//
//
////NOTES
//// compute_fp_correction and bin_size_in_bits from ibf_fpr can help with this. Note that ibf_fpr seems not to be added as a target yet.
//
////UPDATING
//                    //Myrthe Question 14.10: what does this update exactely?
//            //        data.hibf.ibf_vector[ibf_pos] = std::move(ibf);
//            //        data.hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
//            //        data.hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);
//            //        update_user_bins(data, filename_indices, record); //Myrthe Question 14.10: what does this update exactely? I think it does too much, like remaking fiename_indices.
//                        // counting :
//                            // get bin size: size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored);
//                            // cardinality; //!< The size/weight of the bin (either a kmer count or hll sketch estimation) , e.g. using inline void count_kmers
//
//
//
//
//
//    // function that recalls layout
//    // function that recalls hibf construction
//        // new params:
//        // seq_inserts
//        // ub_inserts
//
////    // for new user bins: this file should contain all sequences
////    // for existing user bins: (does the provided file also contain sequences that were added? Then we should not add them again)
