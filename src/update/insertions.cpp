#include <raptor/update/insertions.hpp>
#include <raptor/update/rebuild.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <numeric>
#include <chopper/count/count_kmers.hpp>
#include <chopper/sketch/hyperloglog.hpp>
#include <chopper/count/output.hpp>
#include <chopper/sketch/user_bin_sequence.hpp>

namespace raptor
{

/*!\brief Finds a location in the HIBF for a new UB.
 * \details The algorithm finds a suitable location in the existing HIBF for a new UB, by .
 * 1) looking for a proper IBF with one of the three algorithms
 * 2) calculating the number of TBs needed to insert the UB in this IBF
 * 3) finding the bin_idx in the
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \return ibf_idx, bin_idx, number_of_bins An index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \author Myrthe
 */
std::tuple <uint64_t, uint64_t, uint16_t> get_location(size_t kmer_count, robin_hood::unordered_flat_set<size_t> &kmers,
                                                       raptor_index<index_structure::hibf> & index){
    size_t root_idx = 0;
    size_t ibf_idx = find_ibf_idx_traverse_by_similarity(kmers, index);
    // TODO: change to some input argument, to decide which function to use.  if function_find_ibf == "find_ibf_idx_traverse_by_similarity"
    //size_t ibf_idx = find_ibf_idx_ibf_size(kmers.size(, index); // OR OTHER FIND LOCATION ALGORITHM
    //size_t ibf_idx = find_ibf_idx_traverse_by_fpr(kmer_count, index, root_idx); // OR OTHER FIND LOCATION ALGORITHM
    size_t number_of_bins = 1; // calculate number of user bins needed.
    if (ibf_idx == 0 and index.ibf().ibf_max_kmers(ibf_idx) < (size_t) kmer_count){      // Only if we are at the root we might have to split bins. Delete ibf_idx==0, if using similarity method
        number_of_bins = index.ibf().number_of_bins(ibf_idx, kmer_count);       // calculate among how many bins we should split
    }
    uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, number_of_bins);
    return {ibf_idx, bin_idx, number_of_bins};
    }


/*!\brief Inserts a UB in the assigned TBs and its parent MBs in higher level IBFs
 * \details The algorithm inserts the UB in the leave bins (LBs) and parent merged bins
 * iteratively until reaching the root ibf at the top level. For each level the function
 * `insert_into_ibf` will be called. This function updates the `rebuild_index_tuple`, indicating
 * whether part of the index needs to be rebuild. This is initialized with the ibf_idx being the HIBF size.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index_triple (ibf_idx, bin_idx, number_of_bins) indicating where the UB needs to be inserted.
 * \param[in] index the original HIBF
 * \return rebuild_index_tuple containing (ibf_idx, bin_idx). If the ibf_idx is lower than the total number of IBFs,
 * then the respective IBF needs to be rebuild.
 * \author Myrthe Willemsen
 */
std::tuple <uint64_t, uint64_t> insert_tb_and_parents(robin_hood::unordered_flat_set<size_t> & kmers,
                                                      std::tuple <uint64_t, uint64_t, uint16_t> & index_triple,
                                                      raptor_index<index_structure::hibf> & index){
    std::tuple <uint64_t, uint64_t> rebuild_index_tuple = std::make_tuple(index.ibf().ibf_vector.size(), 0);
    while (true){
        insert_into_ibf(kmers, index_triple, index, rebuild_index_tuple); // add to HIBF. If at some point the IBF will be rebuild, then make sure that kmers stay the same. . The next step will be fine, as you only need the index of the IBF, and the IBF that is rebuild only grows. However, this approach is inefficient if it would trigger again a rebuild on the next level.
        if (std::get<0>(index_triple) == 0) break;
        auto index_tuple = index.ibf().previous_ibf_id[std::get<0>(index_triple)]; // update index triple:
        index_triple = std::make_tuple(std::get<0>(index_tuple), std::get<1>(index_tuple), 1); // number of bins will be 1 for the merged bins (i assume)
        }
    return rebuild_index_tuple;
    }

/*!\brief Insert (multiple) new UBs
 * \details The algorithm inters new UBs At last it computes and stores a sketch for the new UB and stores this.
 * The user is supposed to add the filename to the file containing the bin paths (all_bins_path) himself, if desired.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \param[out] ibf_idx, bin_idx, number_of_bins A index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \author Myrthe Willemsen
 */
void insert_ubs(update_arguments const & update_arguments,
                raptor_index<index_structure::hibf> & index){
    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        if (index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin, inserts it if it does not exist yet.
            std::cout << "The user bin ... that you want to insert does already exist. Please use the --insert-sequences option";
        }else{
            raptor::hibf::compute_kmers(kmers, update_arguments, filename); // or  std::vector<std::string> filename2 = {{filename}};
            size_t kmer_count = kmers.size();
            std::tuple <uint64_t, uint64_t, uint16_t> index_triple = get_location(kmer_count, kmers, index);       //  index_triple; bin_idx, ibf_idx, number_of_bins
            std::tuple <uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers,  index_triple, index);
            if (std::get<0>(rebuild_index_tuple) < index.ibf().ibf_vector.size()) // If the ibf_idx is lower than the total number of IBFs (the initilization value), then the respective IBF needs to be rebuild.
                partial_rebuild(rebuild_index_tuple, index, update_arguments);

            //Compute and store sketches for the new UBs
            chopper::configuration layout_arguments = layout_config(index, update_arguments); // This creates the arguments to be able to compute and write the sketch. The bin path filename is not relevant here.
            chopper::sketch::hyperloglog sketch(layout_arguments.sketch_bits);
            chopper::count::process_sequence_files(filename, layout_arguments, sketch);
            chopper::count::write_sketch_file(std::make_pair("_", filename) , sketch, layout_arguments);
        }
    }
}


/*!\brief Inserts sequences in existing UBs
 * \details The algorithm inserts new sequence content in the technical bins of an existing sample,
 * as well as its parent merged bins.
 * If sketches are used, then also the sketch of this UB is updated.
 * \guideline The user can provide a file specifically containing the part of the sequence that should be added, not the whole sequence.
 * By default this file ends with "_insertsequences", but can be set to a different appendix
 * If this is not possible, the user can provide the
 * I should update hyperloglog sketches, wheres the user: updates sequence files yourself., e.g. using cat.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \author Myrthe Willemsen
 */
void insert_sequences(update_arguments const & update_arguments, raptor_index<index_structure::hibf> & index){
    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize k-mers.
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        if (not index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin. Note that this function inserts it if it does not exist yet.
            std::cout << "The user bin ... that you want to insert to does not exist. If you want to add a new user bin, use the flag -insert-UB";
        }else{
                std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);

                std::string filename_new_sequences = filename[0] + update_arguments.insert_sequence_appendix;
                raptor::hibf::compute_kmers(kmers, update_arguments, std::vector{filename_new_sequences});
                std::tuple <uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers,  index_triple, index);
                if (std::get<0>(rebuild_index_tuple) < index.ibf().ibf_vector.size()) // If the ibf_idx is lower than the total number of IBFs (the initilization value), then the respective IBF needs to be rebuild.
                    partial_rebuild(rebuild_index_tuple, index, update_arguments, 2);

                if (update_arguments.sketch_directory != ""){ // if sketches are used, then update the sketch of this UB.
                    std::vector<chopper::sketch::hyperloglog> sketches;
                    chopper::sketch::user_bin_sequence::read_hll_files_into(update_arguments.sketch_directory, filename, sketches); // instead of hll_dir, use update_arguments.sketch_directory
                    chopper::configuration layout_arguments = layout_config(index, update_arguments); // Create the arguments to run the layout algorithm with.
                    chopper::count::process_sequence_files(filename, layout_arguments, sketches[0]);
                    chopper::count::write_sketch_file(std::make_pair("_", filename) , sketches[0], layout_arguments);
                }
        }
    }
}



/*!\brief Delete sequences from UBs //TODO DOC & TODO CODE
 * \details The algorithm inserts new sequence content in the technical bins of an existing sample,
 * If sketches are used, then also the sketch of this UB is updated.
 * \guideline The user can provide a file specifically containing the part of the sequence that should be added, not the whole sequence.
 * By default this file ends with "_insertsequences", but can be set to a different appendix
 * If this is not possible, the user can provide the
 * I should update hyperloglog sketches, wheres the user: updates sequence files yourself., e.g. using cat.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \author Myrthe Willemsen
 */
void delete_sequences(update_arguments const & arguments,
                  raptor_index<index_structure::hibf> & index){
    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
    for  (auto &filename: arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        if (not index.ibf().user_bins.exists_filename(filename[0])){ //            // check if filename exists in the UBs
            std::cout << "The user bin ... that you want to insert to does not exist. If you want to add a new user bin, use the flag -insert-UB";
        }else{
                std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
                //filename[0].append("_deletesequences");
                raptor::hibf::compute_kmers(kmers, arguments, filename); // or  std::vector<std::string> filename2 = {{filename}};
                //update occupancy table
                // check if rebuild is required --> rebuild the UB rather than splitting it.
                // create new sketch.
            }

    }
}

//DELETE UBS

// deleting the actual sequence file and from all_bins_path is to the user. ?
void delete_ubs(update_arguments const & update_arguments,
                  raptor_index<index_structure::hibf> & index){
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        delete_ub(filename, index); // delete a single user bin from the index.

        if (update_arguments.sketch_directory != ""){ // if sketches are used, then delete the sketch of this UB.
            std::filesystem::path path = update_arguments.sketch_directory / std::filesystem::path(filename[0]).stem();
            path += ".hll";
            std::filesystem::remove(path);
        }
    }
}

/*!\brief Delete a single user bin from the index.
 * \param[in] filename filename of the user bin to be removed.
 * \param[in] index The HIBF.
 * \details Delete a single UB by deleting one ore more leaf bins,
 * with the help of the the `delete_tbs` function.
 * \author Myrthe Willemsen
 */
void delete_ub(std::vector<std::string> const & filename,
                    raptor_index<index_structure::hibf> & index){

    if (not index.ibf().user_bins.exists_filename(filename[0])) // first find index
    {
        std::cout << "Warning: the user bin you want to delete is not present in the HIBF: "; // + filename; //--> if not: return error , make sure by doing this, we dont accidently add the filename..
    }else{
        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
        size_t const ibf_idx = std::get<0>(index_triple); // Create an empty UB
        size_t const start_bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        index.ibf().delete_tbs(ibf_idx, start_bin_idx, number_of_bins);
    }
    index.ibf().user_bins.delete_filename(filename[0]);  // Update filename tables. Even if the UB did not exist, it might have been added through the STL's .find() function.
}





/*!\brief Finds an appropiate IBF for a UB insertion based on FPR.
 * \param[in] kmer_count The number of k-mers to be inserted.
 * \param[in] ibf_idx The IBF index of the current IBF. Should be the root (0) at the start of the traversal.
 * \param[in] index The HIBF.
 * \return the IBF index which is best appropiate based on the FPR, using the FPR table.
 * \details the merged bin with the lowest false positive rate is selected.
 * The algorithm recurses until an IBF's bin size is too small to accommodate the new UB.
 * \author Myrthe Willemsen
 */
size_t find_ibf_idx_traverse_by_fpr(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx =0){ //default is root_idx =0, where we start the search.
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF
    if (index.ibf().ibf_max_kmers(ibf_idx) > kmer_count){ // kmer-capacity of IBF > bin size new UB, go down if possible. Instead of maximal capcity, you can calculate the optimal kmer_ size.
        size_t best_mb_idx = ibf.bin_count(); double best_fpr = 1; // initialize the best idx outside of the ibf, such that we can use this after the loop to check if a MB was found.
         for (size_t bin_idx=0; bin_idx < ibf.bin_count(); ++bin_idx){ //loop over bins to find the bext merged bin
            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
                auto fpr = index.ibf().get_fpr(ibf_idx, bin_idx);
                if (fpr < best_fpr){
                    best_fpr = fpr;
                    best_mb_idx = bin_idx;
                }
            }
         }
         if (best_mb_idx >= ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
             return ibf_idx;
         }else{
             auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //next_ibf_id[ibf_id_high].size()
             return (find_ibf_idx_traverse_by_fpr(kmer_count, index, next_ibf_idx));
         }
    }else{ // kmer-capacity of IBF < bin size new UB, go up if possible
        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
    }
    }

    //!\brief Takes a number and calculates the next multiple of 64
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value) noexcept
{
    return ((value + 63) >> 6) << 6;
}

    /*!\brief Finds empty bins (EBs) within a certain IBF, where the new UB can be inserted.
     * \param[in] ibf_idx the IBF in which EBs need to be found.
     * \param[in] number_of_bins the number of bins needed to store the UB in question in this particular IBF.
     * \param[in] index The HIBF.
     * \return The starting index of the empty TBs where the new bin can be inserted.
     * \details Try inserting in the first EB encountered. Check if there are sufficient adjacent empty bins.
     * Note: this could be improved using empty bin data structure and rank operation.
     * If there is no empty bin, the IBF must be resized.
     * \author Myrthe Willemsen
     */
uint64_t find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t ibf_idx, size_t number_of_bins){
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    size_t bin_idx=0; // The variable is initialized outside the for loop, such that afterwards it can still be used.
    for (; bin_idx + number_of_bins < ibf_bin_count; bin_idx++){ // This could be implemented more efficiently.
        if (std::reduce(&index.ibf().occupancy_table[ibf_idx][bin_idx],
                        &index.ibf().occupancy_table[ibf_idx][bin_idx+number_of_bins])==0){ //empty bin
            break;
        }
    } // If the 'break' has not been triggered, i.e. bin_idx = ibf_bin_count, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    if (bin_idx == ibf_bin_count){  // then the IBF must be resized.
        double EB_percentage = 0.1;
        size_t new_ibf_bin_count = next_multiple_of_64(std::max((size_t) std::round((1+EB_percentage)*ibf_bin_count), ibf_bin_count + number_of_bins));
        index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
    }
    return bin_idx;
}

    /*!\brief Finds an appropiate IBF for a UB insertion based on the number of k-mers and IBF sizes.
     * \param[in] kmer_count The number of k-mers to be inserted.
     * \param[in] index The HIBF.
     * \return the IBF index to insert the new UB in.
     * \details A second algorithm picks the IBF that has the smallest size able to store the UB without it being split.
     * It finds this IBF by doing a binary search on a sorted array with IBF bin sizes in terms of the number of
     * $k$-mers that they can maximally store, and the corresponding IBF indexes. The search is logarithmic in the time of the number of IBFs,
     * \author Myrthe Willemsen
     */
size_t find_ibf_idx_ibf_size(size_t & kmer_count, raptor_index<index_structure::hibf> & index){
    auto & array = index.ibf().ibf_sizes;
        int low = 0;
        int high = array.size()-1;
        while (low <= high) {
            int mid = (low + high) >> 1;
            if (std::get<0>(array[mid]) < kmer_count)
                {low = mid + 1;}
            else if (std::get<0>(array[mid]) > kmer_count)
                {high = mid - 1;}
            else if (std::get<0>(array[mid]) == kmer_count)
                {return std::get<1>(array[mid]);} // exact kmer_count found
        }
        return std::get<1>(array[low]);
    }


    /*!\brief Finds an appropiate IBF for a UB insertion based on sequence similarity
     * \param[in] kmers The kmers to be inserted.
     * \param[in] ibf_idx The IBF index of the current IBF. Should be the root (0) at the start of the traversal.
     * \param[in] index The HIBF.
     * \return ibf_idx of the IBF to insert the new UB in.
     * \details the merged bin with the maximal similarity is selected.
     * The algorithm recurses until an IBF's bin size is too small to accommodate the new UB.
     * The similarities to all TBs in each traversed IBF are approximated by querying a percentage of the k-mers from the new UBs.
     * This can be efficiently done using the query me
     * thod that is already in place.
     * \author Myrthe Willemsen
     */
size_t find_ibf_idx_traverse_by_similarity(robin_hood::unordered_flat_set<size_t> & kmers, raptor_index<index_structure::hibf> & index, size_t ibf_idx){ //default is root_idx =0, where we start the search.
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //select the IBF
    if (index.ibf().ibf_max_kmers(ibf_idx) > kmers.size()){ // kmer-capacity of IBF > bin size new UB, go down if possible. Instead of maximal capcity, you can calculate the optimal kmer_ size.
        auto agent = ibf.template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(kmers); // count occurrences of the kmers in each of the bins in the current IBF.
        size_t best_mb_idx = ibf.bin_count(); int best_similarity = -1; // initialize the best index outside of the ibf, such that we can use this after the loop to check if a MB was found.
        for (size_t bin_idx=0; bin_idx < ibf.bin_count(); ++bin_idx){ // loop over bins to find the next merged bin
            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
                auto similarity = result[bin_idx];
                if (similarity > best_similarity){
                    best_similarity = similarity;
                    best_mb_idx = bin_idx;
                }
            }
        }
        if (best_mb_idx >= ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
            return ibf_idx;
        }else{
            auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //next_ibf_id[ibf_id_high].size()
            return (find_ibf_idx_traverse_by_similarity(kmers, index, next_ibf_idx));
        }
    }else{ // if the kmer-capacity of the IBF is lower than the required bin size for the new UB, then go up if possible
        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
    }
}






// TODO : use update_filename_indices?

} // end namespace
