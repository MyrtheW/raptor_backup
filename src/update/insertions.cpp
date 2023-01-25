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

/*!\brief Finds a new location
 * \details The algorithm //TODO
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \param[out] ibf_idx, bin_idx, number_of_bins A index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \author Myrthe
 */
std::tuple <uint64_t, uint64_t, uint16_t> get_location(size_t kmer_count, robin_hood::unordered_flat_set<size_t> &kmers,
                                                       raptor_index<index_structure::hibf> & index){
    size_t root_idx = 0;
    size_t ibf_idx = find_ibf_idx_traverse_by_similarity(kmers, index);
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
 * \details The algorithm //TODO
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index_triple
 * \param[in] index the original HIBF
 * \param[out] rebuild_index_tuple
 * \author Myrthe
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

/*!\brief Finds a new location
 * \details The algorithm //TODO
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \param[out] ibf_idx, bin_idx, number_of_bins A index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \author Myrthe
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
            if (std::get<0>(rebuild_index_tuple) == index.ibf().ibf_vector.size()) partial_rebuild(rebuild_index_tuple, index, update_arguments);

            //add sketches
            chopper::configuration layout_arguments = layout_config("all_bin_paths.txt", index, update_arguments); // some filename.. // create the arguments to run the layout algorithm with.
            chopper::sketch::hyperloglog sketch(layout_arguments.sketch_bits);
//                if (config.precomputed_files)
//                    process_minimizer_files(cluster_vector[i].second, sketch);
//                else
            chopper::count::process_sequence_files(filename, layout_arguments, sketch);
            // if (!config.disable_sketch_output)
            chopper::count::write_sketch_file(std::make_pair("_", filename) , sketch, layout_arguments);
//            std::filesystem::path path = layout_arguments.sketch_directory /
//                    std::filesystem::u8path(filename[0]).stem();
//            path += ".hll";
//            std::ofstream hll_fout(path, std::ios::binary);
//            sketch.dump(hll_fout);
        }
    }
}


/*!\brief Inserts sequences in existing UBs
 * \details The algorithm //TODO
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index_triple
 * \param[in] index the original HIBF
 * \param[out] rebuild_index_tuple
 * \author Myrthe
 */
void insert_sequences(update_arguments const & update_arguments, raptor_index<index_structure::hibf> & index){
    robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
            // check if filename exists in the UBs
        if (not index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin, inserts it if it does not exist yet.
            std::cout << "The user bin ... that you want to insert to does not exist. If you want to add a new user bin, use the flag -insert-UB";
        }else{
                std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
                // check if filename_insertsequences exits as a file using validate function?
//                const std::basic_string<char> appendix2= "_insertsequences";
//                const std::basic_string<char>& appendix = "_insertsequences";
//                // const basic_string& __str
//                const std::string & filename2 = (const std::string &) filename[0];
//                filename2.append(filename2);
//                "test".append("test");
//                appendix.append("hi");
//                (filename[0]).append((std::basic_string<char>)reinterpret_cast<const char *>('_insertsequences'));
//                (filename[0]).string("_insertsequences");
//                  filename[0].append(appendix);
//https://cplusplus.com/reference/string/basic_string/append/
// Question: how do I do this?

                //filename[0].append("_insertsequences");
                // todo ask Svenja on files.

                raptor::hibf::compute_kmers(kmers, update_arguments, filename); // or  std::vector<std::string> filename2 = {{filename}};
                std::tuple <uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers,  index_triple, index);
                if (std::get<0>(rebuild_index_tuple) == index.ibf().ibf_vector.size())
                    partial_rebuild(rebuild_index_tuple, index, update_arguments);

                // Sketches:
                chopper::configuration layout_arguments = layout_config("all_bin_paths.txt", index, update_arguments); // some filename.. // create the arguments to run the layout algorithm with.
                chopper::sketch::hyperloglog sketch(layout_arguments.sketch_bits);
                // read sketch of original file
                std::vector<size_t> test{kmers.size()}; // todo why doesn't the construction below work?
                auto sketch_toolbox = chopper::sketch::user_bin_sequence{(std::vector<std::string>) filename, test}; //, std::vector{(size_t) kmers.size()}};
                sketch_toolbox.read_hll_files(layout_arguments.sketch_directory);}
                auto sketch  = sketch_toolbox.sketches[0];
                chopper::count::process_sequence_files(filename, layout_arguments, sketch);
                chopper::count::write_sketch_file(std::make_pair("_", filename) , sketch, layout_arguments);
            }

    }
}


//DELETE SEQUENCES
void delete_sequences(update_arguments const & arguments,
                  raptor_index<index_structure::hibf> & index)   //std::move is not correct to use here. https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used
{
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
            }

    }
}

//DELETE UBS
void delete_ubs(update_arguments const & arguments,
                  raptor_index<index_structure::hibf> & index)   //std::move is not correct to use here. https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used
{
    for  (auto &filename: arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
            delete_ub(filename, index);
    } // todo delete sketch
}


void delete_ub(std::vector<std::string> const & filename,
                    raptor_index<index_structure::hibf> & index){

    if (not index.ibf().user_bins.exists_filename(filename[0])) // first find index
    {
        std::cout << "Warning: the user bin you want to delete is not present in the HIBF: "; // + filename; //--> if not: return error , make sure by doing this, we dont accidently add the filename..
    }else{
        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
        //todo: delete hyperloglog sketch. Possibly delete filename from counts. Also detele from bin_paths.
        size_t const ibf_idx = std::get<0>(index_triple); // Create an empty UB
        size_t const start_bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        index.ibf().delete_tbs(ibf_idx, start_bin_idx, number_of_bins);
    }
    index.ibf().user_bins.delete_filename(filename[0]);  // Update filename tables. Even if the UB did not exist, it might have been added through the STL .find() function.
}





//TRAVERSE HIBF
size_t find_ibf_idx_traverse_by_fpr(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx =0){ //default is root_idx =0, where we start the search.
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF , or data.hibf.ibf_vector[] auto test = index.ibf().ibf_max_kmers(ibf_idx); //todo remove
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

[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value) noexcept
{
    return ((value + 63) >> 6) << 6;
}

    /*!\brief Finds empty TBs within a certain IBF, where the new UB can be inserted.
     * \param[in] ibf_idx
     * \param[in] number_of_bins
     * \param[in] index The HIBF.
     * \return The starting index of the empty TBs where the new bin can be inserted.
     * \details
     * // try inserting in the first EB encountered. improve using empty bin datastructure and rank operation.
     * // then there is no empty bin, the IBF must be resized.  TODO
     * \author Myrthe Willemsen
     */
uint64_t find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t ibf_idx, size_t number_of_bins){
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    size_t bin_idx=0; // The variable is initialized outside the for loop, such that afterwards it can still be used.
    for (; bin_idx < ibf_bin_count; bin_idx++){ // This could be implemented  more efficiently.
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
    return bin_idx; //index_tuple
}

    /*!\brief Finds an appropiate IBF for a UB insertion based on the number of k-mers and IBF sizes.
     * \param[in] kmer_count The number of k-mers to be inserted.
     * \param[in] index The HIBF.
     * \return true if value already existed in the bin, i.e. at all hashes the value in de bin was set to 1 already.
     * \details A second algorithm picks the IBF that has the smallest size able to store the UB without it being split.
     * It finds this IBF by doing a binary search on a sorted array with IBF bin sizes in terms of the number of
     * $k$-mers that they can maximally store. The search is logarithmic in the time of the number of IBFs,
     * \author Myrthe Willemsen
     */
size_t find_ibf_idx_ibf_size(size_t & kmer_count, raptor_index<index_structure::hibf> & index){
    auto & array = index.ibf().ibf_sizes;
        int low = 0;//array[0];
        int high = array.size()-1;// array.size()-1];
        while (low <= high) {
            int mid = (low + high) >> 1;
            if (array[mid] < kmer_count)
                {low = mid + 1;}
            else if (array[mid] > kmer_count)
                {high = mid - 1;}
            else if (array[mid] == kmer_count)
                {return mid;} // key found
        }
        return low;  // key not found. or low +1 ?
    }


    /*!\brief Finds an appropiate IBF for a UB insertion based on sequence similarity
     * \param[in] kmers The kmers to be inserted.
     * \param[in] ibf_idx The IBF index of the current IBF. Should be the root (0) at the start of the traversal.
     * \param[in] index The HIBF.
     * \return true if value already existed in the bin, i.e. at all hashes the value in de bin was set to 1 already.
     * \details the merged bin with the maximal similarity is selected.
     * The algorithm recurses until an IBF's bin size is too small to accommodate the new UB.
     * The similarities to all TBs in each traversed IBF are approximated by querying a percentage of the k-mers from the new UBs.
     * This can be efficiently done using the query method that is already in place.
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
    }else{ // kmer-capacity of IBF < bin size new UB, go up if possible
        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
    }
}








} // end namespace
