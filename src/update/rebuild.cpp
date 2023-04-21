//#include <raptor/search/search_single.hpp> // to make sure index_structure::hibf_compresse exists here.
#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/update/load_hibf.hpp>
#include <raptor/update/insertions.hpp>
#include <raptor/build/store_index.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/count/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <raptor/update/rebuild.hpp>

#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp> // when adding this, it recognizes std as raptor::std
#include <raptor/build/hibf/insert_into_ibf.hpp>

namespace raptor
{

// // METHOD 1
//void split_ibf(size_t ibf_idx,
//                  raptor_index<index_structure::hibf> & index, update_arguments const & arguments)
//{   int number_of_splits = 2;
//    std::vector<std::tuple <uint64_t, uint64_t>> index_tuples;  // or initialize with size {number_of_splits};
//    index_tuples[0] = index.ibf().previous_ibf_id[ibf_idx];     // get indices of the current merged bin on the higher level
//    index.ibf().delete_tbs(std::get<0>(index_tuples[0]), std::get<1>(index_tuples[0]));    // empty it
//
////    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level to serve as new merged bins.
////            index_tuples[split] = find_empty_bin_idx(index, ibf_idx);   //import from insertions.
////    }
//    // get empty bin(s) on higher level, for number of splits -1
//    std::vector<std::vector<size_t>> split_indices = find_best_split();
//    // create new IBFs
//
//        // OPTION 1: resize the two new IBFs. Calculate new size.
//        // Delete original IBF
//        for (int split = 0; split < number_of_splits; split++){
//            index_tuples[split] = find_empty_bin_idx(index, ibf_idx);   // get indices of the empty bins on the higher level to serve as new merged bins.. import function from insertions.
//            // create new IBF
//            robin_hood::unordered_flat_set<size_t> parent_kmers{};
//            auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, arguments, false); // for now use current size. perhaps create a new construct_ibf function.
//            // add ibf to ibf_vector
//                //new_ibf_idx  =
//            for (auto bin_idx in split_indices[split];){
//                auto filenames = find_filename(ibf_idx, bin_idx); // get filenames
//                auto number_of_bins = is_split_bin(ibf_idx, bin_idx);
//                kmers  =  // load kmers
//                insert_into_ibf(parent_kmers, kmers, std::make_tuple(new_ibf_idx, bin_idx, number_of_bins), index, false);  // insert in new IBF in
//                //bin_idx ++ number_of_bins
//            }
//            //for each split, also insert 'parent_kmers' in the merged bin at index_tuples[split];
//            insert_into_ibf(parent_kmers, std::make_tuple(std::get<0>(index_tuples[0]), std::get<1>(index_tuples[0]), 1), index);
//        }
//
//
//        // OPTION 2: or only pull them apart?
//        // using bit shifts
//
//
//    //update auxiliiary data structures.
//    //      next_ibf_id
//    //      previous_ibf_id
//    //      occupancy table
//}
//
//std::vector<std::vector<size_t>> find_best_split(ibf_idx, number_of_splits=2){
//    std::vector<std::vector<size_t>> tb_idxs; //array of (number_of_splits) arrays, with bin_idxes
//    // finds best split of an IBF based on kmer counts or
//    // make sure split bins remain together
//    // default/easy: just split in 2.
//    // group/sort based on the occupancy table.
//
//    // can you use a function from chopper for this?
//
//    return (tb_idxs);
//}
//

//// METHOD 2
/*!\brief Does a partial rebuild for the subtree at index ibf_idx
 * \details The algorithm rebuilds a subtree of the hierarchical interleaved bloom filter completely.
 * (1) It computes the set of filenames belonging to all user bins in the subtree, and stores those in a text file
 * (2) It computes a good layout by calling Choppers's layout algorithm.
 * (3) It rebuilds the subtree by calling hierachical build function.
 * (4) It merges the index of the HIBF of the newly obtained subtree with the original index.
 * \param[in] ibf_idx the location of the subtree in the HIBF
 * \param[in] index the original HIBF
 * \param[in] update_arguments the arguments that were passed with the update that was to be done on the HIBF.
 * \author Myrthe
 */
 // Function to convert Vector to Set
std::set<std::string> convertToSet(std::vector<std::string> v)
{
    // Declaring the set
    std::set<std::string> s;

    // Traverse the Vector
    for (std::string x : v) {

        // Insert each element
        // into the Set
        s.insert(x);
    }

    // Return the resultant Set
    return s;
}


void full_rebuild(std::tuple<size_t,size_t> index_tuple,
                  raptor_index<index_structure::hibf> & index,
                  update_arguments const & update_arguments) {

    index.ibf().initialize_ibf_sizes();
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin(); // also update filename_to_idx
    //0) Create layout arguments
    chopper::configuration layout_arguments = layout_config(index, update_arguments); // create the arguments to run the layout algorithm with.
    //write_filenames(layout_arguments.data_file, filenames_subtree); // should be run if chopper count is run.
    //1) Store kmer counts together with the filenames as text file.
    std::vector<std::string> filenames = index.ibf().user_bins.user_bin_filenames; // all filenames
    auto kmer_counts_filenames = get_kmer_counts(index, convertToSet(filenames));
    write_kmer_counts(kmer_counts_filenames, layout_arguments.count_filename);
    //2) call chopper layout on the stored filenames.
    call_layout(layout_arguments);
    //3) call hierarchical build.
    raptor_index<index_structure::hibf> new_index{}; //create the empty HIBF of the subtree.
    build_arguments build_arguments = build_config(layout_arguments.data_file, update_arguments, layout_arguments); // create the arguments to run the build algorithm with.
    call_build(build_arguments, new_index, true);
    index = new_index;
    //4) initialize additional datastructures for the subindex.
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin();
    index.ibf().initialize_previous_ibf_id();
    index.ibf().initialize_ibf_sizes();
    std::filesystem::remove_all("tmp");
}

/*!\brief Does a partial rebuild for the subtree at index ibf_idx
 * \details The algorithm rebuilds a subtree of the hierarchical interleaved bloom filter completely.
 * (1) It computes the set of filenames belonging to all user bins in the subtree, and stores those in a text file
 * (2) It computes a good layout by calling Choppers's layout algorithm.
 * (3) It rebuilds the subtree by calling hierachical build function.
 * (4) It merges the index of the HIBF of the newly obtained subtree with the original index.
 * \param[in] ibf_idx the location of the subtree in the HIBF
 * \param[in] index the original HIBF
 * \param[in] update_arguments the arguments that were passed with the update that was to be done on the HIBF.
 * \author Myrthe
 */
void partial_rebuild(std::tuple<size_t,size_t> index_tuple,
                  raptor_index<index_structure::hibf> & index,
                  update_arguments const & update_arguments,
                  int number_of_splits) // instead of setting a default parameter here, add this to the update arguments
{
    size_t ibf_idx = std::get<0>(index_tuple);
    size_t bin_idx = std::get<1>(index_tuple);
    size_t next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][bin_idx];
    //1.1) Obtain filenames from all lower bins and kmer counts. Perhaps using occupancy table.
    auto filenames_subtree = index.ibf().filenames_children(next_ibf_idx);
    //1.2) Obtain the kmer counts
    auto kmer_counts_filenames = get_kmer_counts(index, filenames_subtree);
    //1.3) Define how to split the IBF, in terms of merged bin indexes and user bin filenames
    auto split_filenames = find_best_split(kmer_counts_filenames, number_of_splits);
    auto split_idxs = split_ibf(index_tuple, index, number_of_splits); // returns the new indices
    //1.4) remove IBFs of the to-be-rebuild subtree from the  index.
    remove_ibfs(index, next_ibf_idx);

    for (int split = 0; split < number_of_splits; split++){ // for each subindex that is rebuild.
        //0) Create layout arguments
        chopper::configuration layout_arguments = layout_config(index, update_arguments); // create the arguments to run the layout algorithm with.
        //write_filenames(layout_arguments.data_file, filenames_subtree); // should be run if chopper count is run.
        //1) Store kmer counts together with the filenames as text file.
        write_kmer_counts(split_filenames[split], layout_arguments.count_filename);
        //2) call chopper layout on the stored filenames.
        call_layout(layout_arguments);
        //3) call hierarchical build.
        raptor_index<index_structure::hibf> subindex{}; //create the empty HIBF of the subtree.
        build_arguments build_arguments = build_config(layout_arguments.data_file, update_arguments, layout_arguments); // create the arguments to run the build algorithm with.
        robin_hood::unordered_flat_set<size_t> root_kmers = call_build(build_arguments, subindex, false); // the last argument sets the 'is_root' parameter to false such that kmers are added to `parent_kmers`. However, this causes that no shuffling takes place reducing the efficiency a bit. This could be further optimized perhaps by adding an extra parameter.
        insert_into_ibf(root_kmers, std::make_tuple(ibf_idx, split_idxs[split], 1), index, std::make_tuple(0,0));       // also fills the (new) MB.
        //4) initialize additional datastructures for the subindex.
        subindex.ibf().initialize_previous_ibf_id();
        subindex.ibf().initialize_ibf_sizes(); // TODO fpr and occupancy table are not initialized, apart from with 0's (bcs fprmax is 0?)? But why is the occupancy table assigned with 0s
        //5) merge the index of the HIBF of the newly obtained subtree with the original index.
        attach_subindex(index, subindex, std::make_tuple(ibf_idx, split_idxs[split]));
    }

    // update datastructures
    index.ibf().initialize_ibf_sizes();
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin(); // this also updates the filename_to_idx datastructure
    // remove temporary files
    std::filesystem::remove_all("tmp");
}

/*!\brief
 * \details splits a merged bin and its child IBF.
 * \param[in] index_tuple the ibf_idx and bin_idx of the merged bin that has reached the FPR limit.
 * \param[in] index the original HIBF
 * \param[in] number_of_splits this may be 1 if it is not split. It is 2 by default.
 * \author Myrthe Willemsen
 */
std::vector<uint64_t> split_ibf(std::tuple<size_t,size_t> index_tuple, //rename to split_mb
               raptor_index<index_structure::hibf> & index,
               int number_of_splits)
{
    std::vector<uint64_t> tb_idxs(number_of_splits);
    index.ibf().delete_tbs(std::get<0>(index_tuple), std::get<1>(index_tuple));    // Empty the merged bin.

    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level IBF to serve as new merged bins.
            tb_idxs[split] = find_empty_bin_idx(index, std::get<0>(index_tuple));       // find an empty bin for a new MB on this IBF or resize. Import from insertions.
    }
    return tb_idxs;
}


/*!\brief Splits the set of user bins in the subtree in `n` similarly sized subsets.
 * \details User bins are devided among `n` subindexes, `n` being th number of splits.
 * The IBF is split on a user bin level, to make sure that no user bin enters both IBFs, which is otherwise difficult because of split bins.
 * Alternatively, one can use the counts for these filenames by calling Chopper count.
 * To take into account sequence similarity, the splitting should be done after rearranging user bins as part of Chopper layout.
 * \param[in] kmer_counts_filenames a vector of tuples with k-mer counts (1) for each filename (2) present in the subtree of the IBF.
 * \param[in] number_of_splits t
 * \author Myrthe
 */
std::vector<std::vector<std::tuple<size_t, std::string>>> find_best_split( //rename to split_filenames
        std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
        size_t number_of_splits){
    size_t sum_kmer_count = 0;
    for (auto& tuple : kmer_counts_filenames) sum_kmer_count += std::get<0>(tuple);
    size_t percentile = sum_kmer_count/number_of_splits; // sum over k-mer counts and find when the threshold of the sum is first exceeded.
    size_t cumulative_sum = 0; size_t filename_idx = 0; size_t split_idx =0;
    std::vector<std::vector<std::tuple<size_t, std::string>>> split_filenames;
    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level to serve as new merged bins.
        while (cumulative_sum < percentile * (split + 1) // while the threshold for the next split not yet exceeded
        && filename_idx < kmer_counts_filenames.size()){ // and while we did not reach the last file
            cumulative_sum += std::get<0>(kmer_counts_filenames[filename_idx]);
            filename_idx += 1;
        }
        split_filenames.push_back(std::vector(        // create a new vector with filename indices.
            std::ranges::next(kmer_counts_filenames.begin(), split_idx, kmer_counts_filenames.end()),  // std::ranges::next(iterator, number, bound) is the same as iterator + number, but bound: it cannot go out of range.
            std::ranges::next(kmer_counts_filenames.begin(), filename_idx + 1, kmer_counts_filenames.end())));
        split_idx = filename_idx;
    }
    return split_filenames;
}


/*!\brief Store kmer or minimizer counts for each specified file stored within the HIBF.
 * \details The algorithm obtains the total kmer count for each file and stores those counts together with the
 * filenames to a text file. Alternatively, one can create/extract counts for these filenames by calling Chopper count
 * \param[in] filenames a set of filenames of which the counts should be obtained
 * \param[in] index the original HIBF
 * \return[out]
 * \author Myrthe
 */
std::vector<std::tuple<size_t, std::string>> get_kmer_counts(raptor_index<index_structure::hibf> const & index,
                                                             std::set<std::string> const & filenames)
{
    std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames{};
    for (std::string const & filename : filenames)
    {
        if (std::filesystem::path filename_as_path{filename}; filename_as_path.extension() != ".empty_bin") // construct filename_as_path such that it can be used for comparison without going out of scope
        {
            int const kmer_count = index.ibf().get_occupancy_file(filename);
            kmer_counts_filenames.push_back(std::make_tuple(kmer_count, filename));
        }
    }
    std::ranges::sort(kmer_counts_filenames);
    return kmer_counts_filenames; // array of tuples with filename and k-mer count, sorted by kmer count.
}

/*!\brief Store kmer or minimizer counts for each specified file stored within the HIBF.
 * \details The algorithm obtains the total kmer count for each file and stores those counts together with the
 * filenames to a text file. Alternatively, one can create/extract counts for these filenames by calling Chopper count
 * prints "filename kmer_count  filename" on each line.
 * \param[in] kmer_counts_filenames a set of filenames of which the counts should be stored.
 * \param[in] count_filename the filename to which the counts should be stored.
 * \author Myrthe
 */
void write_kmer_counts(std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
                       std::filesystem::path count_filename){
    std::ofstream out_stream(count_filename.c_str());
    for (auto const & tuple : kmer_counts_filenames){
                out_stream << std::get<1>(tuple) + '\t' + std::to_string(std::get<0>(tuple)) + '\t' + std::get<1>(tuple)  + '\n' ;
        }
    }


/*!\brief Writes each of the specified files to a text file
 * \param[in] filenames a set of filenames of which the counts should be stored.
 * \param[in] user_bin_filenames the filename to which the counts should be stored.
 * \author Myrthe
 */
void write_filenames(std::string bin_path, std::set<std::string> user_bin_filenames){
    std::ofstream out_stream(bin_path.c_str());
    for (auto const & filename : user_bin_filenames)
        { if (std::filesystem::path(filename).extension() !=".empty_bin"){
                out_stream << filename + '\n';
            }
        }
    }

/*!\brief Creates a configuration object which is passed to chopper's layout algorithm.
 * \param[in] update_arguments the file containing all paths to the user bins for which a layout should be computed
 * \param[in] index the original HIBF
 * \param[in] bin_paths
 * \author Myrthe Willemsen
 */
chopper::configuration layout_config(raptor_index<index_structure::hibf> & index,
                                     update_arguments const & update_arguments){
    std::filesystem::create_directories("tmp");
    chopper::configuration config{};
    config.input_prefix = "tmp/temporary_layout"; // all temporary files ill be stored in the temporary folder. //input_prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = config.input_prefix + ".tsv";
    config.data_file = "tmp/subtree_bin_paths.txt"; // input of bins paths
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory); // here the count filename and output_prefix are set.
    config.sketch_directory = update_arguments.sketch_directory;
    config.rearrange_user_bins = update_arguments.similarity; // indicates whether updates should account for user bin's similarities.
    config.update_ubs = update_arguments.empty_bin_percentage; // percentage of empty bins drawn from distributions //makes sure to use empty bins.
    config.tmax = update_arguments.tmax;
    config.false_positive_rate = index.ibf().fpr_max;
    return config;
}


/*!\brief Calls the layout algorithm from the chopper library
* \param[in] layout_arguments configuration object with parameters required for calling the layout algorithm
* \warning an extra enter/return at the start of the bin file will cause segmentation faults in chopper.
* \author Myrthe
*/
void call_layout(chopper::configuration & layout_arguments){
    int exit_code{};
    try
    {
        exit_code |= chopper::layout::execute(layout_arguments);
    }
    catch (sharg::parser_error const & ext){    // GCOVR_EXCL_START
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
    }
}

/*!\brief Creates a configuration object which is passed to hierarchical build function.
* \param[in] subtree_bin_paths the file containing all paths to the user bins for which a layout should be computed
* \author Myrthe
*/
build_arguments build_config(std::string subtree_bin_paths, update_arguments const & update_arguments, chopper::configuration layout_arguments){
    build_arguments build_arguments{};
    build_arguments.kmer_size = 20;
    build_arguments.window_size = 23;
    build_arguments.fpr = layout_arguments.false_positive_rate; //index.false_positive_rate
    build_arguments.is_hibf = true;
    build_arguments.bin_file = layout_arguments.output_filename; //layout_file
    return build_arguments;
}

/*!\brief Calls the hierarchical build algorithm
* \param[in] build_arguments configuration object with parameters required for calling the building algorithm
* \param[out] index the newly created index
* \return root_kmers the total set of k-mers stored at the root, and that can be used to store in a merged bin if this regards a subindex.
* \author Myrthe
*/
template <seqan3::data_layout data_layout_mode>
robin_hood::unordered_flat_set<size_t> call_build(build_arguments & arguments,
                raptor_index<hierarchical_interleaved_bloom_filter<data_layout_mode>> & index,
                bool is_root){
    hibf::build_data<data_layout_mode> data{};
    robin_hood::unordered_flat_set<size_t> root_kmers = raptor::hibf::create_ibfs_from_chopper_pack(data, arguments, is_root);
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hibf.user_bins.num_user_bins(); ++i)
        bin_path.push_back(std::vector<std::string>{data.hibf.user_bins.filename_of_user_bin(i)});
    index.ibf() = std::move(data.hibf); //instead of creating the index object here.
    return root_kmers;
}
/*!\brief Helper function that removes the given indices from a vector.
* \author Myrthe Willemsen
*/
template <typename T> void remove_indices(std::vector<size_t> indices_to_remove, std::vector<T> & vector) {
    if (indices_to_remove.size() > 1) std::sort(indices_to_remove.rbegin(), indices_to_remove.rend());   // sort descending to prevent wrong indices being removed on the way. Note that this has resulted in an error in the sdsl library. If it still does, use normal sort and a reverse loop
        for (int i : indices_to_remove) {
            vector.erase(vector.begin() + i);
        }
    }

/*!\brief Prunes subtree from the original HIBF
 * \details One should remove the IBFs in the original index which were part of the subtree that had to be rebuild.
 * `indices_map` is used to map IBF indices of the original HIBF to those in the new HIBF. If the IBF at an old index is not in the new HIBF, it will return -1.
 * If using some sort of splitting, then removing only needs to happen once since both new subindexes share the same original ibfs.
 * \param[in|out] index the original HIBF
 * \param[in] ibf_idx the index of the IBF where the subtree needs to be removed, including the ibf_idx itself.
 * \author Myrthe Willemsen
 */

void remove_ibfs(raptor_index<index_structure::hibf> & index, size_t ibf_idx){
    // Store which original indices in the IBF were the subindex that had to be rebuild?
    std::vector<size_t> indices_to_remove = index.ibf().ibf_indices_childeren(ibf_idx);  // Create a map that maps remaining IBF indices of the original HIBF to their new indices
    indices_to_remove.push_back(ibf_idx);
    std::vector<int> indices_map; indices_map.resize(index.ibf().ibf_vector.size()); // The map has the size in number of IBF's of the original index.
    std::ranges::fill(indices_map, -1); // initialize with -1, such that the merged bins pointing to a removed IBF will have a pointer to -1.
    int counter = 0;// Initialize the result vector
    for (int i = 0; i <= index.ibf().ibf_vector.size(); i++) {
        if (std::find(indices_to_remove.begin(), indices_to_remove.end(), i) == indices_to_remove.end()) {  // If the current element is not in indices_to_remove.
            indices_map[i] = counter; // Add it to the result vector, such that indices_map[i] = counter. Do not use .push_back, because that will make the list not long enough
            counter += 1;
        }
    };
    // Remove vectors of indices of subindex datastructures like next_ibf and previous_ibf
    remove_indices(indices_to_remove, index.ibf().ibf_vector);
    remove_indices(indices_to_remove, index.ibf().next_ibf_id);
    remove_indices(indices_to_remove, index.ibf().previous_ibf_id);
    remove_indices(indices_to_remove, index.ibf().fpr_table);
    remove_indices(indices_to_remove, index.ibf().occupancy_table);
    remove_indices(indices_to_remove, index.ibf().user_bins.ibf_bin_to_filename_position);
    index.ibf().ibf_sizes.erase(index.ibf().ibf_sizes.begin());
    for (size_t ibf_idx{0}; ibf_idx < index.ibf().next_ibf_id.size(); ibf_idx++) {     // Replace the (remaining) indices that have to be replaced.
        for (size_t i{0}; i < index.ibf().next_ibf_id[ibf_idx].size(); ++i){
            auto & next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][i];
            if (next_ibf_idx >=0){
                assert(static_cast<size_t>(next_ibf_idx) < indices_map.size()); // this assertion prevents an out of range segmentation fault
                next_ibf_idx = indices_map[next_ibf_idx];
            }
        }
    }
    for (size_t ibf_idx = 0; ibf_idx < index.ibf().previous_ibf_id.size(); ++ibf_idx){
        auto & previous_ibf_idx = std::get<0>(index.ibf().previous_ibf_id[ibf_idx]);
        previous_ibf_idx = indices_map[previous_ibf_idx];
    }
}

/*!\brief Merges the original HIBF with the pruned subtree, with the rebuild subtree.
 * \details One should remove the IBFs in the original index which were part of the subtree that had to be rebuild.
 * When doing some sort of splitting on the merged bins, run this function twice, once for each subindex.
 * \param[in|out] index the original HIBF
 * \param[in] subindex the HIBF subtree that has been rebuild
 * \param[in] ibf_idx the index of the IBF where the subtree needs to be removed. (irrelevant?)
 * \param[in] index_tuple the index of the IBF and the bin index of the MB to which the subtree needs to be attached.
 * \author Myrthe Willemsen
 */
void attach_subindex(raptor_index<index_structure::hibf> & index,
                     raptor_index<index_structure::hibf> & subindex,
                     std::tuple<size_t, size_t> index_tuple){
    // Add new rows representing the subindex.
    size_t ibf_count_before_appending = index.ibf().ibf_count();
    for (size_t ibf_idx{0}; ibf_idx < subindex.ibf().next_ibf_id.size(); ++ibf_idx){ // Add the size of the `index`, in number of IBFs, to all IBF indices in subindex's next_ibf_id.
        for (size_t bin_idx{0}; bin_idx < subindex.ibf().next_ibf_id[ibf_idx].size(); ++bin_idx){
             subindex.ibf().next_ibf_id[ibf_idx][bin_idx] += ibf_count_before_appending;
        }
        std::get<0>(subindex.ibf().previous_ibf_id[ibf_idx]) += ibf_count_before_appending; // Add the size of the `index`, in number of IBFs, to the IBF indices present in previous_ibf_id of the subindex.
    }

    auto append_to_vector = [] (auto & index_vector, auto & subindex_vector){
        index_vector.insert(index_vector.end(), subindex_vector.begin(), subindex_vector.end());
    };
    append_to_vector(index.ibf().ibf_vector, subindex.ibf().ibf_vector); // append the vector of the new index to the original index.
    append_to_vector(index.ibf().next_ibf_id, subindex.ibf().next_ibf_id);
    append_to_vector(index.ibf().previous_ibf_id, subindex.ibf().previous_ibf_id);
    append_to_vector(index.ibf().fpr_table, subindex.ibf().fpr_table);
    append_to_vector(index.ibf().occupancy_table, subindex.ibf().occupancy_table);


    // TODO test
    //  : (1) map filename indices of the new index to those in the original index. Then you do not additionally need to append the filenames.
    // the next line will suffice.
    // use idx_to_filname (subindex)  --> filename --> filename_to_idx (original index) --> idx
                // (2) Alternatively, remove the filenames from the retained index. And add the filename lists. And update the filename positions from the subindex.
            // user_bin_filenames --> use a list of filenames that are in the subindex, and re.



    size_t subindex_n_files = subindex.ibf().user_bins.user_bin_filenames.size();
    std::vector<int> indices_map; indices_map.resize(subindex_n_files); // The map has the size in number of files in the subindex.
    std::ranges::fill(indices_map, -1);
    for (size_t filename_idx{0}; filename_idx < subindex_n_files; filename_idx++){ // TODO correct + doc
        auto filename = subindex.ibf().user_bins.user_bin_filenames[filename_idx];
        subindex.ibf().user_bins.user_bin_filenames[filename_idx];
        indices_map[filename_idx] = (index.ibf().user_bins.filename_to_idx.find(filename))->second;
    }
    auto & bin_to_file = subindex.ibf().user_bins.ibf_bin_to_filename_position;
    for (size_t ibf_idx{0}; ibf_idx < bin_to_file.size(); ++ibf_idx){ // Add the size of the `index`, in number of IBFs, to all IBF indices in subindex's next_ibf_id.
        for (size_t bin_idx{0}; bin_idx < bin_to_file[ibf_idx].size(); ++bin_idx){
                 bin_to_file[ibf_idx][bin_idx] =
                         indices_map[bin_to_file[ibf_idx][bin_idx]];
            }
        }
    append_to_vector(index.ibf().user_bins.ibf_bin_to_filename_position, bin_to_file);

    // Update the indices in one entry of the supporting tables, where are subindex must be attached, such that they refer to our new subindex.
    auto ibf_idx = std::get<0>(index_tuple);
    auto bin_idx = std::get<1>(index_tuple);
    index.ibf().next_ibf_id[ibf_idx][bin_idx] = ibf_count_before_appending;
    index.ibf().previous_ibf_id[ibf_count_before_appending] = std::make_tuple(ibf_idx, bin_idx);
}

} // namespace raptor

