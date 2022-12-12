#include <raptor/search/search_single.hpp> // to make sure index_structure::hibf_compresse exists here.
#include <raptor/upgrade/load_hibf.hpp>
#include <raptor/upgrade/get_fpr.hpp>
#include "raptor/build/store_index.hpp"


namespace raptor
{
 // METHOD 1
void split_ibf(size_t ibf_idx,
                  raptor_index<index_structure::hibf> & index, upgrade_arguments const & arguments)
{   int number_of_splits = 2;
    std::vector<std::tuple <uint64_t, uint64_t>> index_tuples;  // or initialize with size {number_of_splits};
    index_tuples[0] = index.ibf().previous_ibf_id[ibf_idx];     // get indices of the current merged bin on the higher level
    index.ibf().delete_tbs(std::get<0>(index_tuples[0]), std::get<1>(index_tuples[0]));    // empty it

//    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level to serve as new merged bins.
//            index_tuples[split] = find_empty_bin_idx(index, ibf_idx);   //import from insertions.
//    }
    // get empty bin(s) on higher level, for number of splits -1
    std::vector<std::vector<size_t>> split_indices = find_best_split();
    // create new IBFs

        // OPTION 1: resize the two new IBFs. Calculate new size.
        // Delete original IBF
        for (int split = 0; split < number_of_splits; split++){
            index_tuples[split] = find_empty_bin_idx(index, ibf_idx);   // get indices of the empty bins on the higher level to serve as new merged bins.. import function from insertions.
            // create new IBF
            robin_hood::unordered_flat_set<size_t> parent_kmers{};
            auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, arguments, false); // for now use current size. perhaps create a new construct_ibf function.
            // add ibf to ibf_vector
                //new_ibf_idx  =
            for (auto bin_idx in split_indices[split];){
                auto filenames = find_filename(ibf_idx, bin_idx); // get filenames
                auto number_of_bins = is_split_bin(ibf_idx, bin_idx);
                kmers  =  // load kmers
                insert_into_ibf(parent_kmers, kmers, std::make_tuple(new_ibf_idx, bin_idx, number_of_bins), index, false);  // insert in new IBF in
                //bin_idx ++ number_of_bins
            }
            //for each split, also insert 'parent_kmers' in the merged bin at index_tuples[split];
            insert_into_ibf(parent_kmers, std::make_tuple(std::get<0>(index_tuples[0]), std::get<1>(index_tuples[0]), 1), index);
        }


        // OPTION 2: or only pull them apart?
        // using bit shifts


    //update auxiliiary data structures.
    //      next_ibf_id
    //      previous_ibf_id
    //      occupancy table
}

std::vector<std::vector<size_t>> find_best_split(ibf_idx, number_of_splits=2){
    std::vector<std::vector<size_t>> tb_idxs; //array of (number_of_splits) arrays, with bin_idxes
    // finds best split of an IBF based on kmer counts or
    // make sure split bins remain together
    // default/easy: just split in 2.
    // group/sort based on the occupancy table.

    // can you use a function from chopper for this?

    return (tb_idxs);
}


 // METHOD 2
void recall_dp(size_t ibf_idx,
                  raptor_index<index_structure::hibf> & index, upgrade_arguments const & arguments)
{    // first do it without updating hyper loglog (sketch_directory, chopper_sketch_sketches) Perhaps use config from layout file. --> shows the count_filename and sketch_directory.

    //1) obtain filenames from all lower bins and kmer counts. Perhaps using occupancy table.
    // 1.2) and store filenames as text file., with count_filename = "chopper_sketch.count"
   get_kmer_counts(index, arguments);
    //3) call Chopper layout on output,
    //see test file in chopper, which call execute.
}

// store kmer/minimizer counts for each filename
void get_kmer_counts(raptor_index<index_structure::hibf> & index, upgrade_arguments const & arguments){
    // alternatively, create/extract counts for these filenames by calling Chopper count
    // clear the count file first, with count_filename = "chopper_sketch.count"

    for (auto & filename : index.bin_path()){        //for (auto & filename : index.ibf().user_bins){// for each file in index.user_bins, perhaps I have to use .user_bin_filenames, but this is private... or index.filenames
        if (std::filesystem::path(filename).extension() !=".empty_bin"){
                index.ibf().get_occupancy_file(filename);
                write_count_file_line(cluster_vector[i], weight, fout); // function in chopper, to store filenames to the text file.
        }

    }
}

} // namespace raptor

