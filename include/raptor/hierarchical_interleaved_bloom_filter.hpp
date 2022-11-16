// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <ranges>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#ifndef RAPTOR_HIBF_HAS_COUNT
#    define RAPTOR_HIBF_HAS_COUNT 0
#endif

namespace raptor
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [seqan3::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [seqan3::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [seqan3::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [seqan3::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [seqan3::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
 * , or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
 * another IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used
 * to further distinguish between merged bins.
 *
 * In this example, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single technical
 * bin, and another IBF was added for the merged bin.
 * \image html hibf.svg width=90%
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;

#if RAPTOR_HIBF_HAS_COUNT
    // Forward declaration
    template <std::integral value_t>
    class counting_agent_type;
#endif

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of an individual Bloom filter.
    using ibf_t = seqan3::interleaved_bloom_filter<data_layout_mode_>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter() = default;                                              //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter const &) = default;                        //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default;            //!< Defaulted.

    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector<ibf_t> ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

        /*!\brief Stores for each bin in each IBF of the HIBF the ID of the previous IBF.
     * \details
     * Assume we look up an IBF `i`, i.e. `previous_ibf_id[i]`.
     * If `(j =number_ibfs, b =0)` is returned, there is no higher level IBF, and IBF 'i' is thus the root IBF.
     * If `j != number_ibfs` is returned, there is a higher level IBF , with `j` is the ID of the higher
     * level IBF in ibf_vector and bin `b` is the corresponding merged bin in this higher level IBF.
          * \author Myrthe

     */
    std::vector<std::tuple<int64_t, int64_t>> previous_ibf_id;

    //!\brief The table with FPR or occupancy values, [ibf_idx][bin_idx]
    std::vector<std::vector<double>> fpr_table;

    //!\brief The table with FPR or occupancy values, [ibf_idx][bin_idx]
    std::vector<std::vector<size_t>> occupancy_table; // again it may make more sense to have this as part of ibf_vector, interleaved_bloom_filter does not belong to this project.

    //!\brief The underlying user bins.
    user_bins user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent{*this};
    }

    bool is_merged_bin(size_t ibf_idx, size_t bin_idx){
        // Alternative to current implementation: use next_ibf_id.
        // assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`. If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
        auto const current_filename_index = user_bins.filename_index(ibf_idx, bin_idx);
        if (current_filename_index < 0){
            return true;
        }else{ return false;} // merged bin
    }

    void initialize_previous_ibf_id()
    {    //higher_ibf_id[ibf_idx] --> (higher_ibf_idx, bin_idx) . If higher_ibf_idx = std::vector<ibf_t>.size(), it does not exist .
        // size is het zelfde als ibf_vector
                // dit moet ook geinitializeerd worden, bij traversing IBF, en dan naar lagere te gaan.
        // loop over next_ibf_id --> werkt niet, want dan weet je niet waar de laagste vandaan komen. Of kan het wel ? next_ibf_id[ibf_id_high][bin_id_high]=[ibf_id_low]
        // traverse HIBF. [ibf_id
        //auto number_ibfs = ibf_vector.size();
        previous_ibf_id.resize(ibf_vector.size());
        std::fill(previous_ibf_id.begin(), previous_ibf_id.end(), std::make_tuple(0,0)); // instead of number_ibfs, initialize with root_idx (=0). Then if previous_ibf_id[ibf_idx]=ibf_idx, then ibf_idx=root_idx. Previously std::make_tuple(number_of_ibs, 0)
        for (uint64_t ibf_id_high=0; ibf_id_high < next_ibf_id.size(); ibf_id_high++){
            for (uint64_t bin_idx=0;  bin_idx < next_ibf_id[ibf_id_high].size(); bin_idx++){
                uint64_t ibf_id_low = next_ibf_id[ibf_id_high][bin_idx];
                if (ibf_id_low != ibf_id_high){ // it might happen that ibf_id_low == ibf_id_high for LBs, and we should not overwrite the entry.
                     previous_ibf_id[ibf_id_low] = std::make_tuple(ibf_id_high, bin_idx);
                }
            }
        }
    }

    /*!\brief Approximates the false positive rate
     * \details Calculates the approximate probability of a false positive hit of a k-mer for a technical bin, based on the bin size, the number of kmers and the number of hash functions.
     * The second function also uses the number of split bins, and returns the FPR for the user bin that was split among those bin. The probability to obtain a false positive answer for the user bin when querying those s techinical bins is 1 - (1 - fpr_tb)^s.
     * The third function also accounts for the .  number_deleted_kmers/(alphabet_size**k) is the probabibility that a too be searched kmer is the same as a deleted kmer (that is still present in the bloom filter)
     * \ref
     * Find the formulas in the corresponding publication
     * \param[in] m size, in number of bits, of the bloom filter.
     * \param[in] n k-mer count, i.e. the number of kmers stored in the technical bin.
     * \param[in] h number of hash functions used.
     * \param[in] s split bins: number of bins the among which the n k-mers are divided.
     * \param[in] d deleted k-mers: number of k-mers that were virtually deleted from this bin.
     * \author Myrthe
     * \comment It might be computationally faster to work with log2 fpr values, to prevent the power to the h
     */
    double approximate_fpr(size_t & m, size_t & n, size_t & h){// n=#occupied bits / number of kmers, m=length BF, h =#hash functions, deleted_kmers=0
       // when actually using deleted kmers, these need to be stored somewhere, and the value alphabet**k should be pre computed.
        return pow((1-exp(-h*n/m)), h); // - deleted_kmers/(alphabet**k)
    }

    double approximate_fpr(size_t & m, size_t & n, size_t & h, size_t & s){// n=#occupied bits / number of kmers, m=length BF, h =#hash functions, deleted_kmers=0
        return (1-pow(1 - approximate_fpr(m,n,h), s)); // - deleted_kmers/(alphabet**k)
    }


    /*!\brief Update false positive rate of a TB.
     * \details Calculates the approximate false positive rate of a certain technical bin, given its indices in the HIBF, using the approximate_fpr function and updates .
     * \param[in] ibf_idx index of IBF in HIBF
     * \param[in] bin_idx index of TB in IBF
     * \param[in] number_of_bins index of TB in IBF.
     * \author Myrthe
     */
    double update_fpr(size_t ibf_idx, size_t bin_idx, size_t number_of_bins = 1)
    {
        auto& ibf = ibf_vector[ibf_idx]; //  select the IBF
        auto kmer_count = occupancy_table[ibf_idx][bin_idx];
        auto bin_size = ibf.bin_size();
        auto hash_funs = ibf.hash_function_count();
        auto fpr = approximate_fpr(bin_size, kmer_count, hash_funs, number_of_bins); // if a split bin, i think it would be best to store the joint fpr, e.g. of both bins plus multiple testing. We can assume that the occupancy is alomst equal for the tbs among which a UB was split.
        for (size_t offset=0; offset < number_of_bins; offset++){  //loop over split bins and add multiple testing correction
            fpr_table[ibf_idx][bin_idx+offset] = fpr;
        }
        return fpr;
    }

    double get_fpr(size_t ibf_idx, size_t const bin_idx){
        return(fpr_table[ibf_idx][bin_idx]);
    }

    void update_occupancy_table(size_t kmer_count, size_t ibf_idx, size_t const bin_idx, size_t const number_of_bins = 1){
        auto occupancy_per_bin = std::round(kmer_count/number_of_bins);
        for (size_t offset=0; offset < number_of_bins; offset++){ // update FPR table and occupancy=#kmer table. Perhaps do this before inserting.
            occupancy_table[ibf_idx][bin_idx + offset] += occupancy_per_bin;
        }
    }

        /*!\brief Update false positive rate of a TB.
     * \details Calculates the approximate false positive rate of a certain technical bin, given its indices in the HIBF, using the approximate_fpr function and updates .

     * \author Myrthe
     */
//    void resize_ibf_occupancy_table(size_t ibf_idx, size_t new_bin_count){
//        assert(new_bin_count >= occupancy_table[ibf_idx].size()); // check that new bin count is larger then the size.
//        occupancy_table[ibf_idx].resize(new_bin_count);
//    }
//    void resize_ibf_fpr_table(size_t ibf_idx, size_t new_bin_count){
//        assert(new_bin_count >= fpr_table[ibf_idx].size()); // check that new bin count is larger then the size.
//        fpr_table[ibf_idx].resize(new_bin_count);
//    }
//    void resize_hibf_fpr_table(size_t ibf_idx, size_t new_ibf_size){ //add new IBF to last index, with a vector of the size of the new IBF.
//        assert(new_bin_count >= fpr_table[ibf_idx].size()); // check that new bin count is larger then the size.
//        fpr_table.resize(fpr_table.size()+1);
//    }

// function for adding new ibf to the above

    double ibf_max_kmers(size_t ibf_idx, double fpr=0.05){ // get fpr of a TB in a certain IBF. Makes more sense if this is part of interleaved_bloom_filter.hpp, but that is not part of the project. Myrthe. or bin_index const bin_idx
        auto& ibf = ibf_vector[ibf_idx]; //  select the IBF       or hibf_ptr->ibf_vector.[ibf_idx] OR  auto& ibf = index.ibf().ibf_vector[ibf_idx]
        size_t bin_size = ibf.bin_size();
        size_t hash_funs = ibf.hash_function_count();//arguments.hash
        return approximate_kmer_capacity(bin_size, fpr, hash_funs);
    }

    size_t approximate_kmer_capacity(size_t & m, double & fpr, size_t & h) //max number of kmers that fits in a TB, a value specific to an IBF
    // -> (although for split bins also possible, taking into account fpr correction)
    // -> should also possibly take into account fpr rate of merged bin above.
    {
        return std::floor((m*log(1-pow(fpr, 1/h))/h));
    } // can be a property stored per ibf, so it does not need to be recalculated. (Can be 0 for empty IBFs)
        // bin_size_in_bits https://github.com/seqan/raptor/blob/7fe02401bb4f191e2ef4e1454ea9a1c7756816ca/src/build/hibf/bin_size_in_bits.cpp can be used to calculate m from fpr, n and h

    size_t number_of_bins(size_t ibf_idx, size_t kmer_count,  double fpr=0.05){
        // given an ibf, with a certain bin size, how many technical bins are needed to store a certain number of kmers, considering the multiple testing problem?
        auto& ibf = ibf_vector[ibf_idx]; //  select the IBF       or hibf_ptr->ibf_vector.[ibf_idx] OR  auto& ibf = index.ibf().ibf_vector[ibf_idx]
        size_t bin_size = ibf.bin_size();
        size_t hash_funs = ibf.hash_function_count();
        size_t number_of_bins = std::ceil( kmer_count / ibf_max_kmers(ibf_idx));         // first guess
        while (approximate_fpr(bin_size, kmer_count, hash_funs, number_of_bins) > fpr){
            number_of_bins++;
        }
        return number_of_bins;
    }


    // Max_n (fpr, length bf in ibf, update_seqs=false)  calculate maximum
    // if not update_seqs, then it is simply dependend on the fpr_max and bin size.

//    double average_fpr(ibf_idx){
//        // average over the fprs, skipping the empty bins
//    }
//    double min_fpr(ibf_idx){
//        // average over the fprs, skipping the empty bins
//    }

#if RAPTOR_HIBF_HAS_COUNT
    /*!\brief Returns a counting_agent_type to be used for counting.
     * \tparam value_t The type to use for the counters; must model std::integral.
     */
    template <std::integral value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
#endif

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond
};

/*!\brief Bookkeeping for user and technical bins.
 */
template <seqan3::data_layout data_layout_mode>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::user_bins
{
private:
    //!\brief Contains filenames of all user bins.
    std::vector<std::string> user_bin_filenames;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_filename_position[i][b]`.
     * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
     * Otherwise, the returned value `j` can be used to access the corresponding filename `user_bin_filenames[j]`.
     */
    std::vector<std::vector<int64_t>> ibf_bin_to_filename_position{};

    //!\brief Maps filenames to their indices isn the list of filenames
    std::unordered_map<std::string, uint64_t> filename_to_idx;     //filename --> filename_index


    /*!\brief Stores for each filename ID each of its bins and corresponding IBF in the HIBF
     * \details
     * .
     */
    std::vector<std::tuple<uint64_t, uint64_t, uint16_t>> filename_position_to_ibf_bin{}; // filename index --> (ibf_idx, bin_idx, number_of_bins)
public:
    //!\brief Creates
    void initialize_filename_position_to_ibf_bin()
    {
        filename_position_to_ibf_bin.resize(user_bin_filenames.size());
        std::fill(filename_position_to_ibf_bin.begin(), filename_position_to_ibf_bin.end(), std::make_tuple(0,0,0));
        for (uint64_t idx=0; idx < user_bin_filenames.size(); idx++){ // warning: comparison of integer expressions of different signedness: solve this by declaring idx as size_t ipv int ?  â€˜intâ€™ and â€˜std::vector<std::__cxx11::basic_string<char> >::size_typeâ€™ {aka â€˜long unsigned intâ€™} [-Wsign-compare]
            std::string filename = user_bin_filenames[idx];             // Question: create string view from filename?
            filename_to_idx.emplace(filename, idx);
        } // or something similar to parse_user_bin_ids

        for (uint64_t ibf_idx=0; ibf_idx < ibf_bin_to_filename_position.size(); ibf_idx++){ // or should this be uint64_t instead of size_t, here and in the following.
            for (uint64_t bin_idx=0; bin_idx < ibf_bin_to_filename_position[ibf_idx].size(); bin_idx++){
                int64_t filename_position = ibf_bin_to_filename_position[ibf_idx][bin_idx]; // Question: should I use references here?
                if (std::get<2>(filename_position_to_ibf_bin[filename_position])){ //for split bins.
                    ++std::get<2>(filename_position_to_ibf_bin[filename_position]);
                }else{
                    filename_position_to_ibf_bin[filename_position] = std::make_tuple(ibf_idx, bin_idx, 1); // as a user bin can take up multiple bins, this should consist of multiple bin_idx or multiple tuples.
                }
            }
        }

    }

    /*!\brief Returns the bin and IBF indices within the HIBF of a given user bin, specified by filename.
     * \details
     * Should only be used when filename_position_to_ibf_bin has been created, and after checking the filename is present in the filename_to_idx map
     * \author Myrthe
     */
    std::tuple <uint64_t, uint64_t, uint16_t> find_filename(std::string filename)
    {
        return filename_position_to_ibf_bin[filename_to_idx[filename]];
    }

    //!\brief Checks if the filename is already present in the HIBF.
    bool exists_filename(const std::string & filename)
    {
    if (filename_to_idx.find(filename) == filename_to_idx.end())
        return false; // filename/user bin is not yet present in HIBF
    else
        return true; // filename/user bin does already exist in HIBF
    }

    void update_filename_indices(std::string filename, size_t const ibf_idx, size_t const bin_idx, size_t const number_of_bins){
        user_bin_filenames.push_back(filename); // or resize it first to add filename to the end of "filenames"
        filename_to_idx[filename] = user_bin_filenames.size(); // We should not assume filename_to_idx has the same size as user_bin_filenames, but as we do not remove deleted bins from 'user_bin_filenames' nor from 'filename_position_to_ibf_bin', we should use user_bin_filenames.size()
        // for index_pair in index_pairs:
        filename_position_to_ibf_bin[user_bin_filenames.size()] = std::make_tuple(ibf_idx, bin_idx, number_of_bins); // or resize it first?
        ibf_bin_to_filename_position[ibf_idx][bin_idx] = user_bin_filenames.size() ; // possibly resize to update ibf_bin_to_filename_position and filenames
    }

    void delete_filename(const std::string & filename){
        // --> from user_bin_filenames
    // --> from filename_to_idx
    // set to 'empty_bin', is probably cheaper than removing the value? perhaps related to a number?
    // The other 2 ( ) do not have to be changed per se, but make sure that when inserting a new bin, the 'number of bins' is updated.
        // put the original cardinality or ibf size instead of x
        user_bin_filenames[filename_to_idx[filename]] ="x.empty_bin"; //if you delete the item, and decrease the list size, you would have to remake ibf_bin_to_filename_position. Empty bins should also be initiliaze as empty_bin in this list (if building is correct/as should). Check this.
        filename_to_idx.erase(filename);
    }

    //!\brief Returns the number of managed user bins.
    size_t num_user_bins() const noexcept
    {
        return user_bin_filenames.size();
    }

    //!\brief Changes the number of managed IBFs.
    void set_ibf_count(size_t const size)
    {
        ibf_bin_to_filename_position.resize(size);
    }

    //!\brief Changes the number of managed user bins.
    void set_user_bin_count(size_t const size)
    {
        user_bin_filenames.resize(size);
    }

    //!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
    std::vector<int64_t> & bin_indices_of_ibf(size_t const idx)
    {
        return ibf_bin_to_filename_position[idx];
    }

    //!\brief Returns the filename of the `idx`th user bin.
    std::string & filename_of_user_bin(size_t const idx)
    {
        return user_bin_filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return user_bin_filenames[ibf_bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`th IBF.
     *        An empty string is returned for merged bins.
     */
    auto operator[](size_t const ibf_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx]
             | std::views::transform(
                   [this](int64_t i)
                   {
                       if (i == -1)
                           return std::string{};
                       else
                           return user_bin_filenames[i];
                   });
    }

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
    template <typename stream_t>
    void write_filenames(stream_t & out_stream) const
    {   size_t position{};
        std::string line{};
        for (auto const & filename : user_bin_filenames)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        archive(user_bin_filenames);
        archive(ibf_bin_to_filename_position);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the hibf::hierarchical_interleaved_bloom_filter.
 * \see hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
 * \details
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 */
template <seqan3::data_layout data_layout_mode> // TODO: value_t as template?
class hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);

        uint16_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold) // if ibf contains the query.
                    bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold); //recursive
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer.emplace_back(current_filename_index);
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default;                                     //!< Defaulted.
    membership_agent(membership_agent const &) = default;             //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default;                  //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default;      //!< Defaulted.
    ~membership_agent() = default;                                    //!< Defaulted.

    /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent(hibf_t const & hibf) : hibf_ptr(std::addressof(hibf))
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector<int64_t> result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a user bin if there are at least this many hits.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
    {
        assert(hibf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        bulk_contains_impl(values, 0, threshold);

        std::ranges::sort(result_buffer); // TODO: necessary?

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                             size_t const threshold) && noexcept = delete;
    //!\}
};


#if RAPTOR_HIBF_HAS_COUNT
/*!\brief Manages counting ranges of values for the hibf::hierarchical_interleaved_bloom_filter.
 */
template <seqan3::data_layout data_layout_mode>
template <std::integral value_t>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::counting_agent_type
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive bulk counting.
    template <std::ranges::forward_range value_range_t>
    void bulk_count_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<value_t>();
        auto & result = agent.bulk_count(values);

        value_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];
            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_count_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer[current_filename_index] = sum;
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default;                                        //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default;             //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default;                  //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default;      //!< Defaulted.
    ~counting_agent_type() = default;                                       //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit counting_agent_type(hibf_t const & hibf) :
        hibf_ptr(std::addressof(hibf)),
        result_buffer(hibf_ptr->user_bins.num_user_bins())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    seqan3::counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::forward_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     * \param[in] threshold Do not recurse into merged bins with less than this many hits. Default: 1.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) & noexcept
    {
        assert(hibf_ptr != nullptr);
        assert(threshold > 0u);
        assert(result_buffer.size() == hibf_ptr->user_bins.num_user_bins());

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, static_cast<value_t>(0u));

        bulk_count_impl(values, 0, threshold);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) && noexcept = delete;
    //!\}
};
#endif // RAPTOR_HIBF_HAS_COUNT

} // namespace raptor
