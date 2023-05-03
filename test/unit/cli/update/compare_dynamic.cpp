// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include "../../../include/cli_test.hpp"

struct update : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>>
{};

// build indexes, if it doesnt exist already.

// specify parameters.

//    static inline std::filesystem::path const hibf_path(size_t const number_of_empty_bins,
//                                                       size_t const window_size,) noexcept
//    {
//        std::string name{};
//        name += std::to_string(number_of_empty_bins);
//        name += "bins";
//        name += std::to_string(window_size);
//        name += compressed ? "windowc." : "window.";
//        name += "hibf";
//        return cli_test::data(name);
//    }
//
//
//TEST_P(build_hibf, with_file)
//{
//    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
//    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;
//
//    cli_test_result const result = execute_app("raptor",
//                                               "build",
//                                               "--hibf",
//                                               "--kmer 19",
//                                               "--window",
//                                               std::to_string(window_size),
//                                               "--fpr 0.05",
//                                               "--threads",
//                                               run_parallel ? "2" : "1",
//                                               "--output raptor.index",
//                                               pack_path(number_of_repeated_bins));
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//    RAPTOR_ASSERT_ZERO_EXIT(result);
//
//    compare_index<raptor::index_structure::hibf>(
//        ibf_path(number_of_repeated_bins, window_size, is_compressed::no, is_hibf::yes),
//        "raptor.index");
//}
///////////////////////

size_t number_of_lines(std::string filename){
    std::ifstream inFile(filename);
  auto number_of_lines = std::count(std::istreambuf_iterator<char>(inFile),
             std::istreambuf_iterator<char>(), '\n');
  return (size_t) number_of_lines + 1;
}

std::string hibf_ebs = "hibf_ebs.index"; // alternatively use GetParam(); and instantiate testa suite
std::string hibf_no_ebs = "hibf_no_ebs.index";

std::string query_file = "query.fq";
size_t number_of_queries = number_of_lines(query_file);

    static inline std::vector<std::string>  extract_results(
            std::string_view const filename){
        std::ifstream search_result{filename.data()};
        std::string line;
        std::string expected_hits;
        std::vector<std::string> filenames;
        std::vector<std::string> query_results;

        while (true){
            std::getline(search_result, line);
            std::string_view line_view{line};
            if (line == "#QUERY_NAME\tUSER_BINS"){
                break;
            }
            std::string bin_filename{};
            bin_filename += line_view.substr(1, line_view.find('\t'));
            filenames.push_back(bin_filename);
        }            //  skip line

        for (size_t i = 0; i < number_of_queries; ++i){
            std::getline(search_result, line);
            std::string_view line_view{line};
            std::string bin_filename_number{};
            bin_filename_number += line_view.substr(1, line_view.find('\t'));
            std::cout << bin_filename_number << std::flush;
            int bin_number = std::stoi( bin_filename_number ); // convert to integer
            query_results.push_back(filenames[bin_number]);
        }
    return(query_results);
    // TODO : it is more complicated, because I have to map the number of the UB to the filename, and compare those.
    }

    static inline void compare_two_searches(
            std::string_view const filename1,
            std::string_view const filename2){

//        std::ifstream search_result1{filename1.data()};
//        std::ifstream search_result2{filename2.data()};
//        std::string line1;
//        std::string line2;
        std::vector<std::string> query_result1 = extract_results(filename1);
        std::vector<std::string> query_result2 = extract_results(filename2);

        for (size_t i = 0; i < number_of_queries; ++i)
        {
//            ASSERT_TRUE(std::getline(search_result1, line1));
//            ASSERT_TRUE(std::getline(search_result2, line2));
            EXPECT_EQ(query_result1[i], query_result2[i]);
            }

    }




TEST_P(update, ebs_and_no_ebs)
{

    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    // Without any empty bins
    cli_test_result const result_no_ebs =
        execute_app("raptor",
                    "search",
                    "--fpr 0.05",
                    "--output result_no_ebs.out",
                    "--error ",
                    std::to_string(number_of_errors),
                    "--p_max 0.4",
                    "--hibf",
                    "--index ",
                    data(hibf_no_ebs),
                    "--query ",
                    data(query_file));
    EXPECT_EQ(result_no_ebs.out, std::string{});
    EXPECT_EQ(result_no_ebs.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result_no_ebs);

    // With empty bins
    cli_test_result const result_ebs =
            execute_app("raptor",
                        "search",
                        "--fpr 0.05",
                        "--output result_ebs.out",
                        "--error ",
                        std::to_string(number_of_errors),
                        "--p_max 0.4",
                        "--hibf",
                        "--index ",
                        data(hibf_ebs),
                        "--query ",
                        data(query_file));
    EXPECT_EQ(result_ebs.out, std::string{});
    EXPECT_EQ(result_ebs.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result_ebs);


    compare_two_searches("result_ebs.out", "result_no_ebs.out");

}

//
//
//TEST_P(search_hibf, with_threshold)
//{
//    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();
//
//    cli_test_result const result =
//        execute_app("raptor",
//                    "search",
//                    "--fpr 0.05",
//                    "--output search.out",
//                    "--threshold 0.50",
//                    "--hibf",
//                    "--index ",
//                    ibf_path(number_of_repeated_bins, window_size, is_compressed::no, is_hibf::yes),
//                    "--query ",
//                    data("query.fq"));
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//    RAPTOR_ASSERT_ZERO_EXIT(result);
//
//    compare_search(number_of_repeated_bins, 1 /* Always finds everything */, "search.out");
//}
//
//TEST_P(search_hibf, no_hits)
//{
//    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();
//
//    cli_test_result const result =
//        execute_app("raptor",
//                    "search",
//                    "--fpr 0.05",
//                    "--output search.out",
//                    "--error ",
//                    std::to_string(number_of_errors),
//                    "--hibf",
//                    "--tau 0.99",
//                    "--index ",
//                    ibf_path(number_of_repeated_bins, window_size, is_compressed::no, is_hibf::yes),
//                    "--query ",
//                    data("query_empty.fq"));
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//    RAPTOR_ASSERT_ZERO_EXIT(result);
//
//    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::yes);
//}
//
INSTANTIATE_TEST_SUITE_P(update_suite,
                         update,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(0, 1)), // form parameters for different tests.
                         // --kmer 20, -- window 23 => parameters used for building
                         [](testing::TestParamInfo<update::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

//TEST_F(search_hibf, three_levels)
//{
//    cli_test_result const result = execute_app("raptor",
//                                               "search",
//                                               "--fpr 0.05",
//                                               "--output search.out",
//                                               "--error 0",
//                                               "--hibf",
//                                               "--index ",
//                                               data("three_levels.hibf"),
//                                               "--query ",
//                                               data("query.fq"));
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//    RAPTOR_ASSERT_ZERO_EXIT(result);
//
//    compare_search(32, 0, "search.out");
//}
