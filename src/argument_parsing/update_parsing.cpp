// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/update_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/update/update.hpp>

namespace raptor
{

void init_update_parser(sharg::parser & parser, update_arguments & arguments)
{
    init_shared_meta(parser);

    parser.add_option(arguments.bin_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "bins",
                                    .description = "File containing one file per line per bin.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.in_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "The index to upgrade. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(
        arguments.out_file,
        sharg::config{.short_id = '\0', .long_id = "output", .description = "Path to new index.", .required = true});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The original window size.",
                                    .required = true,
                                    .validator = positive_integer_validator{}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The original kmer size.",
                                    .required = true,
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.parts,
                      sharg::config{.short_id = '\0',
                                    .long_id = "parts",
                                    .description = "Original index consisted of this many parts.",
                                    .validator = power_of_two_validator{}});

    parser.add_flag(
        arguments.compressed,
        sharg::config{.short_id = '\0', .long_id = "compressed", .description = "Original index was compressed."});

    parser.add_flag(
        arguments.is_hibf,
        sharg::config{.short_id = '\0', .long_id = "hibf", .description = "Index is an HIBF.", .advanced = true}); // should we also require the fpr as input? it should be deduced from the datastruct?

        // make sure that one of the options, delete UBs, insert UBs, insert sequences or delete sequences is chosen.
    parser.add_flag(
        arguments.delete_ubs,
        sharg::config{.short_id = '\0', .long_id = "delete-UBs", .description = "Delete user bins. Provide filenames of the user bins to be deleted", .advanced = true}); // should we also require the fpr as input? it should be deduced from the datastruct?
    parser.add_flag(
        arguments.insert_ubs,
        sharg::config{.short_id = '\0', .long_id = "insert-UBs", .description = "Insert user bins. Provide filenames of the user bins to be inserted", .advanced = true}); // should we also require the fpr as input? it should be deduced from the datastruct?
    parser.add_flag(
        arguments.insert_sequences,
        sharg::config{.short_id = '\0', .long_id = "insert-sequences", .description = "Insert sequences user bins. Provide filenames of the user bins to be inserted into. Also make sure that in the same location files exists with UB_filename_insertsequence.txt", .advanced = true}); // should we also require the fpr as input? it should be deduced from the datastruct?
    parser.add_flag(
        arguments.delete_sequences,
        sharg::config{.short_id = '\0', .long_id = "delete-sequences", .description = "Insert sequences user bins. Provide filenames of the user bins to be inserted into. Also make sure that in the same location files exists with UB_filename_insertsequence.txt", .advanced = true}); // should we also require the fpr as input? it should be deduced from the datastruct?

     parser.add_option(arguments.sketch_directory,
                      sharg::config{.short_id = '\0',
                                    .long_id = "sketch-directory",
                                    .description = "If sketches should be used during rebuilding then provide their file directory",
                                    .advanced = true,
                                    .validator = sharg::input_directory_validator{}});

     parser.add_flag(arguments.similarity,
                      sharg::config{.short_id = '\0',
                                    .long_id = "sequence-similarity",
                                    .description = "I user bins should be rearranged based on sequence similarity during rebuilding",
                                    .advanced = true});

     parser.add_option(arguments.empty_bin_percentage,
                      sharg::config{.short_id = '\0',
                                    .long_id = "empty-bin-sampling",
                                    .description = "The percentage of empty bins sampled during layout computation. Default: 0.1",
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{0, 1}});

     parser.add_option(arguments.insert_sequence_appendix,
                      sharg::config{.short_id = '\0',
                                    .long_id = "sequence-insertions-appendix",
                                    .description = "User bin file appendix when inserting sequences in existing UBs",
                                    .advanced = true});

}

void update_parsing(sharg::parser & parser)
{
    update_arguments arguments{};
    init_update_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================
    if (arguments.kmer_size > arguments.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

    arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    // GCOVR_EXCL_START
    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
    // GCOVR_EXCL_STOP

    if (arguments.parts == 1)
    {
        sharg::input_file_validator{}(arguments.in_file);
    }
    else
    {
        sharg::input_file_validator validator{};
        for (size_t part{0}; part < arguments.parts; ++part)
        {
            validator(arguments.in_file.string() + std::string{"_"} + std::to_string(part));
        }
    }

    // ==========================================
    // Process bin_path
    // ==========================================
    parse_bin_path(arguments);

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_update(arguments);
}

} // namespace raptor
