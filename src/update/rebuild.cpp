#include <raptor/search/search_single.hpp> // to make sure index_structure::hibf_compresse exists here.
#include <raptor/upgrade/load_hibf.hpp>
#include <raptor/upgrade/get_fpr.hpp>
#include "raptor/build/store_index.hpp"


namespace raptor
{
void split_IBF(upgrade_arguments const & arguments)
// arguements:
//  1 higher level, ibf_idx,
// ibf_idx of merged bin
//
{

}


}


// for each of the bins in bin_paths: (maybe sort? or do this smartly somehow?


// request_ibf_idx
// request_user_bin_idx

// when using such an existing function, make sure the arguments have the right arguments type, by using a template (see store_index.hpp how to) , typename arguments_t