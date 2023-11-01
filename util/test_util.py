# Should probably move const.py into the util folder, though requires refactoring.
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import const

from util import successful_custom_truncation, get_AA_from_nuc, get_grantham_category, get_max_grantham_score


assert(successful_custom_truncation({'Sequence Change': const.SYNTHETIC_BRESEQ_TRUNC_MUT_SEQ_CHANGE, 'Details': 'coding (18/ nt)'})==False)
assert(successful_custom_truncation({'Sequence Change': const.SYNTHETIC_BRESEQ_TRUNC_MUT_SEQ_CHANGE, 'Details': 'coding (19/ nt)'})==True)
assert(successful_custom_truncation({'Sequence Change': const.SYNTHETIC_BRESEQ_TRUNC_MUT_SEQ_CHANGE, 'Details': 'coding (20/ nt)'})==False)
assert(successful_custom_truncation({'Sequence Change': const.SYNTHETIC_BRESEQ_TRUNC_MUT_SEQ_CHANGE, 'Details': 'coding (1/ nt)'})==True)


assert(get_AA_from_nuc(1)==1)
assert(get_AA_from_nuc(2)==1)
assert(get_AA_from_nuc(3)==1)
assert(get_AA_from_nuc(4)==2)
assert(get_AA_from_nuc(6)==2)
assert(get_AA_from_nuc(7)==3)


import numpy as np


assert(get_grantham_category(1)=='conservative')
assert(get_grantham_category(52)=='moderately conservative')
assert(get_grantham_category(120)=='moderately radical')
assert(get_grantham_category(151)=='radical')
assert(get_grantham_category(np.nan)=='')
assert(get_grantham_category('K')=='')


assert(get_max_grantham_score('A', 'A')==0)
assert(get_max_grantham_score('A', 'C')==195)
assert(get_max_grantham_score('AC', 'CI')==198)
assert(get_max_grantham_score('K', '*') is np.nan)