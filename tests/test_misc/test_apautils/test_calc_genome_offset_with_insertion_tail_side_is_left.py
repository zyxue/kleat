import pytest

from kleat.misc.apautils import calc_genome_offset
import kleat.misc.settings as S


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [2, 2],
    [3, 3],
])
def test_for_contig_with_3_base_insertion_with_clv_before_the_insertion(ctg_offset, expected_gnm_offset):
    """
    TT         <-bridge read tail
     └CG       <-bridge read
      |        # blank line to separate the bridge read the insertion
      |AGC     <-inserted sequence
      |456     <-contig offset coord for inserted sequence
      | ┬
    ATCG GT    <-contig
    0123 789   <-contig offset coord
      ^ctg_clv   
    0123 456   <-genome offset coord
      ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [4, 3],
    [5, 3],
    [6, 3],
])
def test_for_contig_with_3_base_insertion_with_clv_inside_the_insertion(ctg_offset, expected_gnm_offset):
    """
            TT         <-bridge read tail
             └AG       <-bridge read
              ||       # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              456      <-contig offset coord for inserted sequence
       ctg_clv^|
              |┬
           ATCG GT    <-contig
           0123 789   <-contig offset coord
           0123 456   <-genome offset coord
              ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    # Here only left-tail case is tested, see
    # ./test_bridge_clv_inside_insertion.py for more comprehensive test cases
    [5, 3],
    [5, 3],
    [6, 3],
])
def test_for_contig_with_3_base_insertion_with_clv_inside_the_insertion_and_clv_is_in_the_middle_of_insertion(ctg_offset, expected_gnm_offset):
    """
             TT        <-bridge read tail
              └G       <-bridge read
               |       # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              456      <-contig offset coord for inserted sequence
       ctg_clv^|
              |┬
           ATCG GT     <-contig
           0123 789    <-contig offset coord
           0123 456    <-genome offset coord
              ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [7, 4],
    [8, 5],
])
def test_for_contig_with_3_base_insertion_with_clv_after_the_insertion(ctg_offset, expected_gnm_offset):
    """
       TT      <-bridge read tail
        └GT    <-bridge read
         |     # blank line to separate the bridge read the insertion
      AGC|     <-inserted sequence
      456|     <-contig offset coord for inserted sequence
       ┬ |
    ATCG GT    <-contig
    0123 78    <-contig offset coord
         ^ctg_clv   
    0123 45    <-genome offset coord
         ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


# Another test case with only one-base case, might be more sensitive to edge
# cases

@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    # before the insertion, see example in the docstring
    [2, 2],
    [3, 3],
])
def test_for_contig_with_1_base_insertion_with_clv_before_insertion(ctg_offset, expected_gnm_offset):
    """
    TT
     └CG       <-bridge read
      |        # blank line to separate the bridge read the insertion
      | G      <-inserted sequence
      | 4      <-contig offset coord for inserted sequence
      | ┬
    ATCG AC    <-contig
    0123 56    <-contig offset coord
      ^ctg_clv   
    0123 45    <-genome offset coord
      ^gnm_offset
    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    # in the insertion
    [4, 3],
])
def test_for_contig_with_1_base_insertion_with_clv_inside_insertion(ctg_offset, expected_gnm_offset):
    """
      TT
       └GC     <-bridge read
        |      # blank line to separate the bridge read the insertion
        G      <-inserted sequence
        4      <-contig offset coord for inserted sequence
        ^ctg_clv   
        ┬
    ATCG AC    <-contig
    0123 56    <-contig offset coord
    0123 45    <-genome offset coord
       ^gnm_offset
    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [5, 4],
    [6, 5],
])
def test_for_contig_with_1_base_insertion_with_clv_after_insertion(ctg_offset, expected_gnm_offset):
    """
       TT
        └AC    <-bridge read
         |     # blank line to separate the bridge read the insertion
        G|     <-inserted sequence
        4|     <-contig offset coord for inserted sequence
        ┬|
    ATCG AC    <-contig
    0123 56    <-contig offset coord
         ^ctg_clv   
    0123 45    <-genome offset coord
         ^gnm_offset
    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset
