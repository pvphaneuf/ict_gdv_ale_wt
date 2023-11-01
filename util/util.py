from aledbmutil.mut import get_DEL_INS_MOB_nuc_start_pos
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import const


# Code to address initial inaccurate representation of truncation predictions.
def successful_custom_truncation(mut):
    is_truncation = False
    if mut['Sequence Change'] == const.SYNTHETIC_BRESEQ_TRUNC_MUT_SEQ_CHANGE:
        nuc_pos = get_DEL_INS_MOB_nuc_start_pos(mut['Details'])
        if nuc_pos % 3 == 1:  # For the custom truncations to be successful, the insertion must happen at position 1 a codon.
            is_truncation = True
    return is_truncation


import math

def get_AA_from_nuc(nuc):
    return math.ceil(nuc/3)


def get_AA_stop_pos(var):
    AA_stop_pos = -1
    if var.INSCRIPTA_CoordinateType=='Amino Acid':
        if var.INSCRIPTA_EditType=='Insertion':
            AA_stop_pos = var['AA start pos']
        if var.INSCRIPTA_EditType=='Deletion':
            AA_stop_pos = var['AA start pos'] + var.INSCRIPTA_NumberOfCoordinatesToDelete-1  # -1 accounting for current pos
        if var.INSCRIPTA_EditType=='Substitution':
            AA_stop_pos = var['AA start pos'] + len(var.INSCRIPTA_InsertionSequence)-1  # -1 accounting for current pos
    if var.INSCRIPTA_CoordinateType=='Nucleotide':
        if var.INSCRIPTA_EditType=='Insertion':
            AA_stop_pos = var['AA start pos']
        if var.INSCRIPTA_EditType=='Deletion':
            AA_stop_pos = get_AA_from_nuc(var.INSCRIPTA_PositionValue + var.INSCRIPTA_NumberOfCoordinatesToDelete-1)  # -1 accounting for current pos
        if var.INSCRIPTA_EditType=='Substitution':
            AA_stop_pos = get_AA_from_nuc(var.INSCRIPTA_PositionValue + len(var.INSCRIPTA_InsertionSequence)-1)  # -1 accounting for current pos
    assert(AA_stop_pos!=-1)
    return AA_stop_pos


# https://www.exeterlaboratory.com/genomics/in-silico-tools-used-for-the-interpretation-of-novel-variants/
def get_grantham_category(grantham_score):
    gc = ''
    if isinstance(grantham_score, (int, float)):
        if grantham_score <= 50:
            gc = 'conservative'
        elif (grantham_score >= 51) & (grantham_score <= 100):
            gc = 'moderately conservative'
        elif (grantham_score >= 101) & (grantham_score <= 150):
            gc = 'moderately radical'
        elif grantham_score > 150:
            gc = 'radical'
        else:
            gc = ''
    return gc


import numpy as np


def get_max_grantham_score(aa_ref_seq, aa_change_seq):
    gs_max = np.nan
    for aa_ref, aa_change in zip(aa_ref_seq, aa_change_seq):
        if aa_ref.isalpha() and aa_change.isalpha():
            if gs_max is np.nan:
                gs_max = int(const.GRANTHAM_SCORE_MAT[(aa_ref, aa_change)])
            else:
                gs_max = int(max(gs_max, const.GRANTHAM_SCORE_MAT[(aa_ref, aa_change)]))
    return gs_max


# combining Grantham score and initial variant effect predictions.
def get_predicted_variant_effect(var):
    pve = ''
    if (var['variant effect']=='synonymous') & (var['Grantham category']==''):
        pve = 'other'
    if var['variant effect'] in ['truncation', 'other']:
        pve = var['variant effect']
    if var['Grantham category'] != '':
        pve = var['Grantham category']
    return pve