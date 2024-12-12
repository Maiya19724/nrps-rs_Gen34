# -*- coding: utf-8 -*-
import os
import subprocess
import logging
from typing import Optional
from collections import namedtuple
from Bio import SearchIO

ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"

# Modify to your actual local path to ensure the hmmer data file exists
ADOMAINS_HMM = "/home/symlab/anaconda3/envs/antismash7.1/lib/python3.10/site-packages/antismash/modules/nrps_pks/data/aa-activating.aroundLys.hmm"
ACTIVE_SITE_PROFILE_NAME = "aa-activating-core.198-334"

# Positions specified according to the original code
HMM_SITE_POSITIONS_10 = [46, 47, 50, 92, 124, 126, 153, 161, 162]
HMM_SITE_POSITIONS_34 = [
    12, 15, 16, 40, 45, 46, 47, 48, 49, 50, 51, 54, 92, 93, 124, 125, 126, 127,
    128, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165,
]

def verify_good_sequence(sequence: str) -> bool:
    """Ensures a sequence is valid"""
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True

def extract_by_reference_positions(query_seq, profile_seq, positions):
    """
    Modified logic:
    Instead of immediately returning None when encountering out-of-bounds and gaps, use '-' as placeholders to ensure a string is always returned.
    """
    extracted = []
    seq_len = len(profile_seq)
    for pos in positions:
        if pos < 0 or pos >= seq_len:
            # Out of alignment range, replace with '-'
            extracted.append('-')
        else:
            # Profile corresponding position
            p_char = profile_seq[pos]
            if p_char == '-':
                # The corresponding profile position is a gap
                extracted.append('-')
            else:
                # Extract the amino acid corresponding to the query
                q_aa = query_seq[pos] if pos < len(query_seq) else '-'
                extracted.append(q_aa)
    return "".join(extracted)

class ModularDomain:
    """
    Simplified ModularDomain class, only for demonstration.
    """
    def __init__(self, translation, name="mock_domain"):
        self.translation = translation
        self._name = name

    def get_name(self):
        return self._name

def get_alignment_against_profile(sequence: str, hmm_file: str, domain_name: str, max_evalue: float = 0.1):
    from io import StringIO

    # Write to temporary FASTA
    query_fasta = "query_tmp.faa"
    with open(query_fasta, "w") as fh:
        fh.write(">query\n")
        fh.write(sequence + "\n")

    # Run hmmsearch
    hmmsearch_cmd = [
        "hmmsearch",
        "--domE", str(max_evalue),
        hmm_file,
        query_fasta
    ]
    result = subprocess.run(hmmsearch_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logging.error("hmmsearch failed: %s", result.stderr)
        os.remove(query_fasta)
        return None

    handle = StringIO(result.stdout)
    try:
        parse_res = list(SearchIO.parse(handle, 'hmmer3-text'))
    except ValueError:
        os.remove(query_fasta)
        return None

    os.remove(query_fasta)

    if len(parse_res) == 0:
        return None

    qresult = parse_res[0]

    # Ensure qresult.id matches the HMM profile name
    if qresult.id != domain_name:
        return None

    for hit in qresult.hits:
        for hsp in hit.hsps:
            query_aln_seq = hsp.hit.seq  # Aligned segment of the actual protein sequence
            profile_aln_seq = hsp.query.seq  # Aligned segment of the HMM profile
            hit_start = hsp.query_start  # Starting position in the profile

            AlnRec = namedtuple("AlnRec", ["seq"])
            HitObj = namedtuple("HitObj", ["hit_start", "aln"])
            hit_obj = HitObj(
                hit_start=hit_start,
                aln=[AlnRec(seq=query_aln_seq),
                     AlnRec(seq=profile_aln_seq)]
            )
            return hit_obj

    return None

def get_a_dom_signatures(domain: ModularDomain, max_evalue: float = 0.1) -> tuple[Optional[str], Optional[str]]:
    """Extract 10 / 34 AA NRPS signatures from A domains"""
    assert verify_good_sequence(domain.translation)

    hit = get_alignment_against_profile(domain.translation, ADOMAINS_HMM,
                                        ACTIVE_SITE_PROFILE_NAME, max_evalue=max_evalue)

    if not hit:
        return None, None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    aa10 = extract_by_reference_positions(query, profile, [p - offset for p in HMM_SITE_POSITIONS_10])
    if not aa10 or set(aa10) == {'-'}:
        logging.debug("extracting the 10 AA profile failed for %s", domain.get_name())
        return None, None
    aa10 += "K"  # since only 9 positions are included and the 10th is generally K

    aa34 = extract_by_reference_positions(query, profile, [p - offset for p in HMM_SITE_POSITIONS_34])
    if not aa34 or set(aa34) == {'-'}:
        logging.debug("extracting the 34 AA profile failed for %s", domain.get_name())
        # If you want to output results even if a complete match is not possible, do not return None here
        # return None, None

    return aa10, aa34

if __name__ == "__main__":
    # Test case (requires a real existing A domain sequence)
    test_sequence = "MINLPQSDSINPLEKQKYYQVSHGQRRLWILQHISESSSAYNIRLAMRINGKLDVSILTAAFQELVNRHEILRTTFTSVGGNIKQVIHEEIPTEQLIIFKDLRAEKDAELLADDLIQESANLPFSLEKLPLIRVLLVQINSEEFIFGLTVHHIIGDAQSLDIIFQEFFTLYSAYTQQKKAELPPLALQYKDYAAWQNQWLESEEVREQHEYWCEVFSEHPPILNLPIDFPRPQVKSSQAASHTYQFSQSLSQQLQTFATQNNTTLFITLLTFFKILLYRHTGQRDLVVGIPISARNHPDLENQIGFYVNTLALRSLLPEGVTFKQVLGEVTNTCLDAYEYRDYPFDKLVSALNLERDLSRNPLFDVMFSLLSKESKTVLKIPGLEHQPYPLKPRMAQFDMTWSFFEDSHNLTLVIEYEPDLFLPETIARMNGHFLQIIKSVIQNPDGKLSEINLLTPEEQHQLLIEWNQTEVAYSQQCLHHLFEEKVRDNSEAIALIFEGKKLTYRELNGRANQVGHYLQEKGVTSEVLVGIFIERSFEMIIGILGIMKAGGAYVPLDPNYPPERLDYMISDSAISLLLTQQSLVQFLPENQAEILCLDTDWSRIANYSQENLTSPVKTENLAYVIYTSGSTGKPKGVMNIHQGICNTLKYNIDNYNLNSEDRILQITPFSFDVSVWEVFSSLTSGATLVVTKPDGYKDIDYLIDLIVQEQVTCFTCVPSILRVFLQHPKSKDCHCLKRVIVGGEALSYELNQRFFQQLNCELYNAYGPTEVAVDATVWCCQPNSQLISIGRPIANVQVYILDSYLQPVPIGVAGELHIGGMGLARGYLNQAELTAEKFIPHPFAEGKLYKTGDLARYLPDGNIEYLGRIDNQVKLRGFRIELGEIQTVLETHPNVEQTVVIMREDTLYNQRLVAYVMRKEPLLTSQDLRRFLQQKLPVYMIPSAFVMLSDFPLNPNGKIDFHKLPIPDETSLVESPYLAPRNQTETILVSLWQQLLQAGKIGVNDNFFELGGHSLQAMNLMALIYEKIAIEIPLSMIYEKPTVAELSDYIIYAQEMNIQPKERPYVVFNKEQEKAVFLFPPALGFAAAYANLAEYITDYALYTFRYISDELMLEKYAELIEDLAQDQDIKLMGHSAGGFLAMLMAQQLESRGRVVSDVILLDTYRGGREAKKADMSEIQEGVDAFLLNPKRQELRGYFLGNQKLRDRTYHQVWEYFNFLWNSDLKNVQIQGTIHLIRAEGNYEARDDWTQATKGERINHYATGIHREMIDPPYLPKNALIINSILNPKQS"
    domain = ModularDomain(test_sequence, name="test_domain")
    aa10, aa34 = get_a_dom_signatures(domain, max_evalue=0.1)
    print("10 AA Signature:", aa10)
    print("34 AA Signature:", aa34)
