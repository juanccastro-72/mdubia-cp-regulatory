"""Offline regression tests for the core pipeline functions (no network / no genomes)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "cpreg"))
import cpreg_pipeline as cp


def test_revcomp():
    assert cp.revcomp("ACGTN") == "NACGT"
    assert cp.revcomp(cp.revcomp("ACGTACGT")) == "ACGTACGT"


def test_gc_and_clean():
    assert abs(cp.gc("GGCC") - 1.0) < 1e-9
    assert abs(cp.gc("ATAT") - 0.0) < 1e-9
    assert cp.clean("acgtx") == "ACGTN"


def test_circ_slice_wraps():
    s = "ABCDEFGHIJ"
    assert cp.circ_slice(s, 8, 12) == "IJAB"      # wraps the origin
    assert cp.circ_slice(s, 2, 5) == "CDE"


def test_dinuc_shuffle_preserves_composition():
    rng = np.random.default_rng(0)
    seq = "ATGGCAATTGACATATAATGGCCGGTTAACCGGATTACAGG" * 4
    sh = cp.dinuc_shuffle(seq, rng)
    assert len(sh) == len(seq)
    assert sorted(sh) == sorted(seq)              # same mononucleotide counts
    def dinucs(x):
        from collections import Counter
        return Counter(x[i:i+2] for i in range(len(x)-1))
    # dinucleotide counts are preserved by the Altschul-Erickson shuffle
    assert dinucs(sh) == dinucs(seq)


def test_bh_fdr_monotone_and_bounds():
    p = np.array([0.001, 0.01, 0.02, 0.5, 1.0])
    q = cp.bh_fdr(p)
    assert np.all(q >= p - 1e-9)
    assert np.all(q <= 1.0 + 1e-9)


def test_empirical_p_range():
    null = np.zeros(99)
    assert abs(cp.empirical_p(1.0, null, greater=True) - 1/100) < 1e-9


def test_pep_model_recovers_consensus():
    bg = np.array([0.34, 0.16, 0.16, 0.34])
    model = cp.PEPModel(bg)
    # planted canonical promoter: -35 TTGACA ... 17 nt ... -10 TATAAT
    prom = "GG" + "TTGACA" + "N"*17 + "TATAAT" + "GGACTGACT"
    prom = prom.replace("N", "C")
    X = cp.encode(prom)[None, :]
    score, pos, sp = model.scan(X)
    assert prom[int(pos[0]):int(pos[0])+6] == "TTGACA"
    assert int(sp[0]) == 17
    assert score[0] > 0


def test_sd_duplex_energy_negative():
    # a strong SD (AGGAGG) should hybridise favourably with the anti-SD
    dG = cp.duplex_dG("AAAAAGGAGGAAAA", "GATCACCTCCTTT")
    assert dG < -4.4


def test_best_window_mfe_and_hairpin():
    seq = "GGGGGGCCCCAAAAAGGGGGGCCCCCC" * 3
    e, struct, off, sub = cp.best_window_mfe(seq, w=40, step=5)
    assert e < 0 and len(sub) <= 40
    helix, loops = cp.hairpin_stats("(((...)))")
    assert helix == 3 and loops == 1


def test_reverse_complement_palindrome_stemloop():
    st, e = cp.mfe("GGGGGAAAAACCCCC")
    assert e < 0
