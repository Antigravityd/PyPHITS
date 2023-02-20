from fortranformat import FortranRecordReader


# Read in an ASCII dump file produced by a PHITS tally
def kf_decode(n):
    simple = {2212: "proton", 2112: "neutron", 211: "pion+", 111: "pion0", -211: "pion-", -13: "muon+", 13: "muon-",
              321: "kaon+", 311: "kaon0", -321: "kaon-", 11: "electron", -11: "positron", 22: "photon",
              1000002: "deuteron", 1000003: "triton", 2000003: "3He", 2000004: "alpha"}
    if n in simple:
        return simple[n]
    elif n > 1000000:
        a = int(str(n)[:5])
        z = (n - a) / 1000000
    else:
        # "other particle"
        return "wacky"


def read_dump(name, columns):
    rd = FortranRecordReader('(30(1p1e24.15))')
    acc = dict.fromkeys(columns)
    with open(name, 'r') as dmp:
        for line in dmp:
            for i, val in enumerate(rd.read(line)):
                if val is not None:
                    acc[col[i]] = val

    return acc
