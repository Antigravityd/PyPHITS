from fortranformat import FortranRecordReader
from pandas import DataFrame
from valspec import elements, kf_decode

def read_dump(name: str, columns: list[str], return_type: str) -> dict:
    """Given a path to a PHITS dump file and names of the record entries in order,
    produce a semantically equivalent dictionary of lists/numpy array/Pandas dataframe of the contents."""
    rd = FortranRecordReader('(30(1p1d24.15))') # PHITS documentation says e, but code says d---latter is consistent with behavior
    acc = dict.fromkeys(columns)
    with open(name, 'r') as dmp:
        for line in dmp:
            for i, val in enumerate(rd.read(line)):
                if val is not None:
                    col = columns[i]
                    if acc[col] is None:
                        acc[col] = [val if col != "particle" else kf_decode(val)]
                    else:
                        acc[col] += val if col != "particle" else kf_decode(val)

    if return_type == "dict":
        return acc
    elif return_type == "numpy":
        return np.fromiter(acc.items())
    elif return_type == "pandas":
        return DataFrame.from_dict(acc)


    return acc
