import os
import argparse
import h5py
import numpy as np
import tables
import collections
import scipy.sparse as sp_sparse


FeatureBCMatrix = collections.namedtuple(
    "FeatureBCMatrix", ["ids", "names", "barcodes", "matrix"]
)


def read_10X_h5(filename):
    """Read 10X HDF5 files, support both gene expression and peaks."""
    with tables.open_file(filename, "r") as f:
        try:
            group = f.get_node(f.root, "matrix")
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        feature_group = getattr(group, "features")
        ids = getattr(feature_group, "id").read()
        names = getattr(feature_group, "name").read()
        barcodes = getattr(group, "barcodes").read()
        data = getattr(group, "data").read()
        indices = getattr(group, "indices").read()
        indptr = getattr(group, "indptr").read()
        shape = getattr(group, "shape").read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return FeatureBCMatrix(ids, names, barcodes, matrix)


def write_10X_h5(
    filename, matrix, features, barcodes, genome="GRCh38", datatype="Peak"
):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, "w")
    if datatype == "Peak":
        M = sp_sparse.csc_matrix(matrix, dtype=np.int32)
    else:
        M = sp_sparse.csc_matrix(matrix, dtype=np.float32)
    B = np.array(barcodes, dtype="|S200")
    P = np.array(features, dtype="|S100")
    GM = np.array([genome] * len(features), dtype="|S10")
    FT = np.array([datatype] * len(features), dtype="|S100")
    AT = np.array(["genome"], dtype="|S10")
    mat = f.create_group("matrix")
    mat.create_dataset("barcodes", data=B)
    mat.create_dataset("data", data=M.data)
    mat.create_dataset("indices", data=M.indices)
    mat.create_dataset("indptr", data=M.indptr)
    mat.create_dataset("shape", data=M.shape)
    fet = mat.create_group("features")
    fet.create_dataset("_all_tag_keys", data=AT)
    fet.create_dataset("feature_type", data=FT)
    fet.create_dataset("genome", data=GM)
    fet.create_dataset("id", data=P)
    fet.create_dataset("name", data=P)
    f.close()


def Filter(rawmatrix, feature, barcode, peak_cutoff, cell_cutoff, outprefix):
    peaks_per_cell = np.asarray((rawmatrix > 0).sum(axis=0))

    passed_cell = peaks_per_cell >= peak_cutoff

    cells_per_peak = np.asarray((rawmatrix > 0).sum(axis=1))
    passed_peak = cells_per_peak >= cell_cutoff
    passed_peak = np.transpose(passed_peak)

    # gene = [True]*rawmatrix.shape[0]
    passed_cell_matrix = rawmatrix[
        np.ix_(passed_peak.tolist()[0], passed_cell.tolist()[0])
    ]

    passed_barcodes = np.array(barcode)[passed_cell.tolist()[0]].tolist()
    passed_peaks = np.array(feature)[passed_peak.tolist()[0]].tolist()

    # passed_barcodes = [bc.decode('utf-8') for bc in passed_barcodes]

    write_10X_h5(
        outprefix + "_filtered_peak_count.h5",
        matrix=passed_cell_matrix,
        features=passed_peaks,
        barcodes=passed_barcodes,
        datatype="Peak",
    )


def scatac_qc(outprefix, peakcount, peak_cutoff, cell_cutoff):
    try:
        os.makedirs("peak")
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    scatac_count = read_10X_h5(peakcount)
    peakmatrix = scatac_count.matrix
    features = scatac_count.names.tolist()
    barcodes = scatac_count.barcodes.tolist()

    if isinstance(features[0], bytes):
        features = [i.decode() for i in features]
    if isinstance(barcodes[0], bytes):
        barcodes = [i.decode() for i in barcodes]

    filename = os.path.join("peak", outprefix)

    Filter(
        rawmatrix=peakmatrix,
        feature=features,
        barcode=barcodes,
        peak_cutoff=peak_cutoff,
        cell_cutoff=cell_cutoff,
        outprefix=filename,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filter peak count matrix")
    parser.add_argument(
        "--peakcount", help="Location of peak count matrix file", required=True
    )
    parser.add_argument("--outprefix", help="Prefix of output files", required=True)
    parser.add_argument(
        "--peak_cutoff",
        type=int,
        help="Minimum number of peaks included in each cell",
        required=True,
    )
    parser.add_argument(
        "--cell_cutoff",
        type=int,
        help="Minimum number of cells covered by each peak",
        required=True,
    )
    args = parser.parse_args()
    scatac_qc(args.outprefix, args.peakcount, args.peak_cutoff, args.cell_cutoff)
