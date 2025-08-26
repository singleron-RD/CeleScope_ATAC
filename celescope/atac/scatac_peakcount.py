import os
import shutil
import itertools
import multiprocessing as mp
import numpy
import h5py
import scipy.sparse as sp_sparse
import random
import string
import gzip
import argparse
from collections import defaultdict
from functools import partial


def randomString(stringLength=10):
    """Generate a random string of fixed length"""
    letters = string.ascii_lowercase
    return "".join(random.choice(letters) for i in range(stringLength))


tmp = randomString()


def is_gzip(filename):
    """Check if the file is gzipped."""
    with gzip.open(filename, "r") as f:
        try:
            f.read(1)
        except OSError:
            return False
    return True


def universal_open(filename, mode):
    if is_gzip(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def bedtools_intersect(barcode, peak_bed):
    """Intersect frag file with peak file to genearate the count output."""
    os.system(
        "bedtools intersect -wa -a "
        + peak_bed
        + " -b "
        + tmp
        + "/"
        + barcode
        + " -c | awk '{if ($NF>0) print $0}' > "
        + tmp
        + "/"
        + barcode
        + ".bed"
    )
    return tmp + "/" + barcode + ".bed"


def filter_fragment_file(barcode_file, frag_file, count_cutoff):
    """Filter fragment file and only keep the valid barcode ones."""

    barcode_list = []
    barcode_out = defaultdict(list)
    barcode_count = defaultdict(lambda: 0)

    if barcode_file == "":
        # no barcode file
        # read fragment file to get barcode
        fhd = universal_open(frag_file, "rt")

        for line in fhd:
            line = line.strip().split("\t")
            barcode_out[line[3]].append(line[0] + "\t" + line[1] + "\t" + line[2])
            barcode_count[line[3]] += int(line[4])

        fhd.close()

        for k in barcode_count:
            if barcode_count[k] >= count_cutoff:
                barcode_list.append(k)

    else:
        # read barcode file
        fhd = universal_open(barcode_file, "rt")

        for line in fhd:
            barcode_list.append(line.strip())
            barcode_out[line.strip()] = []
        fhd.close()

        fhd = universal_open(frag_file, "rt")

        for line in fhd:
            line = line.strip().split("\t")
            barcode_out[line[3]].append(line[0] + "\t" + line[1] + "\t" + line[2])
        fhd.close()

    for k in barcode_list:
        outf = open(tmp + "/" + k, "w")
        for line in sorted(barcode_out[k]):
            print(line, file=outf)
        outf.close()

    return barcode_list


def write_10X_h5(
    filename, matrix, features, barcodes, genome="GRCh38", datatype="Peak"
):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, "w")
    if datatype == "Peak":
        M = sp_sparse.csc_matrix(matrix, dtype=numpy.int32)
    else:
        M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)
    B = numpy.array(barcodes, dtype="|S200")
    P = numpy.array(features, dtype="|S100")
    GM = numpy.array([genome] * len(features), dtype="|S10")
    FT = numpy.array([datatype] * len(features), dtype="|S100")
    AT = numpy.array(["genome"], dtype="|S10")
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


def generate_count_matrix(count_list, peak_list):
    peak_count = {}
    for peak in peak_list:
        peak_count[peak] = sp_sparse.dok_matrix((1, len(count_list)), dtype=numpy.int32)

    barcodes = []
    for i in range(0, len(count_list)):
        barcodes.append(count_list[i].split("/")[-1][:-4])
        for line in open(count_list[i], "r"):
            line = line.strip().split("\t")
            if line[0] + ":" + line[1] + "-" + line[2] in peak_count:
                peak_count[line[0] + ":" + line[1] + "-" + line[2]][0, i] = int(
                    line[-1]
                )

    matrix = []
    for k in peak_list:
        matrix.append(peak_count[k])
    matrix = sp_sparse.vstack(matrix)

    return (matrix, barcodes)


def merge_count_file(peak_file, count_list, count_file, cores, genome="GRCh38"):
    """Merge the intersectBed result into the count table."""

    peak_list = []

    fhd = universal_open(peak_file, "rt")
    for line in fhd:
        line = line.strip().split("\t")
        peak_list.append(line[0] + ":" + line[1] + "-" + line[2])
    fhd.close()
    peak_list = sorted(peak_list)

    count_list_split = []
    i = 0
    while i < len(count_list):
        if i + 1000 < len(count_list):
            count_list_split.append(count_list[i : (i + 1000)])
            i += 1000
        else:
            count_list_split.append(count_list[i : len(count_list)])
            i = len(count_list)

    pool = mp.Pool(processes=int(cores))
    partial_generate_count_matrix = partial(generate_count_matrix, peak_list=peak_list)
    result = pool.map_async(partial_generate_count_matrix, count_list_split)
    pool.close()
    pool.join()

    result_list = result.get()
    matrix_list = [i[0] for i in result_list]
    barcode_list = [i[1] for i in result_list]

    matrix = sp_sparse.hstack(matrix_list)
    barcodes = list(itertools.chain.from_iterable(barcode_list))

    write_10X_h5(
        count_file, matrix, peak_list, barcodes, genome=genome, datatype="Peaks"
    )


def peakcount(barcode, count_cutoff, cores, peak, fragment, outprefix):
    try:
        os.makedirs("peak")
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    os.makedirs(tmp)

    count_file = os.path.join("peak", outprefix + "_peak_count.h5")

    barcode_list = filter_fragment_file(barcode, fragment, count_cutoff)

    pool = mp.Pool(processes=int(cores))
    partial_bedtools_intersect = partial(bedtools_intersect, peak_bed=peak)
    result = pool.map_async(partial_bedtools_intersect, barcode_list)
    pool.close()
    pool.join()

    count_list = result.get()
    merge_count_file(peak, count_list, count_file, cores)
    shutil.rmtree(tmp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate peak count matrix")
    parser.add_argument(
        "--barcode", help="Location of valid cell barcode file", required=True
    )
    parser.add_argument(
        "--count_cutoff",
        help="Cutoff for the number of count in each cell",
        required=True,
    )
    parser.add_argument("--cores", help="Number of cores to use", required=True)
    parser.add_argument("--peak", help="Location of peak file", required=True)
    parser.add_argument(
        "--fragment", help="Location of fragments.tsv file", required=True
    )
    parser.add_argument("--outprefix", help="Prefix of output files", required=True)
    args = parser.parse_args()
    peakcount(
        args.barcode,
        args.count_cutoff,
        args.cores,
        args.peak,
        args.fragment,
        args.outprefix,
    )
