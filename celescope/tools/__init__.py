import os
import celescope

# barcode
PATTERN_DICT = {
    "auto": None,
    "atac1": "C6L4C8L4C6",
    "atac2": "C4L4C7L4C5",
    "atac3": "C4L6C4L6C4",
    "customized": None,
}

TOOLS_DIR = os.path.dirname(celescope.__file__) + "/tools"
DATA_DIR = os.path.dirname(celescope.__file__) + "/data"

# count
OUTS_DIR = "outs"
RAW_MATRIX_DIR_SUFFIX = "raw"
FILTERED_MATRIX_DIR_SUFFIX = "filtered"
MATRIX_FILE_NAME = "matrix.mtx.gz"
FEATURE_FILE_NAME = "features.tsv.gz"
BARCODE_FILE_NAME = "barcodes.tsv.gz"
