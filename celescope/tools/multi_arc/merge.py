import io
import os
import glob
import json
import argparse
from jinja2 import Environment, FileSystemLoader, select_autoescape


def render_html(rna_path, atac_path, outdir):
    # jinja env
    env = Environment(
        loader=FileSystemLoader(os.path.dirname(__file__) + "/templates/"),
        autoescape=select_autoescape(["html", "xml"]),
    )

    rna_json = glob.glob(f"{rna_path}/.data.json")[0]
    atac_json = glob.glob(f"{atac_path}/.data.json")[0]

    with open(rna_json, "r") as f:
        rna_data = json.load(f)
    with open(atac_json, "r") as f:
        atac_data = json.load(f)
    rna_data = {f"rna_{k}": v for k, v in rna_data.items()}
    atac_data = {f"atac_{k}": v for k, v in atac_data.items()}
    data = rna_data.copy()
    data.update(atac_data)

    template = env.get_template("html/merge/base.html")
    with io.open(f"{outdir}/multi_arc_report.html", "w", encoding="utf8") as f:
        html = template.render(data)
        f.write(html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="merge atac rna report")
    parser.add_argument("--rna", help="rna_path", required=True)
    parser.add_argument("--atac", help="atac_path", required=True)
    parser.add_argument("--outdir", help="out_path", required=True)
    args = parser.parse_args()
    render_html(args.rna, args.atac, args.outdir)
