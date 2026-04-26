#!/usr/bin/env python
import argparse
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


VELOCITY_SUFFIX_RE = re.compile(r"__(exon|intron|transcript)$")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ASTROutils",
        description="Small utility commands for ASTRO outputs.",
    )
    subparsers = parser.add_subparsers(dest="command")

    velocity_parser = subparsers.add_parser(
        "velocity",
        help="Run RNA velocity from a loom file or 10x-style exon/intron matrix.",
    )
    velocity_parser.add_argument(
        "--10xinput",
        dest="tenxinput",
        required=True,
        help="10x output folder, Cell Ranger outs folder, or matrix directory.",
    )
    velocity_parser.add_argument(
        "--output",
        required=True,
        help="Output directory for the velocity h5ad, tables, and plots.",
    )
    velocity_parser.add_argument(
        "--loom",
        help="Optional loom file with spliced and unspliced layers.",
    )
    velocity_parser.add_argument(
        "--gtf",
        help=(
            "ASTRO modified GTF used to count exon/intron features from a 10x "
            "BAM when the 10x matrix does not already contain velocity layers."
        ),
    )
    velocity_parser.add_argument(
        "--bam",
        help="Optional 10x BAM path. Defaults to possorted_genome_bam.bam under --10xinput.",
    )
    velocity_parser.add_argument(
        "--positions",
        help="Optional Visium tissue positions file. Defaults to tissue_positions*.csv under --10xinput.",
    )
    velocity_parser.add_argument(
        "--exon-matrix",
        help="Optional TSV/CSV matrix with cells as rows and genes as columns.",
    )
    velocity_parser.add_argument(
        "--intron-matrix",
        help="Optional TSV/CSV matrix with cells as rows and genes as columns.",
    )
    velocity_parser.add_argument(
        "--umap",
        help="Optional CSV/TSV with cell id, UMAP_1, and UMAP_2 columns.",
    )
    velocity_parser.add_argument(
        "--clusters",
        help="Optional CSV/TSV with cell id and cluster columns.",
    )
    velocity_parser.add_argument(
        "--cluster-key",
        default="clusters",
        help="obs column name for clusters when --clusters is provided.",
    )
    velocity_parser.add_argument(
        "--cluster-column",
        help="Column name in --clusters to use as cluster labels.",
    )
    velocity_parser.add_argument(
        "--basis",
        default="umap",
        help="Embedding basis for velocity plots, for example umap or space.",
    )
    velocity_parser.add_argument(
        "--color",
        help="obs column or variable used to color plots.",
    )
    velocity_parser.add_argument(
        "--mode",
        default="steady_state",
        choices=["steady_state", "deterministic", "stochastic", "dynamical"],
        help="scVelo velocity mode. steady_state is treated as deterministic.",
    )
    velocity_parser.add_argument(
        "--min-shared-counts",
        type=int,
        default=20,
        help="Minimum shared counts for scv.pp.filter_and_normalize.",
    )
    velocity_parser.add_argument(
        "--n-top-genes",
        type=int,
        help="Optional number of top genes for scv.pp.filter_and_normalize.",
    )
    velocity_parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCs for scv.pp.moments.",
    )
    velocity_parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of neighbors for scv.pp.moments.",
    )
    velocity_parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of jobs for scVelo steps that support parallel execution.",
    )
    velocity_parser.add_argument(
        "--qualityfilter",
        default="25:0.75",
        help="ASTRO-style feature filter, for example 25:0.75 or 0:0.",
    )
    velocity_parser.add_argument(
        "--barcode-tag",
        default="CB",
        help="BAM tag containing the corrected cell or spot barcode.",
    )
    velocity_parser.add_argument(
        "--umi-tag",
        default="UB",
        help="BAM tag containing the corrected UMI.",
    )
    velocity_parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep ASTRO intermediate BED and read map files.",
    )
    velocity_parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip PNG plot generation.",
    )
    velocity_parser.set_defaults(func=run_velocity)

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 1
    args.func(args)
    return 0


def run_velocity(args):
    libs = _load_velocity_libraries()
    np, pd, sc, scv, ad, sparse, plt = libs

    input_dir = Path(args.tenxinput).expanduser().resolve()
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    adata, source = _read_velocity_input(args, input_dir, output_dir, libs)
    adata = _add_optional_metadata(adata, args, libs)

    if "X_space" not in adata.obsm:
        _add_space_from_obs_names(adata, np)

    scv.settings.verbosity = 2
    filter_kwargs = {"min_shared_counts": args.min_shared_counts}
    if args.n_top_genes is not None:
        filter_kwargs["n_top_genes"] = args.n_top_genes
    scv.pp.filter_and_normalize(adata, **filter_kwargs)
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    basis_name = _basis_name(args.basis)
    basis_key = _basis_key(args.basis)
    if basis_name == "umap" and basis_key not in adata.obsm:
        sc.tl.umap(adata)

    scvelo_mode = "deterministic" if args.mode == "steady_state" else args.mode
    if scvelo_mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
    scv.tl.velocity(adata, mode=scvelo_mode)
    scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)

    for step_name, step_func in (
        ("latent_time", scv.tl.latent_time),
        ("velocity_pseudotime", scv.tl.velocity_pseudotime),
    ):
        try:
            step_func(adata)
        except Exception as exc:
            _warn(f"Skipping {step_name}: {exc}")

    outputs = []
    h5ad_path = output_dir / "rna_velocity.h5ad"
    adata.write(h5ad_path)
    outputs.append(str(h5ad_path))

    obs_path = _write_obs_table(adata, args, output_dir, pd)
    if obs_path is not None:
        outputs.append(str(obs_path))

    rank_path = _write_ranked_genes(adata, args, output_dir, pd, scv)
    if rank_path is not None:
        outputs.append(str(rank_path))

    if not args.no_plots:
        outputs.extend(_write_velocity_plots(adata, args, output_dir, scv, plt))

    summary = {
        "source": source,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "mode": args.mode,
        "scvelo_mode": scvelo_mode,
        "basis": basis_name,
        "outputs": outputs,
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    outputs.append(str(summary_path))

    print(f"Wrote RNA velocity results to {output_dir}")


def _load_velocity_libraries():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import scvelo as scv
        from scipy import sparse
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(
            "ASTROutils velocity requires the velocity extra dependencies. "
            "Install them with: pip install '.[velocity]'. "
            f"Dependency import failed with: {exc}"
        ) from exc
    return np, pd, sc, scv, ad, sparse, plt


def _read_velocity_input(args, input_dir, output_dir, libs):
    np, pd, sc, scv, ad, sparse, _plt = libs

    loom_path = Path(args.loom).expanduser().resolve() if args.loom else _find_loom(input_dir)
    if loom_path is not None:
        adata = scv.read(str(loom_path), cache=False)
        _require_velocity_layers(adata, f"loom file {loom_path}")
        return adata, str(loom_path)

    exon_matrix, intron_matrix = _find_exon_intron_tables(args, input_dir)
    if exon_matrix is not None and intron_matrix is not None:
        adata = _read_exon_intron_tables(exon_matrix, intron_matrix, libs)
        return adata, f"{exon_matrix} and {intron_matrix}"

    bam_path = _find_10x_bam(args, input_dir)
    if args.gtf and bam_path is not None:
        return _read_10x_bam_with_astro(args, input_dir, output_dir, bam_path, libs)

    matrix_dir = _find_10x_matrix_dir(input_dir)
    if matrix_dir is None:
        raise SystemExit(
            "Could not find a loom file, exon/intron TSV files, 10x BAM, or a 10x matrix "
            f"under {input_dir}."
        )
    return _read_10x_velocity_matrix(matrix_dir, libs)


def _find_loom(input_dir):
    loom_files = sorted(input_dir.rglob("*.loom")) if input_dir.exists() else []
    if len(loom_files) > 1:
        raise SystemExit("Multiple loom files found. Please select one with --loom.")
    return loom_files[0] if loom_files else None


def _find_exon_intron_tables(args, input_dir):
    if args.exon_matrix or args.intron_matrix:
        if not args.exon_matrix or not args.intron_matrix:
            raise SystemExit("--exon-matrix and --intron-matrix must be provided together.")
        return (
            Path(args.exon_matrix).expanduser().resolve(),
            Path(args.intron_matrix).expanduser().resolve(),
        )

    pairs = (
        ("fixed_exon.tsv", "fixed_intron.tsv"),
        ("exon.tsv", "intron.tsv"),
        ("spliced.tsv", "unspliced.tsv"),
        ("spliced.csv", "unspliced.csv"),
    )
    folders = (input_dir, input_dir / "outs")
    for folder in folders:
        for exon_name, intron_name in pairs:
            exon_path = folder / exon_name
            intron_path = folder / intron_name
            if exon_path.exists() and intron_path.exists():
                return exon_path, intron_path
    return None, None


def _read_exon_intron_tables(exon_matrix, intron_matrix, libs):
    np, pd, _sc, ad, sparse = libs[0], libs[1], libs[2], libs[4], libs[5]
    exon = pd.read_csv(exon_matrix, sep=_table_sep(exon_matrix), index_col=0)
    intron = pd.read_csv(intron_matrix, sep=_table_sep(intron_matrix), index_col=0)

    common_cells = exon.index.intersection(intron.index)
    common_genes = exon.columns.intersection(intron.columns)
    if common_cells.empty or common_genes.empty:
        raise SystemExit("The exon and intron matrices do not share cells and genes.")

    exon = exon.loc[common_cells, common_genes]
    intron = intron.loc[common_cells, common_genes]

    spliced = sparse.csr_matrix(exon.to_numpy())
    unspliced = sparse.csr_matrix(intron.to_numpy())
    adata = ad.AnnData(X=spliced)
    adata.obs_names = common_cells.astype(str)
    adata.var_names = common_genes.astype(str)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    return adata


def _find_10x_matrix_dir(input_dir):
    candidates = (
        input_dir,
        input_dir / "filtered_feature_bc_matrix",
        input_dir / "raw_feature_bc_matrix",
        input_dir / "outs" / "filtered_feature_bc_matrix",
        input_dir / "outs" / "raw_feature_bc_matrix",
        input_dir / "STAR",
    )
    for candidate in candidates:
        if _is_10x_matrix_dir(candidate):
            return candidate

    if input_dir.exists():
        for matrix_path in sorted(input_dir.rglob("matrix.mtx*")):
            candidate = matrix_path.parent
            if _is_10x_matrix_dir(candidate):
                return candidate
    return None


def _find_10x_bam(args, input_dir):
    if args.bam:
        bam_path = Path(args.bam).expanduser().resolve()
        if not bam_path.exists():
            raise SystemExit(f"BAM file does not exist: {bam_path}")
        return bam_path

    candidates = (
        input_dir / "possorted_genome_bam.bam",
        input_dir / "outs" / "possorted_genome_bam.bam",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _read_10x_bam_with_astro(args, input_dir, output_dir, bam_path, libs):
    np, pd, _sc, _scv, ad, sparse, _plt = libs

    gtf_path = Path(args.gtf).expanduser().resolve()
    if not gtf_path.exists():
        raise SystemExit(f"GTF file does not exist: {gtf_path}")

    positions_path = _find_positions_file(args, input_dir)
    if positions_path is None:
        raise SystemExit(
            "A 10x BAM needs spot coordinates. Provide --positions or put "
            "tissue_positions_list.csv(.gz), tissue_positions.csv(.gz), or "
            "spot_meta.tsv under --10xinput."
        )

    count_dir = output_dir / "astro_velocity_count"
    count_dir.mkdir(parents=True, exist_ok=True)

    barcode_to_xy, ordered_locs, loc_to_xy = _load_10x_spot_coordinates(positions_path, libs)
    barcodes_file = count_dir / "astro_barcodes.tsv"
    _write_astro_barcodes(barcodes_file, barcode_to_xy)

    astro_bed = count_dir / "10x_astro.bed"
    read_map = count_dir / "read_map.tsv"
    stats = _write_astro_bed_from_10x_bam(
        bam_path=bam_path,
        output_bed=astro_bed,
        read_map=read_map,
        barcode_to_xy=barcode_to_xy,
        args=args,
    )
    if stats["written_reads"] == 0:
        raise SystemExit(
            "No mapped BAM records could be matched to tissue positions. Check "
            f"the BAM barcode tag ({args.barcode_tag}) and the positions file."
        )

    from .countfeature import inter_bed2geneFile

    expmat_bed = count_dir / "expmat.bed"
    inter_bed2geneFile(str(astro_bed), str(expmat_bed), str(gtf_path))

    adata, matrix_stats = _adata_from_astro_expmat(
        expmat_bed=expmat_bed,
        read_map=read_map,
        output_dir=output_dir,
        ordered_locs=ordered_locs,
        loc_to_xy=loc_to_xy,
        args=args,
        libs=libs,
    )
    adata.uns["astro_velocity_count"] = {
        "bam": str(bam_path),
        "gtf": str(gtf_path),
        "positions": str(positions_path),
        "bam_records": stats,
        "matrix": matrix_stats,
    }

    if not args.keep_temp:
        for path in (astro_bed, read_map):
            try:
                path.unlink()
            except OSError:
                pass

    return adata, f"ASTRO count from {bam_path} with {gtf_path}"


def _find_positions_file(args, input_dir):
    if args.positions:
        positions_path = Path(args.positions).expanduser().resolve()
        if not positions_path.exists():
            raise SystemExit(f"Positions file does not exist: {positions_path}")
        return positions_path

    candidates = (
        input_dir / "tissue_positions_list.csv.gz",
        input_dir / "tissue_positions_list.csv",
        input_dir / "tissue_positions.csv.gz",
        input_dir / "tissue_positions.csv",
        input_dir / "spatial" / "tissue_positions_list.csv.gz",
        input_dir / "spatial" / "tissue_positions_list.csv",
        input_dir / "spatial" / "tissue_positions.csv.gz",
        input_dir / "spatial" / "tissue_positions.csv",
        input_dir / "spot_meta.tsv",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _load_10x_spot_coordinates(positions_path, libs):
    _np, pd = libs[0], libs[1]
    sep = "\t" if str(positions_path).endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    df = pd.read_csv(positions_path, sep=sep)
    if "barcode" not in df.columns and df.shape[1] >= 6:
        df = pd.read_csv(
            positions_path,
            sep=sep,
            header=None,
            names=[
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ],
        )

    if "barcode" not in df.columns:
        df = df.rename(columns={df.columns[0]: "barcode"})

    if {"array_col", "array_row"}.issubset(df.columns):
        x_col, y_col = "array_col", "array_row"
    elif {"pxl_col_in_fullres", "pxl_row_in_fullres"}.issubset(df.columns):
        x_col, y_col = "pxl_col_in_fullres", "pxl_row_in_fullres"
    elif df.shape[1] >= 3:
        x_col, y_col = df.columns[1], df.columns[2]
    else:
        raise SystemExit(f"Could not find coordinate columns in {positions_path}.")

    if "in_tissue" in df.columns:
        df = df[df["in_tissue"].astype(str).isin(("1", "True", "true", "TRUE"))]

    barcode_to_xy = {}
    ordered_locs = []
    loc_to_xy = {}
    for _index, row in df.iterrows():
        barcode = str(row["barcode"])
        x = str(row[x_col])
        y = str(row[y_col])
        loc = f"{x}_{y}"
        barcode_to_xy[barcode] = (x, y)
        if loc not in loc_to_xy:
            ordered_locs.append(loc)
            loc_to_xy[loc] = (x, y)

    if not barcode_to_xy:
        raise SystemExit(f"No tissue positions were loaded from {positions_path}.")
    return barcode_to_xy, ordered_locs, loc_to_xy


def _write_astro_barcodes(barcodes_file, barcode_to_xy):
    with barcodes_file.open("w") as fout:
        for barcode, (x, y) in barcode_to_xy.items():
            fout.write(f"{barcode}\t{x}\t{y}\n")


def _write_astro_bed_from_10x_bam(bam_path, output_bed, read_map, barcode_to_xy, args):
    stats = {
        "sam_records": 0,
        "written_reads": 0,
        "written_bed_blocks": 0,
        "missing_barcode": 0,
        "barcode_not_in_positions": 0,
        "missing_umi": 0,
        "missing_blocks": 0,
    }
    cmd = ["samtools", "view", "-F", "2308", str(bam_path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    read_id = 0
    with output_bed.open("w") as bed_out, read_map.open("w") as map_out:
        for line in proc.stdout:
            stats["sam_records"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            tags = _parse_sam_tags(fields[11:])
            barcode = tags.get(args.barcode_tag) or tags.get("CB") or tags.get("CR")
            if not barcode:
                stats["missing_barcode"] += 1
                continue
            barcode = _normalize_barcode(barcode, barcode_to_xy)
            if barcode is None:
                stats["barcode_not_in_positions"] += 1
                continue

            umi = tags.get(args.umi_tag) or tags.get("UB") or tags.get("UR")
            if not umi:
                stats["missing_umi"] += 1
                umi = fields[0]

            chrom = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            blocks = _cigar_to_bed_blocks(pos, cigar)
            if not blocks:
                stats["missing_blocks"] += 1
                continue

            read_id += 1
            x, y = barcode_to_xy[barcode]
            loc = f"{x}_{y}"
            flag = int(fields[1])
            strand = "-" if flag & 16 else "+"
            score = tags.get("AS", ".")
            read_length = _read_length(fields[9], cigar)
            read_name = f"{read_id}_{x}_{y}:{barcode}:{score}:{read_length}"

            map_out.write(f"{read_id}\t{loc}\t{barcode}\t{umi}\n")
            for start, end in blocks:
                bed_out.write(f"{chrom}\t{start}\t{end}\t{read_name}\t0\t{strand}\n")
                stats["written_bed_blocks"] += 1
            stats["written_reads"] += 1

    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"samtools view failed on {bam_path}")
    return stats


def _parse_sam_tags(tag_fields):
    tags = {}
    for field in tag_fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def _normalize_barcode(barcode, barcode_to_xy):
    if barcode in barcode_to_xy:
        return barcode
    stripped = barcode.split("-", 1)[0]
    if stripped in barcode_to_xy:
        return stripped
    with_suffix = f"{stripped}-1"
    if with_suffix in barcode_to_xy:
        return with_suffix
    return None


def _cigar_to_bed_blocks(pos, cigar):
    if cigar == "*":
        return []

    blocks = []
    ref_pos = pos - 1
    block_start = None
    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in ("M", "D", "=", "X"):
            if block_start is None:
                block_start = ref_pos
            ref_pos += length
        elif op == "N":
            if block_start is not None and ref_pos > block_start:
                blocks.append((block_start, ref_pos))
            ref_pos += length
            block_start = None
        elif op in ("I", "S", "H", "P"):
            continue

    if block_start is not None and ref_pos > block_start:
        blocks.append((block_start, ref_pos))
    return blocks


def _read_length(sequence, cigar):
    if sequence and sequence != "*":
        return len(sequence)

    length = 0
    for length_text, op in CIGAR_RE.findall(cigar):
        if op in ("M", "I", "S", "=", "X"):
            length += int(length_text)
    return length or "."


def _adata_from_astro_expmat(expmat_bed, read_map, output_dir, ordered_locs, loc_to_xy, args, libs):
    np, pd, _sc, _scv, ad, sparse, _plt = libs
    read_info = _load_read_map(read_map)
    filter_settings = _parse_qualityfilter(args.qualityfilter)

    counts = {
        "spliced": defaultdict(int),
        "unspliced": defaultdict(int),
    }
    seen_umis = set()
    stats = {
        "feature_lines": 0,
        "counted_umis": 0,
        "missing_read_map": 0,
        "missing_velocity_suffix": 0,
        "quality_filtered": 0,
    }

    with expmat_bed.open() as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            stats["feature_lines"] += 1
            parsed = _parse_astro_feature_line(line, filter_settings)
            if parsed is None:
                stats["quality_filtered"] += 1
                continue
            read_id, gene_feature = parsed
            if read_id not in read_info:
                stats["missing_read_map"] += 1
                continue

            gene_name, layer = _split_astro_velocity_feature(gene_feature)
            if layer is None:
                stats["missing_velocity_suffix"] += 1
                continue

            loc, barcode, umi = read_info[read_id]
            umi_key = (layer, gene_name, loc, barcode, umi)
            if umi_key in seen_umis:
                continue
            seen_umis.add(umi_key)
            counts[layer][(loc, gene_name)] += 1
            stats["counted_umis"] += 1

    genes = sorted(
        {gene for _loc, gene in counts["spliced"]}
        | {gene for _loc, gene in counts["unspliced"]}
    )
    used_locs = [
        loc
        for loc in ordered_locs
        if any((loc, gene) in counts["spliced"] or (loc, gene) in counts["unspliced"] for gene in genes)
    ]
    if not used_locs or not genes:
        raise SystemExit("ASTRO overlap completed, but no exon/intron UMI counts were produced.")

    obs_names = [loc.replace("_", "x") for loc in used_locs]
    exon = pd.DataFrame(0, index=obs_names, columns=genes, dtype=int)
    intron = pd.DataFrame(0, index=obs_names, columns=genes, dtype=int)
    loc_to_obs = dict(zip(used_locs, obs_names))
    for (loc, gene), value in counts["spliced"].items():
        if loc in loc_to_obs:
            exon.at[loc_to_obs[loc], gene] = value
    for (loc, gene), value in counts["unspliced"].items():
        if loc in loc_to_obs:
            intron.at[loc_to_obs[loc], gene] = value

    exon_path = output_dir / "fixed_exon.tsv"
    intron_path = output_dir / "fixed_intron.tsv"
    exon.to_csv(exon_path, sep="\t")
    intron.to_csv(intron_path, sep="\t")

    spliced = sparse.csr_matrix(exon.to_numpy())
    unspliced = sparse.csr_matrix(intron.to_numpy())
    adata = ad.AnnData(X=spliced)
    adata.obs_names = exon.index.astype(str)
    adata.var_names = exon.columns.astype(str)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.obsm["X_space"] = np.asarray(
        [[float(loc_to_xy[loc][0]), float(loc_to_xy[loc][1])] for loc in used_locs],
        dtype=float,
    )
    stats["cells"] = len(used_locs)
    stats["genes"] = len(genes)
    stats["fixed_exon"] = str(exon_path)
    stats["fixed_intron"] = str(intron_path)
    return adata, stats


def _load_read_map(read_map):
    read_info = {}
    with read_map.open() as fin:
        for line in fin:
            read_id, loc, barcode, umi = line.rstrip("\n").split("\t", 3)
            read_info[read_id] = (loc, barcode, umi)
    return read_info


def _parse_qualityfilter(filter_str):
    if not filter_str or filter_str == "NA":
        return None
    parts = filter_str.split(":")
    if len(parts) != 2:
        return None
    try:
        as_threshold = float(parts[0])
        ratio_threshold = float(parts[1])
    except ValueError:
        return None
    if as_threshold == 0 and ratio_threshold == 0:
        return None
    return as_threshold, ratio_threshold


def _parse_astro_feature_line(line, filter_settings):
    fields = line.split("\t")
    if len(fields) < 9:
        return None
    read_id = fields[8]
    fields.pop(8)

    i = 0
    while i + 7 < len(fields):
        gene_feature = fields[i + 2]
        as_score = fields[i + 3]
        gene_length = fields[i + 6]
        i += 8

        if not _passes_quality(as_score, gene_length, filter_settings):
            continue
        return read_id, gene_feature
    return None


def _passes_quality(as_score, gene_length, filter_settings):
    if filter_settings is None or as_score == ".":
        return True
    try:
        as_value = float(as_score)
        gene_length_value = float(gene_length)
    except ValueError:
        return True
    as_threshold, ratio_threshold = filter_settings
    return not (
        as_value <= as_threshold
        and as_value <= ratio_threshold * gene_length_value
    )


def _split_astro_velocity_feature(feature_name):
    feature_name = str(feature_name)
    match = VELOCITY_SUFFIX_RE.search(feature_name)
    if match is None:
        return feature_name, None

    gene_name = feature_name.split("__", 1)[0]
    if match.group(1) == "exon":
        return gene_name, "spliced"
    return gene_name, "unspliced"


def _is_10x_matrix_dir(folder):
    return (
        folder.is_dir()
        and _has_any(folder, ("matrix.mtx", "matrix.mtx.gz"))
        and _has_any(folder, ("barcodes.tsv", "barcodes.tsv.gz"))
        and _has_any(folder, ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))
    )


def _has_any(folder, names):
    return any((folder / name).exists() for name in names)


def _read_10x_velocity_matrix(matrix_dir, libs):
    np, _pd, sc, _scv, ad, sparse, _plt = libs

    try:
        source = sc.read_10x_mtx(
            str(matrix_dir),
            var_names="gene_symbols",
            make_unique=False,
            gex_only=False,
        )
    except TypeError:
        source = sc.read_10x_mtx(
            str(matrix_dir),
            var_names="gene_symbols",
            make_unique=False,
        )

    feature_names = _feature_names_with_velocity_suffixes(source)
    layer_names = []
    gene_names = []
    for feature_name in feature_names:
        gene_name, layer_name = _split_velocity_feature(feature_name)
        gene_names.append(gene_name)
        layer_names.append(layer_name)

    layer_names = np.asarray(layer_names, dtype=object)
    gene_names = np.asarray(gene_names, dtype=object)
    spliced_mask = layer_names == "spliced"
    unspliced_mask = layer_names == "unspliced"

    if not spliced_mask.any() or not unspliced_mask.any():
        raise SystemExit(
            "No exon/intron features were found in the 10x matrix. Standard "
            "Cell Ranger matrices contain total counts only and cannot produce "
            "RNA velocity directly. Provide --gtf so ASTROutils can count "
            "exon/intron features from possorted_genome_bam.bam, pass "
            "--exon-matrix and --intron-matrix, or use a matrix with __exon "
            "and __intron/__transcript feature suffixes."
        )

    matrix = sparse.csr_matrix(source.X)
    spliced, spliced_genes = _collapse_columns(matrix, gene_names, spliced_mask, np, sparse)
    unspliced, unspliced_genes = _collapse_columns(matrix, gene_names, unspliced_mask, np, sparse)

    common_genes, spliced_index, unspliced_index = np.intersect1d(
        spliced_genes,
        unspliced_genes,
        return_indices=True,
    )
    if common_genes.size == 0:
        raise SystemExit("The exon and intron features do not share gene names.")

    spliced = spliced[:, spliced_index].tocsr()
    unspliced = unspliced[:, unspliced_index].tocsr()

    adata = ad.AnnData(X=spliced, obs=source.obs.copy())
    adata.obs_names = source.obs_names.astype(str)
    adata.var_names = common_genes.astype(str)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    return adata, str(matrix_dir)


def _feature_names_with_velocity_suffixes(adata):
    candidates = [adata.var_names.astype(str)]
    if "gene_ids" in adata.var:
        candidates.append(adata.var["gene_ids"].astype(str).to_numpy())

    for candidate in candidates:
        if any(VELOCITY_SUFFIX_RE.search(str(name)) for name in candidate):
            return candidate
    return candidates[0]


def _split_velocity_feature(feature_name):
    match = VELOCITY_SUFFIX_RE.search(str(feature_name))
    if match is None:
        return str(feature_name), None

    gene_name = str(feature_name)[: match.start()]
    if match.group(1) == "exon":
        return gene_name, "spliced"
    return gene_name, "unspliced"


def _collapse_columns(matrix, gene_names, mask, np, sparse):
    selected_genes = gene_names[mask]
    selected_matrix = matrix[:, mask]
    unique_genes, inverse = np.unique(selected_genes, return_inverse=True)
    groups = sparse.csr_matrix(
        (
            np.ones(inverse.size, dtype=selected_matrix.dtype),
            (inverse, np.arange(inverse.size)),
        ),
        shape=(unique_genes.size, inverse.size),
    )
    return selected_matrix @ groups.T, unique_genes


def _require_velocity_layers(adata, source):
    missing = [layer for layer in ("spliced", "unspliced") if layer not in adata.layers]
    if missing:
        raise SystemExit(f"{source} is missing required layers: {', '.join(missing)}")


def _add_optional_metadata(adata, args, libs):
    np, pd = libs[0], libs[1]

    if args.umap:
        umap = pd.read_csv(args.umap, sep=_table_sep(args.umap))
        if umap.shape[1] < 3:
            raise SystemExit("--umap must contain a cell id column and two coordinate columns.")
        cell_col = umap.columns[0]
        umap[cell_col] = umap[cell_col].astype(str)
        coord_cols = ["UMAP_1", "UMAP_2"] if {"UMAP_1", "UMAP_2"}.issubset(umap.columns) else list(umap.columns[1:3])
        adata = _subset_to_metadata(adata, umap, cell_col)
        umap = umap.set_index(cell_col).reindex(adata.obs_names)
        adata.obsm["X_umap"] = umap[coord_cols].to_numpy(dtype=float)

    if args.clusters:
        clusters = pd.read_csv(args.clusters, sep=_table_sep(args.clusters))
        if clusters.shape[1] < 2:
            raise SystemExit("--clusters must contain a cell id column and a cluster column.")
        cell_col = clusters.columns[0]
        clusters[cell_col] = clusters[cell_col].astype(str)
        cluster_col = args.cluster_column or clusters.columns[1]
        adata = _subset_to_metadata(adata, clusters, cell_col)
        clusters = clusters.set_index(cell_col).reindex(adata.obs_names)
        adata.obs[args.cluster_key] = clusters[cluster_col].astype(str).to_numpy()

    return adata


def _subset_to_metadata(adata, metadata, cell_col):
    metadata_cells = set(metadata[cell_col].astype(str))
    common_cells = [cell for cell in adata.obs_names if cell in metadata_cells]
    if not common_cells:
        raise SystemExit(f"No cells overlap with metadata column {cell_col}.")
    return adata[common_cells].copy()


def _add_space_from_obs_names(adata, np):
    coords = []
    for obs_name in adata.obs_names:
        match = re.search(
            r"^(?:x)?(-?\d+(?:\.\d+)?)x(-?\d+(?:\.\d+)?)$",
            str(obs_name),
        )
        if match is None:
            return
        coords.append((float(match.group(1)), float(match.group(2))))
    adata.obsm["X_space"] = np.asarray(coords, dtype=float)


def _table_sep(path):
    path = str(path)
    if path.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "\t"
    return ","


def _basis_name(basis):
    return basis[2:] if basis.startswith("X_") else basis


def _basis_key(basis):
    return basis if basis.startswith("X_") else f"X_{basis}"


def _write_obs_table(adata, args, output_dir, pd):
    columns = []
    for column in (args.cluster_key, "velocity_pseudotime", "latent_time"):
        if column in adata.obs:
            columns.append(column)
    if not columns:
        return None

    output_path = output_dir / "velocity_obs.csv"
    adata.obs[columns].to_csv(output_path)
    return output_path


def _write_ranked_genes(adata, args, output_dir, pd, scv):
    if args.cluster_key not in adata.obs:
        return None

    try:
        scv.tl.rank_velocity_genes(adata, groupby=args.cluster_key)
    except Exception as exc:
        _warn(f"Skipping rank_velocity_genes: {exc}")
        return None

    output_path = output_dir / "rank_velocity_genes.csv"
    pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]).to_csv(output_path, index=False)
    return output_path


def _write_velocity_plots(adata, args, output_dir, scv, plt):
    basis_name = _basis_name(args.basis)
    basis_key = _basis_key(args.basis)
    if basis_key not in adata.obsm:
        _warn(f"Skipping plots because {basis_key} is not available.")
        return []

    color = args.color
    if color is None:
        color = args.cluster_key if args.cluster_key in adata.obs else "velocity_pseudotime"

    plot_dir = output_dir / "figures"
    plot_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    plot_specs = (
        (
            "velocity_stream",
            lambda: scv.pl.velocity_embedding_stream(
                adata,
                basis=basis_name,
                color=color,
                show=False,
                title="",
            ),
        ),
        (
            "velocity_grid",
            lambda: scv.pl.velocity_embedding_grid(
                adata,
                basis=basis_name,
                color=color,
                show=False,
                title="",
            ),
        ),
    )
    if "velocity_pseudotime" in adata.obs:
        plot_specs = plot_specs + (
            (
                "velocity_pseudotime",
                lambda: scv.pl.scatter(
                    adata,
                    basis=basis_name,
                    color="velocity_pseudotime",
                    show=False,
                    title="",
                ),
            ),
        )

    for name, plot_func in plot_specs:
        try:
            plot_func()
            fig = plt.gcf()
            output_path = plot_dir / f"{name}_{basis_name}.png"
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            outputs.append(str(output_path))
        except Exception as exc:
            plt.close("all")
            _warn(f"Skipping {name}: {exc}")

    return outputs


def _warn(message):
    print(f"Warning: {message}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
