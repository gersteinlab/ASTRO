#!/usr/bin/env python
import argparse
import gzip
import json
import multiprocessing
import os
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
        help="Run or resume RNA velocity analysis.",
    )
    velocity_parser.add_argument(
        "--10xinput",
        dest="tenxinput",
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
            "or Cell Ranger BAM when the 10x matrix does not already contain "
            "velocity layers."
        ),
    )
    velocity_parser.add_argument(
        "--bam",
        help=(
            "Optional 10x/Cell Ranger BAM path. ASTROutils reads barcode and "
            "UMI tags from this BAM and converts it to ASTRO counting BED "
            "internally. Defaults to possorted_genome_bam.bam under --10xinput."
        ),
    )
    velocity_parser.add_argument(
        "--bam-format",
        choices=["10x", "astro"],
        default="10x",
        help=(
            "Input BAM format. Use 10x for Cell Ranger BAMs with barcode/UMI "
            "tags. ASTRO BAM input is not converted by this velocity command."
        ),
    )
    velocity_parser.add_argument(
        "--positions",
        help="Optional Visium tissue positions file. Defaults to tissue_positions*.csv under --10xinput.",
    )
    velocity_parser.add_argument(
        "--barcode-whitelist",
        help=(
            "Optional filtered barcode whitelist. Supports one-column barcode "
            "files or two-column species,barcode CSV files."
        ),
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
        help="Number of jobs for parallel overlap and scVelo steps that support it.",
    )
    velocity_parser.add_argument(
        "--overlap-chunk-lines",
        type=int,
        default=2000000,
        help="Target BED lines per parallel overlap chunk when --n-jobs is greater than 1.",
    )
    velocity_parser.add_argument(
        "--qualityfilter",
        default="25:0.75",
        help="ASTRO-style feature filter, for example 25:0.75 or 0:0.",
    )
    velocity_parser.add_argument(
        "--progress-interval",
        type=int,
        default=1000000,
        help="Print progress every N records while reading large intermediate files.",
    )
    velocity_parser.add_argument(
        "--dedup-buckets",
        type=int,
        default=512,
        help="Number of on-disk hash buckets used for memory-bounded UMI deduplication.",
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
    velocity_parser.add_argument(
        "--no-notebook",
        action="store_true",
        help="Skip writing the output inspection notebook.",
    )
    velocity_parser.add_argument(
        "--no-force-exit",
        action="store_true",
        help="Do not force the CLI process to exit after all outputs are written.",
    )
    velocity_parser.set_defaults(func=run_velocity)

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 1
    args.func(args)
    if args.command == "velocity" and not getattr(args, "no_force_exit", False):
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(0)
    return 0


def run_velocity(args):
    input_dir = Path(args.tenxinput).expanduser().resolve() if args.tenxinput else Path(".").resolve()
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if _resume_velocity_from_output(args, input_dir, output_dir, required=False):
        return

    if not args.tenxinput and not (args.bam or args.loom or (args.exon_matrix and args.intron_matrix)):
        raise SystemExit(
            "Provide --10xinput, or use --bam with --gtf, --loom, or "
            "--exon-matrix together with --intron-matrix."
        )

    libs = _load_velocity_libraries()
    adata, source = _read_velocity_input(args, input_dir, output_dir, libs)
    _run_velocity_analysis(adata, args, source, output_dir, libs)


def _resume_velocity_from_output(args, input_dir, output_dir, required=False):
    h5ad_path = output_dir / "rna_velocity.h5ad"
    command_name = getattr(args, "command", "velocity")

    if h5ad_path.exists():
        if args.barcode_whitelist:
            _warn(
                "Existing rna_velocity.h5ad found; --barcode-whitelist is ignored "
                "unless the h5ad is removed and the command rebuilds from expmat.bed."
            )
        source = "existing RNA velocity output"
        summary_path = output_dir / "run_summary.json"
        if summary_path.exists():
            try:
                summary = json.loads(summary_path.read_text())
                source = summary.get("source", source)
            except json.JSONDecodeError:
                pass
        notebook_path = _write_velocity_notebook(output_dir, source, args)
        print(f"Wrote RNA velocity report notebook to {notebook_path}")
        return True

    expmat_bed = _existing_path_with_gz(output_dir / "expmat.bed")
    expmat_done = output_dir / "expmat.bed.done"
    read_map = _existing_path_with_gz(output_dir / "read_map.tsv")

    if expmat_bed is None and read_map is None:
        if required:
            raise SystemExit(
                f"{command_name} requires an existing RNA velocity h5ad or "
                f"ASTRO feature BED at {output_dir / 'expmat.bed'}"
            )
        return False
    if expmat_bed is None:
        if not required:
            return False
        raise SystemExit(
            f"{command_name} resume requires an existing ASTRO feature BED at "
            f"{output_dir / 'expmat.bed'}"
        )
    if not expmat_done.exists():
        _warn(f"Found {expmat_bed} without {expmat_done.name}; rerunning overlap.")
        return False
    if read_map is None:
        raise SystemExit(
            f"{command_name} resume requires the read-to-UMI map from velocity counting at "
            f"{output_dir / 'read_map.tsv'}"
        )

    _log_status(f"{command_name}: loading scVelo/Scanpy libraries")
    libs = _load_velocity_libraries()
    _log_status(f"{command_name}: loading spot coordinates")
    positions_path = _find_positions_file(args, input_dir)
    if positions_path is not None:
        _barcode_to_xy, ordered_locs, loc_to_xy = _load_10x_spot_coordinates(positions_path, libs)
    else:
        barcodes_file = _existing_path_with_gz(output_dir / "astro_barcodes.tsv")
        if barcodes_file is None:
            raise SystemExit(
                f"{command_name} needs spot coordinates. Provide --positions, put a "
                "10x tissue_positions file under --10xinput, or keep "
                f"{output_dir / 'astro_barcodes.tsv'} from the first velocity run."
        )
        ordered_locs, loc_to_xy = _load_astro_barcode_locations(barcodes_file)

    _log_status(f"{command_name}: building sparse AnnData from {expmat_bed}")
    adata, matrix_stats = _sparse_adata_from_astro_expmat(
        expmat_bed=expmat_bed,
        read_map=read_map,
        output_dir=output_dir,
        ordered_locs=ordered_locs,
        loc_to_xy=loc_to_xy,
        args=args,
        libs=libs,
    )
    _log_status(
        f"{command_name}: sparse matrix built "
        f"({matrix_stats['cells']} cells, {matrix_stats['genes']} genes, "
        f"{matrix_stats['counted_umis']} counted UMIs)"
    )
    _log_status(f"{command_name}: running scVelo analysis")
    _run_velocity_analysis(
        adata,
        args,
        f"sparse ASTRO matrix from {expmat_bed}",
        output_dir,
        libs,
        summary_extra={"sparse_matrix": matrix_stats},
    )
    return True


def _existing_path_with_gz(path):
    if path.exists():
        return path
    gz_path = Path(f"{path}.gz")
    if gz_path.exists():
        return gz_path
    return None


def _open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def _run_velocity_analysis(adata, args, source, output_dir, libs, summary_extra=None):
    np, pd, sc, scv, ad, sparse, plt = libs
    adata = _add_optional_metadata(adata, args, libs)

    if "X_space" not in adata.obsm:
        _add_space_from_obs_names(adata, np)

    scv.settings.verbosity = 2
    filter_kwargs = {"min_shared_counts": args.min_shared_counts}
    if args.n_top_genes is not None:
        filter_kwargs["n_top_genes"] = args.n_top_genes
    _log_status("scVelo: filter_and_normalize")
    scv.pp.filter_and_normalize(adata, **filter_kwargs)
    _log_status("scVelo: moments")
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    basis_name = _basis_name(args.basis)
    basis_key = _basis_key(args.basis)
    if basis_name == "umap" and basis_key not in adata.obsm:
        sc.tl.umap(adata)

    scvelo_mode = "deterministic" if args.mode == "steady_state" else args.mode
    if scvelo_mode == "dynamical":
        _log_status("scVelo: recover_dynamics")
        scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
    _log_status(f"scVelo: velocity mode={scvelo_mode}")
    scv.tl.velocity(adata, mode=scvelo_mode)
    _log_status("scVelo: velocity_graph")
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

    if not getattr(args, "no_notebook", False):
        notebook_path = _write_velocity_notebook(output_dir, source, args)
        outputs.append(str(notebook_path))

    summary = {
        "source": source,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "mode": args.mode,
        "scvelo_mode": scvelo_mode,
        "basis": basis_name,
        "outputs": outputs,
    }
    if summary_extra:
        summary.update(summary_extra)
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
            "Install them in the same environment with: "
            "python -m pip install -e './python[velocity]'. "
            f"Dependency import failed with: {exc}"
        ) from exc
    return np, pd, sc, scv, ad, sparse, plt


def _read_velocity_input(args, input_dir, output_dir, libs):
    np, pd, sc, scv, ad, sparse, _plt = libs

    loom_path = Path(args.loom).expanduser().resolve() if args.loom else None
    if loom_path is not None:
        adata = scv.read(str(loom_path), cache=False)
        _require_velocity_layers(adata, f"loom file {loom_path}")
        return adata, str(loom_path)

    exon_matrix, intron_matrix = _find_exon_intron_tables(args, input_dir) if (args.exon_matrix or args.intron_matrix) else (None, None)
    if exon_matrix is not None and intron_matrix is not None:
        adata = _read_exon_intron_tables(exon_matrix, intron_matrix, libs)
        return adata, f"{exon_matrix} and {intron_matrix}"

    if args.bam:
        if not args.gtf:
            raise SystemExit("--bam input requires --gtf for ASTRO exon/intron feature counting.")
        if args.bam_format != "10x":
            raise SystemExit(
                "ASTROutils velocity --bam currently supports --bam-format 10x "
                "only. Use --bam-format 10x for Cell Ranger BAMs, or resume "
                "from existing ASTRO velocity count intermediates with "
                "ASTROutils velocity."
            )
        return _read_10x_bam_with_astro(args, input_dir, output_dir, _find_10x_bam(args, input_dir), libs)

    loom_path = _find_loom(input_dir)
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
        if args.bam_format != "10x":
            raise SystemExit(
                "The BAM found under --10xinput is handled as a 10x/Cell "
                "Ranger BAM. Use --bam-format 10x for this path."
            )
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

    _log_status(
        "velocity: treating BAM as 10x/Cell Ranger BAM and converting "
        "barcode/UMI tags to ASTRO counting BED"
    )
    barcode_whitelist = _load_barcode_whitelist(args.barcode_whitelist)
    if barcode_whitelist is not None:
        _log_status(f"velocity: loaded {len(barcode_whitelist)} whitelist barcodes")

    gtf_path = Path(args.gtf).expanduser().resolve()
    if not gtf_path.exists():
        raise SystemExit(f"GTF file does not exist: {gtf_path}")

    positions_path = _find_positions_file(args, input_dir)
    use_barcode_locs = positions_path is None
    if positions_path is None and not args.bam:
        raise SystemExit(
            "A 10x BAM needs spot coordinates. Provide --positions or put "
            "tissue_positions_list.csv(.gz), tissue_positions.csv(.gz), or "
            "spot_meta.tsv under --10xinput."
        )
    if positions_path is None:
        _warn("No spot positions found; treating BAM barcodes as single-cell observations.")

    gtf_chroms = _load_gtf_chroms(gtf_path)
    if positions_path is not None:
        barcode_to_xy, ordered_locs, loc_to_xy = _load_10x_spot_coordinates(positions_path, libs)
    else:
        barcode_to_xy, ordered_locs, loc_to_xy = {}, [], {}
    existing_barcodes_file = _existing_path_with_gz(output_dir / "astro_barcodes.tsv")
    barcodes_file = existing_barcodes_file or output_dir / "astro_barcodes.tsv"

    existing_astro_bed = _existing_path_with_gz(output_dir / "10x_astro.bed")
    existing_read_map = _existing_path_with_gz(output_dir / "read_map.tsv")
    astro_bed = existing_astro_bed or output_dir / "10x_astro.bed"
    gtf_chrom_bed = output_dir / "10x_astro.gtfchrom.bed"
    read_map = existing_read_map or output_dir / "read_map.tsv"
    reused_intermediate = existing_astro_bed is not None and existing_read_map is not None
    if reused_intermediate and use_barcode_locs:
        if existing_barcodes_file is not None:
            ordered_locs, loc_to_xy = _load_astro_barcode_locations(
                barcodes_file,
                use_barcode_as_loc=True,
            )
        else:
            reused_intermediate = False
    if reused_intermediate:
        _log_status(f"velocity: reusing existing ASTRO BED and read map under {output_dir}")
        stats = {
            "reused_intermediate": True,
            "bed": str(astro_bed),
            "read_map": str(read_map),
        }
    else:
        _log_status(f"velocity: converting BAM to ASTRO BED at {astro_bed}")
        stats = _write_astro_bed_from_10x_bam(
            bam_path=bam_path,
            output_bed=astro_bed,
            read_map=read_map,
            barcode_to_xy=barcode_to_xy,
            gtf_chroms=gtf_chroms,
            args=args,
            use_barcode_locs=use_barcode_locs,
            ordered_locs=ordered_locs,
            loc_to_xy=loc_to_xy,
            barcode_whitelist=barcode_whitelist,
        )
        _log_status(
            "velocity: BAM conversion completed "
            f"({stats['written_reads']} reads, {stats['written_bed_blocks']} BED blocks)"
        )
    if not reused_intermediate and stats["written_reads"] == 0:
        raise SystemExit(
            "No mapped BAM records could be used for velocity counting. Check "
            f"the BAM barcode tag ({args.barcode_tag})"
            + ("." if use_barcode_locs else " and the positions file.")
        )
    if existing_barcodes_file is None or not reused_intermediate:
        barcodes_file = output_dir / "astro_barcodes.tsv"
        _write_astro_barcodes(barcodes_file, barcode_to_xy)
    overlap_bed, chrom_stats = _prepare_astro_bed_for_gtf(
        astro_bed,
        gtf_chrom_bed,
        gtf_chroms,
    )
    stats["chromosome_fix"] = chrom_stats

    expmat_bed = output_dir / "expmat.bed"
    _log_status(f"velocity: running ASTRO feature overlap into {expmat_bed}")
    _run_astro_feature_overlap(overlap_bed, expmat_bed, gtf_path, output_dir, args)
    _log_status("velocity: ASTRO feature overlap completed")

    _log_status("velocity: building sparse AnnData from ASTRO feature BED")
    adata, matrix_stats = _sparse_adata_from_astro_expmat(
        expmat_bed=expmat_bed,
        read_map=read_map,
        output_dir=output_dir,
        ordered_locs=ordered_locs,
        loc_to_xy=loc_to_xy,
        args=args,
        libs=libs,
    )
    adata.uns["astro_velocity"] = {
        "bam": str(bam_path),
        "bam_format": "10x_cellranger",
        "gtf": str(gtf_path),
        "positions": str(positions_path),
        "barcode_whitelist": str(Path(args.barcode_whitelist).expanduser()) if args.barcode_whitelist else None,
        "bam_records": stats,
        "matrix": matrix_stats,
    }

    if not args.keep_temp and not reused_intermediate:
        for path in (gtf_chrom_bed,):
            try:
                path.unlink()
            except OSError:
                pass

    return adata, f"ASTRO count from 10x/Cell Ranger BAM {bam_path} with {gtf_path}"


def _run_astro_feature_overlap(input_bed, expmat_bed, gtf_path, output_dir, args):
    done_path = expmat_bed.with_suffix(expmat_bed.suffix + ".done")
    tmp_output = expmat_bed.with_suffix(expmat_bed.suffix + ".tmp")
    for stale_path in (tmp_output, done_path):
        try:
            stale_path.unlink()
        except OSError:
            pass

    n_jobs = max(1, int(getattr(args, "n_jobs", 1) or 1))
    if n_jobs == 1 and not str(input_bed).endswith(".gz"):
        _run_overlap_chunk_worker((str(input_bed), str(tmp_output), str(gtf_path)))
        os.replace(tmp_output, expmat_bed)
        done_path.write_text("ok\n")
        return

    chunk_lines = max(1, int(getattr(args, "overlap_chunk_lines", 2000000) or 2000000))
    chunk_paths = _split_bed_for_parallel_overlap(input_bed, output_dir, chunk_lines)
    if len(chunk_paths) == 1:
        _log_status("velocity: running ASTRO overlap on one decompressed BED chunk")
        _run_overlap_chunk_worker((str(chunk_paths[0]), str(tmp_output), str(gtf_path)))
        os.replace(tmp_output, expmat_bed)
        done_path.write_text("ok\n")
        try:
            chunk_paths[0].unlink()
        except OSError:
            pass
        return

    _log_status(
        "velocity: running parallel ASTRO overlap "
        f"({len(chunk_paths)} chunks, {n_jobs} jobs)"
    )
    chunk_outputs = [
        output_dir / f".astro_overlap_chunk_{index:05d}.expmat.bed"
        for index in range(len(chunk_paths))
    ]
    tasks = [
        (str(chunk_path), str(chunk_output), str(gtf_path))
        for chunk_path, chunk_output in zip(chunk_paths, chunk_outputs)
    ]
    with multiprocessing.Pool(processes=n_jobs) as pool:
        pool.map(_run_overlap_chunk_worker, tasks)

    with tmp_output.open("w") as fout:
        for chunk_output in chunk_outputs:
            with chunk_output.open() as fin:
                for line in fin:
                    fout.write(line)

    os.replace(tmp_output, expmat_bed)
    done_path.write_text("ok\n")
    for path in chunk_paths + chunk_outputs:
        try:
            path.unlink()
        except OSError:
            pass


def _split_bed_for_parallel_overlap(input_bed, output_dir, chunk_lines):
    for stale_path in output_dir.glob(".astro_overlap_chunk_*.bed"):
        try:
            stale_path.unlink()
        except OSError:
            pass

    chunk_paths = []
    chunk_handle = None
    chunk_index = -1
    lines_in_chunk = 0
    previous_read_name = None

    def open_next_chunk():
        nonlocal chunk_handle, chunk_index, lines_in_chunk
        if chunk_handle is not None:
            chunk_handle.close()
        chunk_index += 1
        chunk_path = output_dir / f".astro_overlap_chunk_{chunk_index:05d}.bed"
        chunk_paths.append(chunk_path)
        chunk_handle = chunk_path.open("w")
        lines_in_chunk = 0

    with _open_text(input_bed) as fin:
        for line in fin:
            fields = line.split("\t", 4)
            read_name = fields[3] if len(fields) > 3 else ""
            if chunk_handle is None:
                open_next_chunk()
            elif lines_in_chunk >= chunk_lines and read_name != previous_read_name:
                open_next_chunk()
            chunk_handle.write(line)
            lines_in_chunk += 1
            previous_read_name = read_name

    if chunk_handle is not None:
        chunk_handle.close()
    if not chunk_paths:
        raise SystemExit(f"No BED records were found in {input_bed}")
    _log_status(f"velocity: split ASTRO BED into {len(chunk_paths)} overlap chunks")
    return chunk_paths


def _run_overlap_chunk_worker(task):
    input_bed, output_bed, gtf_path = task
    from .countfeature import inter_bed2geneFile

    inter_bed2geneFile(input_bed, output_bed, gtf_path)
    return output_bed


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


def _load_barcode_whitelist(whitelist_path):
    if not whitelist_path:
        return None

    path = Path(whitelist_path).expanduser().resolve()
    if not path.exists():
        raise SystemExit(f"Barcode whitelist file does not exist: {path}")

    whitelist = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if "," in line:
                barcode = line.split(",")[-1].strip()
            else:
                barcode = line.split()[-1].strip()
            if barcode.lower() in {"barcode", "barcodes"}:
                continue
            whitelist.add(barcode)

    if not whitelist:
        raise SystemExit(f"No barcodes were loaded from whitelist: {path}")
    return whitelist


def _barcode_in_whitelist(barcode, whitelist):
    if barcode in whitelist:
        return True
    stripped = barcode.split("-", 1)[0]
    return stripped in whitelist or f"{stripped}-1" in whitelist


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


def _load_astro_barcode_locations(barcodes_file, use_barcode_as_loc=False):
    ordered_locs = []
    loc_to_xy = {}
    with _open_text(barcodes_file) as fin:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            barcode, x, y = parts[0], parts[1], parts[2]
            loc = barcode if use_barcode_as_loc else f"{x}_{y}"
            if loc not in loc_to_xy:
                ordered_locs.append(loc)
                loc_to_xy[loc] = (x, y)
    if not ordered_locs:
        raise SystemExit(f"No barcode coordinates were loaded from {barcodes_file}.")
    return ordered_locs, loc_to_xy


def _load_gtf_chroms(gtf_path):
    chroms = set()
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0]:
                chroms.add(fields[0])
    if not chroms:
        raise SystemExit(f"No chromosome names were found in GTF: {gtf_path}")
    return chroms


def _chrom_for_gtf(chrom, gtf_chroms):
    for candidate in _chrom_candidates(chrom):
        if candidate in gtf_chroms:
            return candidate
    return None


def _chrom_candidates(chrom):
    candidates = [chrom]
    if chrom.startswith("chr"):
        stripped = chrom[3:]
        candidates.append(stripped)
        if stripped == "M":
            candidates.extend(["MT", "M"])
    else:
        candidates.append(f"chr{chrom}")
        if chrom in ("M", "MT"):
            candidates.extend(["chrM", "chrMT"])
    return candidates


def _prepare_astro_bed_for_gtf(input_bed, output_bed, gtf_chroms):
    stats = {
        "input_bed": str(input_bed),
        "bed_for_overlap": str(input_bed),
        "records": 0,
        "converted": 0,
        "skipped_chrom_not_in_gtf": 0,
    }
    needs_output = False
    with _open_text(input_bed) as fin:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = _chrom_for_gtf(fields[0], gtf_chroms) if fields else None
            if chrom is None:
                needs_output = True
            elif chrom != fields[0]:
                needs_output = True
            break

    if not needs_output:
        return input_bed, stats

    with _open_text(input_bed) as fin, output_bed.open("w") as fout:
        for line in fin:
            if not line.strip():
                continue
            stats["records"] += 1
            fields = line.rstrip("\n").split("\t")
            chrom = _chrom_for_gtf(fields[0], gtf_chroms)
            if chrom is None:
                stats["skipped_chrom_not_in_gtf"] += 1
                continue
            if chrom != fields[0]:
                stats["converted"] += 1
                fields[0] = chrom
            fout.write("\t".join(fields) + "\n")

    if stats["records"] == stats["skipped_chrom_not_in_gtf"]:
        raise SystemExit(
            "All BED records were skipped because their chromosomes are absent "
            "from the GTF. Check that BAM and GTF use the same reference."
        )
    stats["bed_for_overlap"] = str(output_bed)
    return output_bed, stats


def _write_astro_bed_from_10x_bam(
    bam_path,
    output_bed,
    read_map,
    barcode_to_xy,
    gtf_chroms,
    args,
    use_barcode_locs=False,
    ordered_locs=None,
    loc_to_xy=None,
    barcode_whitelist=None,
):
    stats = {
        "sam_records": 0,
        "written_reads": 0,
        "written_bed_blocks": 0,
        "missing_barcode": 0,
        "barcode_not_in_whitelist": 0,
        "barcode_not_in_positions": 0,
        "missing_umi": 0,
        "missing_blocks": 0,
        "chrom_name_converted": 0,
        "chrom_not_in_gtf": 0,
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
            raw_barcode = tags.get(args.barcode_tag) or tags.get("CB") or tags.get("CR")
            if not raw_barcode:
                stats["missing_barcode"] += 1
                continue
            if barcode_whitelist is not None and not _barcode_in_whitelist(raw_barcode, barcode_whitelist):
                stats["barcode_not_in_whitelist"] += 1
                continue
            barcode = _normalize_barcode(raw_barcode, barcode_to_xy)
            if barcode is None:
                if not use_barcode_locs:
                    stats["barcode_not_in_positions"] += 1
                    continue
                barcode = raw_barcode
                cell_number = len(barcode_to_xy) + 1
                barcode_to_xy[barcode] = (str(cell_number), "0")
                if ordered_locs is not None:
                    ordered_locs.append(barcode)
                if loc_to_xy is not None:
                    loc_to_xy[barcode] = barcode_to_xy[barcode]

            umi = tags.get(args.umi_tag) or tags.get("UB") or tags.get("UR")
            if not umi:
                stats["missing_umi"] += 1
                umi = fields[0]

            bam_chrom = fields[2]
            chrom = _chrom_for_gtf(bam_chrom, gtf_chroms)
            if chrom is None:
                stats["chrom_not_in_gtf"] += 1
                continue
            if chrom != bam_chrom:
                stats["chrom_name_converted"] += 1
            pos = int(fields[3])
            cigar = fields[5]
            blocks = _cigar_to_bed_blocks(pos, cigar)
            if not blocks:
                stats["missing_blocks"] += 1
                continue

            read_id += 1
            x, y = barcode_to_xy[barcode]
            loc = barcode if use_barcode_locs else f"{x}_{y}"
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


def _sparse_adata_from_astro_expmat(expmat_bed, read_map, output_dir, ordered_locs, loc_to_xy, args, libs):
    np, pd, _sc, _scv, ad, sparse, _plt = libs
    command_name = getattr(args, "command", "velocity")
    barcode_whitelist = _load_barcode_whitelist(getattr(args, "barcode_whitelist", None))
    if barcode_whitelist is not None:
        _log_status(f"{command_name}: loaded {len(barcode_whitelist)} whitelist barcodes")
    filter_settings = _parse_qualityfilter(args.qualityfilter)

    loc_to_index = {loc: idx for idx, loc in enumerate(ordered_locs)}
    gene_to_index = {}
    genes = []
    counts = {
        "spliced": defaultdict(int),
        "unspliced": defaultdict(int),
    }
    used_loc_indices = set()
    stats = {
        "feature_lines": 0,
        "candidate_umis": 0,
        "counted_umis": 0,
        "deduplicated_umis": 0,
        "missing_read_map": 0,
        "read_map_records_scanned": 0,
        "missing_velocity_suffix": 0,
        "quality_filtered": 0,
    }
    bucket_count = max(1, int(getattr(args, "dedup_buckets", 512) or 512))
    bucket_dir = output_dir / ".velocity_umi_buckets"
    bucket_dir.mkdir(exist_ok=True)
    for stale_path in bucket_dir.glob("bucket_*.tsv"):
        try:
            stale_path.unlink()
        except OSError:
            pass
    bucket_paths = [bucket_dir / f"bucket_{idx:04d}.tsv" for idx in range(bucket_count)]

    read_map_fh = _open_text(read_map)
    current_read_id = None
    current_read_record = None

    def advance_read_map(target_read_id):
        nonlocal current_read_id, current_read_record
        target = int(target_read_id)
        while current_read_id is None or current_read_id < target:
            map_line = read_map_fh.readline()
            if not map_line:
                current_read_id = None
                current_read_record = None
                return None
            stats["read_map_records_scanned"] += 1
            if args.progress_interval and stats["read_map_records_scanned"] % args.progress_interval == 0:
                _log_status(
                    f"{command_name}: scanned "
                    f"{stats['read_map_records_scanned']} read_map records, "
                    f"{stats['candidate_umis']} candidate UMIs"
                )
            read_id_text, loc, barcode, umi = map_line.rstrip("\n").split("\t", 3)
            current_read_id = int(read_id_text)
            if barcode_whitelist is not None and not _barcode_in_whitelist(barcode, barcode_whitelist):
                current_read_record = None
            else:
                current_read_record = (loc, barcode, umi)
        if current_read_id == target:
            return current_read_record
        return None

    bucket_handles = []
    try:
        try:
            bucket_handles = [path.open("w", buffering=1024 * 1024) for path in bucket_paths]
        except OSError as exc:
            raise SystemExit(
                f"Could not open {bucket_count} UMI dedup bucket files. "
                "Use a smaller --dedup-buckets value."
            ) from exc

        with read_map_fh, _open_text(expmat_bed) as fin:
            for line in fin:
                line = line.strip()
                if not line:
                    continue
                stats["feature_lines"] += 1
                if args.progress_interval and stats["feature_lines"] % args.progress_interval == 0:
                    _log_status(
                        f"{command_name}: parsed "
                        f"{stats['feature_lines']} expmat.bed lines, "
                        f"{stats['candidate_umis']} candidate UMIs"
                    )
                parsed = _parse_astro_feature_line(line, filter_settings)
                if parsed is None:
                    stats["quality_filtered"] += 1
                    continue
                read_id, gene_feature = parsed
                try:
                    read_record = advance_read_map(read_id)
                except ValueError:
                    stats["missing_read_map"] += 1
                    continue
                if read_record is None:
                    stats["missing_read_map"] += 1
                    continue

                gene_name, layer = _split_astro_velocity_feature(gene_feature)
                if layer is None:
                    stats["missing_velocity_suffix"] += 1
                    continue

                loc, _barcode, umi = read_record
                loc_idx = loc_to_index.get(loc)
                if loc_idx is None:
                    loc_idx = len(ordered_locs)
                    loc_to_index[loc] = loc_idx
                    ordered_locs.append(loc)
                    loc_to_xy[loc] = _xy_from_loc(loc)

                gene_idx = gene_to_index.get(gene_name)
                if gene_idx is None:
                    gene_idx = len(genes)
                    gene_to_index[gene_name] = gene_idx
                    genes.append(gene_name)

                layer_id = 0 if layer == "spliced" else 1
                bucket_index = hash((layer_id, loc_idx, gene_idx, umi)) % bucket_count
                bucket_handles[bucket_index].write(f"{layer_id}\t{loc_idx}\t{gene_idx}\t{umi}\n")
                stats["candidate_umis"] += 1
    finally:
        for handle in bucket_handles:
            try:
                handle.close()
            except OSError:
                pass

    if stats["candidate_umis"] == 0 or not genes:
        raise SystemExit("Existing expmat.bed did not produce exon/intron UMI counts.")

    _log_status(
        f"{command_name}: aggregating UMI hash buckets "
        f"({bucket_count} buckets, {stats['candidate_umis']} candidate UMIs)"
    )
    bucket_log_interval = max(1, bucket_count // 10)
    for bucket_index, bucket_path in enumerate(bucket_paths):
        if args.progress_interval and bucket_index % bucket_log_interval == 0:
            _log_status(
                f"{command_name}: aggregating UMI bucket "
                f"{bucket_index + 1}/{bucket_count}"
            )
        if not bucket_path.exists() or bucket_path.stat().st_size == 0:
            continue
        seen_bucket = set()
        with bucket_path.open() as bucket_in:
            for line in bucket_in:
                key = line.rstrip("\n")
                if key in seen_bucket:
                    continue
                seen_bucket.add(key)
                fields = key.split("\t", 3)
                if len(fields) != 4:
                    continue
                layer_id_text, loc_idx_text, gene_idx_text, _umi = fields
                try:
                    layer_id = int(layer_id_text)
                    loc_idx = int(loc_idx_text)
                    gene_idx = int(gene_idx_text)
                except ValueError:
                    continue
                layer = "spliced" if layer_id == 0 else "unspliced"
                counts[layer][(loc_idx, gene_idx)] += 1
                used_loc_indices.add(loc_idx)
                stats["counted_umis"] += 1
                if args.progress_interval and stats["counted_umis"] % args.progress_interval == 0:
                    _log_status(
                        f"{command_name}: aggregated "
                        f"{stats['counted_umis']} unique UMIs"
                    )

    stats["deduplicated_umis"] = stats["candidate_umis"] - stats["counted_umis"]

    if not used_loc_indices:
        raise SystemExit("Existing expmat.bed did not produce exon/intron UMI counts.")

    used_loc_indices = sorted(used_loc_indices)
    row_remap = {old_idx: new_idx for new_idx, old_idx in enumerate(used_loc_indices)}
    used_locs = [ordered_locs[idx] for idx in used_loc_indices]

    spliced = _counts_to_sparse_matrix(counts["spliced"], row_remap, len(used_locs), len(genes), sparse)
    unspliced = _counts_to_sparse_matrix(counts["unspliced"], row_remap, len(used_locs), len(genes), sparse)
    obs_names = [loc.replace("_", "x") for loc in used_locs]

    obs = pd.DataFrame(index=obs_names)
    obs["x"] = [float(loc_to_xy[loc][0]) for loc in used_locs]
    obs["y"] = [float(loc_to_xy[loc][1]) for loc in used_locs]
    var = pd.DataFrame(index=pd.Index(genes, dtype=str))

    adata = ad.AnnData(X=spliced, obs=obs, var=var)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.obsm["X_space"] = obs[["x", "y"]].to_numpy(dtype=float)

    spliced_path = output_dir / "spliced_counts.npz"
    unspliced_path = output_dir / "unspliced_counts.npz"
    obs_path = output_dir / "sparse_obs_names.tsv"
    var_path = output_dir / "sparse_var_names.tsv"
    sparse.save_npz(spliced_path, spliced)
    sparse.save_npz(unspliced_path, unspliced)
    pd.Series(obs_names, name="obs_name").to_csv(obs_path, sep="\t", index=False)
    pd.Series(genes, name="gene").to_csv(var_path, sep="\t", index=False)

    stats["cells"] = len(used_locs)
    stats["genes"] = len(genes)
    stats["spliced_nnz"] = int(spliced.nnz)
    stats["unspliced_nnz"] = int(unspliced.nnz)
    stats["spliced_counts"] = str(spliced_path)
    stats["unspliced_counts"] = str(unspliced_path)
    stats["obs_names"] = str(obs_path)
    stats["var_names"] = str(var_path)
    stats["umi_deduplication"] = "hash_buckets"
    stats["umi_dedup_buckets"] = bucket_count
    if getattr(args, "keep_temp", False):
        stats["umi_bucket_dir"] = str(bucket_dir)
    else:
        for path in bucket_paths:
            try:
                path.unlink()
            except OSError:
                pass
        try:
            bucket_dir.rmdir()
        except OSError:
            pass
    return adata, stats


def _counts_to_sparse_matrix(counts, row_remap, n_rows, n_cols, sparse):
    rows = []
    cols = []
    data = []
    for (old_row, col), value in counts.items():
        new_row = row_remap.get(old_row)
        if new_row is None:
            continue
        rows.append(new_row)
        cols.append(col)
        data.append(value)
    return sparse.coo_matrix((data, (rows, cols)), shape=(n_rows, n_cols)).tocsr()


def _xy_from_loc(loc):
    parts = loc.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return "0", "0"


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


def _write_velocity_notebook(output_dir, source, args):
    notebook_path = output_dir / "rna_velocity_report.ipynb"
    figures = sorted((output_dir / "figures").glob("*.png"))
    command_name = getattr(args, "command", "velocity")

    cells = [
        _markdown_cell(
            [
                "# RNA velocity report\n",
                "\n",
                f"Generated by `ASTROutils {command_name}`.\n",
                "\n",
                f"Source: `{source}`\n",
            ]
        ),
        _code_cell(
            [
                "import json\n",
                "from pathlib import Path\n",
                "import pandas as pd\n",
                "import scvelo as scv\n",
                "\n",
                "out = Path('.')\n",
                "adata = scv.read(out / 'rna_velocity.h5ad')\n",
                "adata\n",
            ]
        ),
        _code_cell(
            [
                "with open(out / 'run_summary.json') as fh:\n",
                "    summary = json.load(fh)\n",
                "summary\n",
            ]
        ),
        _code_cell(
            [
                "obs_path = out / 'velocity_obs.csv'\n",
                "pd.read_csv(obs_path, index_col=0).head() if obs_path.exists() else adata.obs.head()\n",
            ]
        ),
        _code_cell(
            [
                "adata.layers.keys(), adata.obsm.keys()\n",
            ]
        ),
    ]

    if figures:
        cells.append(_markdown_cell(["## Generated figures\n"]))
        for figure in figures:
            rel_path = figure.relative_to(output_dir)
            cells.append(
                _markdown_cell(
                    [
                        f"### {figure.stem}\n",
                        "\n",
                        f"![{figure.stem}]({rel_path.as_posix()})\n",
                    ]
                )
            )

    cells.append(
        _code_cell(
            [
                "# Example: make additional scVelo plots from the saved AnnData object.\n",
                "# scv.pl.velocity_embedding_stream(adata, basis='umap', color='velocity_pseudotime')\n",
            ]
        )
    )

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "pygments_lexer": "ipython3",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    notebook_path.write_text(json.dumps(notebook, indent=1) + "\n")
    return notebook_path


def _markdown_cell(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source,
    }


def _code_cell(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source,
    }


def _log_status(message):
    print(message, file=sys.stderr, flush=True)


def _warn(message):
    print(f"Warning: {message}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
