import os
import subprocess


def run_old_pipeline(
    R1,
    R2,
    barcode_file,
    outputfolder,
    starref="data/reference/human",
    gtffile="data/reference/RNAcentral_hsa",
    PrimerStructure1="AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150}",
    StructureUMI="b_ATCCACGTGCTTGAGAGGCCAGAGCATTCG_b_GTGGCCGATGTTTCGCATCGGCGTACGACT_10",
    StructureBarcode="b_ATCCACGTGCTTGAGAGGCCAGAGCATTCG_b_GTGGCCGATGTTTCGCATCGGCGTACGACT_10",
    scriptFolder="scripts/",
    barcodeposition="NANA",
    barcodelengthrange="NANA",
    threadnum=16,
):

    scriptFolder = scriptFolder.rstrip("/") + "/"

    outputfolder = outputfolder.rstrip("/")
    os.makedirs(outputfolder, exist_ok=True)
    os.makedirs(os.path.join(outputfolder, "temps"), exist_ok=True)

    logf = os.path.join(outputfolder, "log.txt")
    with open(logf, "w") as f:
        f.write(
            f"{R1}\n{R2}\n{barcode_file}\n{outputfolder}\n{starref}\n{gtffile}\n{PrimerStructure1}\n"
        )
        f.write(
            f"{StructureUMI}\n{StructureBarcode}\n{scriptFolder}\n{barcodeposition}\n{barcodelengthrange}\n"
        )

    Cleanr1Fq1 = os.path.join(outputfolder, "temps", "cleanr1fq1.fq")
    Cleanr1Fq2 = os.path.join(outputfolder, "temps", "cleanr1fq2.fq")
    CombineFq = os.path.join(outputfolder, "combine.fq")

    barcode_db_fa = os.path.join(outputfolder, "temps", "barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps", "barcode_db")

    index_fq = os.path.join(outputfolder, "temps", "index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps", "UMI.fastq")
    index0_out = os.path.join(outputfolder, "temps", "index0.out")
    index_out = os.path.join(outputfolder, "temps", "index.out")
    unmapfq = os.path.join(outputfolder, "temps", "unmap.fastq")

    STARoutput = os.path.join(outputfolder, "STAR", "temp")
    STARbam = f"{STARoutput}Aligned.sortedByCoord.out.bam"
    pureSTARbam = f"{STARoutput}filtered.bam"

    expmatbed = os.path.join(outputfolder, "expmat.bed")
    expmattsv = os.path.join(outputfolder, "expmat.tsv")

    if R1 != "NA":

        if "_" in PrimerStructure1:
            prefixread1 = PrimerStructure1.split("_", 1)[0]
            suffixread1 = PrimerStructure1.rsplit("_", 1)[1]
        else:
            prefixread1 = PrimerStructure1
            suffixread1 = ""

        cmd_cutadapt = [
            "cutadapt",
            "-e",
            "0.25",
            "-a",
            suffixread1,
            "--times",
            "4",
            "-g",
            prefixread1,
            "-j",
            str(threadnum),
            "-o",
            Cleanr1Fq1,
            "-p",
            Cleanr1Fq2,
            R1,
            R2,
        ]
        subprocess.run(cmd_cutadapt, check=True)

        cmd_umi = [
            "python3",
            f"{scriptFolder}singleCutadapt.py",
            "-i",
            Cleanr1Fq2,
            "-o",
            UMI_fq,
            "-b",
            StructureUMI,
            "-t",
            str(threadnum),
        ]
        subprocess.run(cmd_umi, check=True)

        cmd_idx = [
            "python3",
            f"{scriptFolder}singleCutadapt.py",
            "-i",
            Cleanr1Fq2,
            "-o",
            index_fq,
            "-b",
            StructureBarcode,
            "-t",
            str(threadnum),
        ]
        subprocess.run(cmd_idx, check=True)

        with open(barcode_file, "r") as fin, open(barcode_db_fa, "w") as fout:
            for line in fin:
                line = line.strip()
                if not line:
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    continue
                seq = cols[0]
                xval = cols[1]
                yval = cols[2]
                fout.write(f">{xval}_{yval}\n{seq}\n")

        cmd_bowtie2build = ["bowtie2-build", barcode_db_fa, barcode_db_path]
        subprocess.run(cmd_bowtie2build, check=True)

        if barcodeposition == "NANA":
            cmd_runbowtie = [
                "python3",
                f"{scriptFolder}runBowtie.py",
                "-i",
                index_fq,
                "-o",
                index_out,
                "-r",
                barcode_db_path,
                "-t",
                str(threadnum),
            ]
            subprocess.run(cmd_runbowtie, check=True)

            cmd_barcodefq = [
                "perl",
                f"{scriptFolder}barcodedFq.pl",
                "-r",
                os.path.join(outputfolder, "temps", "index.out"),
                "-u",
                UMI_fq,
                "-o",
                CombineFq,
                "-i",
                Cleanr1Fq1,
            ]
            subprocess.run(cmd_barcodefq, check=True)

        else:
            cmd_ezbar = [
                "perl",
                f"{scriptFolder}ez_barcode.pl",
                "-u",
                unmapfq,
                "-o",
                index0_out,
                "-r",
                barcode_file,
                "-i",
                index_fq,
                "-b",
                barcodeposition,
                "-l",
                barcodelengthrange,
            ]
            subprocess.run(cmd_ezbar, check=True)

            cmd_runbowtie2 = [
                "python3",
                f"{scriptFolder}runBowtie.py",
                "-i",
                unmapfq,
                "-o",
                index_out,
                "-r",
                barcode_db_path,
                "-t",
                str(threadnum),
            ]
            subprocess.run(cmd_runbowtie2, check=True)

            merged_temp = index_out + ".temp"
            with open(merged_temp, "w") as fout, open(index0_out, "r") as f1, open(
                index_out, "r"
            ) as f2:
                for line in f1:
                    fout.write(line)
                for line in f2:
                    fout.write(line)
            os.replace(merged_temp, index_out)

            cmd_barcodefq2 = [
                "perl",
                f"{scriptFolder}barcodedFq.pl",
                "-r",
                index_out,
                "-u",
                UMI_fq,
                "-o",
                CombineFq,
                "-i",
                Cleanr1Fq1,
            ]
            subprocess.run(cmd_barcodefq2, check=True)

    params = (
        f" --runThreadN {threadnum}"
        " --outSAMattributes NH HI AS nM NM"
        " --genomeLoad NoSharedMemory"
        " --limitOutSAMoneReadBytes 200000000"
        " --outFilterMultimapNmax -1"
        " --outFilterMultimapScoreRange 0"
        " --readMatesLengthsIn NotEqual"
        " --limitBAMsortRAM 0"
        " --outMultimapperOrder Random"
        " --outSAMtype BAM SortedByCoordinate"
        " --outSAMunmapped Within"
        " --outSAMorder Paired"
        " --outSAMprimaryFlag AllBestScore"
        " --outSAMmultNmax -1"
        " --outFilterType Normal"
        " --outFilterScoreMinOverLread 0"
        " --alignSJDBoverhangMin 30"
        " --outFilterMatchNmin 15"
        " --outFilterMatchNminOverLread 0"
        " --outFilterMismatchNoverLmax 0.1"
        " --outFilterMismatchNoverReadLmax 0.15"
        " --alignIntronMin 20"
        " --alignIntronMax 1000000"
        " --alignEndsType Local"
        " --twopassMode Basic"
    )

    if starref != "NA":
        cmd_star = [
            "STAR",
            "--genomeDir",
            starref,
            "--readFilesIn",
            CombineFq,
            "--outFileNamePrefix",
            STARoutput,
            "--sjdbGTFfile",
            gtffile,
        ]

        cmd_star.extend(params.split())
        subprocess.run(cmd_star, check=True)

    cmd_redup = [
        "perl",
        f"{scriptFolder}redupBAM.pl",
        "-i",
        STARbam,
        "-o",
        pureSTARbam,
        "-p",
        str(threadnum),
    ]
    subprocess.run(cmd_redup, check=True)

    cmd_interBed = [
        "python3",
        f"{scriptFolder}interBED2GeneFile.py",
        "-g",
        gtffile,
        "-b",
        pureSTARbam,
        "-o",
        expmatbed,
    ]
    subprocess.run(cmd_interBed, check=True)

    cmd_genemat2 = [
        "perl",
        f"{scriptFolder}genemat2tsv.pl",
        "-i",
        expmatbed,
        "-o",
        expmattsv,
        "-d",
        barcode_file,
    ]
    subprocess.run(cmd_genemat2, check=True)

    print("[run_old_pipeline] Done. The output should match Fq2Mat_old.sh results.")
