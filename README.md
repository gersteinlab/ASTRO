The usage instructions for the ASTRO Python package are as followings.  

1.Functional Overview  
Demultiplexing: Adapter trimming, UMI, and Barcode splitting.  
Genome Mapping: Uses STAR to align reads to the genome and optionally removes duplicate reads using either samtools markdup or a custom deduplication module.  
Feature Counting: Calculates gene expression from alignment results based on GTF annotation files and outputs an expression matrix.  
Feature Filtering: Filters out low-quality or abnormal genes/barcodes based on user-defined thresholds.  

2.Installation Guide  
External dependencies: ASTRO requires the following external tools to function- STAR, bedtools, samtools, and cutadapt. If you have Python 3.6+ installed locally, follow these steps to set it up.  
Install from source repository or compressed package:  
2.1 Clone the repository:  
git clone git@github.com:gersteinlab/ASTRO.git  
2.2 Enter the directory named "python":  
cd python  
2.3 Install dependencies and build/install:  
pip install -e .  
2.4 Check if the installation was successful (If you see the help documentation, the installation is complete.):  
ASTRO --help  
filtmatbyrt --help  
2.5 After installation, the following executable scripts will be available in the command line:  
ASTRO: Main pipeline entry point.  
filtmatbyrt: Independent row-column filtering/merging tool.  

Alternative: the Docker image is hosted on Docker Hub and can be downloaded using the following command

```bash
docker pull zhiyuanchu/astro:v1.0
```

3.Parameter Description
The ASTRO script accepts parameters via the command line or a JSON file. The main parameters are listed below. Certain parameters are only required for specific steps. If these steps are executed (--steps control) but their parameters are missing, the program will perform a runtime check and exit with an error.
<table>
  <thead>
    <tr>
      <th>Parameter Name</th>
      <th>Required</th>
      <th>Relevant Steps</th>
      <th>Default</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>R1</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>None</td>
      <td>FASTQ file 1 containing main RNA sequences</td>
    </tr>
    <tr>
      <td>R2</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>None</td>
      <td>FASTQ file 2 containing spatial barcode and UMI information</td>
    </tr>
    <tr>
      <td>barcode_file</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>None</td>
      <td>Text file with barcode information; must have at least three columns: barcode, tX, tY</td>
    </tr>
    <tr>
      <td>outputfolder</td>
      <td>Yes</td>
      <td>All Steps</td>
      <td>-</td>
      <td>Output directory for results</td>
    </tr>
    <tr>
      <td>starref</td>
      <td>Yes</td>
      <td>Step 2</td>
      <td>-</td>
      <td>STAR genome index directory</td>
    </tr>
    <tr>
      <td>gtffile</td>
      <td>Yes</td>
      <td>Steps 2+4</td>
      <td>None</td>
      <td>Path to the GTF file</td>
    </tr>
    <tr>
      <td>PrimerStructure1</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>-</td>
      <td>Primer structure for R1, e.g. <code>AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A&#123;10&#125;N&#123;150&#125;</code></td>
    </tr>
    <tr>
      <td>StructureUMI</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>-</td>
      <td>UMI structure definition, e.g. <code>CAAGCGTTGGCTTCTCGCATCT_10</code></td>
    </tr>
    <tr>
      <td>StructureBarcode</td>
      <td>Yes</td>
      <td>Step 1</td>
      <td>-</td>
      <td>Spatial barcode structure, e.g., 25_ATCCACGTGCTTGAGAGGCCAAGATCG: ATCCACGTGCTTGAGAGGCCAAGATCG... GTGGCCGATGTTTCGCATCGGCGTAGACT</td>
    </tr>
    <tr>
      <td>threadnum</td>
      <td>No</td>
      <td></td>
      <td>16</td>
      <td>Number of threads for parallel tasks like cutadapt, STAR, samtools, etc.</td>
    </tr>
    <tr>
      <td>options</td>
      <td>No</td>
      <td></td>
      <td>""</td>
      <td>
        Enables extra modes, including:<br>
        H: Hard mode<br>
        M: Use samtools markup for deduplication
      </td>
    </tr>
    <tr>
      <td>steps</td>
      <td>No</td>
      <td></td>
      <td>7</td>
      <td>
        Specifies steps to execute; uses bitwise integers:<br>
        1: Demultiplexing<br>
        2: Genome Mapping<br>
        4: Feature Counting and Filtering<br>
        For example, 7 = 1+2+4
      </td>
    </tr>
    <tr>
      <td>STARparamfile4genome</td>
      <td>No</td>
      <td></td>
      <td>NA</td>
      <td>File containing extra STAR parameters; specify the file path for custom STAR parameters</td>
    </tr>
    <tr>
      <td>qualityfilter</td>
      <td>No</td>
      <td></td>
      <td>"25:0.75"</td>
      <td>Quality filter threshold; default "25:0.75" means filtering if AS &le; 25 and &le; 0.75 of gene length. Set to 0:0 or NA to disable</td>
    </tr>
    <tr>
      <td>removeByDim</td>
      <td>No</td>
      <td></td>
      <td>True</td>
      <td>Whether to remove genes/barcodes with abnormal row/column variance. Set to False to skip</td>
    </tr>
    <tr>
      <td>filterlogratio</td>
      <td>No</td>
      <td></td>
      <td>2</td>
      <td>Filters based on log2 variance differences; default is log2 &gt; 2 (4x difference) for removal</td>
    </tr>
    <tr>
      <td>workflow</td>
      <td>No</td>
      <td></td>
      <td>new</td>
      <td>
      Specifies which pipeline to execute. If set to
      <code>"old"</code>, the older workflow will run.
      If not specified or left as default, it will run the latest
      workflow.
      </td>
    </tr>
    <tr>
      <td>barcodemode</td>
      <td>No</td>
      <td>Step 1</td>
      <td>"spatial"</td>
      <td>
      When set to "singlecell", enables single-cell mode. If kept as "spatial", the pipeline runs in the conventional spatial mode.
      </td>
    </tr>
    <tr>
      <td>barcode_threshold</td>
      <td>No</td>
      <td>Step 1</td>
      <td>100</td>
      <td>
      In single-cell mode, when extracting barcodes automatically, any barcode whose occurrence is below this threshold will be discarded.
      </td>
    </tr>
    <tr>
      <td>barcodelength</td>
      <td>No</td>
      <td>Step 1</td>
      <td>0</td>
      <td>
      In single-cell mode, if greater than 0, only barcodes of length barcodelength will be retained. If set to 0, no length filtering is applied.
      </td>
    </tr>
    <tr>
      <td>barcode_file</td>
      <td>No (for single-cell mode)<br/>Yes (for spatial mode)</td>
      <td>Step 1</td>
      <td>None or "notavailable"</td>
      <td>
      In spatial mode, this parameter must be provided by the user and must include barcode coordinate information. In single-cell mode, if a ready-made file is not available, set this to "notavailable", so ASTRO can generate a three-column barcode file automatically.
      </td>
    </tr>
    <tr>
      <td>genes2check</td>
      <td>No</td>
      <td>Step 4</td>
      <td>False</td>
      <td>
      If provided (a text file listing suspect genes, e.g. piRNAs/miRNAs), ASTRO runs an advanced check before feature counting. Genes failing this check are removed from the GTF annotation, effectively excluding them from the final expression matrix.
      </td>
    </tr>
  </tbody>
</table>

3.1 Parameter Priority  
ASTRO can receive parameters from the command line or a JSON file. JSON parameters can be specified as positional arguments (json_file_path1) or via --json_file_path. Priority is as follows:

1st: Command Line Parameters: Overrides JSON values if explicitly provided.  
E.G.:ASTRO --R1 R1.fq --R2 R2.fq \
--barcode_file spatial_barcodes.txt \
--gtffile hsa.no_piRNA.gtf --starref StarIndex/ \
--PrimerStructure1 AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150} \
--StructureUMI CAAGCGTTGGCTTCTCGCATCT_10 \
--StructureBarcode 20_ATCCACGTGCTTGAGAGGCCAGAGCATTCG:... \
--threadnum 16 \
--steps 7 \
--outputfolder output/
To enable hard mode (H) and samtools markdup deduplication (M), add:
--options HM

2nd:JSON File  
--json_file_path: If specified, the program reads this file first.  
E.g.:ASTRO --json_file_path myparams.json  
json_file_path1: If --json_file_path is not provided, the program reads the positional JSON file. E.g.:ASTRO parameter.json

3rd:Default Values  
Used if neither command-line parameters nor JSON specify the value.

3.2 Key Modes Explanation  
(1) hard mode (H)  
If options include H, it means that when a single read aligns to multiple genes, it will no longer only take the first gene but will record all aligned genes as a multi-gene form, separated by a hyphen.  
(2) M decides which way is used for remove depulicate reads.
If options include M, samtools markdup will be used to mark and remove duplicate reads.  If M is not included, the built-in ASTRO deduplication logic will be used: This logic relies on UMIs and barcodes to determine duplicates and uses alignment score for filtering.

3.3 Quality Filtering and Additional Tool: filtmatbyrt  
This is a standalone script. If the removeByDim parameter is set to True, it will be automatically called in step 4 to perform row and column variance filtering, and the output will be saved as finalexpmat.tsv. If removeByDim is set to False, the script can also be run independently to further filter or merge an already generated expression matrix.  
Command example:  
filtmatbyrt good_expmat.tsv bad_expmat.tsv final_expmat.tsv 2  
Here, good_expmat.tsv represents a "high-quality" expression matrix, bad_expmat.tsv represents a "suspicious or to-be-merged" expression matrix, and final_expmat.tsv is the filtered and merged output matrix. 2 indicates --filterlogratio=2.  
It can also be used with parameter syntax:  
filtmatbyrt --expmatgood good_expmat.tsv --expmatbad bad_expmat.tsv --finalexpmat final_expmat.tsv --filterlogratio 2

4.A Simple Example

4.1 Assume the following files are prepared:  
R1.fq, R2.fq: Input sequencing reads.  
spatial_barcodes.txt: Records coordinates and barcode sequences.  
StarIndex/: STAR genome index directory.  
hsa.no_piRNA.gtf: Gene annotation file.  

4.2 Create the following JSON file (parameter.json):  
{  
  "R1": "R1.fq",  
  "R2": "R2.fq",  
  "barcode_file": "spatial_barcodes.txt",  
  "PrimerStructure1": "AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150}",
  "StructureUMI": "CAAGCGTTGGCTTCTCGCATCT_10",  
  "StructureBarcode": "20_ATCCACGTGCTTGAGAGGCCAGAGCATTCG:...GTGGCCGATGTTTCGCATCGGCGTACGACT",  
  "threadnum": 16,  
  "steps": 7,  
  "outputfolder": "output/",  
  "gtffile": "hsa.no_piRNA.gtf",  
  "starref": "StarIndex/"  
}

4.3 Run the command:
ASTRO parameter.json

5. Single-Cell Mode

In the conventional spatial transcriptomics (spatial) mode, ASTRO requires a `barcode_file` that contains at least three columns (barcode sequence, X coordinate, and Y coordinate).

However, in single-cell mode (`barcodemode="singlecell"`), if the user does not provide an existing barcode file (or sets it to `"notavailable"`), ASTRO will, during **Step 1 (Demultiplexing)**, automatically enumerate all possible barcodes from R2 (which typically contains cell barcode information) based on **`structurebarcode`**, and then filter them according to **`barcode_threshold`** (the minimum count of occurrences required for each barcode) and **`barcodelength`** (if greater than 0, only barcodes of that specific length are retained). It then generates an “auto-generated” three-column barcode file (in the form of `barcode, i, i`) for subsequent steps.

In the feature counting output, ASTRO will replace the temporary column names (`ixi`) with the actual single-cell barcodes, allowing these barcodes to be used directly as column names in the expression matrix. This mode is particularly suitable for single-cell sequencing scenarios where the precise barcodes are not known in advance, or when they need to be automatically extracted and filtered from raw sequence data.

If a barcode file is provided in single-cell mode, ASTRO will treat that file as a “whitelist”—after automatically enumerating barcodes, it retains only those that are both present in the whitelist and meet the threshold/length requirements. A three-column barcode file is then generated for the downstream workflow.

In single-cell mode, Step 4 (Feature Filtering) by default only creates the basic expression matrix and does not perform further row/column variance filtering. If needed, you can run that manually afterward or configure the relevant parameters.

6.1 Symbol meanings in StructureBarcode / StructureUMI

A colon (:) means concatenating multiple segments. 

“20_CGTTGGCTTCT”Means that, on the 3' end, we recognize the fixed adapter CGTTGGCTTCT.Then we take the 20 nucleotides to the left of that fixed sequence.In other words, find CGTTGGCTTCT at the read’s 3' end and extract the preceding 20 bases.

“CGTTGGCTTCT_20”Means that, on the 5' end, we recognize the fixed adapter CGTTGGCTTCT.Then we take the 20 nucleotides to the right of that fixed sequence.In other words, find CGTTGGCTTCT at the read’s 5' end and keep the next 20 bases.

“TTCTCGCATCT...ATCCACGTGCTTGA”Means we take what lies between two stable (known) sequences.Here, the left boundary is TTCTCGCATCT and the right boundary is ATCCACGTGCTTGA, so we extract whatever is in between them.

6.2 Two ways to write (or “point to”) barcode/UMI locations in R2

**Method A**: Use explicit numeric positions

When you already know the barcode is strictly located at positions (23–30, 61–68, 99–106) and the UMI is at (137–146), you can specify:

"StructureBarcode": "22_8:60_8:98_8"

"StructureUMI": "136_10"

For 22_8, it means: skip the first 22 bases, then keep the next 8 bases (covering R2 positions 23–30).
For 60_8, it means: skip the first 60 bases, then keep the next 8 (positions 61–68).
For 98_8, it means skip the first 98 bases, keep the next 8 (positions 99–106).

Similarly, 136_10 means skip 136 bases, then keep the next 10 (positions 137–146) for the UMI.

Caution: If the read does not strictly match these exact positions (for instance, if some reads have shorter length or shifted inserts), then using hard-coded positions can cause errors or incorrect trimming.

**Method B**: Rely on stable flanking sequences

In many protocols, the positions can shift a bit, but each barcode or UMI region is still bounded by stable/fixed adapter sequences. In that scenario, you can specify the left and right flanking adapters that sandwich the barcode or UMI.

For example, assume we know the read has the following layout of R2 (showing only some parts for illustration):

XX AGCGTTGGCTTCTCGCATCT BBBBBBBB ATCCACGTGCTTGAGCGCGCTGCATACTTG BBBBBBBB CCCATGATCGTCCGAAGGCCAGAGCATTCG BBBBBBBB GTGGCCGATGTTTCGCATCGGCGTACGACT UUUUUUUUUU XXXXX

X stands for arbitrary nucleotides (non-barcode).B is the actual barcode.U stands for the UMI region.The bold segments are stable, known sequences that mark the boundaries before and after each barcode/UMI.

By observation:The segment at positions 23–30 is bounded by AGCGTTGGCTTCTCGCATCT (left) and ATCCACGTGCTTGAGCGCGCTGCATACTTG (right). The 8 bp sandwiched between these two adapters is the barcode.Similarly, for positions 61–68, the stable adapters around that 8 bp region are ATCCACGTGCTTGAGCGCGCTGCATACTTG (left) and CCCATGATCGTCCGAAGGCCAGAGCATTCG (right).For positions 99–106, the bounding adapters are CCCATGATCGTCCGAAGGCCAGAGCATTCG and GTGGCCGATGTTTCGCATCGGCGTACGACT.

Hence, you can write:
"StructureBarcode": "AGCGTTGGCTTCTCGCATCT...ATCCACGTGCTTGAGCGCGCTGCATACTTG : ATCCACGTGCTTGAGCGCGCTGCATACTTG...CCCATGATCGTCCGAAGGCCAGAGCATTCG : CCCATGATCGTCCGAAGGCCAGAGCATTCG...GTGGCCGATGTTTCGCATCGGCGTACGACT"

(where each colon : indicates concatenating those three sub-barcodes in order).

If the UMI is at positions 137–146, and it always follows the stable prefix GTGGCCGATGTTTCGCATCGGCGTACGACT, then you can define:

"StructureUMI": "GTGGCCGATGTTTCGCATCGGCGTACGACT_10"

meaning we look for that fixed adapter on the 5' end, and once found, we keep the next 10 bases as the UMI.

This approach gives flexibility in cases where the read length or positions vary slightly, as long as the bounding sequences remain identifiable.

7. Genes2Check Advanced Filtering

If your experiment includes some suspicious genes (for example, piRNA or miRNA), you can list their gene IDs or gene names in a text file (e.g., genes2check.txt) and specify this file in the ASTRO parameters using --genes2check genes2check.txt (or by adding "genes2check": "genes2check.txt" in the JSON). When running Step 4 (Feature Counting), ASTRO will call its built-in advanced detection logic (getvalidedgtf_parallel) to determine whether to discard these suspicious genes based on their over-enrichment in corresponding control regions. If anomalous enrichment is detected, those entries will be removed before formal counting, thereby reducing false positives and yielding a more accurate gene expression matrix.

This feature relies on the BAM file (STAR/tempfiltered.bam) and its index (which will be created automatically if absent) generated in the previous step. It then checks each gene interval in a multithreaded manner. The detection algorithm uses statistical tests (such as Poisson tests), and if it concludes that a gene interval is significantly higher than the background, it deems the interval likely to be junk or a spurious mapping. You can customize the list of genes in genes2check.txt according to your needs (e.g., including lncRNA, miRNA, or piRNA).
