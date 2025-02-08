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
git clone https://github.com/dingyaozhang/WholeST.git
2.2 Enter the directory named "python":
cd python
2.3 Install dependencies and build/install:
pip install -e .
2.4 Check if the installation was successful:
ASTRO --help
filtmatbyrt --help
2.5 If you see the help documentation, the installation is complete.
After installation, the following executable scripts will be available in the command line:
ASTRO: Main pipeline entry point.
filtmatbyrt: Independent row-column filtering/merging tool.

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
(2) samtools markdup (M)
If options include M, samtools markdup will be used to mark and remove duplicate reads.
Markdup is suitable for standard RNA-seq analysis. It also supports additional checks when specific tags (Barcode/UMI) are present.
If M is not included, the built-in ASTRO deduplication logic will be used:
This logic relies on UMIs and barcodes to determine duplicates and uses alignment score for filtering.

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
