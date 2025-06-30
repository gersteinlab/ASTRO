import subprocess
import re
import argparse


ap = argparse.ArgumentParser()
ap.add_argument("-r", "--referencefile", required=True, help="Path to input BAM file")
ap.add_argument("-i", "--inputfa", required=True, help="Path to input BAM file")
ap.add_argument("-o", "--outputfile", required=True, help="Path to input BAM file")
ap.add_argument(
    "-t", "--thread", required=False, default="8", type=str, help="Number of threads"
)
args = vars(ap.parse_args())

outputfile = args["outputfile"]
inputfa = args["inputfa"]
referencefile = args["referencefile"]
thread = args["thread"]


cmd = [
    "bowtie2",
    "-p",
    thread,
    "-x",
    referencefile,
    "-U",
    inputfa,
    "-D",
    "7",
    "-R",
    "2",
    "-N",
    "1",
    "-L",
    "6",
    "-i",
    "C,7",
    "--score-min",
    "C,10",
    "--no-hd",
    "--no-sq",
    "-a",
    "--ma",
    "1",
    "--rdg",
    "0,1",
    "--rfg",
    "0,1",
    "--mp",
    "1,1",
    "--gbar",
    "1",
    "--local",
]


process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

with open(outputfile, "w") as output_file:
    major_read = None
    major_map = None
    major_value = None
    see1 = 0
    for line in process.stdout:
        parts = line.strip().split()
        current_read = parts[0]
        match = re.match(r"^AS:i:(\d+)$", parts[11])
        current_map = parts[2]
        if current_map != "*":
            if match:
                current_value = int(match.group(1))
            else:
                print("Error in line:", line)
                raise ValueError("Tag string does not match the required format AS:i:x")
        else:
            output_file.write(f"{major_read}\t*" + "\n")
            continue
        if current_read != major_read:
            if major_read is not None:
                output_file.write(f"{major_read}\t{major_map}" + "\n")
            see1 = 1
            major_read = current_read
            major_value = current_value
            major_map = current_map
        else:
            if see1 == 1:
                if major_value <= current_value:
                    major_map = "*"
                see1 = 0
    if major_read is not None:
        output_file.write(f"{major_read}\t{major_map}" + "\n")


process.wait()
