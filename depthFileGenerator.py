import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-v","--vcf_file", required=True)

args = parser.parse_args()

vcfFile = open(args.vcf_file)
vcfFileLines = vcfFile.readlines()

vcfLines = len(vcfFileLines)

data = ['Position', 'Depth', 'RefCount', 'Altcount', 'AltPercent']
with open('depth.txt', 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerow(data)

for y in range(0, vcfLines):
    for d in range(0,12):
        if ((vcfFileLines[y+1].split()[7].split(";")[d][:3]) == "DP4"):
            dp4 = (((((vcfFileLines[y+1].split())[7]).split(";")[d])[4:]).split(","))
            ref = (int(dp4[0])) + (int(dp4[1]))
            alt = int(dp4[2]) + int(dp4[3])
            pos = int(vcfFileLines[y].split()[1])
            depth = ref + alt
            altPct = float(alt / depth)
            line = [pos, depth, ref, alt, altPct]
            with open('depth.txt', 'a', newline='') as file_output:
                tsv_output = csv.writer(file_output, delimiter='\t')
                tsv_output.writerow(line)
            break
