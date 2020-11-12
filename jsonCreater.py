import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s","--samstats", required=True)
parser.add_argument("-o","--host_samstats", required=True)
parser.add_argument("-q","--qc_stats", required=True)
parser.add_argument("-k","--kraken2_report", required=True)
parser.add_argument("-m","--mean_coverage", required=True)
parser.add_argument("-v","--vcf_file", required=True)




args = parser.parse_args()
data = {}

samstatsHostMapped = open(args.host_samstats)
samstatsLinesHostMapped = samstatsHostMapped.readlines()

samstatsPrimerTrimmed = open(args.samstats)
samstatsLinesPrimerTrimmed = samstatsPrimerTrimmed.readlines()

meanCoverage = open(args.mean_coverage)
meanCoverageLines = meanCoverage.readlines()

qcStats = open(args.qc_stats)
qcStatsLines = qcStats.readlines()

kraken2Report = open(args.kraken2_report)
kraken2ReportLines = kraken2Report.readlines()

vcfFile = open(args.vcf_file)
vcfFileLines = vcfFile.readlines()



qcStatsHeadings = (qcStatsLines[0].split(","))
qcStatsHeadings[7] = qcStatsHeadings[7][:-1]
qcStatsValues = (qcStatsLines[1].split(","))
qcStatsValues[7] = qcStatsValues[7][:-1]

totalReads = (samstatsLinesHostMapped[1].split())
hostMappedReads = (samstatsLinesHostMapped[7].split())

mappedTrimmedReads = (samstatsLinesPrimerTrimmed[7].split())
meanLengthRead = (samstatsLinesPrimerTrimmed[25].split())
meanInsertSize = (samstatsLinesPrimerTrimmed[32].split())
insertSizeSD = (samstatsLinesPrimerTrimmed[33].split())

totalKraken = 0
unclassifiedReads = 0
unclassifiedPct=0
sarsPct=0
sars2Pct=0
humanPct=0

sum =0

for n in range(0, len(meanCoverageLines)):
    sum = sum + ((int(meanCoverageLines[n].split()[1]))*(int(meanCoverageLines[n].split()[2])))

size = int((meanCoverageLines[1].split()[2]))
meanCvg = float(sum/size)

for x in range(0, len(kraken2ReportLines)):
    totalKraken = totalKraken + int(kraken2ReportLines[x].split()[2])
    if (kraken2ReportLines[x].split()[3] == 'U'):
        unclassifiedReads = int(kraken2ReportLines[x].split()[2])
        unclassifiedPct = float(kraken2ReportLines[x].split()[0])
    elif (kraken2ReportLines[x].split()[3] == 'S') and (kraken2ReportLines[x].split()[5] == 'Severe'):
        sarsPct = float(kraken2ReportLines[x].split()[0])
    elif (kraken2ReportLines[x].split()[3] == 'S1') and (kraken2ReportLines[x].split()[5] == 'Severe'):
        sars2Pct = float(kraken2ReportLines[x].split()[0])
    elif (kraken2ReportLines[x].split()[3] == 'S') and (kraken2ReportLines[x].split()[5] == 'Homo'):
        humanPct = float(kraken2ReportLines[x].split()[0])

data["QCStats"] = {}
for x in range(0,8):
    data["QCStats"][qcStatsHeadings[x]] = qcStatsValues[x]
    if (x == 1) or (x == 2):
        data["QCStats"][qcStatsHeadings[x]] = float(qcStatsValues[x])
    elif (x==3) or (x==4):
        data["QCStats"][qcStatsHeadings[x]] = int(qcStatsValues[x])

vcfLines = len(vcfFileLines)

data["TaxonomicClassification"] = {}
data["TaxonomicClassification"]["Reads"] = totalKraken
data["TaxonomicClassification"]["Unclassified reads"] = unclassifiedReads
data["TaxonomicClassification"]["Unclassified percent"] = unclassifiedPct
data["TaxonomicClassification"]["Severe acute respiratory syndrome-related coronavirus"] = sarsPct
data["TaxonomicClassification"]["Severe acute respiratory syndrome coronavirus 2"] = sars2Pct
data["TaxonomicClassification"]["Homo sapiens"] = humanPct


data["VariantDetection"] = {}

for y in range(0, vcfLines):
    for d in range(0,12):
        if ((vcfFileLines[y].split()[7].split(";")[d][:3]) == "DP4"):
            dp4 = (((((vcfFileLines[y].split())[7]).split(";")[d])[4:]).split(","))
            ref = (int(dp4[0])) + (int(dp4[1]))
            alt = int(dp4[2]) + int(dp4[3])
            break

    data["VariantDetection"]["Variant_" + str(y+1)] = {}
    data["VariantDetection"]["Variant_" + str(y+1)]["POS"] = int(vcfFileLines[y].split()[1])
    data["VariantDetection"]["Variant_" + str(y + 1)]["REF"] = vcfFileLines[y].split()[3]
    data["VariantDetection"]["Variant_" + str(y + 1)]["ALT"] = vcfFileLines[y].split()[4]
    data["VariantDetection"]["Variant_" + str(y + 1)]["DEPTH"] = ref + alt
    data["VariantDetection"]["Variant_" + str(y + 1)]["REFCOUNT"] = ref
    data["VariantDetection"]["Variant_" + str(y + 1)]["ALTCOUNT"] = alt
    data["VariantDetection"]["Variant_" + str(y + 1)]["VAF"] = round(float(alt / (ref + alt)), 2)




data["SequencingStats"] = {}
data["SequencingStats"]["TotalReads"] = int(totalReads[3])
data["SequencingStats"]["HostDepletedReads"] = (int(totalReads[3]) - int(hostMappedReads[2]))
data["SequencingStats"]["HostMappedReads"] = int(hostMappedReads[2])
data["SequencingStats"]["MappedTrimmedReads"] = int(mappedTrimmedReads[2])
data["SequencingStats"]["MeanCoverage"] = round(meanCvg, 2)
data["SequencingStats"]["MeanReadLength"] = float(meanLengthRead[2])
data["SequencingStats"]["MeanInsertSize"] = float(meanInsertSize[3])
data["SequencingStats"]["InsertSizeSD"] = float(insertSizeSD[4])






with open('metrics.json', 'w') as f:
    json.dump(data, f, indent=4)
