import csv
from operator import itemgetter
import pandas as pd
import time
import bisect
from bisect import bisect_left
import sys, os

class Stack:
    def __init__(self):
        self.items = []

    def size(self):
        return len(self.items)

    def isEmpty(self):
        return not self.items

    def push(self, val):
        self.items.append(val)

    def top(self):
        return self.items[-1] if self.items else None

    def pop(self):
        return self.items.pop() if self.items else None

def bi_contains(lst, item):
    return bisect_left(lst, item)

def MergeIntervals(inputlist):
    if not inputlist:
        return []

    inputlist.sort(key=itemgetter(0))

    merged = []
    current_start, current_end = inputlist[0]

    for start, end in inputlist[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return merged

def CountTotalReadCount(chrom, exList, bam_list, position_row):
    totalCount = 0
    positions = pd.Series(position_row)
    reads = pd.Series([int(x[1]) for x in bam_list])

    for start, end in exList:
        pos1 = bisect_left(position_row, start)
        pos2 = bisect_left(position_row, end)

        if pos1 < len(position_row) and pos2 < len(position_row):
            if int(bam_list[pos2][0]) != end:
                pos2 -= 1

            totalCount += reads[pos1:pos2 + 1].sum()

    return totalCount



def writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list):
    targetRC = CountTotalReadCount(chrom, [(start, end)], bam_list, position_row)
    targetLength = end - start + 1

    averageTargetRC = targetRC / targetLength if targetLength else 0
    averageRCothers = (RC - targetRC) / (mergedExListLength - targetLength) if mergedExListLength != targetLength else 0

    writer_list.append((chrom, gene, start, end, targetRC, targetLength, RC, mergedExListLength, RC - targetRC, mergedExListLength - targetLength, averageTargetRC, averageRCothers))
    return writer_list

def Find_Novel_splicing_events(ChromDict_merged, ChromDict, chromosomes, AS, input_dir, species_folder, sample, output_dir):
    tt = time.time()
    AS_flag = []
    as_df = pd.read_csv(os.path.join(species, AS+'.csv'), delimiter='\t')
    writer_list = []
    output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'splicedExonReadCount(rc)', 'splicedExonLength(sl)', 'othersExonsReadCount(RC)', 'othersExonsLength(L)', 'RC - rc', 'L - sl', 'splicedExonAverageReadCoverage(n)', 'otherExonsAverageReadCoverage(N)']

    for chrom in chromosomes:
        # print("Starting:",chrom)
        tts = time.time()
        GeneDict = ChromDict[chrom]
        GeneDict_merged = ChromDict_merged[chrom]
        if os.path.getsize(os.path.join(input_dir, sample, chrom+".txt")) > 0:
            bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom+".txt"), delimiter='\t')
            position_row = bam_df.iloc[:, 0].tolist()
            bam_list = bam_df.values.tolist()

            for gene in GeneDict.keys():
                exonList = list(set(GeneDict[gene.upper()]))
                mergedExList = GeneDict_merged[gene.upper()]

                mergedExListLength = sum(end - start + 1 for start, end in mergedExList)
                RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

                for ex in exonList:
                    start, end = int(ex[0]), int(ex[1])

                    if (chrom, gene, start, end) not in AS_flag:
                        writer_list = writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list)
                        AS_flag.append((chrom, gene, start, end))

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, sample+"_"+AS+".csv"), sep='\t', index=False)

    print("Elapsed time: ", round(((time.time()-tt)/60), 2), "minutes")

def Find_splicing_events(ChromDict_merged, chromosomes, AS, input_dir, species, sample, output_dir):
    tt = time.time()
    AS_flag = []
    as_df = pd.read_csv(os.path.join(species, AS+'.csv'), delimiter='\t')
    writer_list = []
    output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'splicedExonReadCount(rc)', 'splicedExonLength(sl)', 'othersExonsReadCount(RC)', 'othersExonsLength(L)', 'RC - rc', 'L - sl', 'splicedExonAverageReadCoverage(n)', 'otherExonsAverageReadCoverage(N)']

    for chrom in chromosomes:
        # print("Starting:",chrom)
        tts = time.time()
        GeneDict = ChromDict_merged[chrom]
        if os.path.getsize(os.path.join(input_dir, sample, chrom+".txt")) > 0:
            bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom+".txt"), delimiter='\t')
            position_row = bam_df.iloc[:, 0].tolist()
            bam_list = bam_df.values.tolist()

            as_chr_rows = as_df[as_df['chrom'] == chrom]
            for ind1, t_row in as_chr_rows.iterrows():
                gene = t_row['gene'].strip().upper()
                mergedExList = GeneDict[gene]

                mergedExListLength = sum(end - start + 1 for start, end in mergedExList)
                RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

                if AS in ['SE', 'RI']:
                    exonStart, exonEnd = t_row['exonStart'], t_row['exonEnd']
                    if (chrom, gene, exonStart, exonEnd) not in AS_flag:
                        writer_list = writeResult(chrom, gene, exonStart, exonEnd, bam_list, position_row, RC, mergedExListLength, writer_list)
                        AS_flag.append((chrom, gene, exonStart, exonEnd))

                elif AS == 'MXE':
                    exon1Start, exon1End = t_row['exon1Start'], t_row['exon1End']
                    exon2Start, exon2End = t_row['exon2Start'], t_row['exon2End']
                    if (chrom, gene, exon1Start, exon1End) not in AS_flag:
                        writer_list = writeResult(chrom, gene, exon1Start, exon1End, bam_list, position_row, RC, mergedExListLength, writer_list)
                        AS_flag.append((chrom, gene, exon1Start, exon1End))

                    if (chrom, gene, exon2Start, exon2End) not in AS_flag:
                        writer_list = writeResult(chrom, gene, exon2Start, exon2End, bam_list, position_row, RC, mergedExListLength, writer_list)
                        AS_flag.append((chrom, gene, exon2Start, exon2End))

                else:
                    longExonStart, longExonEnd, shortExonStart, shortExonEnd, strand = t_row['longExonStart'], t_row['longExonEnd'], t_row['shortExonStart'], t_row['shortExonEnd'], t_row['strand']
                    start, end = 0, 0
                    if AS == 'A5SS':
                        if strand == '+':
                            start, end = longExonEnd + 1, shortExonEnd
                        else:
                            start, end = shortExonStart, longExonStart - 1

                    elif AS == 'A3SS':
                        if strand == '+':
                            start, end = shortExonStart, longExonStart - 1
                        else:
                            start, end = longExonEnd + 1, shortExonEnd

                    if (chrom, gene, start, end) not in AS_flag:
                        writer_list = writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list)
                        AS_flag.append((chrom, gene, start, end))

    df_out = pd.DataFrame(writer_list, columns=output_columns)
    df_out.to_csv(os.path.join(output_dir, sample+"_"+AS+".csv"), sep='\t', index=False)

    print("Elapsed time: ", round(((time.time()-tt)/60), 2), "minutes")


def MakeFullDictionary(ann_df, chromosomes):
    ChromDict = {}
    for chrom in chromosomes:
        chr_rows = ann_df[ann_df['chrom'] == chrom]
        gene_list = chr_rows['gene'].unique()
        GeneDict = {}
        for gene in gene_list:
            gene_rows = chr_rows[chr_rows['gene'] == gene]
            exon_starts = gene_rows['exonStarts'].str.split(',').explode()
            exon_ends = gene_rows['exonEnds'].str.split(',').explode()
            
            # Convert to integers, coercing errors to NaN and then dropping them
            exon_starts = pd.to_numeric(exon_starts, errors='coerce').dropna().astype(int)
            exon_ends = pd.to_numeric(exon_ends, errors='coerce').dropna().astype(int)
            
            exList = list(set(zip(exon_starts, exon_ends)))
            GeneDict[gene.strip().upper()] = exList
        ChromDict[chrom] = GeneDict
    return ChromDict


def merge_ChromDict(ChromDict, chromosomes):
    ChromDict_merged = {}
    for chrom in chromosomes:
        GeneDict_merged = {}
        GeneDict = ChromDict[chrom]
        for gene, exonList in GeneDict.items():
            mergedExonList = methods.MergeIntervals(exonList)
            GeneDict_merged[gene] = mergedExonList
        ChromDict_merged[chrom] = GeneDict_merged
    return ChromDict_merged

