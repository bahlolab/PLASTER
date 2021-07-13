#!/usr/bin/env python3

from pysam import AlignmentFile
import argparse
from threading import Thread
from queue import Queue
import json
import operator
import Levenshtein as Lev
import csv

rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
cigar_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'EQ': 7, 'X': 8, 'B': 9}


class BamWriteThread(Thread):
    def __init__(self, af, q):
        Thread.__init__(self, daemon=True)
        self.q = q
        self.af = af

    def run(self):
        while True:
            read = self.q.get()
            self.af.write(read)
            self.q.task_done()


def parse_manifest(manifest):
    smbc_amps = {}
    all_amps = set()
    with open(manifest, 'rt') as handler:
        reader = csv.reader(handler, dialect=csv.excel_tab)
        assert next(reader) == ['sample', 'barcode', 'amplicons']
        for sample, barcode, amplicons in reader:
            amps = amplicons.split(';')
            smbc_amps[(sample, barcode)] = amps
            all_amps.update(amps)
    return smbc_amps, all_amps


def match_pattern(subject, pattern, max_dist=1, refine=True):
    """
    returns the start and end position leftmost best match less within levenstein distance of max_dist, or None
    """
    ds = [Lev.distance(s, pattern) for s in [subject[i:i+len(pattern)]
                                             for i in range(max(0, len(subject) - len(pattern)))]]
    idx, min_dist = min(enumerate(ds), key=operator.itemgetter(1))
    end = min(idx + len(pattern), len(subject))
    if refine and min_dist > 0 and min_dist - max_dist <= max_dist:
        se = [(idx+i, end+j) for i in range(-max_dist, max_dist+1) for j in range(-max_dist, max_dist+1)
              if abs(i) + abs(j) <= max_dist and idx+i >= 0 and end+j <= len(subject)]
        ds = [Lev.distance(subject[s:e], pattern) for s, e in se]
        i, min_dist = min(enumerate(ds), key=operator.itemgetter(1))
        idx, end = se[i]
    if min_dist <= max_dist:
        return idx, end
    else:
        return None


def rev_comp(dna):
    rc = ''.join([rc_dict[dna[i]] for i in range(len(dna)-1, -1, -1)])
    return rc


def parse_amplicons_json(json_fn):
    in_dict = json.load(open(json_fn, 'rt'))
    amp_dict = {}
    for amp, data in in_dict.items():
        assert {'chrom', 'start', 'end', 'fwd_primer', 'rvs_primer', 'strand'}.issubset(data.keys())
        amp_dict[amp] = {'amplicon': amp, 'chrom': data['chrom'], 'start': data['start'], 'end': data['end'],
                         'len':  data['end'] - data['start'] + 1, 'fwd': data['fwd_primer'], 'rev': data['rvs_primer'],
                         'is_fwd': data['strand'] == '+'}
    return amp_dict


def overlap(r_start, r_end, q_start, q_end):
    ol_len = min(q_end, r_end) - max(q_start, r_start) + 1
    tot_len = max(q_end, r_end) - min(q_start, r_start) + 1
    return ol_len / tot_len


def ident(cigar_stats):
    return cigar_stats[cigar_map['EQ']] / (
            cigar_stats[cigar_map['EQ']] + cigar_stats[cigar_map['X']] + cigar_stats[cigar_map['I']] +
            cigar_stats[cigar_map['D']] + cigar_stats[cigar_map['N']])


def is_proper_pair(match, is_fwd, length):
    if is_fwd and {'fwd', 'rev_rc'}.issubset(match):
        return match['fwd'][1] <= length - match['rev_rc'][1]
    elif not is_fwd and {'rev', 'fwd_rc'}.issubset(match):
        return match['rev'][1] <= length - match['fwd_rc'][1]
    return False


def trim_primers(read, match, is_fwd, inclusive=True):
    """trim primers from reads and revert bam fields to unmapped"""
    i = 0 if inclusive else 1
    if is_fwd and not read.is_reverse:
        trim = (match['fwd'][i], read.query_length - match['rev_rc'][i])
    elif is_fwd and read.is_reverse:
        trim = (match['rev_rc'][i], read.query_length - match['fwd'][i])
    elif not is_fwd and not read.is_reverse:
        trim = (match['rev'][i], read.query_length - match['fwd_rc'][i])
    else:  # not is_fwd and read.is_reverse:
        trim = (match['fwd_rc'][i], read.query_length - match['rev'][i])
    quals = read.get_forward_qualities()[trim[0]:trim[1]]
    read.query_sequence = read.get_forward_sequence()[trim[0]:trim[1]]
    read.query_qualities = quals
    read.cigar = []
    read.flag = 4
    read.reference_name = '*'
    read.reference_start = 0
    return read


def main(in_bam, amp_json, manifest, window, max_dist, out):
    amp_dict = parse_amplicons_json(amp_json)
    smbc_amps, all_amps = parse_manifest(manifest)
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af:
        q = Queue(5)
        af_out = AlignmentFile(out, 'wb', check_sq=False, header=af.header)
        BamWriteThread(af_out, q).start()
        for read in af:
            if read.is_unmapped:
                read.set_tag('AM', 'unmapped')
            else:
                if read.has_tag('SM') and read.has_tag('BC'):
                    amp_set = smbc_amps[(read.get_tag('SM'), read.get_tag('BC'))]
                else:
                    amp_set = all_amps
                matches = []
                for amp in [amp_dict[a] for a in amp_set]:
                    if read.reference_name != amp['chrom']:
                        continue
                    ol = overlap(read.reference_start, read.reference_end, amp['start'], amp['end'])
                    if ol <= 0:
                        continue
                    match = {'amplicon': amp['amplicon'], 'overlap': ol}

                    # match left end of read
                    primer = 'fwd' if amp['is_fwd'] else 'rev'
                    mp = match_pattern(read.seq[:window], amp[primer], max_dist=max_dist, refine=True)
                    match.update({primer: mp}) if mp else None

                    # match rev_comp right end of read
                    primer = 'rev' if amp['is_fwd'] else 'fwd'
                    mp = match_pattern(rev_comp(read.seq[-window:]), amp[primer], max_dist=max_dist, refine=True)
                    primer = 'rev_rc' if amp['is_fwd'] else 'fwd_rc'
                    match.update({primer: mp}) if mp else None

                    match['pp'] = is_proper_pair(match, amp['is_fwd'], read.query_length)
                    ori = ','.join([x for x in ['fwd', 'rev', 'fwd_rc', 'rev_rc'] if x in match])
                    match['ori'] = ori if len(ori) else 'NA'
                    matches.append(match)
                if matches:
                    best = 0
                    if len(matches) > 1:
                        for i in range(1, len(matches)):
                            if matches[i]['pp'] is matches[best]['pp']:
                                if matches[i]['overlap'] > matches[best]['overlap']:
                                    best = i
                            elif matches[i]['pp']:
                                best = i
                    best_match = matches[best]
                    cigar_stats = read.get_cigar_stats()[0]
                    read.set_tag('AM', best_match['amplicon'])  # amplicon
                    read.set_tag('PP', '1' if best_match['pp'] else '0')  # proper pair
                    read.set_tag('PO', best_match['ori'])  # pair orientation
                    read.set_tag('OL', best_match['overlap'], 'f')  # Jaccard overlap
                    read.set_tag('ID', ident(cigar_stats), 'f')  # percent ident aligned
                    read.set_tag('PM', cigar_stats[cigar_map['EQ']] / amp_dict[best_match['amplicon']]['len'], 'f') # coverage - i.e. percent amp based covered
                    if best_match['pp']:
                        read = trim_primers(read, best_match, amp_dict[best_match['amplicon']]['is_fwd'])
                else:
                    read.set_tag('AM', 'off-target')
            q.put(read)
        q.join()
        af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split aligned bam files based on overlap with genomic targets. '
                                                 'percent matched of target is recorded in PM tag in bam file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('--out', required=True, help='name of output bam file', )
    parser.add_argument('--amplicons', required=True, help='json file describing input amplicons')
    parser.add_argument('--window', type=int, help='window at end of reads to search for primers', default=500)
    parser.add_argument('--max-dist', type=int, help='maximum Levenshtein distance to primer sequence for match',
                        default=2)
    parser.add_argument('--manifest', required=True,
                        help='Tab delimited file with sample names in first column, barcode in second. '
                             'First row is header.')
    args = parser.parse_args()
    main(args.in_bam, args.amplicons, args.manifest, args.window, args.max_dist, args.out)
