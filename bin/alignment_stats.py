#!/usr/bin/env python3

from pysam import AlignmentFile
import argparse
import sys
import gzip

cigar_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'EQ': 7, 'X': 8, 'B': 9}


def get_tag_or_na(read, tag):
    try:
        return str(read.get_tag(tag))
    except KeyError:
        return 'NA'


def main(in_bam, out, as_zmw, tags, ext):
    if ext is not None:
        key, value = ext.split(':', 1)
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af_in:
        header = ['query_name', 'query_length', 'query_alignment_start', 'query_alignment_end', 'reference_name',
                  'reference_length', 'reference_start', 'reference_end', 'match', 'mismatch', 'ins', 'del']
        if tags:
            header.extend(tags)
        if ext is not None:
            header.append(key)
        out.write('\t'.join(header) + '\n')
        for read in af_in:
            cigar_stats = read.get_cigar_stats()[0]
            if as_zmw:
                qname = read.query_name.split('/')[1]
            else:
                qname = read.query_name
            fields = [str(qname), str(read.query_length), str(read.query_alignment_start), str(read.query_alignment_end),
                      str(read.reference_name), str(read.reference_length), str(read.reference_start),
                      str(read.reference_end), str(cigar_stats[cigar_map['EQ']]), str(cigar_stats[cigar_map['X']]),
                      str(cigar_stats[cigar_map['I']]),
                      str(cigar_stats[cigar_map['D']] + cigar_stats[cigar_map['N']])]
            if tags:
                fields.extend([get_tag_or_na(read, tag) for tag in tags])
            if ext is not None:
                fields.append(value)
            out.write('\t'.join(fields) + '\n')
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input bam file, \'-\' for stdin')
    parser.add_argument('-o', '--out', default='-', help='output .tsv(.gz) file, \'-\' for stdout')
    parser.add_argument('--zmw', action='store_true', help='report ZMW hole number instead of qname')
    parser.add_argument('--ext', help='additional metadata to include, in form key:value')
    parser.add_argument('--tag', action='append', help='additional tags to include')
    args = parser.parse_args()

    if args.out is '-':
        out = sys.stdout
    elif args.out.endswith('.gz') or args.out.endswith('.gzip'):
        out = gzip.open(args.out, 'wt')
    else:
        out = open(args.out, 'wt')

    main(args.input, out, args.zmw, args.tag, args.ext)
