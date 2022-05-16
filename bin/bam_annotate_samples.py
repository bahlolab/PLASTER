#!/usr/bin/env python3
import sys
from pysam import AlignmentFile, FastxFile
import argparse
from threading import Thread
from queue import Queue


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


def main(in_bam, barcodes, sample, out_bam):
    if barcodes:
        barcodes = [(r.name, r.sequence) for r in FastxFile(barcodes)]
    else:
        sm = sample
        bc = 'N'
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af:
        q = Queue(5)
        af_out = AlignmentFile(out_bam, 'wb', check_sq=False, header=af.header)
        BamWriteThread(af_out, q).start()
        for read in af:
            if barcodes:
                sm, bc = barcodes[tuple(read.get_tag('bc'))[0]]
            read.set_tag('BC', bc)
            read.set_tag('SM', sm)
            q.put(read)
        q.join()
        af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split bam files based on pacbio barocde tag.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('--barcodes', help='Fasta file containing barcodes')
    parser.add_argument('--sample', help='Sample ID to set if not using barcodes')
    parser.add_argument('--out', default='out.bam', required=True, help='output bam filename')
    args = parser.parse_args()
    if not (args.barcodes or args.sample):
        sys.exit(1)
    main(args.in_bam, args.barcodes, args.sample, args.out)
