#!/usr/bin/env python3

from pysam import AlignmentFile
import argparse
from threading import Thread
from queue import Queue
import csv


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


def bam_filename(sm_bc, prefix, suffix):
    return '{}.BC-{}.SM-{}.{}'.format(prefix, sm_bc[1], sm_bc[0], suffix)


def get_sample_map(manifest, order):
    bc_lines = open(order, 'rt').readlines()
    bc_i_map = {bc.strip(): i for i, bc in enumerate(bc_lines)}
    i_sm_bc_map = {}
    with open(manifest, 'rt') as handler:
        reader = csv.reader(handler, dialect=csv.excel_tab)
        assert next(reader) == ['sample', 'barcode', 'amplicons']
        for sample, barcode, amplicons in reader:
            i_sm_bc_map[bc_i_map[barcode]] = (sample, barcode)
    return i_sm_bc_map


def main(in_bam, manifest, order, out_bam):
    sm_bc_map = get_sample_map(manifest, order)
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af:
        q = Queue(5)
        af_out = AlignmentFile(out_bam, 'wb', check_sq=False, header=af.header)
        BamWriteThread(af_out, q).start()
        for read in af:
            sm, bc = sm_bc_map[tuple(read.get_tag('bc'))[0]]
            read.set_tag('BC', bc)
            read.set_tag('SM', sm)
            q.put(read)
        q.join()
        af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split bam files based on pacbio barocde tag.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('--manifest', required=True,
                        help='Tab delimited file with sample names in first column, barcode in second. '
                             'First row is header.')
    parser.add_argument('--bc-order', required=True,
                        help='Text file with order of barcodes in fasta file used by lima')
    parser.add_argument('--out', default='out.bam', required=True, help='output bam filename')

    args = parser.parse_args()
    main(args.in_bam, args.manifest, args.bc_order, args.out)
