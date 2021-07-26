#!/usr/bin/env python3

from pysam import AlignmentFile, AlignmentHeader
import argparse
from threading import Thread
from queue import Queue


class BamWriteThread(Thread):
    def __init__(self, q):
        Thread.__init__(self, daemon=True)
        self.q = q

    def run(self):
        while True:
            af, read = self.q.get()
            af.write(read)
            self.q.task_done()


def bam_file(prefix, lb_tag, sample, amplcion, header_dict):
    header_dict['RG'] = header_dict['RG'][0:1]
    header_dict['RG'][0].update({'ID': '1', 'SM': sample, 'LB': lb_tag})
    fn = '{}LB-{}.SM-{}.AM-{}.bam'.format(prefix, lb_tag, sample, amplcion)
    return AlignmentFile(fn, 'wb', check_sq=False, header=AlignmentHeader.from_dict(header_dict))


def get_tag(read, tag):
    try:
        return str(read.get_tag(tag))
    except KeyError:
        return None


def main(in_bam, lb_tag, prefix):
    out_bams = {}
    counts = {}
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af:
        header_dict = af.header.to_dict()
        q = Queue(5)
        BamWriteThread(q).start()
        for read in af:
            if get_tag(read, 'PP') == '1':
                read.set_tag('RG', '1')
                sm_am = (read.get_tag('SM'), read.get_tag('AM'))
                if sm_am not in out_bams:
                    out_bams[sm_am] = bam_file(prefix, lb_tag, *sm_am, header_dict)
                    counts[sm_am] = 0
                q.put((out_bams[sm_am], read))
                counts[sm_am] += 1
        q.join()
    for k, bam in out_bams.items():
        bam.close()
        print(counts[k], file=open(bam.filename.decode() + '.count', 'wt'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split bam file by SM tag',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('--lb-tag', help='library tag', required=True)
    parser.add_argument('--prefix', help='prefix for output bam file', default='')

    args = parser.parse_args()
    main(args.in_bam, args.lb_tag, args.prefix)
