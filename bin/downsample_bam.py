#!/usr/bin/env python3
from pysam import AlignmentFile
import argparse
from queue import Queue
from threading import Thread
from random import sample, seed
from sys import stderr


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


def main(in_bam, out_bam, fraction, count, seed_val, uncompressed):
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af_in:
        total_count = af_in.count()
        if fraction:
            count = round(total_count * fraction)
        if count > total_count:
            stderr.write("Warning - selecting all reads\n")
            count = total_count
        seed(seed_val)
        at = sample(range(total_count), count)
        at.sort()

        mode = 'wbu' if uncompressed else 'wb'
        af_out = AlignmentFile(out_bam, mode, header=af_in.header)
        q = Queue(3)
        BamWriteThread(af_out, q).start()

        af_in.reset()
        read_idx = 0
        at_idx = 0

        for read in af_in:
            if read_idx == at[at_idx]:
                q.put(read)
                at_idx += 1
                if at_idx >= count:
                    break
            read_idx += 1
        q.join()

    af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Take a random subsample of reads from a bam file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='input bam file')
    parser.add_argument('-o', '--out', help='output bam filename, "-" for stdout', default='-')
    parser.add_argument('--fraction', type=float, help='fraction of reed to keep')
    parser.add_argument('--count', type=int, help='number of reads to keep')
    parser.add_argument('--seed', type=int, default=0, help='integer seed for random')
    parser.add_argument('-u', action='store_true', help='enable uncompressed bam output')
    args = parser.parse_args()
    if not (args.fraction or args.count):
        print('Error: one of arguments --count or --fraction is required.')
        exit(1)
    main(args.input, args.out, args.fraction, args.count, args.seed, args.u)

