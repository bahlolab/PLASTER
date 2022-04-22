#!/usr/bin/env python3
from pysam import AlignmentFile
import argparse
from queue import Queue
from threading import Thread


class BamWriteThread(Thread):
    def __init__(self, q, af):
        Thread.__init__(self, daemon=True)
        self.q = q
        self.af = af

    def run(self):
        while True:
            read = self.q.get()
            self.af.write(read)
            self.q.task_done()


def main(super_bam, sub_bam, output_bam, as_intersect, uncompressed, count):
    with AlignmentFile(super_bam, 'rb', check_sq=False) as af_super:
        mode = 'wbu' if uncompressed else 'wb'
        af_out = AlignmentFile(output_bam, mode, header=af_super.header)
        q = Queue(3)
        BamWriteThread(q, af_out).start()
        n = 0
        with AlignmentFile(sub_bam, 'rb', check_sq=False) as af_sub:
            sub_zmw = int(next(af_sub).query_name.split('/')[1])
            for read in af_super:
                zmw = int(read.query_name.split('/')[1])
                while sub_zmw < zmw:
                    try:
                        sub_zmw = int(next(af_sub).query_name.split('/')[1])
                    except StopIteration:
                        break
                if zmw == sub_zmw:
                    if as_intersect:
                        q.put(read)
                elif not as_intersect:
                    q.put(read)
                    n += 1
        if not as_intersect:
            for read in af_super:
                q.put(read)
                n += 1
        q.join()
        af_out.close()
        if count:
            print(n, file=open(count, 'wt'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('super_bam', help='input superset pb bam file, sorted by zmw')
    parser.add_argument('sub_bam', help='intput subset pb bam, sorted by zmw')
    parser.add_argument('--intersect', help='output set intersection instead of set difference', action='store_true')
    parser.add_argument('-o', '--out', help='output bam filename, "-" for stdout', default='-')
    parser.add_argument('-u', action='store_true', help='enable uncompressed bam output')
    parser.add_argument('--count', help='output count filename', default=None)
    args = parser.parse_args()
    main(args.super_bam, args.sub_bam, args.out, args.intersect, args.u, args.count)
