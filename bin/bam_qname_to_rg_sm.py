#!/usr/bin/env python3
from pysam import AlignmentFile, AlignmentHeader
import argparse
from queue import Queue
from threading import Thread


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


def get_rg_list(rg_dict, qnames):
    rg_list = []
    for qn in qnames:
        rg_list.append(rg_dict.copy())
        rg_list[-1].update({'ID': qn, 'SM': qn})
    return rg_list


def main(in_bam, out_bam, as_zmw):
    with AlignmentFile(in_bam, 'rb', check_sq=False) as af_in:
        if as_zmw:
            qnames = [read.qname.split('/')[1] for read in af_in]
        else:
            qnames = [read.qname for read in af_in]
        af_in.reset()

        header_dict = af_in.header.to_dict()
        header_dict['RG'] = get_rg_list(header_dict['RG'][0], qnames)

        af_out = AlignmentFile(out_bam, 'wb', header=AlignmentHeader.from_dict(header_dict))
        q = Queue(3)
        BamWriteThread(af_out, q).start()

        for read in af_in:
            if as_zmw:
                read.set_tag('RG', read.qname.split('/')[1])
            else:
                read.set_tag('RG', read.qname)
            q.put(read)
        q.join()
        af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Set RG ID to qname for each read and append qname to SM.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('--zmw', action='store_true', help='report ZMW hole number instead of qname')
    parser.add_argument('-o', '--out', help='output bam filename, "-" for stdout', default='-')
    args = parser.parse_args()
    main(args.in_bam, args.out, args.zmw)
