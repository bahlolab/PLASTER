#!/usr/bin/env python3
from pysam import AlignmentFile, AlignmentHeader
import argparse
from queue import Queue
from threading import Thread
import gzip
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


def get_rg_list(rg_dict, phases):
    rg_list = []
    for ps in phases:
        rg_list.append(rg_dict.copy())
        rg_list[-1].update({'ID': ps, 'SM': '{}_ph{}'.format(rg_dict['SM'], ps)})
    return rg_list


def main(in_bam, read_phases, as_zmw, update_rg, out_bam):
    qn_ps = {}
    ps_set = set()
    with gzip.open(read_phases, 'rt') if read_phases.endswith('.gz') else open(read_phases, 'rt') as handle:
        reader = csv.DictReader(handle, dialect='excel-tab')
        assert {'qname', 'phase'}.issubset(reader.fieldnames)
        for row in reader:
            if row['phase'] != 'NA':
                qn_ps[row['qname']] = row['phase']
                ps_set.add(row['phase'])

    with AlignmentFile(in_bam, 'rb', check_sq=False) as af_in:
        header_dict = af_in.header.to_dict()
        if update_rg:
            header_dict['RG'] = get_rg_list(header_dict['RG'][0], ps_set)
        af_out = AlignmentFile(out_bam, 'wb', header=AlignmentHeader.from_dict(header_dict))

        q = Queue(3)
        BamWriteThread(af_out, q).start()

        for read in af_in:
            qn = read.qname.split('/')[1] if as_zmw else read.qname
            if qn in qn_ps:
                if update_rg:
                    read.set_tag('RG', qn_ps[qn])
                q.put(read)
        q.join()
        af_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add reads to separate read groups based on phase set',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_bam', help='input bam file')
    parser.add_argument('read_phases', help='input tsv file with read phases')
    parser.add_argument('--zmw', action='store_true', help='report ZMW hole number instead of qname')
    parser.add_argument('--update-rg', action='store_true', help='add phase set to readgroup')
    parser.add_argument('-o', '--out', help='output bam filename, "-" for stdout', default='-')
    args = parser.parse_args()
    main(args.in_bam, args.read_phases, args.zmw, args.update_rg, args.out)
