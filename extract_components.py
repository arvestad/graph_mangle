#! /usr/bin/env python3

import sys
import json
import argparse


def main():
    minsize=0
    maxsize=1000000

    parser = argparse.ArgumentParser()
    parser.add_argument("component_file", help="Input file containing components in JSON format.")
    parser.add_argument("-m", "--minsize", type=int, help="Minimum component size")
    parser.add_argument("-M", "--maxsize", type=int, help="Maximum component size")
    parser.add_argument("-l", "--lst", action="store_true", help="Use the simple format of one contig id per line, instead of the default JSON output. Output goes to files named 'component_<num>.lst', where <num> is an integer.")
    args=parser.parse_args()

    if args.minsize:
        minsize = args.minsize

    if args.maxsize:
        maxsize = args.maxsize
    
    selected_components = []
    with open(args.component_file) as ih:
        data = json.load(ih)
        components = data['components']
        for c in components:
            sz = len(c['nodes'])
            if sz >= minsize and sz <= maxsize:
                selected_components.append(c)

    if args.lst:
        i = 0
        for c in selected_components:
            filename = 'component_' + str(i) + '.lst'
            with open(filename, "w") as oh:
                for contig in c['nodes']:
                    oh.write(contig)
                    oh.write("\n")
            i += 1
    else:
        output = { "graph_file" : "selected from " + args.component_file,
                   "components" : selected_components
        }
        print(json.dumps(output, indent=4))

if __name__ == '__main__':
    main()
