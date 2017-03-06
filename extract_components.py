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
            sz = len(c)
            if sz >= minsize and sz <= maxsize:
                selected_components.append(c)

    output = { "graph_file" : "selected from " + args.component_file,
               "components" : selected_components
    }

    print(json.dumps(output, indent=4))

if __name__ == '__main__':
    main()
