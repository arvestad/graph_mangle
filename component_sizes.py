#! /usr/bin/env python3

import sys
import json
import argparse
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("component_file", help="Input file containing components in JSON format.")
    args=parser.parse_args()

    with open(args.component_file) as ih:
        data = json.load(ih)
        components = data['components']
        sizes = map(len, components)
        for c in sizes:
            print(c)
        plt.hist(sizes, 50)
        plt.xlabel("Component sizes")
        plt.savefig("component_size_histo.pdf")

if __name__ == '__main__':
    main()
