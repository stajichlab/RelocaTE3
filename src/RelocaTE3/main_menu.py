#!/usr/bin/env python3
"""RelocaTE3 Interface and CLI logic."""

import argparse
import sys


def main():
    """Primary RelocaTE3 CLI logic."""
    parser = argparse.ArgumentParser(description='EXAMPLE')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        default=sys.stdin,
                        help='Input file (or by stdin)')
