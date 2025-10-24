#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command line interface for PymHM
"""

import argparse
import sys
from pathlib import Path


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description="PymHM - Python Mesoscale Hydrological Model",
        prog="pymhm"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0"
    )
    
    parser.add_argument(
        "--info",
        action="store_true",
        help="Show package information"
    )
    
    args = parser.parse_args()
    
    if args.info:
        print("PymHM - Python Mesoscale Hydrological Model")
        print("Version: 0.1.0")
        print("Author: Sanjeev Bashyal")
        print("Email: sanjeev.bashyal01@gmail.com")
        print("Description: Python package for mesoscale Hydrological Model")
        print("Homepage: https://github.com/SanjeevBashyal/pymhm")
        return 0
    
    print("PymHM CLI - Use --help for more options")
    return 0


if __name__ == "__main__":
    sys.exit(main())
