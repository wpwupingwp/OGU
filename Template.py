#!/usr/bin/python3

from time import process_time
import argparse

def main():
    """docstring
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='./', 
                        help='target path, default is "./"')
    parser.print_help()
    arg = parser.parse_args()
    #start here
    
    #end
    print('Cost {:.3f}s.\n'.format(process_time()))

if __name__ == '__main__':
    main()
