#the script compares the majority and minority variants tabulars and removes the lines from minority tabular that are in the majority tabular
import sys
import argparse

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--input_max_tab', dest='input_max_tab', help='input majority tab file')
    parser.add_argument('-2', '--input_min_tab', dest='input_min_tab', help='input minority tab file')
    parser.add_argument('-o', '--output_min_tab', dest='output_min_tab', help='output minority tab file')
    args = parser.parse_args()

    file_max = open(args.input_max_tab)
    file_min = open(args.input_min_tab)
    read_file_max = file_max.readlines()
    read_file_min = file_min.readlines()
    out_file = open(args.output_min_tab,'w')

    for element in read_file_min[1:]:
    	if element in read_file_max:
    		read_file_min.remove(element)
    
    for line in read_file_min:
    	out_file.write(line)
    out_file.close()
    file_max.close()
    file_min.close()


if __name__ == "__main__":
    __main__()
