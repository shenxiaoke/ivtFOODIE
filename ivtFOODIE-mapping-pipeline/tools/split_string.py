import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-N","--Num",type=int)
parser.add_argument("-n","--name",type=str)
parser.add_argument("-R","--read",type=str)
parser.add_argument("-i","--input_path",type=str)
parser.add_argument("-o","--output_path",type=str)

args = parser.parse_args()

name = args.name
N = args.Num
read = args.read
input_path = args.input_path
ouput_path = args.output_path

split_string = "fastqsplitter -i " + input_path + "/" + name + "_" + read + ".fastq.gz"

for i in range(N):
	split_string += " -o " + ouput_path + "/" + name + "-part" + str(i) + "/" + name + "-part" + str(i) + "_" + read + ".fastq.gz"

print(split_string)
