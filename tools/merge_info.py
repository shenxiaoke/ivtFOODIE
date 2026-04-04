import sys, os

in_file_path, out_file_path, out_file_name = sys.stdin.readline().strip().split(" ")

conversion_array = [0]*32

with open(out_file_path+ "/" + out_file_name, "w") as fto:
	in_files = os.listdir(in_file_path)
	for in_file in in_files:
		with open(in_file_path + "/" + in_file) as fin:
			try:
				# _ = next(fin)
				line = next(fin).strip().split("\t")
				for i in range(32):
					conversion_array[i] += int(line[i])
			except:
				pass
	fto.write("\t".join([str(conversion) for conversion in conversion_array]))
