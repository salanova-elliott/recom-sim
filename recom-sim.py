from __future__ import print_function
import argparse, random, re, math, sys

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("import_file", type=str, help="filepath to genepop file")
parser.add_argument("introgress", type=int, choices=[1,2,3], help="level of introgression to simulate")
parser.add_argument("--num-offs", type=int, default=100, help="number of offspring (def = 100)")
parser.add_argument("--p1name", type=str, default="POP1", help="name of population 1 (def = POP1)")
parser.add_argument("--p2name", type=str, default="POP2", help="name of population 2 (def = POP2)")
parser.add_argument("--exclude", action="store_true", help="exclude parentals from output")
parser.add_argument("--out", type=str, default="out", help="name of output file (def = out.txt)")
args = parser.parse_args()

#Parses genepop file collecting loci, samples, and data
def read_genepop(genepop_file):

	names_pop1 = []
	names_pop2 = []
	allele_list_pop1 = []
	allele_list_pop2 = []
	loci_list = []

	sector_loci = True
	current_pop = 1
	pop_count = 0

	pop_catch = '(?:Pop|pop|POP)\s'
	sample_catch = '(\w+.*)\,(?:((?:\s*\d{6})+)|((?:\s*\d{4})+))'
	comma_loci_catch = '((?:\w+\,)+\w+)'

	#Skips header
	genepop_file.readline()

	#Main loop
	for line in genepop_file:

		pop_search = re.search(pop_catch, line)
		sample_search = re.search(sample_catch, line)

		#If comma separated format
		if sector_loci and "," in line:
			print(line)
			loci_search = re.search(comma_loci_catch, line)
			loci_list = loci_search.group(1).split(",")
			sector_loci = False

		#Elif each loci on new line
		elif not pop_search and "," not in line:

			loci_list.append(line.rstrip())

		#Sample capture
		elif sample_search:

			if current_pop == 1:
				names_pop1.append(sample_search.group(1).rstrip())
				allele_list_pop1.append(sample_search.group(2).lstrip().split(" "))
			if current_pop == 2:
				names_pop2.append(sample_search.group(1).rstrip())
				allele_list_pop2.append(sample_search.group(2).lstrip().split(" "))

		#POP line
		if pop_search:
			pop_count += 1
			sector_loci = False

			if pop_count == 2:
				current_pop = 2

			#If more than two reference populations
			if pop_count == 3:
				print("Woah there. Looks like there might be more than two reference pops in your file.")
				return

	return loci_list, allele_list_pop1, allele_list_pop2, names_pop1, names_pop2

#Allele frequency matrix constructor
def build_matrix(raw_alleles, locus_list):

	allele_freqs = []

	for i in range(len(raw_alleles[0])):

		#Progress display
		print("Reading locus: " + locus_list[i], end="")
		sys.stdout.flush()
		print("\r", end = "")

		allele_a = ""
		allele_b = ""

		count_a = 0
		count_b = 0

		for sample in raw_alleles:

			a = sample[i][0:len(sample[i])/2]
			b = sample[i][len(sample[i])/2:len(sample[i])]

			allele_a, allele_b, count_a, count_b = allele_check(a, allele_a, allele_b, count_a, count_b)
			allele_a, allele_b, count_a, count_b = allele_check(b, allele_a, allele_b, count_a, count_b)

		allele_freqs.append([count_a, count_b, allele_a, allele_b])

	return allele_freqs

#Helper for build_matrix
def allele_check(exam, allele_a, allele_b, count_a, count_b):

	if exam == allele_a:
		count_a += 1
		return allele_a, allele_b, count_a, count_b
	elif exam == allele_b:
		count_b += 1
		return allele_a, allele_b, count_a, count_b
	elif exam != "00" and exam != "000" and not allele_a:
		allele_a = exam
		count_a += 1
		return allele_a, allele_b, count_a, count_b
	elif exam != "00" and exam != "000" and not allele_b:
		allele_b = exam
		count_b += 1
		return allele_a, allele_b, count_a, count_b

	#In case of "00" or "000"
	return allele_a, allele_b, count_a, count_b

#Generates offspring
def offspring_gen(pop1_freq, pop2_freq):
	offspring_m = []
	for i in range(args.num_offs):

		#Progress display
		print("Simulating offspring: " + str(i + 1), end="")
		sys.stdout.flush()
		print("\r", end = "")

		offspring_i = []

		for j in range(len(pop1_freq)):
			locus = ""
			locus += allele_select(pop1_freq, j)
			locus += allele_select(pop2_freq, j)

			offspring_i.append(locus)

		offspring_m.append(offspring_i)

	return offspring_m

#Helper for offspring_gen. Actually selects alleles
def allele_select(pop_freq, j):
	allele_string = ""

	if not pop_freq[j][3]:
		allele_string = pop_freq[j][2]
	elif random.randrange(1,pop_freq[j][0]+pop_freq[j][1]+1) <= pop_freq[j][0]:
		allele_string = pop_freq[j][2]
	else:
		allele_string = pop_freq[j][3]

	return allele_string

#Output
def output_file(loc_list, pop1_raw, pop2_raw, pop1_names, pop2_names, f1_off):
	file_out = open(args.out + ".txt", "w")
	file_out.write("Simulated with recom-sim: https://github.com/salanova-elliott/recom-sim\n")

	for name in loc_list:
		file_out.write(name + "\n")

	#Parentals
	if not args.exclude:

		file_out.write("POP\n")
		for i in range(len(pop1_names)):
			file_out.write(pop1_names[i] + ",  ")
			for locus in pop1_raw[i]:
				file_out.write(locus + " ")
			file_out.write("\n")

		file_out.write("POP\n")
		for i in range(len(pop2_names)):
			file_out.write(pop2_names[i] + ",  ")
			for locus in pop2_raw[i]:
				file_out.write(locus + " ")
			file_out.write("\n")

	#Hybs
	for i in range(len(int_off)):
		#Filters out labels from int_off
		if i % 2 == 0:
			prefix = int_off[i]
			continue

		file_out.write("POP\n")
		for j in range(len(int_off[i])):
			if j < 9:
				zeroes = "00"
			elif j < 99:
				zeroes = "0"
			else:
				zeroes = ""

			file_out.write(prefix + zeroes + str(j + 1) + ",  ")
			for locus in int_off[i][j]:
				file_out.write(locus + " ")
			file_out.write("\n")

#Execution
file_object = open(args.import_file)

loc_list, pop1_raw, pop2_raw, pop1_names, pop2_names = read_genepop(file_object)
pop1_freq = build_matrix(pop1_raw, loc_list)
pop2_freq = build_matrix(pop2_raw, loc_list)

int_off = []
int_off.append("F1HYB_")
int_off.append(offspring_gen(pop1_freq, pop2_freq))

if args.introgress == 2 or args.introgress == 3:
	f1_freq = build_matrix(int_off[1], loc_list)

	int_off.append("B1" + args.p1name + "_")
	int_off.append(offspring_gen(pop1_freq, f1_freq))
	int_off.append("B1" + args.p2name + "_")
	int_off.append(offspring_gen(pop2_freq, f1_freq))

if args.introgress == 3:
	b1pop1_freq = build_matrix(int_off[3], loc_list)
	b1pop2_freq = build_matrix(int_off[5], loc_list)

	int_off.append("F2HYB_")
	int_off.append(offspring_gen(f1_freq, f1_freq))
	int_off.append("B2" + args.p1name + "_")
	int_off.append(offspring_gen(pop1_freq, b1pop1_freq))
	int_off.append("B2" + args.p2name + "_")
	int_off.append(offspring_gen(pop2_freq, b1pop2_freq))

output_file(loc_list, pop1_raw, pop2_raw, pop1_names, pop2_names, int_off)
