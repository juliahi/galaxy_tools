#!/usr/bin/python

import csv, argparse, sys

if __name__ == '__main__': 
     
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files",help = "Files with results od quality controll ", nargs="+",  default = '' )
    parser.add_argument("-o", "--output",help="The name of output file", default='' )
    parser.add_argument("-n", "--names",help="The labels of input files", default='', nargs="*" )
    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        print parser.format_help() #prints help if no arguments
        sys.exit(1)
        
    inputfiles = args.files
    if len(args.names) < len(inputfiles):
        names = inputfiles
    else:
        names = args.names
    
    
    all_stat = False
    
    probs = ["Stat_name"]
    data = {}
    for plik in inputfiles:
        probs.append(names.pop(0))
        inp_file = open(plik)
        for line in inp_file:
            if line == "\n": continue
            elif "Left reads" in line:
                left = True
                continue
            elif "Right reads" in line:
                left = False
                continue
            elif "Input" in line:
                if left:
                    if "Input_left" not in data.keys():
                        data["Input_left"] = [line.split()[1]]
                    else: data["Input_left"].append(line.split()[1])
                else: 
                    if "Input_right" not in data.keys():
                        data["Input_right"] = [line.split()[1]]
                    else: data["Input_right"].append(line.split()[1])
            elif "Mapped" in line:
                if left:
                    if "Mapped_left" not in data.keys():
                        data["Mapped_left"] = [line.split("(")[1].split("%")[0].strip()]
                    else: data["Mapped_left"].append(line.split("(")[1].split("%")[0].strip())
                else: 
                    if "Mapped_right" not in data.keys():
                        data["Mapped_right"] = [line.split("(")[1].split("%")[0].strip()]
                    else: data["Mapped_right"].append(line.split("(")[1].split("%")[0].strip())
            elif "of these" in line:
                if all_stat:
                    if "Pairs_mult_align" not in data.keys():
                        data["Pairs_mult_align"] = [line.split("(")[1].split("%")[0].strip()]
                    else: data["Pairs_mult_align"].append(line.split("(")[1].split("%")[0].strip())
                elif left: 
                    if "Left_mult_align" not in data.keys():
                        data["Left_mult_align"] = [line.split("(")[1].split("%")[0].strip()]
                    else: data["Left_mult_align"].append(line.split("(")[1].split("%")[0].strip())
                else:
                    if "Right_mult_align" not in data.keys():
                        data["Right_mult_align"] = [line.split("(")[1].split("%")[0].strip()]
                    else: data["Right_mult_align"].append(line.split("(")[1].split("%")[0].strip())
            elif "overall" in line:
                all_stat = True
                if "Overall_reads" not in data.keys():
                    data["Overall_reads"] = [line.split("%")[0]]
                else: data["Overall_reads"].append(line.split("%")[0])
            elif "Aligned pairs" in line:
                if "Aligned_pairs" not in data.keys():
                    data["Aligned_pairs"] = [line.split()[2]]
                else: data["Aligned_pairs"].append(line.split()[2])
            elif "discordant" in line:
                if "Discordant" not in data.keys():
                    data["Discordant"] = [line.split("(")[1].split("%")[0].strip()]
                else: data["Discordant"].append(line.split("(")[1].split("%")[0].strip())
            elif "concordant" in line:
                all_stat = False
                if "Concordant" not in data.keys():
                    data["Concordant"] = [line.split("%")[0].strip()]
                else: data["Concordant"].append(line.split("%")[0].strip())
            else: print "Nieznana linia", line    
            
    names = ["Input_left", "Mapped_left", "Left_mult_align", "Input_right", "Mapped_right", "Right_mult_align", "Overall_reads", "Aligned_pairs", "Pairs_mult_align", "Discordant", "Concordant"]
    
    if args.output != '':
        output = args.output
    else: output = "QC_statistics.csv"
    with open(output, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(probs)
        for name in names:
            spamwriter.writerow([name] + data[name])
        #spamwriter.writerow(["Mapped_left"] + data["Mapped_left"])
                
