#!/usr/bin/env python2.7
# This is a script to generate the correct command for an analysis with RAxML
# It saves you remembering what a lot of the common commands do
# Additionally, you can choose the option to run the analysis immediately after creating the input
# Usage: ./raxml_script.py

import os

test_version = raw_input("I'm going to assume you're using 'raxmlHPC-PTHREADS-AVX'. Is this correct? Y or N?\n> ")

if test_version == "Y" or test_version == "y":
    version = "raxmlHPC-PTHREADS-AVX"

elif test_version == "N" or test_version == "n":
    version = raw_input("Well, what version are you using then?\n> ")

elif not test_version == "Y" or test_version == "y" or test_version == "N" or test_version == "n":
    print "Why didn't you answer correctly? I'll have to quit now."
    quit()

threads = raw_input("How many threads do you wanna use?\n")

print "I am going to envoke the rapid bootstrapping algorithm, if you don't like this, you can write your own script.\n*****WARNING: this will probably only work with the latest version(s) of RAxML.*****\nThe options, seeds and model are:\n"
first_option = "-f a -x 19828 -p 1291 -m GTRGAMMA"
print first_option
print "Note: you can change these seeds afterwards if you want.\n"
test_outgroup = raw_input("\nDo you want to use an outgroup? Y or N?\n> ")

if test_outgroup == "Y" or test_outgroup == "y":
    outgroup = raw_input("What is your outgroup? Enter the names of the taxa followed by a comma. If you don't have one and still (strangely) chose this option, leave this blank)\n> ")

elif test_outgroup == "N" or test_outgroup == "n":
    outgroup = ""

elif not test_outgroup == "Y" or test_outgroup == "y" or test_outgroup == "N" or test_outgroup == "n":
    print "You didn't answer correctly, so I'm going to assume you meant no...\n"
    outgroup = ""
    
partitions = raw_input("What is the name of your partitions file?\n> ")

alignment = raw_input("What is the name of your alignment file?\n> ")

analysis_name = raw_input("What do you want to call the analysis?\n> ")

boots = raw_input("How many bootstrap replicates?\n> ")

script_name = raw_input("What should we call the script?\n> ")

print "OK, so your analysis will be run like so:\n"

full_analysis = "./%s -T %s %s -o %s -q %s -s %s -n %s -N %s -k" % (version, threads, first_option, outgroup, partitions, alignment, analysis_name, boots)
print full_analysis

print "\nAnd your script will be called: %s" % script_name

checkpoint = raw_input("Is this OK? Y or N?\n> ")

if checkpoint == "Y" or checkpoint == "y":
    open(script_name,'w').writelines(full_analysis) 
    print "OK, I've written your command to file."
    os.system("say I love you") 

elif checkpoint == "N" or checkpoint == "n": 
    print "Well, you must have made a mistake. Start again! \a\a\a\a\a\a\a"

run_analysis_query = raw_input("Do you want to run the analysis now? Y or N?\n> ")

def run_analysis():
    run_script = "sh %s" % script_name          
    os.system(run_script)
    
if run_analysis_query == "Y" or run_analysis_query == "y":
    run_analysis()

elif not run_analysis_query == "Y" or run_analysis_query == "y":
    quit()
