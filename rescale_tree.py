from sys import argv
from os import system

def multiply_brlens(tree):

	tree_as_list = []

	for thing in tree:
		tree_as_list.append(thing)

	i = 0

	new_tree = []

	while i<len(tree_as_list):
	
		if tree_as_list[i] == ':' or tree_as_list[i] == '{' or tree_as_list[i] == ',':
			new_tree.append(tree_as_list[i])
			i+=1

			num = []

			if tree_as_list[i].isdigit() or tree_as_list[i] =='.':
				while tree_as_list[i].isdigit() or tree_as_list[i] =='.':
					num.append(tree_as_list[i])
					i+=1

				num = float(''.join(num))
				num = str(float(num*100))
				new_tree.append(num)

			else:
				new_tree.append(tree_as_list[i])
				i+=1

		else:
			new_tree.append(tree_as_list[i])
			i+=1

	new_tree = ''.join(new_tree)

	return new_tree


f = open (argv[1],"r")
tree = f.read()


remove_nexus = '#NEXUS\nBEGIN TREES;\n\n\tUTREE 1 = '
remove_end = '\nEND;'

tree2 = tree.replace(remove_nexus,'')
tree3 = tree2.replace(remove_end,'')
tree4 = tree3.replace(' ','')

rescaled_tree = multiply_brlens(tree4)

final_tree = tree.replace(tree3,rescaled_tree)

new_name = 'scaled_' + argv[1]

outfile = open(new_name, 'w')
outfile.write(final_tree)
outfile.close()
