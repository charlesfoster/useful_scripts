#!/bin/sh
# GNU All-Permissive License

SDQUO=$(echo -ne '\u2018\u2019')
RDQUO=$(echo -ne '\u201C\u201D')

HORIZONTAL="===================================================="

FIX_QUOTES () {
	echo -e "\n${HORIZONTAL}"
	echo -e "Paste in the text with dodgy quote symbols, then press enter\n"
	read -p "> " DODGY
	FIXED=$(echo $DODGY | sed -e "s|[${SDQUO}]|\'|g" -e "s|[${RDQUO}]|\"|g")
	echo -e "\nFixed:\n"
	echo $FIXED
	echo -e "${HORIZONTAL}\n"
}

export -f FIX_QUOTES

FIX_QUOTES

exit 0
