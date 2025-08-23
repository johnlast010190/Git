#! /bin/bash

entries=(
#bad       #good
"[ \t]*$"  "" #trailing whitespaces
"forAll (" "forAll("
"Info <<"  "Info<<" 
"Pout <<"  "Pout<<" 
" if("     " if ("
" for("    " for ("
" while("  " while ("
"^if("     "if ("
"^for("    "for ("
"^while("  "while ("
"> > > >"  ">>>>"
"> > >"    ">>>"
"> >"      ">>"
)

for ((i=0; i<${#entries[@]}; i+=2)); do
  bad=${entries[$i]}
  good=${entries[$((i+1))]}
  #echo bad: "$bad"
  #echo good: "$good"

  sed -i "s/$bad/$good/g" $1
done

# remove spaces inside parentheses
array=(
" if ("
" for ("
" while ("
"^if ("
"^for ("
"^while ("
)

for i in "${array[@]}"; do
  #echo "$i"
  sed -i -e "/$i/ s/( /(/; /$i/ s/ )/)/" $1
done


