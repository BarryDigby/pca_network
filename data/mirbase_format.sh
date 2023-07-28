#!/bin/bash

input_file="mirbase_aliases.txt"
output_file="output.txt"

prev_key=""
prev_line=""

while IFS=$'\t' read -r key values; do
    IFS=';' read -ra arr <<< "$values"
    length=${#arr[@]}

    if [[ "$prev_key" != "$key" ]]; then
        # New grouping key found, update previous line with the correct number of tabs
        if [[ -n "$prev_line" ]]; then
            printf "%s\n" "$prev_line" >> "$output_file"
        fi
        prev_key="$key"
        prev_line="$key\t${arr[0]}"
    else
        # Same grouping key, add tabs to the previous line for the additional values
        prev_line+="\t${arr[0]}"
    fi

    for (( i = 1; i < length; i++ )); do
        # For each additional value, add a new line with the grouping key and the value
        printf "%s\t\t%s\n" "$key" "${arr[i]}" >> "$output_file"
    done
done < "$input_file"

# Add the last processed line after the loop ends
if [[ -n "$prev_line" ]]; then
    printf "%s\n" "$prev_line" >> "$output_file"
fi
