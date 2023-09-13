import os
import numpy as np


def checkTSD(regions, output, path):
    reg_left = []
    reg_right = []

    with open(regions, "r") as file:
        for line in file:
            parts = line.split("||")

            # Get the genomoic region before the first "||" and after the one after the second "||"
            l = parts[0]
            r = parts[2]
            reg_left.append(l)
            reg_right.append(r)
    
    reg_right = [string.replace('\n', '') for string in reg_right]

    tsd = []
    tsd_len = []
    for _ in range(len(reg_left)-1):
        x = len(reg_left[_])
        left = reg_left[_]
        right = reg_right[_]
        while x > 0:
            if left == right:
                tsd.append(left)
                tsd_len.append(len(left))
                x = 0
            elif x == 1:
                tsd.append("No TSD")
                tsd_len.append(0)
            
            if x > 1:
                left = left[1:]
                right = right[:-1]
            x = x - 1

    # Save the tsd of each insertion in a file
    tsd_path = os.path.join(path, output +'_tsd.txt')
    with open(tsd_path, "w") as file:
        for item in tsd:
            file.write(str(item) + "\n")

    tsd_len_path = os.path.join(path, output +'_tsd_len.txt')
    with open(tsd_len_path, "w") as file:
        for item in tsd_len:
            file.write(str(item) + "\n")


    print(f"Target site duplication for each insertion saved to {tsd_path}")
    print(f"Target site duplication length for each insertion saved to {tsd_len_path}")
