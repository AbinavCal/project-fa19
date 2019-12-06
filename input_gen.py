import string
import random
import re

def main():
    file1 = "\\26TSP.txt"
    file2 = "\\26TSPsol.txt"
    homes = in_writer(file1, 26)
    out_writer(file2, 26, homes)

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def in_writer(FILE_NAME, n):
    read_file = open("test" + FILE_NAME, 'r')
    write_file = open("test\\input" + "\\100.in", 'w')
    
    homes = []
    for i in range(2 * n):
        homes.append(randomString(6))
    homeString = " ".join(homes)
    locString = " ".join([str(num) for num in range(1, n + 1)]) + " " + homeString

    write_file.write(str(3 * n) + "\n")
    write_file.write(str(2 * n) + "\n")
    write_file.write(locString  + "\n")
    write_file.write(homeString + "\n")
    write_file.write("1" + "\n")
    
    rands = []
    for i in range(n):    
        rnum = str(random.randint(20, 60))
        rands.append(rnum)
        line = read_file.readline()
        line = " ".join(re.sub('\\b0\\b', 'x', line).split())
        aux = [rnum if 2*i == num or 2*i+1 == num else "x" for num in range(2 * n)]
        auxString = " ".join(aux)
        lineFin = line + " " + auxString
        write_file.write(lineFin + "\n")
    
    for i in range(2*n):
        aux2 = [rands[num] if i//2 == num  else "x" for num in range(n)]
        aux3 = ["x" for num in range(2 * n)]
        aux2.extend(aux3)
        line = " ".join(aux2)
        write_file.write(line + "\n")

    read_file.close()
    write_file.close()
    return homes

def out_writer(FILE_NAME, n, HOME_LOCS):
    read_file = open("test" + FILE_NAME, 'r')
    write_file = open("test\\output" + "\\100.out", 'w')

    lines = read_file.readlines()
    line = " ".join(" ".join(lines).split())
    write_file.write(line + "\n")
    write_file.write(str(n) + "\n")
    for i in range(n):
        line = str(i + 1) + " " + HOME_LOCS[2 * i] + " " + HOME_LOCS[2 * i + 1]
        write_file.write(line + "\n")

    read_file.close()
    write_file.close()

main()

