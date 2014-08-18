import fileinput
import re

def main():
    state = "searching"
    for line in fileinput.input():
        if state == "searching" and re.search("Dependencies Found", line):
            state = "listing"
        elif state == "listing" and len(line.strip()) == 0:
            state = "done"
            break
        elif state == "listing":
            print line.strip()

if __name__ == '__main__':
    main()
