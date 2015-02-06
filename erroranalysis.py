results = []
with open('out.txt') as inputfile:
    for line in inputfile:
        results.append(line.strip().split(','))

print results[0]
