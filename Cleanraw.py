

with open('raw.txt', 'r') as f:
  lines = f.readlines()
  for line in lines:
    if line[:-1].count('=') > 1:
      No = line[:-1].count('=') 
      while line[-1] in ['\n',',', ' ']:
        line = line[:-1]
      for each in line.split(', '):
        print(each + '|\\')  
    elif len(line) > 1:
      while line[-1] in ['\n',',', ' ']:
        line = line[:-1]
      print(line + '|\\') 
     
  
  
