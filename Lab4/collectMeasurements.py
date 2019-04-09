from datetime import datetime, timedelta
from urllib.request import urlopen
import math
import scipy.io

OUTPUT_MAT = '/Users/mohiuddi/Desktop/svn/class-perfeval/trunk/hw/l4/' #<--------------------------------- Here specify the folder where you want to save .mat file


''' Process a line of the txt file and return what to save '''
def processLine(line):
  labels = ["ts", "delay", "frequency", "rocof", "PAVM", "PAVA", "PACM", "PACA", "PBVM", "PBVA", "PBCM", "PBCA", "PCVM", "PCVA", "PCCM", "PCCA"]
  tokens = line.split(",")
  
  meas = {}
  
  for i, token in enumerate(tokens):
    meas[labels[i]] = float(token)
    
  Pa = 1e-3*meas["PAVM"]*meas["PACM"]*math.cos(meas["PAVA"] - meas["PACA"])
  Pb = 1e-3*meas["PBVM"]*meas["PBCM"]*math.cos(meas["PBVA"] - meas["PBCA"])
  Pc = 1e-3*meas["PCVM"]*meas["PCCM"]*math.cos(meas["PCVA"] - meas["PCCA"])
  Ptot = Pa + Pb + Pc
  return Ptot # Already multiplied with 1e-3 so the result is given back in KW


''' Read an url and retrieve the content '''
def readurl (url):
  try:
    lines = urlopen(url).readlines()
    print(url, 'successfully read.')
  except:
    print('WARNING: Could not retrieve the following file: '+url)
    lines = []
  return lines

# url example: http://nanotera-stg2.epfl.ch/data/2015/Jan/01/PMU_ID2_DATA-01Jan2015-00h.dat.asc
def composeurl (date):
  url = "http://nanotera-stg2.epfl.ch/data/"
  digits = (date.strftime('%Y'), date.strftime('%b'), date.strftime('%d'), date.strftime('%H'))
  filename = 'PMU_ID2_DATA-' + digits[2]+digits[1]+digits[0]+'-'+digits[3] + 'h.dat.asc'
  url += digits[0] + '/' + digits[1] + '/' + digits[2] + '/' + filename
  return url
  

def main():
  startTime = datetime(2015, 6, 1, 0) #<--------------------------------- Here is to specify the start time (the first time to fetch), here: 1st of March 00:00:00 to 23:59:59.98 
  print('The start time is', startTime)#   To check whether you have a good starting point


  endTime = datetime.now()
  endTime = datetime(2015, 8, 1, 0)#<--------------------------------- Here is to specify the end time (the first time NOT to fetch)
  print('The end time is', endTime)#   To check whether you have a good end point

  
  timeDiff = (endTime - startTime)
  diffHours = timeDiff.days*24 + timeDiff.seconds/3600  #   Check the output (integer number of hours?)
  print('I should grab ' + str(diffHours) + ' files!\n')

  output=[]
  currentTime = startTime
  i = 0
  while currentTime < endTime:
    i += 1
    url = composeurl(currentTime)
    currentTime = currentTime + timedelta(hours=1)
    lines = readurl(url)
    print("[{}/{}] Processing {}\n".format(i, int(diffHours), url))
    #print('Number of lines is ', len(lines))
    if len(lines)<(90000):#  Number of entries is 60minutes*60seconds*50entriesPerSecond = 180000 + additional 19 lines, hence, 180019 in total.
      print('WARNING: The file that corresponds to ', url, ' has small number of lines equal to ', len(lines),'. Putting zeros instead.\n')
      output.append(0)
      continue
    powers = []
    for line in lines:
      dataline = None
      line = line.decode('utf-8')
      if '[data]' in line or '[d]' in line:
        dataline = line.replace('\t\r\r\n', '').replace('[data]\t', '').replace('\t', ',')
      if '[d]' in line:
        dataline = line.replace('\r\r\n', '').replace('\t', ',').replace('[d],', '')
      if dataline is not None:
        datalineCrunched = processLine(dataline) #   Here we get activePower [KW] for one data point
        powers.append(datalineCrunched)
    
    #print('There are ', len(powers), ' values for powers.')
    output.append(sum(powers)/(len(lines)/180019)) # Here we sum up powers [KW], and if number of measurements is between 90000 and 18000 we do the scaling.
  output = [x / 180000 for x in output] # We multiply with 1/180000 h = 20 ms to get [KWh]

    
  scipy.io.savemat(OUTPUT_MAT+'latestActivePowers.mat', {'latestActivePowers': output})
  

if __name__ == "__main__":
  main()
