from time import sleep
import requests

def down(url, logfile):
    r = requests.get(url)
    logfile.write(r.text)
    intervel = 5
    sleep(intervel)

def main():
    url = 'http://202.117.80.219/switches/switches.php'
    handle = open('f:\\log.txt', 'a')
    while 1:
        down(url, handle)
    handle.close()
