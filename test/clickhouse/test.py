import signal
import time
from multiprocessing import Process
import sys
class GracefulKiller:
    kill_now = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self,signum, frame):
        self.kill_now = True

def _j(a):
    killer = GracefulKiller()
    for z in range(0,10):
        for y in range(0,2):
            for x in range(0, 100):
                time.sleep(0.01)
                print(a,':',y, ':',x)
            if killer.kill_now:
                sys.exit()
if __name__ == '__main__':
    jobs = []
    for x in range (0, 10):
        p = Process(target = _j, args = (x,))
        p.start()
        jobs.append(p)

    try:
        for p in jobs:
            p.join()

    except:
        for p in jobs:
            p.terminate()
        sys.exit(10)

    print ("End of the program. I was killed gracefully :)")
