# logging function library and definition
import logging

def set_logging():
    # procedure for initialising logging to a log file (astroh.log)
    # and to the console
    logging.basicConfig(filename='astroh.log', filemode='w', \
                        format='%(levelname)s %(asctime)s,%(msecs)d %(name)s %(message)s', \
                        datefmt='%H:%M:%S', level=logging.DEBUG)
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
def remove_logging():
    h_list = list(logging.getLogger('').handlers)
#    print("Loggers are currently :", h_list)
    for h in h_list:    
        if h.__class__.__name__ == "FileHandler":
#            print("Dealing with FileHandler")
            h.flush()
            h.close()
            logging.getLogger('').removeHandler(h)
        else:
#            print("Dealing with Stream and other Handlers")
            h.flush()
            logging.getLogger('').removeHandler(h)
  
    h_list = list(logging.getLogger('').handlers)
#    print("Updated view of loggers is :", h_list)
    print("***************Logging Handlers Removed! *******************")
    
    
def handle_exception(code):
    print("Error Raised!  CODE: <{}>".format(code))

# class that defines the data structures to store coefficients for a particular band (wavelength)
# and methods for getting K0, K1 and K

# reading the data function and csv library
import csv

class Reddening:
    
    def __init__(self):
        self.coeffs = self.readData()
        self.kBP = self.createkValue('kBP')
        self.kRP = self.createkValue('kRP')
        self.kG = self.createkValue('kG')
    
    def readData(self):
        coeffs = {}
        try:
            csvfile = open('GBPRP_extcoefs_nov17.csv',newline='')
        except Exception as e:
            handle_exception(e)
            
        c_reader = csv.reader(csvfile, dialect='excel')

        next(c_reader) # this is just stepping past the header row

        kRP    = next(c_reader)
        kBP    = next(c_reader)
        kG     = next(c_reader)

        coeffs['kRP'] = [float(kRP[1]), float(kRP[2]), float(kRP[3]), float(kRP[4]), float(kRP[5]), float(kRP[6]), float(kRP[7])]
        coeffs['kBP'] = [float(kBP[1]), float(kBP[2]), float(kBP[3]), float(kBP[4]), float(kBP[5]), float(kBP[6]), float(kBP[7])]
        coeffs['kG']  = [float(kG[1]), float(kG[2]), float(kG[3]), float(kG[4]), float(kG[5]), float(kG[6]), float(kG[7])]

        csvfile.close()
        logging.info("Coefficients read from file!!!")
        return (coeffs)

    def createkValue(self,band):
        return Reddening.kValue(self,band)

    class kValue:
        def __init__(self, Reddening, band):
            self.band = band
            self.Reddening = Reddening
            self.c_1, self.c_BPRP, self.c_BPRP2, self.c_BPRP3, \
                self.c_A0, self.c_A02, self.c_A0BPRP = \
                Reddening.coeffs[band]
            logging.info("Constant values include: c_1 {}, c_BPRP {}, c_A0 {}"\
                .format(self.c_1, self.c_BPRP, self.c_A0))

        def getK0(self, A0):

            return self.c_1 + A0*(self.c_A0)

        def getK1(self, A0, BPRP):

            return self.c_1 + BPRP * self.c_BPRP \
                + A0*(self.c_A0 + A0*self.c_A02)

        def getK(self, A0, BPRP):

            return self.c_1 + BPRP * (self.c_BPRP \
                + BPRP * (self.c_BPRP2 + BPRP * self.c_BPRP3)) \
                + A0 * (self.c_A0 + A0 * self.c_A02 \
                + BPRP * self.c_A0BPRP)

        def getKBPRP(self, BPRP):

            return self.c_1 + BPRP * (self.c_BPRP \
                + BPRP * (self.c_BPRP2 + BPRP * self.c_BPRP3))

    def getReddening(self, A0, BPRP):
        BPRP_0 = BPRP
        
        KB = Reddening.kValue.getK0(self.kBP,A0)
        KR = Reddening.kValue.getK0(self.kRP,A0)

        BPRP -= (KB-KR)*A0

        KB = Reddening.kValue.getK1(self.kBP, A0, BPRP)
        KR = Reddening.kValue.getK1(self.kRP, A0, BPRP)

        BPRP = (BPRP_0 - (KB-KR)*A0)
        BPRP_S = BPRP_0
        icount = 0
        logging.debug("Prior to loop delta is: {}".format(abs(BPRP-BPRP_S)))
        while (abs(BPRP - BPRP_S) > 0.01):
            
            icount +=1
            logging.debug("IN LOOP ({}), delta is: {}".format(icount,abs(BPRP-BPRP_S)))
            BPRP_S = BPRP
            KB = Reddening.kValue.getK(self.kBP, A0, BPRP)
            KR = Reddening.kValue.getK(self.kRP, A0, BPRP)
            BPRP = (float) (BPRP_0 - (KB-KR)*A0)
        logging.debug("PAST LOOP ({}), delta is: {}".format(icount,abs(BPRP-BPRP_S)))

        return (Reddening.kValue.getK(self.kG, BPRP, A0), Reddening.kValue.getK(self.kBP, BPRP, A0),\
            Reddening.kValue.getK(self.kRP, BPRP, A0),icount)


def main():

    set_logging()
    try:
        HelenRed = Reddening()
    except Exception as e:
        print('Fatal error <{}> whilst trying to create HelenRed object - check that the coefficients file was available'.format(e))
        
    else:
        A0 = 0.3
        BPRP = 0.4
        rG, rBP, rRP, ic = Reddening.getReddening(HelenRed, A0, BPRP)
        print("Reddening values calculated are for \nG:{}\nBP:{}\nRP:{}".format(rG,rBP,rRP))
        print("Iterations within 'while loop' was: {}".format(ic))

        del(HelenRed.kBP)
        del(HelenRed.kRP)
        del(HelenRed.kG)
        del(HelenRed)

    remove_logging()
    
if __name__ == "__main__":

    main()